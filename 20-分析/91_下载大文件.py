"""
91_下载大文件.py — 断点续传 / 并发分段下载器
================================================
为什么单独写一个：
    1. 需要断点续传（网络断了几次也不至于从头再来）
    2. 需要进度日志（10 GB 不可能盯着终端看）
    3. 本机设了 http_proxy，拉大文件必须绕过代理
    4. 单连接在 S3 上只有 1–3 MB/s，开 8 个连接能到 5 MB/s

用法
    # 单连接（默认，最省事）
    python 91_下载大文件.py --url <直链> --out <目标路径>

    # 并发 8 路分段（大文件推荐）
    python 91_下载大文件.py --url <直链> --out <目标路径> --connections 8

    # 断点续传：直接把同样的命令再跑一遍

行为
    - 单连接：写 <out>.part，按 32 MB 一块向后追加，中断后从 .part 长度续传
    - 并发：把剩余区间切成 N 段，各写 <out>.seg0…segN-1，全部完成后按序拼接，
      删掉分段文件，把 .part 改名为正式文件名
    - 进度同时打印到 stdout

并发模式的分段起点记录在 <out>.meta 里（JSON）。
参数和上次不一样（比如换了 --connections）会自动作废分段重来，不会拼错。
"""

import argparse
import json
import os
import shutil
import sys
import threading
import time

import requests
from requests.adapters import HTTPAdapter, Retry

_UA = {"User-Agent": "Mozilla/5.0"}

# 每线程一个 Session，保持 keep-alive 连接。
# urllib 每次请求都重做 TLS 握手，实测单请求延迟约 1.7 秒，
# 换成 Session 后大文件下载速度提升 3–5 倍。
_tls = threading.local()


def _session():
    s = getattr(_tls, "session", None)
    if s is None:
        s = requests.Session()
        s.trust_env = False                 # 绕开本机 http_proxy
        retry = Retry(total=3, backoff_factor=0.5,
                      status_forcelist=(500, 502, 503, 504))
        s.mount("https://", HTTPAdapter(pool_connections=2, pool_maxsize=2,
                                        max_retries=retry))
        s.mount("http://", HTTPAdapter(pool_connections=2, pool_maxsize=2,
                                       max_retries=retry))
        _tls.session = s
    return s


def remote_size(url):
    r = _session().head(url, timeout=(10, 60), headers=_UA)
    if r.status_code >= 400:
        raise RuntimeError(f"HEAD 返回 {r.status_code}")
    return int(r.headers["Content-Length"])


def fetch(url, start, end, retries=5, timeout=90):
    """
    取 [start, end] 闭区间的一段字节。

    timeout 不要太长：S3 在一次持续大流量之后会把连接挂住，
    设成 90 秒可以让失败的请求尽快重连，而不是傻等 5 分钟。
    """
    for attempt in range(retries):
        try:
            r = _session().get(url, timeout=(10, timeout),
                               headers={**_UA,
                                        "Range": f"bytes={start}-{end}"})
            if r.status_code not in (200, 206):
                raise RuntimeError(f"HTTP {r.status_code}")
            data = r.content
            if len(data) != end - start + 1:
                raise RuntimeError(f"长度不对：要 {end-start+1}，得 {len(data)}")
            return data
        except Exception as e:
            if attempt == retries - 1:
                raise
            print(f"[警告] Range {start}-{end} 第 {attempt+1} 次失败："
                  f"{type(e).__name__}，重试", flush=True)
            time.sleep(2 * (attempt + 1))


# ---------- 单连接 ----------
def download_single(url, out, chunk):
    part = out + ".part"
    total = remote_size(url)
    have = os.path.getsize(part) if os.path.exists(part) else 0
    print(f"[信息] 总大小 {total/1e9:.2f} GB，已有 {have/1e9:.2f} GB")

    if have >= total:
        print("[信息] 已完整")
        os.replace(part, out)
        return

    t0, base = time.time(), have
    with open(part, "ab" if have else "wb") as fh:
        while have < total:
            end = min(have + chunk, total) - 1
            data = fetch(url, have, end)
            if not data:
                print("[错误] 收到空数据，停止；.part 已保留")
                sys.exit(1)
            fh.write(data)
            fh.flush()
            have += len(data)
            dt = time.time() - t0
            speed = (have - base) / dt / 1e6 if dt else 0
            eta = (total - have) / (speed * 1e6) if speed else 0
            print(f"[进度] {have/total*100:5.2f}%  {have/1e9:6.2f}/{total/1e9:.2f} GB  "
                  f"{speed:5.1f} MB/s  剩余约 {eta/60:.1f} 分钟")

    os.replace(part, out)
    print(f"[完成] {out}  {os.path.getsize(out)/1e9:.2f} GB")


# ---------- 并发分段 ----------
def download_parallel(url, out, nconn):
    part = out + ".part"
    meta_path = out + ".meta"
    total = remote_size(url)
    base = os.path.getsize(part) if os.path.exists(part) else 0

    if base >= total:
        print("[信息] 已完整")
        os.replace(part, out)
        return

    state = {"url": url, "total": total, "base": base, "n": nconn}
    if os.path.exists(meta_path):
        old = json.load(open(meta_path))
        if old != state:
            print("[信息] 参数或文件变了，作废旧分段重来")
            for f in _seg_files(out):
                os.remove(f)
    json.dump(state, open(meta_path, "w"))

    remain = total - base
    seg_len = remain // nconn
    segs = []
    s = base
    for i in range(nconn):
        e = total - 1 if i == nconn - 1 else s + seg_len - 1
        segs.append((i, s, e))
        s = e + 1

    print(f"[信息] 总 {total/1e9:.2f} GB，已有 {base/1e9:.2f} GB，"
          f"剩 {remain/1e9:.2f} GB 切成 {nconn} 段")

    done = [0] * nconn
    errors = []

    def worker(i, s_, e_):
        seg = f"{out}.seg{i}"
        have = os.path.getsize(seg) if os.path.exists(seg) else 0
        done[i] = have
        need = e_ - s_ + 1
        if have >= need:
            return
        with open(seg, "ab") as fh:
            cur = have
            step = 8 << 20
            while cur < need:
                a = s_ + cur
                b = min(a + step, e_) 
                try:
                    data = fetch(url, a, b)
                except Exception as ex:
                    errors.append(f"段{i}: {type(ex).__name__}")
                    return
                if not data:
                    errors.append(f"段{i}: 空数据")
                    return
                fh.write(data)
                fh.flush()
                cur += len(data)
                done[i] = cur

    ths = [threading.Thread(target=worker, args=a, daemon=True) for a in segs]
    t0 = time.time()
    last = 0
    for t in ths:
        t.start()
    while any(t.is_alive() for t in ths):
        time.sleep(3)
        cur = sum(done)
        dt = time.time() - t0
        speed = (cur - last) / dt / 1e6 if dt else 0
        got = base + cur
        eta = (total - got) / (speed * 1e6) if speed > 0 else 0
        print(f"[进度] {got/total*100:5.2f}%  {got/1e9:6.2f}/{total/1e9:.2f} GB  "
              f"{speed:5.1f} MB/s  剩余约 {eta/60:.1f} 分钟")
        last = cur
    for t in ths:
        t.join()

    if errors:
        print(f"[错误] 有分段失败：{errors}；分段文件已保留，重跑可续传")
        sys.exit(1)

    # 校验每段长度
    for i, s_, e_ in segs:
        seg = f"{out}.seg{i}"
        if os.path.getsize(seg) != e_ - s_ + 1:
            print(f"[错误] 段 {i} 长度不对，保留文件以便重跑")
            sys.exit(1)

    # 拼接
    print("[信息] 拼接分段……")
    with open(part, "ab") as dst:
        for i, s_, e_ in segs:
            with open(f"{out}.seg{i}", "rb") as src:
                shutil.copyfileobj(src, dst, 1 << 24)

    if os.path.getsize(part) != total:
        print(f"[错误] 拼接后大小 {os.path.getsize(part)} ≠ {total}")
        sys.exit(1)

    for f in _seg_files(out):
        os.remove(f)
    os.remove(meta_path)
    os.replace(part, out)
    print(f"[完成] {out}  {os.path.getsize(out)/1e9:.2f} GB  "
          f"用时 {(time.time()-t0)/60:.1f} 分钟")


def _seg_files(out):
    d = os.path.dirname(os.path.abspath(out))
    stem = os.path.basename(out) + ".seg"
    return [os.path.join(d, f) for f in os.listdir(d) if f.startswith(stem)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--url", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--chunk-mb", type=int, default=32,
                    help="单连接模式每块大小（MB）")
    ap.add_argument("--connections", type=int, default=1,
                    help="并发连接数；>1 时用分段模式（大文件推荐 8）")
    a = ap.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    if a.connections > 1:
        download_parallel(a.url, a.out, a.connections)
    else:
        download_single(a.url, a.out, a.chunk_mb << 20)


if __name__ == "__main__":
    main()
