"""
92_远程切片.py — 不下载整个 h5ad，只取需要的那部分细胞
===============================================================
为什么需要它
    Tabula Sapiens 全量 h5ad 是 10.02 GB。实测网络：
      - 单连接约 0.24–1.2 MB/s，8 并发反而更慢（0.3 MB/s，S3 惩罚并发）
      - 持续大流量后 S3 会直接挂住连接（两次分别卡在 28% 和 62%）
    整包下载不可行。

    但本项目不需要全部 483,152 个细胞 —— 算的是「每个细胞类型的平均
    基因年龄谱」，每个类型 150–300 个细胞就够（TAI 是加权均值，细胞间
    SD 约 0.15 时，n=200 的标准误只有 0.011，远小于类型间差异）。

关键：raw.X 是 chunked + gzip 的 CSR
    实测 raw.X/data 的布局：
        shape (1,270,410,953,)  dtype float32
        chunks (77540,)        compression gzip    shuffle False
    也就是每 77,540 个非零元一个 chunk，整个数组 16,385 个 chunk。
    选 28,017 个细胞（占 5.8%）只需要其中约 1,000 个 chunk，
    折合下载量约 250 MB —— 只有整包的 2.5%。

做法
    1. 远程读 obs（细胞类型 / 供者 / 组织 / 区室）与 raw.X 的 indptr
    2. 每种细胞类型取最多 N 个细胞（不足则全取），优先从同类型连续段取，
       使选中的行落在少数连续区间里
    3. 把选中行换算成非零元区间，再换算成 chunk 段；
       每段只读一次（h5py 负责 gzip 解压），并对后续字节块做后台预取
    4. 落盘为 ts_slice.npz + ts_slice.obs.tsv + ts_slice.genes.txt

接口与完整 h5ad 完全一致（ts_io.SliceSource），下游 01/03/04/05/06 不用改。

用法
    python 92_远程切片.py --scan                    # 只扫代价
    python 92_远程切片.py --per-type 200 --dry-run  # 只规划
    python 92_远程切片.py --per-type 200            # 真下载

注意
    第一次运行要远程读 obs 与 indptr（约 3 分钟），之后走本地缓存
    10-数据/_ts_remote_cache.npz，秒开。
"""

import argparse
import json
import os
import re
import sys
import threading
import time
import urllib.request
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import requests
from requests.adapters import HTTPAdapter, Retry

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(os.path.dirname(HERE), "10-数据")

DEFAULT_URL = ("https://s3.us-west-2.amazonaws.com/cellxgene-data-public"
               "/cell-census/2023-05-15/h5ads"
               "/53d208b0-2cfd-4366-9866-c3c6114081bc.h5ad")

# ----------------------------------------------------------------- 会话池
# 为什么必须用 requests 而不是 urllib
#     切片要发成百上千个 Range 请求。urllib 每次都新建 TCP + TLS 连接，
#     实测单段延迟约 1.7 秒。requests 的 Session 复用 keep-alive 连接，
#     握手只做一次。每线程一个 Session（Session 不是线程安全的）。
_tls = threading.local()


def _session():
    s = getattr(_tls, "session", None)
    if s is None:
        s = requests.Session()
        s.trust_env = False                 # 别走本机 http_proxy
        retry = Retry(total=2, backoff_factor=0.3,
                      status_forcelist=(500, 502, 503, 504))
        adapter = HTTPAdapter(pool_connections=2, pool_maxsize=2,
                              max_retries=retry)
        s.mount("https://", adapter)
        s.mount("http://", adapter)
        _tls.session = s
    return s


# ----------------------------------------------------------------- HTTP
class RangeFile:
    """
    把 HTTP 上的文件包装成 h5py 能用的只读文件对象。

    特点
      - 块缓存：h5py 会 seek 到很多分散位置，按块缓存避免重复下载
      - 后台预取：请求第 i 块时，同时把 i+1…i+prefetch 块拉起来，
        串行的 h5py 因此能吃满网络带宽
      - 带重试与退避
    """

    def __init__(self, url, block=1 << 20, timeout=45, retries=6,
                 prefetch=24, cache_blocks=3000, verbose=False):
        self.url = url
        self.block = block
        self.timeout = timeout
        self.retries = retries
        self.prefetch = prefetch
        self.cache_blocks = cache_blocks
        self.verbose = verbose
        self.pos = 0
        self.cache = OrderedDict()
        self._inflight = set()
        self._lock = threading.Lock()
        self.requests = 0
        self.bytes_read = 0
        self.size = self._head()
        self._pool = ThreadPoolExecutor(max_workers=6)

    # ---------- 底层 ----------
    def _head(self):
        last = None
        for k in range(self.retries):
            try:
                r = _session().head(self.url, timeout=(10, self.timeout),
                                    headers={"User-Agent": "Mozilla/5.0"})
                if r.status_code >= 400:
                    raise RuntimeError(f"HTTP {r.status_code}")
                n = r.headers.get("Content-Length")
                if n:
                    return int(n)
                raise RuntimeError("没有 Content-Length")
            except Exception as e:            # noqa: BLE001
                last = e
                time.sleep(1.5 * (k + 1))
        raise RuntimeError(f"HEAD 失败：{last}")

    def _fetch(self, start, end):
        last = None
        for k in range(self.retries):
            try:
                r = _session().get(
                    self.url,
                    headers={"User-Agent": "Mozilla/5.0",
                             "Range": f"bytes={start}-{end}"},
                    timeout=(10, self.timeout))
                if r.status_code not in (200, 206):
                    raise RuntimeError(f"HTTP {r.status_code}")
                data = r.content
                if len(data) != end - start + 1:
                    raise RuntimeError(
                        f"长度不对：要 {end-start+1}，得 {len(data)}")
                with self._lock:
                    self.requests += 1
                    self.bytes_read += len(data)
                if self.verbose:
                    print(f"      [net] {start}-{end} "
                          f"{len(data)/1e6:.2f} MB "
                          f"（累计 {self.bytes_read/1e6:.0f} MB，"
                          f"{self.requests} 次）", flush=True)
                return data
            except Exception as e:            # noqa: BLE001
                last = e
                time.sleep(1.0 * (k + 1))
        raise RuntimeError(f"Range {start}-{end} 重试 {self.retries} 次仍失败：{last}")

    def read_range(self, start, end):
        """直接读 [start, end) 字节，不走块缓存。"""
        if end <= start:
            return b""
        return self._fetch(start, end - 1)

    # ---------- 块缓存 + 预取 ----------
    def _block(self, i):
        with self._lock:
            if i in self.cache:
                self.cache.move_to_end(i)
                return self.cache[i]
        for j in range(i + 1, i + 1 + self.prefetch):
            self._prefetch(j)
        data = self._fetch(i * self.block,
                           min((i + 1) * self.block, self.size) - 1)
        with self._lock:
            self.cache[i] = data
            while len(self.cache) > self.cache_blocks:
                self.cache.popitem(last=False)
        return data

    def _prefetch(self, j):
        with self._lock:
            if j in self.cache or j in self._inflight:
                return
            if j * self.block >= self.size:
                return
            self._inflight.add(j)
        self._pool.submit(self._do_prefetch, j)

    def _do_prefetch(self, j):
        try:
            data = self._fetch(j * self.block,
                               min((j + 1) * self.block, self.size) - 1)
            with self._lock:
                self.cache[j] = data
                while len(self.cache) > self.cache_blocks:
                    self.cache.popitem(last=False)
        except Exception:                       # noqa: BLE001
            pass
        finally:
            with self._lock:
                self._inflight.discard(j)

    # ---------- h5py 需要的接口 ----------
    def read(self, n=-1):
        if n is None or n < 0:
            n = self.size - self.pos
        n = min(n, self.size - self.pos)
        out = bytearray()
        while n > 0:
            i, off = divmod(self.pos, self.block)
            blk = self._block(i)
            take = min(n, len(blk) - off)
            out += blk[off:off + take]
            self.pos += take
            n -= take
        return bytes(out)

    def seek(self, off, whence=0):
        if whence == 0:
            self.pos = off
        elif whence == 1:
            self.pos += off
        else:
            self.pos = self.size + off
        self.pos = max(0, min(self.pos, self.size))
        return self.pos

    def tell(self):
        return self.pos

    def seekable(self):
        return True

    def readable(self):
        return True

    def writable(self):
        return False

    def __len__(self):
        return self.size


# ----------------------------------------------------------------- 规划
def find_runs(codes):
    """把 codes 切成「同值连续段」。返回 (起始行, 结束行(不含), 类型)。"""
    n = len(codes)
    if n == 0:
        return (np.array([], dtype=np.int64),) * 3
    chg = np.flatnonzero(np.diff(codes) != 0) + 1
    starts = np.r_[0, chg]
    ends = np.r_[chg, n]
    return starts, ends, codes[starts]


def plan_rows(codes, n_types, per_type):
    """
    为每个类型挑 per_type 个细胞（不足则全取），并让选中的行尽量紧凑。

    做法：把一个类型的细胞按「连续段长度」从长到短排序，优先从长段里
    取，且只取段的前缀。这样选中的行天然落在少数连续区间里。

    返回 (selected_rows_sorted, cover)
    """
    starts, ends, rtypes = find_runs(codes)
    rlen = ends - starts

    tot = np.bincount(codes[codes >= 0], minlength=n_types)[:n_types]
    target = np.minimum(per_type, tot)

    order = np.argsort(-rlen, kind="stable")      # 长段优先
    by_type = [[] for _ in range(n_types)]
    for ri in order:
        t = int(rtypes[ri])
        if 0 <= t < n_types:
            by_type[t].append(int(ri))

    picked = []
    cover = np.zeros(n_types, dtype=np.int64)
    for t in range(n_types):
        need = int(target[t])
        if need <= 0:
            continue
        got = 0
        for ri in by_type[t]:
            if got >= need:
                break
            s, e = int(starts[ri]), int(ends[ri])
            take = min(e - s, need - got)
            picked.append((s, s + take))
            got += take
        cover[t] = got

    if not picked:
        return np.array([], dtype=np.int64), cover
    rows = np.concatenate([np.arange(s, e) for s, e in picked])
    return np.unique(rows), cover


def row_ranges(indptr, rows):
    """选中行 → 非零元区间列表 [(a, b), ...]（已按 a 排序，跳过全零行）。"""
    rows = np.asarray(rows, dtype=np.int64)
    a = indptr[rows]
    b = indptr[rows + 1]
    keep = b > a
    return list(zip(a[keep].tolist(), b[keep].tolist()))


def chunk_plan(ranges, chunk_len, n_elem, max_group=40):
    """
    把非零元区间映射成「chunk 段」读取计划。

    返回 [(chunk_c0, chunk_c1, [(range_index, ...), ...]), ...]
      一个 chunk 段 = 一次 dset[a:b] 读取，涵盖 chunk_c0..chunk_c1
    每个 chunk 只会被读一次（合并相邻需求）。
    """
    groups = []
    i = 0
    n = len(ranges)
    while i < n:
        c0 = ranges[i][0] // chunk_len
        c1 = (ranges[i][1] - 1) // chunk_len
        j = i + 1
        while j < n:
            cj = (ranges[j][1] - 1) // chunk_len
            if cj - c0 > max_group:          # 一组最多 max_group 个 chunk
                break
            c1 = max(c1, cj)
            j += 1
        groups.append((c0, c1, list(range(i, j))))
        i = j
    return groups


def count_chunks(ranges, chunk_len):
    """需要的不同 chunk 数（用于估算下载量）。"""
    s = set()
    for a, b in ranges:
        s.update(range(a // chunk_len, (b - 1) // chunk_len + 1))
    return len(s)


# ----------------------------------------------------------------- 远程元数据
def remote_meta(url, group_key, extra_cols, block=1 << 20, verbose=False):
    """
    读 obs（codes + categories）与 raw.X 的 indptr / dtype / 形状 / chunk 布局。
    结果缓存到 10-数据/_ts_remote_cache.npz，第二次起秒开。
    """
    import h5py

    os.makedirs(DATA_DIR, exist_ok=True)
    cache_path = os.path.join(DATA_DIR, "_ts_remote_cache.npz")

    want = list(dict.fromkeys([group_key] + extra_cols))
    if os.path.exists(cache_path):
        try:
            z = np.load(cache_path, allow_pickle=True)
            if str(z["url"]) == url and str(z["group_key"]) == group_key \
                    and list(z["cols"]) == want:
                print(f"[缓存] {cache_path}")
                return dict(n_obs=int(z["n_obs"]), n_vars=int(z["n_vars"]),
                            codes=z["codes"], cats=list(z["cats"]),
                            indptr=z["indptr"],
                            data_dtype=np.dtype(str(z["data_dtype"])),
                            idx_dtype=np.dtype(str(z["idx_dtype"])),
                            chunk_len=int(z["chunk_len"]),
                            obs={c: z["obs_" + c] for c in want},
                            symbols=z["symbols"], rf=None,
                            cache_hit=True)
        except Exception as e:                     # noqa: BLE001
            print(f"[缓存] 读不了（{type(e).__name__}），重新远程读")

    t0 = time.time()
    rf = RangeFile(url, block=block, verbose=verbose)
    print(f"远程文件 {rf.size/1e9:.2f} GB")

    with h5py.File(rf, "r") as h:
        r = h["raw"]["X"]
        n_obs, n_vars = (int(x) for x in r.attrs["shape"])
        print(f"矩阵 {n_obs:,} × {n_vars:,}  [{time.time()-t0:.0f}s]", flush=True)

        indptr = r["indptr"][:].astype(np.int64)
        print(f"indptr 非零元 {indptr[-1]:,}  [{time.time()-t0:.0f}s]", flush=True)

        data_dtype = r["data"].dtype
        idx_dtype = r["indices"].dtype
        chunk_len = int(r["data"].chunks[0])
        print(f"raw.X: data {data_dtype} chunks {r['data'].chunks} "
              f"comp {r['data'].compression}；indices {idx_dtype} "
              f"chunks {r['indices'].chunks}  [{time.time()-t0:.0f}s]", flush=True)

        obs = {}
        for c in want:
            node = h["obs"][c]
            if isinstance(node, h5py.Group):
                cs = [x.decode() if isinstance(x, bytes) else str(x)
                      for x in node["categories"][:]]
                cc = node["codes"][:].astype(np.int64)
                obs[c] = np.array([cs[i] if 0 <= i < len(cs) else "unknown"
                                   for i in cc], dtype=object)
            else:
                vals = node[:]
                obs[c] = np.array([v.decode() if isinstance(v, bytes) else v
                                   for v in vals], dtype=object)
        print(f"obs {len(want)} 列  [{time.time()-t0:.0f}s]", flush=True)

        node = h["obs"][group_key]
        cats = ([x.decode() if isinstance(x, bytes) else str(x)
                 for x in node["categories"][:]]
                if isinstance(node, h5py.Group) else [])
        codes = (node["codes"][:].astype(np.int64)
                 if isinstance(node, h5py.Group)
                 else np.arange(n_obs, dtype=np.int64))

        # ---- var：结构不标准，必须自适应 ----
        # 不要用 list(h["var"].keys()) 遍历 —— 这个 h5ad 的 var 条目极多，
        # 逐个读 header 会把远程请求拖到十几分钟。只按候选名定点读。
        symbols = None
        for cand in ("feature_name", "gene_symbols", "symbol",
                     "feature_symbol", "gene_symbol"):
            try:
                o = h["var"][cand]
            except KeyError:
                continue
            print(f"  var/{cand}：{type(o).__name__}", flush=True)
            try:
                if isinstance(o, h5py.Dataset):
                    symbols = np.array(
                        [re.sub(r"_ENSG\d+$", "",
                                x.decode() if isinstance(x, bytes) else str(x))
                         for x in o[:]], dtype=object)
                elif isinstance(o, h5py.Group):
                    if "categories" in o and "codes" in o:
                        cs = [x.decode() if isinstance(x, bytes) else str(x)
                              for x in o["categories"][:]]
                        cc = o["codes"][:]
                        symbols = np.array(
                            [re.sub(r"_ENSG\d+$", "",
                                    cs[i] if 0 <= i < len(cs) else "")
                             for i in cc], dtype=object)
                    else:
                        for sub in list(o.keys()):
                            so = o[sub]
                            if isinstance(so, h5py.Dataset) \
                                    and so.size == n_vars:
                                symbols = np.array(
                                    [re.sub(r"_ENSG\d+$", "",
                                            x.decode() if isinstance(x, bytes)
                                            else str(x))
                                     for x in so[:]], dtype=object)
                                break
                if symbols is not None and len(symbols) == n_vars:
                    print(f"  → 基因符号取自 var/{cand}，共 {len(symbols):,}",
                          flush=True)
                    break
                symbols = None
            except Exception as ex:                 # noqa: BLE001
                print(f"  → 读取失败：{type(ex).__name__}: {ex}", flush=True)
                symbols = None
        if symbols is None:
            raise RuntimeError("在 var 里找不到基因符号（试了 feature_name 等）")

    print(f"远程读完，共 {rf.bytes_read/1e6:.1f} MB / {rf.requests} 次请求 "
          f"/ {time.time()-t0:.0f}s")

    np.savez_compressed(
        cache_path, url=np.array(url), group_key=np.array(group_key),
        cols=np.array(want), n_obs=n_obs, n_vars=n_vars,
        codes=codes, cats=np.array(cats, dtype=object), indptr=indptr,
        data_dtype=np.array(str(data_dtype)), idx_dtype=np.array(str(idx_dtype)),
        chunk_len=chunk_len, symbols=symbols,
        **{f"obs_{c}": obs[c] for c in want})
    print(f"已缓存 → {cache_path}")

    return dict(n_obs=n_obs, n_vars=n_vars, codes=codes, cats=cats,
                indptr=indptr, data_dtype=data_dtype, idx_dtype=idx_dtype,
                chunk_len=chunk_len, obs=obs, symbols=symbols, rf=rf,
                cache_hit=False)


# ----------------------------------------------------------------- 主流程
def fetch_by_chunk(hf, dset, ranges, dtype, t0, label):
    """
    按 chunk 粒度读出 ranges 覆盖的元素，顺序与 ranges 一致。

    返回一个数组，长度 = sum(b-a)。
    """
    CH = int(dset.chunks[0])
    N = dset.shape[0]
    total = sum(b - a for a, b in ranges)
    out = np.empty(total, dtype=dtype)

    groups = chunk_plan(ranges, CH, N)
    need_chunks = set()
    for c0, c1, _ in groups:
        need_chunks.update(range(c0, c1 + 1))
    print(f"  {label}: {len(ranges):,} 段 → {len(groups):,} 次读、"
          f"{len(need_chunks):,} 个 chunk，{total:,} 元素", flush=True)

    pos = 0
    for gi, (c0, c1, members) in enumerate(groups):
        a0 = c0 * CH
        a1 = min((c1 + 1) * CH, N)
        arr = dset[a0:a1]
        for k in members:
            a, b = ranges[k]
            out[pos:pos + (b - a)] = arr[a - a0: b - a0]
            pos += b - a
        if (gi + 1) % 20 == 0 or gi + 1 == len(groups):
            el = time.time() - t0
            print(f"    {gi+1:,}/{len(groups):,} 读  {pos:,}/{total:,} 元素  "
                  f"{el:.0f}s", flush=True)
    if pos != total:
        raise RuntimeError(f"填充不完整：{pos} != {total}")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--url", default=DEFAULT_URL)
    ap.add_argument("--out", default=os.path.join(DATA_DIR, "ts_slice"))
    ap.add_argument("--group-key", default="cell_ontology_class")
    ap.add_argument("--per-type", type=int, default=200,
                    help="每个细胞类型目标细胞数（不足则全取）")
    ap.add_argument("--max-group-chunks", type=int, default=40,
                    help="一次读取最多跨多少个 chunk")
    ap.add_argument("--scan", action="store_true",
                    help="只扫不同 per_type 的下载代价，不下载")
    ap.add_argument("--dry-run", action="store_true", help="只规划，不下载")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    t_start = time.time()
    print("=" * 70)
    print("远程切片  Tabula Sapiens")
    print("=" * 70)

    M = remote_meta(args.url, args.group_key,
                    ["donor_id", "tissue", "compartment", "assay", "sex"],
                    verbose=args.verbose)
    codes, indptr = M["codes"], M["indptr"]
    n_types = len(M["cats"])
    CH = M["chunk_len"]
    print(f"分组列 {args.group_key}：{n_types} 类；chunk {CH:,} 元素；"
          f"data {M['data_dtype']} + indices {M['idx_dtype']}")

    tot_per_type = np.bincount(codes[codes >= 0], minlength=n_types)[:n_types]
    short = int((tot_per_type < args.per_type).sum())
    print(f"类型规模：{short} 类全库不足 {args.per_type} 个细胞（将全取），"
          f"其余按 {args.per_type} 抽样")

    if args.scan:
        print("\n" + "-" * 78)
        print(f"{'per_type':>9} {'细胞数':>9} {'非零元':>13} {'chunk数':>8} "
              f"{'读取组数':>9} {'估计GB':>8}")
        print("-" * 78)
        for pt in (100, 150, 200, 300):
            rr, _ = plan_rows(codes, n_types, pt)
            rg = row_ranges(indptr, rr)
            nch = count_chunks(rg, CH)
            grp = chunk_plan(rg, CH, indptr[-1])
            gb = nch * CH * 8 / 1e9
            print(f"{pt:>9} {len(rr):>9,} {sum(b-a for a,b in rg):>13,} "
                  f"{nch:>8,} {len(grp):>9,} {gb:>8.2f}")
        print("\n注：估计GB 按「解压后原始大小 × 8 字节」算，"
              "实际网络流量因 gzip 压缩约为其 1/3。")
        return

    rows, cover = plan_rows(codes, n_types, args.per_type)
    ranges = row_ranges(indptr, rows)
    nch = count_chunks(ranges, CH)
    groups = chunk_plan(ranges, CH, indptr[-1])

    print("-" * 70)
    print(f"选中 {len(rows):,} 细胞（{len(rows)/M['n_obs']*100:.1f}%），"
          f"{len(ranges):,} 个非零元区间")
    print(f"需要 {nch:,} 个 chunk（全部 {int(indptr[-1])//CH:,} 个的 "
          f"{nch/(int(indptr[-1])//CH)*100:.1f}%），{len(groups):,} 次读")

    if args.dry_run:
        print("\n[dry-run] 到此为止")
        return

    import h5py
    from scipy import sparse

    rf = M["rf"] or RangeFile(args.url, verbose=args.verbose)
    with h5py.File(rf, "r") as h:
        r = h["raw"]["X"]
        print("\n== raw.X.data ==")
        data_buf = fetch_by_chunk(h, r["data"], ranges, M["data_dtype"],
                                  t_start, "data")
        print("\n== raw.X.indices ==")
        idx_buf = fetch_by_chunk(h, r["indices"], ranges, M["idx_dtype"],
                                 t_start, "indices")

    # ---------- 组装本地稀疏矩阵 ----------
    counts = indptr[rows + 1] - indptr[rows]
    new_indptr = np.zeros(len(rows) + 1, dtype=np.int64)
    new_indptr[1:] = np.cumsum(counts)
    Xs = sparse.csr_matrix(
        (data_buf.astype(np.float32), idx_buf.astype(np.int32), new_indptr),
        shape=(len(rows), M["n_vars"]))
    print(f"\n本地矩阵 {Xs.shape}  非零 {Xs.nnz:,}  "
          f"密度 {Xs.nnz/(Xs.shape[0]*Xs.shape[1])*100:.2f}%  "
          f"max={Xs.data.max():.0f}  min={Xs.data.min():.0f}")
    if not np.isfinite(Xs.data).all():
        die("矩阵里出现 NaN/Inf，说明字节偏移或解压有问题")

    # ---------- 落盘 ----------
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    sparse.save_npz(args.out + ".npz", Xs)

    obs_cols = list(M["obs"].keys())
    with open(args.out + ".obs.tsv", "w", encoding="utf-8") as f:
        f.write("cell_index\t" + "\t".join(obs_cols) + "\n")
        for i, ri in enumerate(rows):
            f.write(f"{ri}\t" + "\t".join(str(M["obs"][c][ri])
                                          for c in obs_cols) + "\n")

    with open(args.out + ".genes.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(M["symbols"]))

    info = dict(url=args.url, group_key=args.group_key,
                per_type=args.per_type, chunk_len=CH,
                n_cells=int(len(rows)), n_genes=int(M["n_vars"]),
                nnz=int(Xs.nnz),
                frac_of_cells=round(len(rows) / M["n_obs"], 4),
                frac_of_nnz=round(Xs.nnz / int(indptr[-1]), 4),
                n_chunks_needed=int(nch),
                bytes_downloaded=int(rf.bytes_read),
                elapsed_s=round(time.time() - t_start, 1))
    with open(args.out + ".json", "w", encoding="utf-8") as f:
        json.dump(info, f, ensure_ascii=False, indent=2)

    print("\n" + "=" * 70)
    print(f"完成：{time.time()-t_start:.0f}s，远程共读 "
          f"{rf.bytes_read/1e6:.0f} MB（{rf.requests} 次请求）")
    for suf in (".npz", ".obs.tsv", ".genes.txt", ".json"):
        p = args.out + suf
        print(f"  {os.path.basename(p):<24} {os.path.getsize(p)/1e6:>8.1f} MB")


if __name__ == "__main__":
    main()
