"""
90_远程探测h5ad.py — 不下载整个文件，直接读远程 h5ad 的元数据
================================================================
用途
    某些数据源（如 CELLxGENE 的数据 CDN）在校园网/国内网络下会被重置，
    但同一批文件在 S3 上是可访问的。这个脚本用 HTTP Range 按块读取
    HDF5 结构，只取几 MB 就能读出：
        - 文件里的 global attrs（schema_version / title / 组学类型）
        - obs 的列名（有没有 cell_type / donor_id / tissue）
        - obs 里某个列有多少种类别、类别名是什么

为什么可行
    HDF5 是"随机访问"格式：超级块 + B-tree + 对象头，不需要顺序读全文件。
    h5py 接受任何带 read/seek/tell 的文件对象，所以套一个走 HTTP Range
    的只读文件对象就能远程打开。

用法
    python 90_远程探测h5ad.py --url <h5ad 直链> [--col cell_type] [--max-cats 200]

    # 批量：给一个 URL 列表文件（每行一个）
    python 90_远程探测h5ad.py --list urls.txt

注意
    - 只适合 HDF5/h5ad，不适合别的格式
    - 每次 seek 未命中缓存都会发一次请求，所以块要开大一点（默认 4 MB）
"""

import argparse
import os
import sys
import urllib.request

BLOCK = 4 << 20          # 4 MB 一块
_opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))


class HttpRangeFile:
    """把 HTTP 上的文件包装成 h5py 能用的只读文件对象。"""

    def __init__(self, url, block=BLOCK, timeout=60):
        self.url = url
        self.block = block
        self.timeout = timeout
        self.pos = 0
        self.cache = {}          # block_index -> bytes
        self.requests = 0
        self.bytes_read = 0
        self.size = self._head()

    def _head(self):
        req = urllib.request.Request(self.url, method="HEAD",
                                     headers={"User-Agent": "Mozilla/5.0"})
        with _opener.open(req, timeout=self.timeout) as r:
            n = r.headers.get("Content-Length")
        if not n:
            raise RuntimeError("服务器没返回 Content-Length，无法做 Range 读取")
        return int(n)

    def _block(self, i):
        if i in self.cache:
            return self.cache[i]
        start = i * self.block
        end = min(start + self.block, self.size) - 1
        req = urllib.request.Request(
            self.url,
            headers={"User-Agent": "Mozilla/5.0",
                     "Range": f"bytes={start}-{end}"},
        )
        with _opener.open(req, timeout=self.timeout) as r:
            data = r.read()
        if len(data) != end - start + 1:
            raise RuntimeError(f"Range 返回长度不对：要 {end-start+1}，得 {len(data)}")
        self.cache[i] = data
        self.requests += 1
        self.bytes_read += len(data)
        return data

    # ---- h5py 需要的接口 ----
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


def _decode(x):
    return x.decode("utf-8", "ignore") if isinstance(x, bytes) else str(x)


def get_categories(h, col):
    """
    取 obs 里某一列的类别名。

    anndata 写 h5ad 时，分类变量有两种存法：
      新版：obs/__categories/<col>
      旧版：obs/<col> 上的 attrs['categories']
    两种都试一下。
    """
    cats = []
    try:
        if "__categories" in h["obs"] and col in h["obs"]["__categories"]:
            cats = [_decode(x) for x in h["obs"]["__categories"][col][:]]
    except Exception:
        cats = []
    if not cats:
        try:
            raw = h["obs"][col].attrs.get("categories", [])
            cats = [_decode(x) for x in raw]
        except Exception:
            cats = []
    return cats


def describe(url, col=None, max_cats=200, brief=False):
    import h5py

    print("=" * 66)
    print(url)
    f = HttpRangeFile(url)
    print(f"  文件大小 {f.size/1e9:.2f} GB")
    with h5py.File(f, "r") as h:
        # 规模
        shape = None
        for key in ("X", "raw"):
            if key in h and hasattr(h[key], "shape") and len(h[key].shape) == 2:
                shape = h[key].shape
                break
        if shape:
            print(f"  矩阵：{shape[0]:,} 细胞 × {shape[1]:,} 基因")

        if brief:
            n_ct = len(get_categories(h, col)) if (col and "obs" in h and col in h["obs"]) else 0
            n_dn = len(get_categories(h, "donor_id")) if "obs" in h and "donor_id" in h["obs"] else 0
            print(f"  cell_type 类别 {n_ct} 个；donor 类别 {n_dn} 个")
            print(f"  [远程请求 {f.requests} 次，共读 {f.bytes_read/1e6:.1f} MB]")
            return dict(url=url, cells=shape[0] if shape else 0,
                        genes=shape[1] if shape else 1, n_cell_type=n_ct,
                        n_donor=n_dn, size=f.size)

        keys = list(h.keys())
        print(f"  顶层结构：{keys}")
        for k in ("schema_version", "title", "organism", "assay",
                  "cell_count", "default_embedding"):
            if k in h.attrs:
                print(f"    {k} = {_decode(h.attrs[k])}")
        if "obs" in h:
            obs_keys = list(h["obs"].keys())
            obs_keys = [k for k in obs_keys if k != "__categories"]
            print(f"  obs 列（{len(obs_keys)}）：{obs_keys}")
            marked = [k for k in ("cell_type", "donor_id", "tissue",
                                  "tissue_ontology_term_id", "assay",
                                  "development_stage", "sex", "disease")
                      if k in obs_keys]
            print(f"    关键列命中：{marked}")
        if col and "obs" in h and col in h["obs"]:
            cats = get_categories(h, col)
            print(f"  obs['{col}'] 共 {len(cats)} 个类别")
            for c in cats[:max_cats]:
                print("     -", c)
    print(f"  [远程请求 {f.requests} 次，共读 {f.bytes_read/1e6:.1f} MB]")
    return dict(url=url, cells=shape[0] if shape else 0,
                genes=shape[1] if shape else 0,
                n_cell_type=len(get_categories(h, col)) if col else 0,
                size=f.size)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--url")
    ap.add_argument("--list")
    ap.add_argument("--col", default="cell_type")
    ap.add_argument("--max-cats", type=int, default=200)
    ap.add_argument("--brief", action="store_true",
                    help="只打印规模统计，不打印列名和类别名（批量扫描用）")
    args = ap.parse_args()

    if args.list:
        urls = [l.strip() for l in open(args.list, encoding="utf-8")
                if l.strip() and not l.startswith("#")]
    elif args.url:
        urls = [args.url]
    else:
        ap.error("要么给 --url，要么给 --list")

    results = []
    for u in urls:
        try:
            r = describe(u, args.col, args.max_cats, args.brief)
            if r:
                results.append(r)
        except Exception as e:
            print(f"  探测失败：{type(e).__name__}: {e}")

    if args.brief and results:
        print("\n" + "=" * 66)
        print(f"{'细胞数':>10} {'基因数':>7} {'cell_type':>9} {'文件GB':>7}  文件")
        for r in sorted(results, key=lambda x: -x["cells"]):
            print(f"{r['cells']:>10,} {r['genes']:>7,} {r['n_cell_type']:>9} "
                  f"{r['size']/1e9:>7.1f}  {r['url'].split('/')[-1]}")


if __name__ == "__main__":
    main()
