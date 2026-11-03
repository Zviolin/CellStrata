"""
ts_io.py — 直接读写 Tabula Sapiens h5ad 的底层接口
====================================================
为什么不直接用 anndata / scanpy
    1. 这份 10 GB 的 h5ad 如果 read_h5ad 全量载入，稀疏矩阵约 1270 M 个非零元，
       光 X 就 ~10 GB，再叠一份 raw 就 20 GB，普通笔记本直接爆掉。
    2. Tabula Sapiens 的 var index 是 Ensembl ID（ENSG...），基因符号在
       var['feature_name'] 里，而且 2.2% 的符号带 _ENSG 后缀。anndata 不会
       帮你处理，scanpy 的 QC / normalize 也会默认按 var_names 走，容易静默出错。
    3. X 已经是 log 归一化过的，raw.X 才是原始计数。用 anndata 时很容易
       不小心对已经归一化的 X 再来一次 normalize_total + log1p。

    所以干脆自己开 h5py：读 obs 不碰矩阵，读矩阵按块 Range 读，内存可控。

本文件只解决"把 Tabula Sapiens 这个特定 h5ad 读出来"，不做任何分析。

关键的存储约定（已实测确认）
    f['X']            csr，attrs['shape'] = (483152, 58604)，值是 log 归一化后的
    f['raw']['X']     csr，同样形状，值是原始 UMI 计数（整数）
    f['obs'][col]     普通 Dataset（数值/布尔/字符串）或 Group（categories + codes）
    f['var'].attrs['_index'] = 'ensemblid'
    f['var']['feature_name']  分类变量，值是基因符号
    f['obsm']['X_umap']       (n_obs, 2)
"""

import re

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

# 形如 "U4_ENSG00000274197" 的重复符号：把尾巴去掉还原成 "U4"
_ENSG_SUFFIX = re.compile(r"_ENSG\d+$")


def _decode(arr):
    """把 h5py 读出来的 bytes 数组转成 str 数组。"""
    a = np.asarray(arr)
    if a.dtype.kind in ("S", "O"):
        return np.array([x.decode("utf-8", "ignore") if isinstance(x, bytes)
                         else str(x) for x in a], dtype=object)
    return a


class TabulaSapiens:
    """Tabula Sapiens h5ad 的只读视图。"""

    def __init__(self, path):
        self.path = path
        self.f = h5py.File(path, "r")
        self.n_obs, self.n_vars = (int(x) for x in self.f["X"].attrs["shape"])
        self._symbols = None
        self._obs_index = None

    # ---------- 生命周期 ----------
    def close(self):
        self.f.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    # ---------- 元数据 ----------
    @property
    def title(self):
        try:
            t = self.f["uns"]["title"][()]
            return t.decode("utf-8", "ignore") if isinstance(t, bytes) else str(t)
        except Exception:
            return "(无标题)"

    def obs_columns(self):
        return [k for k in self.f["obs"].keys() if k != "__categories"]

    def obs_index(self):
        """细胞名（barcode）数组。"""
        if self._obs_index is None:
            name = self.f["obs"].attrs.get("_index", "_index")
            if isinstance(name, bytes):
                name = name.decode()
            node = self.f["obs"][name]
            self._obs_index = _decode(node[:])
        return self._obs_index

    def obs(self, col):
        """
        取 obs 的一列。
        分类变量返回 str 数组；数值列返回原 dtype。
        """
        if col not in self.f["obs"]:
            raise KeyError(f"obs 里没有 '{col}'。现有列：{self.obs_columns()}")
        node = self.f["obs"][col]

        if isinstance(node, h5py.Group):
            # 分类变量：categories + codes
            cats = _decode(node["categories"][:])
            codes = node["codes"][:]
            out = np.empty(len(codes), dtype=object)
            valid = codes >= 0
            out[valid] = cats[codes[valid]]
            # codes == -1 表示缺失。填成字符串而不是 None，
            # 否则后面转成 pd.Categorical 会得到 -1 的组号，索引直接错位。
            out[~valid] = "unknown"
            return out

        arr = node[:]
        if arr.dtype.kind in ("S", "O"):
            return _decode(arr)
        return arr

    def obs_frame(self, cols):
        """一次取多列，拼成 DataFrame。"""
        return pd.DataFrame({c: self.obs(c) for c in cols}, index=self.obs_index())

    def symbols(self):
        """
        基因符号数组（长度 = n_vars）。

        去掉了重复符号的 _ENSG 后缀，因此可能有重名 —— 这是故意的，
        对年龄分层求和没有影响。
        """
        if self._symbols is None:
            node = self.f["var"]["feature_name"]
            if isinstance(node, h5py.Group):
                cats = _decode(node["categories"][:])
                codes = node["codes"][:]
                raw = np.empty(len(codes), dtype=object)
                valid = codes >= 0
                raw[valid] = cats[codes[valid]]
                raw[~valid] = ""
            else:
                raw = _decode(node[:])
            self._symbols = np.array([_ENSG_SUFFIX.sub("", str(s)) for s in raw],
                                     dtype=object)
        return self._symbols

    def ensembl_ids(self):
        return _decode(self.f["var"]["ensemblid"][:])

    def var_flags(self, name):
        """取 var 里的一个普通列（如 feature_is_filtered / highly_variable）。"""
        return self.f["var"][name][:]

    def umap(self):
        return self.f["obsm"]["X_umap"][:]

    # ---------- 表达矩阵 ----------
    def raw_block(self, start, stop):
        """
        取 raw.X 的第 [start, stop) 行，返回 csr_matrix（原始计数）。
        只读这一段对应的 data/indices，不碰其余部分。
        """
        r = self.f["raw"]["X"]
        ip = r["indptr"][start:stop + 1].astype(np.int64)
        a, b = int(ip[0]), int(ip[-1])
        data = r["data"][a:b]
        indices = r["indices"][a:b]
        return sparse.csr_matrix((data, indices, ip - a),
                                 shape=(stop - start, self.n_vars))

    def raw_rows(self, rows):
        """取任意若干行（不要求连续），按给定顺序返回。"""
        rows = np.asarray(rows, dtype=np.int64)
        order = np.argsort(rows)
        srows = rows[order]
        block = self.raw_block(int(srows[0]), int(srows[-1]) + 1)
        picked = block[srows - int(srows[0])]
        inv = np.empty(len(rows), dtype=np.int64)
        inv[order] = np.arange(len(rows))
        return picked[inv]

    def n_genes_per_cell(self, start, stop):
        """每个细胞检测到的基因数。"""
        ip = self.f["raw"]["X"]["indptr"][start:stop + 1]
        return np.diff(ip)

    def mito_mask(self):
        """线粒体基因（MT- 开头）。"""
        return np.array([str(s).upper().startswith("MT-")
                         for s in self.symbols()], dtype=bool)


class SliceSource:
    """
    远程切片得到的本地子集，接口与 TabulaSapiens 完全一致。

    背景
        全量 h5ad 10.02 GB，在实测网络（单连接 0.24 MB/s、8 并发 1.4 MB/s、
        且持续大流量后会被限流挂断）下不可行。92_远程切片.py 用 HTTP Range
        只取需要的细胞，存成 .npz + .obs.tsv + .genes.txt。

    三个产物文件
        <prefix>.npz       scipy 稀疏矩阵（原始 UMI 计数，csr）
        <prefix>.obs.tsv   每行一个细胞，含 cell_index（原始行号）+ 若干 obs 列
        <prefix>.genes.txt 基因符号，每行一个
    """

    def __init__(self, prefix):
        import os

        self.prefix = prefix
        self.X = sparse.load_npz(prefix + ".npz")
        self.n_obs, self.n_vars = self.X.shape

        obs = pd.read_csv(prefix + ".obs.tsv", sep="\t", dtype=str,
                          keep_default_na=False)
        self._obs = obs
        self._symbols = np.array(
            open(prefix + ".genes.txt", encoding="utf-8").read().split("\n"),
            dtype=object)
        if len(self._symbols) != self.n_vars:
            raise ValueError(f"基因名数量 {len(self._symbols)} "
                             f"与矩阵列数 {self.n_vars} 不符")

    # ---------- 生命周期 ----------
    def close(self):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    # ---------- 元数据 ----------
    def obs_columns(self):
        return [c for c in self._obs.columns if c != "cell_index"]

    def obs(self, col):
        if col not in self._obs.columns:
            raise KeyError(f"obs 里没有 '{col}'。现有列：{self.obs_columns()}")
        return self._obs[col].to_numpy()

    def obs_frame(self, cols):
        return self._obs[[c for c in cols if c != "cell_index"]]

    def symbols(self):
        return self._symbols

    def obs_index(self):
        return self._obs["cell_index"].to_numpy()

    # ---------- 表达矩阵 ----------
    def raw_block(self, start, stop):
        return self.X[start:stop]

    def raw_rows(self, rows):
        return self.X[np.asarray(rows, dtype=np.int64)]

    def n_genes_per_cell(self, start, stop):
        return np.diff(self.X.indptr[start:stop + 1])

    def mito_mask(self):
        return np.array([str(s).upper().startswith("MT-")
                         for s in self.symbols()], dtype=bool)
