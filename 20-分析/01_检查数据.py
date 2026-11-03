"""
01_检查数据.py — 开工前的体检
================================
不加载矩阵，只读元数据，几秒钟出结果。任何一项不对，后面的脚本都会算错，
所以**新拿到数据先跑这个**。

检查内容：
    1. 文件身份（标题、细胞数、基因数）
    2. 表达矩阵是原始计数还是已归一化（决定 03 要不要自己归一化）
    3. obs 里有哪些列，关键列（细胞类型 / 供者 / 组织）齐不齐
    4. 基因名是什么格式，能不能对上年龄表，匹配率多少
    5. 数据规模与分块策略的内存估算

用法：
    python 01_检查数据.py
"""

import numpy as np
import pandas as pd

from common import (H5AD_PATH, GENE_AGE_CLEAN,
                    die, info, ok, require_file, map_genes, STRATUM_NAMES,
                    open_source)
from ts_io import TabulaSapiens


def check_identity(ts):
    print("─" * 60)
    print("1. 文件身份")
    print("─" * 60)
    if hasattr(ts, "path"):
        print(f"  标题：{ts.title}")
        print(f"  规模：{ts.n_obs:,} 细胞 × {ts.n_vars:,} 基因")
        print(f"  大小：{__import__('os').path.getsize(ts.path)/1e9:.2f} GB")
    else:
        print(f"  数据源：远程切片（{ts.n_obs:,} 细胞 × {ts.n_vars:,} 基因）")
        import os
        p = ts.prefix + ".npz"
        print(f"  大小：{os.path.getsize(p)/1e6:.1f} MB")


def check_matrix(ts):
    """判断 X 是原始计数还是已归一化。"""
    print("\n" + "─" * 60)
    print("2. 表达矩阵性质")
    print("─" * 60)

    if not hasattr(ts, "f"):
        d = np.asarray(ts.X.data[:5000], dtype=np.float64)
        is_int = bool(np.all(d == np.round(d)))
        print(f"  切片矩阵 dtype={ts.X.data.dtype}  最大={d.max():.2f}  整数={is_int}")
        print("  → 切片只存 raw.X（原始 UMI 计数），03 脚本自己做 QC + 归一化。")
        return

    def describe(node_name, getter):
        d = getter()[:5000]
        d = np.asarray(d, dtype=np.float64)
        is_int = bool(np.all(d == np.round(d)))
        print(f"  {node_name:10s} dtype={getter().dtype}  "
              f"最大={d.max():.2f}  整数={is_int}")
        return is_int

    try:
        x_int = describe("X", lambda: ts.f["X"]["data"])
    except Exception as e:
        print(f"  X 读取失败：{e}")
        x_int = None

    if "raw" in ts.f:
        try:
            r_int = describe("raw.X", lambda: ts.f["raw"]["X"]["data"])
        except Exception as e:
            print(f"  raw.X 读取失败：{e}")
            r_int = None
    else:
        r_int = None
        print("  raw.X    不存在")

    print()
    if r_int:
        print("  → raw.X 是原始 UMI 计数。03 脚本用 raw.X 自己做 QC + 归一化。")
    elif x_int:
        print("  → X 本身就是原始计数。03 脚本直接用 X。")
    else:
        print("  → X 已经是归一化/对数化的值，03 脚本不能再归一化一遍。")
        print("    ⚠ 如果 raw 层不存在，QC（线粒体比例、基因数）会失真，要在报告里说明。")


def check_obs(ts):
    print("\n" + "─" * 60)
    print("3. obs 元数据")
    print("─" * 60)
    cols = ts.obs_columns()
    print(f"  共 {len(cols)} 列：{cols}")

    key_specs = [
        ("cell_type", "细胞类型（分组变量）"),
        ("cell_ontology_class", "细胞类型（CL 本体注释）"),
        ("donor_id", "供者（划分训练/测试用）"),
        ("tissue", "组织"),
        ("compartment", "区室（谱系代理）"),
        ("assay", "测序平台"),
        ("disease", "疾病状态"),
    ]
    print("\n  关键列：")
    missing = []
    for c, why in key_specs:
        if c in cols:
            vals = pd.Series(ts.obs(c))
            n = vals.nunique()
            sample = ", ".join(map(str, vals.dropna().unique()[:4]))
            print(f"    ✓ {c:22s} {n:>4} 类   {why}")
            print(f"      {sample}{' …' if n > 4 else ''}")
        else:
            missing.append(c)
            print(f"    ✗ {c:22s} 缺失   {why}")

    if "cell_type" in missing and "cell_ontology_class" not in missing:
        print("\n  → 用 cell_ontology_class 当分组变量（--group-key cell_ontology_class）")
    return cols, missing


def check_genes(ts, age_df):
    print("\n" + "─" * 60)
    print("4. 基因名与年龄表匹配")
    print("─" * 60)
    syms = ts.symbols()

    if hasattr(ts, "ensembl_ids"):
        ens = ts.ensembl_ids()
        print(f"  基因标识样例：{list(ens[:3])}  （var index 是 Ensembl ID）")
    else:
        print("  （切片产物只保留基因符号，未保留 Ensembl ID）")
    print(f"  基因符号样例：{list(syms[:8])}")
    n_dup = len(syms) - int(pd.Index(syms).nunique())
    print(f"  去后缀后重名：{n_dup} 个（同名不同 Ensembl ID，不影响按层求和）")

    matched, stratum = map_genes(syms, age_df)
    print(f"  年龄表覆盖的层：{sorted(set(int(s) for s in np.unique(stratum[matched])))}")

    # 验收标准是「年龄表基因的覆盖率」，不是「表达矩阵基因的命中率」
    age_up = set(age_df["symbol"].astype(str).str.upper())
    present = sum(1 for g in pd.Index(syms).astype(str).str.upper() if g in age_up)
    cov = present / len(age_df)
    print(f"  年龄表覆盖率：{present:,} / {len(age_df):,} = {cov:.1%}")
    if cov < 0.8:
        die(f"年龄表覆盖率 {cov:.1%} 低于 80% 的验收标准。检查基因名格式是否一致。")
    return matched, stratum


def check_memory(ts):
    print("\n" + "─" * 60)
    print("5. 规模与内存")
    print("─" * 60)
    if hasattr(ts, "f"):
        nnz = (ts.f["raw"]["X"]["data"].shape[0] if "raw" in ts.f
               else ts.f["X"]["data"].shape[0])
    else:
        nnz = ts.X.nnz
    density = nnz / (ts.n_obs * ts.n_vars)
    per_cell = nnz / ts.n_obs
    print(f"  非零元：{nnz:,}（密度 {density:.2%}，平均每细胞 {per_cell:,.0f} 个基因）")
    print(f"  稀疏存储约：{nnz*8/1e9:.1f} GB（float32 + int32）")
    print(f"  稠密存储约：{ts.n_obs*ts.n_vars*4/1e9:.0f} GB  ← 绝对不能用 dense")

    for chunk in (10000, 20000, 30000):
        peak = chunk * per_cell * 8 / 1e9
        print(f"  分块 {chunk:>6,} 细胞/批：峰值约 {peak:.2f} GB，"
              f"共 {int(np.ceil(ts.n_obs/chunk))} 批")


def main():
    require_file(GENE_AGE_CLEAN, "请先运行 02_构建年龄表.py")

    age_df = pd.read_csv(GENE_AGE_CLEAN, sep="\t")
    info(f"年龄表：{len(age_df):,} 个基因，{age_df['stratum'].nunique()} 层")

    with open_source() as ts:
        check_identity(ts)
        check_matrix(ts)
        cols, missing = check_obs(ts)
        matched, stratum = check_genes(ts, age_df)
        check_memory(ts)

    print("\n" + "=" * 60)
    print("结论")
    print("=" * 60)
    group_key = "cell_type" if "cell_type" in cols else "cell_ontology_class"
    print(f"  ✓ 可以开始分析了")
    print(f"    分组变量：{group_key}")
    print(f"    供者列：{'donor_id' if 'donor_id' in cols else '缺失（05 无法按供者划分）'}")
    print(f"    可用年龄层：{len(STRATUM_NAMES)} 层（Trigos 2017）")
    print(f"\n  下一步：python 03_细胞类型年龄谱.py --group-key {group_key}")
    print("=" * 60)


if __name__ == "__main__":
    main()
