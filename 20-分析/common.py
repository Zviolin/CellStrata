"""
common.py — 公共配置与工具函数
================================
所有脚本共用。改路径只改这里。

依赖：numpy / pandas / scipy。
**不需要 scanpy 和 anndata** —— 本项目直接读写 h5ad 的底层 HDF5 结构
（见 ts_io.py），这样既省内存，也不受 anndata 版本兼容问题影响。
"""

import os
import sys

import numpy as np
import pandas as pd

# ---------- 路径 ----------
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(ROOT, "10-数据")
ANALYSIS_DIR = os.path.join(ROOT, "20-分析")
FIG_DIR = os.path.join(ROOT, "30-可视化", "figures")

for d in (DATA_DIR, ANALYSIS_DIR, FIG_DIR):
    os.makedirs(d, exist_ok=True)

# 单细胞表达矩阵（Tabula Sapiens 1.0，483,152 细胞 × 58,604 基因）
H5AD_PATH = os.path.join(DATA_DIR, "tabula_sapiens.h5ad")

# 远程切片产物（92_远程切片.py 生成）。
# 完整 h5ad 10.02 GB 在实测网络下不可行，所以用 HTTP Range 只取需要的细胞，
# 存成 .npz + .obs.tsv + .genes.txt。接口与 TabulaSapiens 完全一致，
# 下游脚本（03/04/05/06）不用区分数据源。
SLICE_PREFIX = os.path.join(DATA_DIR, "ts_slice")

# 基因年龄表：主用 Trigos et al. 2017（PNAS）的 16 层划分
GENE_AGE_RAW = os.path.join(DATA_DIR, "gene_age_原始_TPPG17.csv")
# 交叉验证用：Liebeskind et al. 2016 的 8 个粗类别（13 个算法共识）
GENE_AGE_LMM16 = os.path.join(DATA_DIR, "gene_age_交叉验证_LMM16.csv")
# 置信度指标：Liebeskind 原始共识表，带每个基因的算法支持数与熵
GENE_AGE_META = os.path.join(DATA_DIR, "gene_age_Liebeskind2016_共识与误差指标.csv")

# 产物
GENE_AGE_CLEAN = os.path.join(DATA_DIR, "gene_age.tsv")
AGGR_PATH = os.path.join(ANALYSIS_DIR, "celltype_mean_expr.tsv")
PROFILE_PATH = os.path.join(ANALYSIS_DIR, "celltype_age_profile.tsv")
LINEAGE_MAP_PATH = os.path.join(ANALYSIS_DIR, "lineage_map.tsv")
DONOR_TABLE_PATH = os.path.join(ANALYSIS_DIR, "celltype_donor_table.tsv")

# ---------- 常量 ----------
RANDOM_STATE = 42

# 年龄层的名字。
# 来源：JMDplots 包 R/evdevH2O.R 里画图例的那段代码，逐条对应 Trigos et al.
# 2017 (PNAS) Dataset S1 的 1–16 号 phylostratum。不是网上抄的二手说法。
#
# 注意：Trigos 的方案只到 16 层，没有 Domazet-Lošo 19 层方案里的
# Deuterostomia / Vertebrata / Tetrapoda 这几级，别混用。
STRATUM_NAMES = {
    1: "Cellular organisms 细胞生物",
    2: "Eukaryota 真核生物",
    3: "Opisthokonta 后鞭毛类",
    4: "Metazoa 后生动物",
    5: "Eumetazoa 真后生动物",
    6: "Bilateria 两侧对称动物",
    7: "Chordata 脊索动物",
    8: "Euteleostomi 硬骨脊椎动物",
    9: "Amniota 羊膜动物",
    10: "Mammalia 哺乳纲",
    11: "Theria 兽亚纲",
    12: "Eutheria 真兽下纲",
    13: "Euarchontoglires 灵长总目",
    14: "Catarrhini 狭鼻小目",
    15: "Homininae 人亚科",
    16: "Homo sapiens 智人",
}
N_STRATA = len(STRATUM_NAMES)


# ---------- 工具 ----------
def die(msg):
    """带格式的致命错误，方便在日志里定位。"""
    print("\n" + "=" * 60)
    print("✗ 中止：" + msg)
    print("=" * 60 + "\n")
    sys.exit(1)


def info(msg):
    print(f"[信息] {msg}")


def ok(msg):
    print(f"[完成] {msg}")


def warn(msg):
    print(f"[注意] {msg}")


def require_file(path, hint=""):
    if not os.path.exists(path):
        die(f"找不到文件：{path}\n  {hint}")
    return path


# ---------- 数据源 ----------
def open_source(prefer=None):
    """
    打开单细胞数据源，返回接口统一的只读视图。

    prefer: None=自动（有切片用切片） / True=只用切片 / False=只用完整 h5ad

    两种源的接口完全一致（obs / symbols / raw_block / raw_rows / mito_mask），
    下游脚本不需要知道底下是哪种。
    """
    from ts_io import TabulaSapiens, SliceSource

    slice_ok = os.path.exists(SLICE_PREFIX + ".npz")
    h5ad_ok = os.path.exists(H5AD_PATH)

    if prefer is True and not slice_ok:
        die(f"没有切片文件：{SLICE_PREFIX}.npz\n"
            f"  先跑：python 92_远程切片.py --per-type 200")
    if prefer is False and not h5ad_ok:
        die(f"没有完整 h5ad：{H5AD_PATH}")

    if slice_ok:
        info(f"数据源：远程切片 {os.path.basename(SLICE_PREFIX)}.*")
        return SliceSource(SLICE_PREFIX)
    if h5ad_ok:
        info(f"数据源：完整 h5ad {os.path.basename(H5AD_PATH)}")
        return TabulaSapiens(H5AD_PATH)
    die("找不到单细胞数据。\n"
        f"  要么跑 92_远程切片.py 生成切片，要么把完整 h5ad 放到：\n"
        f"  {H5AD_PATH}")


# ---------- 基因年龄 ----------
def load_gene_age():
    """读取清洗后的基因年龄表，返回 DataFrame[symbol, stratum]。"""
    require_file(GENE_AGE_CLEAN, "请先运行 02_构建年龄表.py 生成")
    df = pd.read_csv(GENE_AGE_CLEAN, sep="\t")
    need = {"symbol", "stratum"}
    if not need.issubset(df.columns):
        die(f"基因年龄表缺列。现有列：{df.columns.tolist()}，需要 {need}")
    return df


def map_genes(var_names, age_df):
    """
    把表达矩阵的基因名映射到年龄表。
    返回 (matched_mask, stratum_array)。

    matched_mask: 长度 = len(var_names) 的布尔数组
    stratum_array: 对齐后的年龄层，未匹配处为 NaN

    注意：Tabula Sapiens 的基因名有少量重复（同名不同 Ensembl ID），
    这里不做去重 —— 两个同名基因会拿到同一个 age，表达量在后面按层
    求和时会自然合并，这正是我们要的。

    两个方向的比例含义不同，都要报：
      - 「表达矩阵命中率」：表达矩阵里多少基因有年龄标注。
        这个数天然不高（58,604 个基因里只有约 1.9 万是蛋白编码，
        而 Trigos 表只覆盖蛋白编码基因），所以它**不是**验收指标。
      - 「年龄表覆盖率」：年龄表里的基因有多少能在表达矩阵中找到。
        这个才是验收指标（目标 ≥ 80%）。
    """
    age_map = dict(zip(age_df["symbol"].astype(str).str.upper(),
                       age_df["stratum"]))
    names_upper = pd.Index(var_names).astype(str).str.upper()

    stratum = np.array([age_map.get(g, np.nan) for g in names_upper], dtype=float)
    matched = ~np.isnan(stratum)

    rate = matched.mean() if len(var_names) else 0.0
    info(f"表达矩阵基因命中：{matched.sum():,} / {len(var_names):,} = {rate:.1%}")

    # 反向覆盖：年龄表里的基因有多少能在表达矩阵中找到
    hit = sum(1 for g in names_upper if g in age_map)
    cov = hit / len(age_df) if len(age_df) else 0.0
    info(f"年龄表基因覆盖率：{hit:,} / {len(age_df):,} = {cov:.1%}   ← 验收指标")

    # 按层看命中率：某层若大量丢失，TAI 会被系统性拉偏，必须查
    name_set = set(names_upper)
    per_layer = []
    for s in sorted(set(int(x) for x in age_df["stratum"])):
        sub = age_df[age_df["stratum"] == s]["symbol"].astype(str).str.upper()
        n_hit = sum(1 for g in sub if g in name_set)
        per_layer.append((s, len(sub), n_hit, n_hit / len(sub) if len(sub) else 0))
    worst = min(per_layer, key=lambda x: x[3])
    info(f"按层命中率：最低 S{worst[0]}（{worst[1]} 个基因，命中 {worst[2]}，"
         f"{worst[3]:.1%}）")
    info("  " + "  ".join(f"S{s}:{h/n:.0%}" for s, n, h, r in per_layer
                          if n))

    if cov < 0.8:
        print(f"\n  ⚠️  年龄表覆盖率 {cov:.1%} 低于验收标准 80%。")
        print("     常见原因：")
        print("     - 表达矩阵用 Ensembl ID（ENSG...），年龄表用 symbol")
        print("     - 年龄表 symbol 大小写不统一")
        print("     - 年龄表是旧版基因符号")
        print("  处理：先跑 01_检查数据.py 看基因名格式。\n")
    elif cov < 0.9:
        print(f"  ⚠️  年龄表覆盖率 {cov:.1%}，未达 90% 但已过 80% 验收线。\n")

    return matched, stratum


def compute_age_profile(mean_expr, stratum):
    """
    由细胞类型平均表达 + 基因年龄层，算年龄构成与 TAI。

    mean_expr: DataFrame (cell_type x gene)
    stratum:   与 mean_expr.columns 对齐的年龄层数组（1=最老）

    返回 DataFrame: index=cell_type, columns=[S1..S16, TAI, total_expr]
    按 TAI 升序排列（越小越古老）。
    """
    valid = ~np.isnan(stratum)
    E = mean_expr.values[:, valid]            # celltypes x valid_genes
    S = stratum[valid]

    strata = np.sort(np.unique(S))
    profile = np.zeros((E.shape[0], len(strata)), dtype=np.float64)

    for i, s in enumerate(strata):
        cols = (S == s)
        profile[:, i] = E[:, cols].sum(axis=1)

    total_raw = profile.sum(axis=1, keepdims=True)
    total_safe = np.where(total_raw == 0, 1.0, total_raw)
    prop = profile / total_safe               # 占比

    tai = (prop * strata[None, :]).sum(axis=1)

    col_names = [f"S{int(s)}" for s in strata]
    out = pd.DataFrame(prop, index=mean_expr.index, columns=col_names)
    out["TAI"] = tai
    out["total_expr"] = total_raw.ravel()
    return out.sort_values("TAI")
