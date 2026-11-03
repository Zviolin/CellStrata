"""
02_构建年龄表.py — 把下载来的原始年龄表清洗成统一格式
========================================================
输入：
    10-数据/gene_age_原始_TPPG17.csv
        Trigos et al. 2017 (PNAS) Dataset S1，列：GeneID, Entrez, Entry, Phylostrata
        17,318 个蛋白编码基因，Phylostrata 取 1–16（1=最古老）
    10-数据/gene_age_交叉验证_LMM16.csv
        Liebeskind et al. 2016 (GBE)，列：UniProt, modeAge, PS
        作为**第二个独立来源**，用来交叉验证年龄调用是否一致
    10-数据/gene_age_Liebeskind2016_共识与误差指标.csv
        Liebeskind 的原始共识表，带每个基因的误差指标
        （13 个算法里有多少个支持、Shannon 熵、节点误差……）

输出：
    10-数据/gene_age.tsv
        列：symbol, stratum, stratum_name, n_dbs, entropy
        n_dbs / entropy 是"这个基因的年龄判得可不可靠"的代理指标，
        缺失表示该基因不在 Liebeskind 的表里。

用法：
    python 02_构建年龄表.py

为什么要做交叉验证
    phylostratigraphy 有一个被写进文献的方法学争议（Moyers & Zhang 2015/2017/2018）：
    它会系统性低估基因年龄，而且误差非随机。光声明"我有局限性"没用，
    必须给出两个独立来源一致性的定量证据。这个脚本算的就是这个数。
"""

import numpy as np
import pandas as pd

from common import (GENE_AGE_RAW, GENE_AGE_LMM16, GENE_AGE_META, GENE_AGE_CLEAN,
                    STRATUM_NAMES, N_STRATA, die, info, ok, warn, require_file)


def load_main():
    """读主年龄表（Trigos 2017）。"""
    require_file(GENE_AGE_RAW, "请先按 10-数据/00-数据源清单.md 下载")
    info(f"读取主年龄表：{GENE_AGE_RAW}")
    df = pd.read_csv(GENE_AGE_RAW, encoding="utf-8")
    info(f"原始：{df.shape[0]:,} 行 × {df.shape[1]} 列，列名 {df.columns.tolist()}")

    need = {"GeneID", "Phylostrata"}
    if not need.issubset(df.columns):
        die(f"主年龄表缺列，需要 {need}，实际 {df.columns.tolist()}")

    out = pd.DataFrame({
        "symbol": df["GeneID"].astype(str).str.strip().str.upper(),
        "stratum": pd.to_numeric(df["Phylostrata"], errors="coerce"),
        "uniprot": df["Entry"].astype(str).str.strip() if "Entry" in df else "",
    })
    out = out[out["symbol"].str.len() > 0].dropna(subset=["stratum"])
    out["stratum"] = out["stratum"].astype(int)
    return out


def attach_confidence(df):
    """
    把 Liebeskind 的误差指标并进来。

    两个文件都按 UniProt 索引，可以 join。

    ⚠ **重要口径修正（2026-11-03 审核发现）**：
    Liebeskind 的 `NumDBsContributing`（本项目重命名为 `n_dbs`）**不是**
    「有几个算法支持这个年龄调用」，而是**参与投票的算法数（分母）**。
    原始表里 8 个粗类别的 share 之和恒为 1，且 `share_c × n_dbs` 恰好是整数
    （实测最大残差 6e-12），说明 share 就是「各算法投给各粗类别的比例」，
    n_dbs 是这些票的总数。

    因此：
      真正「支持 modeAge 那一类」的算法数 = round(max_share × n_dbs)
      共识强度 = max_share（∈[0,1]），或等价地用 entropy（越低越一致）

    更反直觉的是：n_dbs 与共识度**负相关**（实测 corr = −0.128）——
    n_dbs ≤ 4 的基因 mode_share 中位 1.000（全体一致），
    n_dbs ≥ 10 的基因 mode_share 中位仅 0.692。
    所以「按 n_dbs 过滤」筛的是「投票算法多的基因」，**方向与「高置信」相反**。
    下游的置信度过滤必须用 mode_share / mode_support / entropy。
    """
    import os
    if not os.path.exists(GENE_AGE_META):
        info(f"没有 {os.path.basename(GENE_AGE_META)}，跳过置信度列")
        for c in ("n_dbs", "mode_support", "mode_share", "entropy"):
            df[c] = np.nan
        return df

    info(f"并入置信度指标：{os.path.basename(GENE_AGE_META)}")
    meta = pd.read_csv(GENE_AGE_META, index_col=0, encoding="utf-8")
    meta.index = meta.index.astype(str).str.strip()

    # 8 个粗类别（与 cross_validate 里的 order 一致）
    cat_cols = [c for c in ("Cellular_organisms", "Euk_Archaea", "Euk+Bac",
                             "Eukaryota", "Opisthokonta", "Eumetazoa",
                             "Vertebrata", "Mammalia") if c in meta.columns]
    share = meta[cat_cols].astype(float)
    n_dbs = (meta["NumDBsContributing"].astype(float)
             if "NumDBsContributing" in meta.columns
             else pd.Series(np.nan, index=meta.index))
    ent = (meta["entropy"].astype(float)
           if "entropy" in meta.columns
           else pd.Series(np.nan, index=meta.index))

    # 真正的支持数与共识份额
    S = share.to_numpy()
    mode_share = np.nanmax(S, axis=1)
    with np.errstate(invalid="ignore"):
        mode_support = np.round(mode_share * n_dbs.to_numpy())
    # share×n_dbs 必须是整数票数——若不是，说明文件结构变了，要报警
    votes = S * n_dbs.to_numpy()[:, None]
    resid = np.nanmax(np.abs(votes - np.round(votes)))
    if np.isfinite(resid) and resid > 1e-6:
        warn(f"share×n_dbs 不是整数（最大残差 {resid:.2e}），"
             f"Liebeskind 表结构可能已变，请复核共识度定义")

    meta = meta.assign(
        n_dbs=n_dbs.to_numpy(),
        mode_support=mode_support,
        mode_share=mode_share,
        entropy=ent.to_numpy(),
    )
    meta = meta[["n_dbs", "mode_support", "mode_share", "entropy"]]

    joined = df.merge(meta, left_on="uniprot", right_index=True, how="left")
    hit = joined["n_dbs"].notna().mean()
    info(f"  命中 {joined['n_dbs'].notna().sum():,} / {len(joined):,} = {hit:.1%}")
    ok_g = joined["mode_share"].notna()
    if ok_g.sum():
        info(f"  共识度：mode_share 中位 {joined.loc[ok_g, 'mode_share'].median():.3f}"
             f"，全体一致（mode_share=1）的基因 "
             f"{int((joined.loc[ok_g, 'mode_share'] >= 0.999).sum()):,}")
        info(f"  mode_support 中位 {joined.loc[ok_g, 'mode_support'].median():.0f}"
             f"（= 真正支持该层号的算法数）")
    return joined


def cross_validate(df):
    """
    与第二个独立来源（Liebeskind 2016，8 个粗类别）比对。
    两个来源方法完全不同（13 算法共识 vs 单算法系谱分配），
    如果它们对同一批基因给出的年龄**排序**一致，说明年龄调用不是噪声。
    """
    import os
    if not os.path.exists(GENE_AGE_LMM16):
        info("没有交叉验证表，跳过这一步")
        return

    info(f"交叉验证：{os.path.basename(GENE_AGE_LMM16)}")
    lmm = pd.read_csv(GENE_AGE_LMM16, encoding="utf-8")
    lmm["UniProt"] = lmm["UniProt"].astype(str).str.strip()
    lmm = lmm.drop_duplicates("UniProt")

    m = df.merge(lmm, left_on="uniprot", right_on="UniProt", how="inner")
    n = len(m)
    if n < 1000:
        info(f"  可比对的基因只有 {n} 个，太少，跳过")
        return

    rho = m["stratum"].corr(m["PS"], method="spearman")
    tau = m["stratum"].corr(m["PS"], method="kendall")

    print(f"  可比对基因：{n:,}")
    print(f"  Spearman ρ = {rho:.3f}   Kendall τ = {tau:.3f}   "
          f"（两个来源的年龄排序一致性）")

    print("\n  Liebeskind 的 8 个粗类别 → 落在 TPPG 各层的分布（中位数 / 四分位）")
    print(f"  {'粗类别':<20}{'基因数':>7}{'中位层':>8}{'Q1':>6}{'Q3':>6}")
    print("  " + "-" * 50)
    order = ["Cellular_organisms", "Euk_Archaea", "Euk+Bac", "Eukaryota",
             "Opisthokonta", "Eumetazoa", "Vertebrata", "Mammalia"]
    meds = []
    for cat in order:
        sub = m.loc[m["modeAge"] == cat, "stratum"]
        if len(sub) == 0:
            continue
        q1, med, q3 = np.percentile(sub, [25, 50, 75])
        meds.append(med)
        print(f"  {cat:<20}{len(sub):>7}{med:>8.0f}{q1:>6.0f}{q3:>6.0f}")

    # 单调性：粗类别从老到新，中位层号也应该单调不减
    if len(meds) >= 3:
        violations = sum(1 for a, b in zip(meds, meds[1:]) if b < a)
        print(f"\n  单调性检查：{len(meds)-1} 个相邻对中，违反单调的有 {violations} 个")
        if violations == 0:
            print("  ✓ 两个来源的年龄顺序完全一致 —— 年龄调用不是随机噪声")
        else:
            print("  ⚠ 存在逆序，说明部分基因在两个来源间分歧较大（写进局限性）")

    if rho >= 0.7:
        print(f"  ✓ ρ={rho:.2f} ≥ 0.7，一致性良好")
    else:
        print(f"  ⚠ ρ={rho:.2f} < 0.7。注意：LMM16 只有 8 个取值，"
              f"大量基因并列在同一类别，\n"
              f"    Spearman 在这种强并列情况下会被系统性压低，"
              f"所以以「中位层号单调不逆序」为主证据。")


def main():
    df = load_main()
    n0 = len(df)

    # 同一符号多条：保留最古老（层号最小）的那条，最保守
    dup = df["symbol"].duplicated().sum()
    if dup:
        info(f"发现 {dup} 个重复符号，保留最古老的那条")
        df = df.sort_values("stratum").drop_duplicates("symbol", keep="first")

    out = attach_confidence(df)
    out["stratum_name"] = out["stratum"].map(STRATUM_NAMES).fillna("Unknown")

    # ---------- 自检 ----------
    print("\n" + "─" * 56)
    print("自检")
    print("─" * 56)

    print(f"  行数：{n0:,} → {len(out):,}")
    print(f"  层数：{out['stratum'].nunique()}")
    print(f"  层范围：{out['stratum'].min()} – {out['stratum'].max()}")

    counts = out["stratum"].value_counts().sort_index()
    print("\n  各层基因数：")
    for s, c in counts.items():
        bar = "█" * max(1, int(c / counts.max() * 34))
        print(f"    S{s:>2}  {c:>5}  {bar}")

    problems = []
    if not (15000 <= len(out) <= 23000):
        problems.append(f"基因数 {len(out)} 不在 15000–23000 区间")
    if not (12 <= out["stratum"].nunique() <= 20):
        problems.append(f"层数 {out['stratum'].nunique()} 不在 12–20 区间")
    if out["stratum"].min() != 1:
        problems.append("最老层不是 1，层号可能没有从 1 开始")
    # 最老层应该明显比最新层大（基因越老共享的物种越多）
    oldest, youngest = counts.iloc[0], counts.iloc[-1]
    if oldest < youngest * 3:
        problems.append(f"最老层 {oldest} 个、最新层 {youngest} 个，不符合递减规律")

    print("\n" + "─" * 56)
    if problems:
        print("✗ 自检未通过：")
        for p in problems:
            print("   - " + p)
        print("\n  可能拿错了表。请回到 00-数据源清单.md 换备选来源。")
        print("─" * 56)
        return

    print("✓ 形态自检通过")
    print("─" * 56)

    # ---------- 交叉验证 ----------
    print("\n" + "─" * 56)
    print("交叉验证（两个独立来源）")
    print("─" * 56)
    cross_validate(out)

    # ---------- 落盘 ----------
    keep = ["symbol", "stratum", "stratum_name", "n_dbs", "mode_support",
            "mode_share", "entropy"]
    missing = [c for c in keep if c not in out.columns]
    if missing:
        warn(f"置信度列缺失（{missing}），将不写入")
        keep = [c for c in keep if c in out.columns]
    out[keep].to_csv(GENE_AGE_CLEAN, sep="\t", index=False)
    print()
    ok(f"已写出：{GENE_AGE_CLEAN}  ({len(out):,} 个基因 × {out['stratum'].nunique()} 层)")
    print(f"\n  下一步：python 01_检查数据.py  （确认表达矩阵能对上这张表）")


if __name__ == "__main__":
    main()
