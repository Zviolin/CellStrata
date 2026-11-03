"""
06_画图.py — 生成报告的九张图
==============================
用法：
    python 06_画图.py            # 画 1–5、7–9
    python 06_画图.py --all      # 含 UMAP（需要重新读 h5ad，较慢）

输出：30-可视化/figures/fig1_strata.png ... fig9_shortprot_premise.png
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm

# 这个脚本在 30-可视化/ 下，而 common.py / ts_io.py 在 20-分析/。
# 不把 20-分析 加进 sys.path，从本目录直接跑会 ModuleNotFoundError。
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "20-分析"))

from common import (PROFILE_PATH, ANALYSIS_DIR, FIG_DIR,  # noqa: E402
                    die, info, ok, require_file)

# 中文字体（Windows 常见）
for f in ["Microsoft YaHei", "SimHei", "PingFang SC", "DejaVu Sans"]:
    if f in [x.name for x in matplotlib.font_manager.fontManager.ttflist]:
        plt.rcParams["font.sans-serif"] = [f]
        break
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["figure.dpi"] = 150
plt.rcParams["savefig.bbox"] = "tight"


def color_for_stratum(i, n):
    """从深红（最老）到亮蓝（最新）。"""
    return cm.coolwarm_r(i / max(1, n - 1))


# ---------- 图 1：细胞地层柱 ----------
def fig1(prof):
    s_cols = [c for c in prof.columns if c.startswith("S") and c[1:].isdigit()]
    s_cols.sort(key=lambda x: int(x[1:]))
    n = len(s_cols)

    fig, ax = plt.subplots(figsize=(15, 7))
    x = np.arange(len(prof))
    bottom = np.zeros(len(prof))

    for i, c in enumerate(s_cols):
        vals = prof[c].values
        ax.bar(x, vals, bottom=bottom, width=0.85,
               color=color_for_stratum(i, n), edgecolor="none")
        bottom += vals

    # 叠加 TAI 曲线（右轴）：一眼看出年龄层构成与综合年龄指数是否同向
    ax2 = ax.twinx()
    ax2.plot(x, prof["TAI"].values, color="#111111", lw=1.1, alpha=0.75,
             zorder=5)
    ax2.set_ylabel("TAI（越小越古老）")
    ax2.set_ylim(prof["TAI"].min() - 0.1, prof["TAI"].max() + 0.1)

    ax.set_xlim(-1, len(prof))
    ax.set_ylim(0, 1)
    ax.set_ylabel("表达占比")
    ax.set_xlabel("细胞类型（按 TAI 排序：左=古老，右=年轻）")
    # ⚠ 图内不再放标题：期刊体例下图题由 LaTeX 图注承担，
    #    图内再写一遍会与图注重复，且长标题会把画布撑宽、正文里被压成微缩字。
    ax.set_title("")

    # 每隔若干个标一个类型名，否则 177 个标签完全重叠
    step = max(1, len(prof) // 22)
    ticks = list(range(0, len(prof), step))
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(prof.index[i])[:26] for i in ticks],
                       rotation=90, fontsize=5.5, ha="center")
    ax.tick_params(axis="x", length=0, pad=1)

    sm = cm.ScalarMappable(cmap="coolwarm_r",
                           norm=plt.Normalize(1, n))
    cbar = fig.colorbar(sm, ax=ax, pad=0.01, aspect=30)
    cbar.set_label("基因年龄层（1=最古老 → %d=最新）" % n)

    out = os.path.join(FIG_DIR, "fig1_strata.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    ok(out)


# ---------- 图 2：PCA 散点 ----------
def fig2():
    path = os.path.join(ANALYSIS_DIR, "unsupervised_result.tsv")
    if not os.path.exists(path):
        info("跳过图2：先跑 04_无监督分析.py")
        return
    df = pd.read_csv(path, sep="\t", index_col=0)

    fig, ax = plt.subplots(figsize=(8, 7))
    color_key = "lineage" if "lineage" in df.columns and df["lineage"].notna().any() else "cluster"
    if color_key == "lineage":
        cats = df["lineage"].astype("category")
        for i, c in enumerate(cats.cat.categories):
            m = cats == c
            ax.scatter(df.loc[m, "PC1"], df.loc[m, "PC2"],
                       s=45, alpha=0.8, label=str(c), edgecolor="w", linewidth=0.4)
        ax.legend(fontsize=8, framealpha=0.9, loc="best")
    else:
        sc = ax.scatter(df["PC1"], df["PC2"], c=df["cluster"],
                        cmap="tab10", s=45, alpha=0.85,
                        edgecolor="w", linewidth=0.4)
        ax.legend(*sc.legend_elements(), title="簇", fontsize=8)

    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("")
    ax.grid(alpha=0.2)
    out = os.path.join(FIG_DIR, "fig2_pca.png")
    fig.savefig(out)
    plt.close(fig)
    ok(out)


# ---------- 图 3：热图 ----------
def fig3(prof):
    s_cols = [c for c in prof.columns if c.startswith("S") and c[1:].isdigit()]
    s_cols.sort(key=lambda x: int(x[1:]))

    fig, ax = plt.subplots(figsize=(7, 12))
    im = ax.imshow(prof[s_cols].values, aspect="auto",
                   cmap="RdYlBu_r", interpolation="nearest")
    ax.set_xticks(range(len(s_cols)))
    ax.set_xticklabels([f"S{s[1:]}" for s in s_cols], fontsize=7, rotation=90)
    ax.set_yticks(range(len(prof)))
    ax.set_yticklabels(prof.index, fontsize=3.5)
    # 配色语义必须写清楚，否则读者会误读：
    # RdYlBu_r 是「红→黄→蓝」，红端=数值大=占比高。
    # 由于 S1 最古老、S16 最新，**红色端对应古老层**（不是"年轻"）。
    ax.set_title("")
    fig.colorbar(im, ax=ax, label="表达占比", pad=0.01, aspect=30)
    out = os.path.join(FIG_DIR, "fig3_heatmap.png")
    fig.savefig(out)
    plt.close(fig)
    ok(out)


# ---------- 图 4：重要性 × 年龄 ----------
def fig4():
    imp_path = os.path.join(ANALYSIS_DIR, "feature_importance.tsv")
    enr_path = os.path.join(ANALYSIS_DIR, "age_enrichment.tsv")
    if not os.path.exists(enr_path):
        info("跳过图4：先跑 05_监督分析.py")
        return

    enr = pd.read_csv(enr_path, sep="\t")
    # 背景基因数太少的层（如 S16 只有 1 个基因）统计上不可靠，
    # 且其倍数会拉高 y 轴掩盖主体层次 → 下panel 排除，但上panel 保留以示真实背景。
    MIN_BG = 5
    ok_mask = enr["n_background"] >= MIN_BG

    fig, (a1, a2) = plt.subplots(2, 1, figsize=(10, 7.5), sharex=True,
                                 gridspec_kw={"height_ratios": [1, 1.5]})

    a1.bar(enr["stratum"], enr["n_background"], color="#b0b0b0",
           label="背景（全部 HVG）", width=0.8)
    a1.set_ylabel("基因数")
    a1.set_yscale("log")
    a1.set_title("")
    a1.legend(fontsize=8)

    sub = enr[ok_mask].copy()
    # 显著性来自双单尾：富集看右尾 q，耗竭看左尾 q，各自配对。
    # ⚠ 判定必须用 **BH 校正后的 q**（正式口径），不是原始 p。
    #   本项目每次运行有 16 层 × 2 方向 = 32 次检验，跨三次运行共 96 次；
    #   用未校正 p<0.05 会把 S7/S12/S13/S15 也标成显著，与正文表格矛盾。
    has_dep = "p_depletion" in enr.columns
    p_up = sub["p_value"].values
    p_dn = sub["p_depletion"].values if has_dep else np.full(len(sub), 1.0)
    q_up = (sub["q_value"].values if "q_value" in sub.columns
            else p_up)
    q_dn = (sub["q_depletion"].values
            if "q_depletion" in sub.columns else p_dn)
    fold = sub["fold"].fillna(1).values

    colors = []
    for f, qu, qd in zip(fold, q_up, q_dn):
        if qu < 0.05 and f > 1:
            colors.append("#d62728")
        elif qd < 0.05 and f < 1:
            colors.append("#1f77b4")
        else:
            colors.append("#999999")
    # 富集层标红并加★；耗竭层标蓝并加▼
    marks = []
    for f, qu, qd in zip(fold, q_up, q_dn):
        marks.append("★" if (qu < 0.05 and f > 1)
                     else ("▼" if (qd < 0.05 and f < 1) else ""))

    a2.bar(sub["stratum"], fold, color=colors, width=0.8)
    a2.axhline(1, color="k", lw=0.8, ls="--")
    a2.set_ylabel("富集倍数")
    a2.set_xlabel("基因年龄层（1 = 最古老）")
    a2.set_ylim(0, max(1.15, float(np.nanmax(fold)) * 1.25))

    # 倍数 + 显著性标记（只标显著的，避免文字堆叠）
    for x, f, mk in zip(sub["stratum"], fold, marks):
        if mk:
            a2.text(x, f + 0.10, f"{f:.2f}{mk}", ha="center", fontsize=8,
                    color="#d62728" if mk == "★" else "#1f77b4")

    # 被排除的层（背景基因过少）在 x 轴上用灰底注明
    for _, r in enr[~ok_mask].iterrows():
        a2.axvspan(r["stratum"] - 0.45, r["stratum"] + 0.45,
                   color="#eeeeee", zorder=0)
        a2.text(r["stratum"], 0.08, f"背景仅\n{r['n_background']}基因\n不参与检验",
                ha="center", va="bottom", fontsize=6.5, color="#666666")

    fig.tight_layout()
    out = os.path.join(FIG_DIR, "fig4_importance_age.png")
    fig.savefig(out)
    plt.close(fig)
    ok(out)


# ---------- 图 5：各细胞类型的 TAI 排序 ----------
def fig5(prof):
    fig, ax = plt.subplots(figsize=(9, max(6, len(prof) * 0.09)))
    y = np.arange(len(prof))
    norm = plt.Normalize(prof["TAI"].min(), prof["TAI"].max())
    colors = cm.coolwarm_r(norm(prof["TAI"].values))

    # 细胞数太少的类型（低置信）用浅色标出，避免读者误当结论
    if "low_confidence" in prof.columns:
        low = prof["low_confidence"].to_numpy(dtype=bool)
        colors = np.array(colors, dtype=float).copy()
        colors[low] = np.array([0.85, 0.85, 0.85, 1.0])

    # 误差棒：**供者级为主口径**（正文 3.1 的信号/噪声 4.3 倍就是用它算的），
    # 细胞级作为对照细棒一并画出，让读者看到"按细胞算是多乐观"。
    #
    # ⚠ NaN 的处理：
    #   · 供者级 SE 对 81 个单供者类型是 NaN（跨供者不确定性无法估计）
    #   · 细胞级 SE 对 1 个类型（QC 后仅 1 个细胞）是 NaN
    #   NaN **不能填 0**——matplotlib 会画出一根**零长度误差棒**，
    #   看起来像"误差极小"，比不画更误导。
    #   `Axes.errorbar` 内部会过滤非有限值，NaN 对应的棒自动不画。
    labels = list(prof.index)
    se_cell = (prof["TAI_se"].values.astype(float)
               if "TAI_se" in prof.columns else None)
    se_donor = (prof["TAI_se_donor"].values.astype(float)
                if "TAI_se_donor" in prof.columns else None)

    n_no_cell_se = 0
    if se_cell is not None:
        n_no_cell_se = int(np.isnan(se_cell).sum())
        for i in np.where(np.isnan(se_cell))[0]:
            labels[i] = str(labels[i]) + " *"

    ax.barh(y, prof["TAI"].values, color=colors, edgecolor="none", zorder=1)

    # 细胞级（对照）：浅色细棒
    if se_cell is not None:
        ax.errorbar(prof["TAI"].values, y, xerr=se_cell, fmt="none",
                    ecolor="#8fb8d4", elinewidth=0.6, capsize=1.2,
                    alpha=0.95, zorder=2)
    # 供者级（主口径）：深色粗棒
    if se_donor is not None:
        ax.errorbar(prof["TAI"].values, y, xerr=se_donor, fmt="none",
                    ecolor="#333333", elinewidth=1.0, capsize=2.4, zorder=3)

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=3.2)
    ax.invert_yaxis()
    ax.set_xlabel("TAI（越小越古老）")

    note = ""
    if se_donor is not None and se_cell is not None:
        note += "　（深色棒=供者级标准误（主口径），浅色棒=细胞级标准误（对照）"
    elif se_cell is not None:
        note += "　（横线为细胞级标准误"
    if "low_confidence" in prof.columns:
        note += "；灰色为细胞数不足 20 的低置信类型"
    if n_no_cell_se:
        note += (f"；标 * 的 {n_no_cell_se} 个类型细胞数过少，"
                 f"细胞级标准误无法估计")
    if se_donor is not None and se_cell is not None:
        n_single = int((prof["n_donors"] < 2).sum()) if "n_donors" in prof.columns else 0
        med_d = np.nanmedian(se_donor)
        med_c = np.nanmedian(se_cell)
        note += (f"；供者级中位 {med_d:.4f}，是细胞级 {med_c:.4f} 的 "
                 f"{med_d / med_c:.2f} 倍")
        if n_single:
            note += (f"，故 {n_single} 类仅 1 名供者、供者级标准误不可估，"
                     f"只画了浅色对照棒")
    if note:
        note += "）"
    ax.set_title("")
    out = os.path.join(FIG_DIR, "fig5_tai_rank.png")
    fig.savefig(out)
    plt.close(fig)
    ok(out)


# ---------- 图 6：UMAP 按 TAI 上色 ----------
def fig6(prof, group_key="cell_type", max_points=80000):
    """
    在细胞水平上色，看"年龄"是否在细胞图谱上成片分布。

    注意：TAI 是按细胞类型算的，所以图上同一个细胞类型的所有点颜色相同。
    这张图要表达的不是"单个细胞有多老"，而是"古老基因主导的细胞类型
    在图谱上是否聚集"。
    """
    from common import open_source
    from scipy import sparse

    with open_source() as ts:
        cols = ts.obs_columns()
        key = group_key if group_key in cols else "cell_ontology_class"
        labels = np.asarray(ts.obs(key)).astype(str)

        if hasattr(ts, "umap"):
            # 完整 h5ad 里有官方算好的 UMAP 坐标
            emb = ts.umap()
            ax_names = ("UMAP1", "UMAP2")
            emb_note = "UMAP"
        else:
            # 远程切片没有 obsm，改用稀疏 SVD 的前两个主成分当"细胞图谱"
            from sklearn.decomposition import TruncatedSVD
            X = ts.X.copy()
            tot = np.asarray(X.sum(axis=1)).ravel()
            scale = np.where(tot > 0, 1e4 / np.maximum(tot, 1e-12), 0.0)
            X = sparse.diags(scale.astype(np.float32)) @ X
            X = X.tocsr()
            X.data = np.log1p(X.data)
            info("切片没有 UMAP 坐标，改用稀疏 SVD 的 PC1/PC2")
            svd = TruncatedSVD(n_components=2, random_state=42)
            emb = svd.fit_transform(X)
            ax_names = ("PC1", "PC2")
            emb_note = "PC"

    tai_map = prof["TAI"].to_dict()
    tai_cell = np.array([tai_map.get(l, np.nan) for l in labels])

    valid = ~np.isnan(tai_cell)
    if valid.sum() == 0:
        info("跳过图6：细胞类型名对不上（检查 --group-key）")
        return

    idx = np.where(valid)[0]
    if len(idx) > max_points:
        rng = np.random.default_rng(42)
        idx = rng.choice(idx, max_points, replace=False)

    fig, ax = plt.subplots(figsize=(8.5, 7.5))
    sc = ax.scatter(emb[idx, 0], emb[idx, 1], c=tai_cell[idx],
                    cmap="coolwarm_r", s=2.0, alpha=0.6, linewidths=0, rasterized=True)
    cb = fig.colorbar(sc, ax=ax, pad=0.01, aspect=30)
    cb.set_label("TAI（越小越古老）")
    ax.set_xlabel(ax_names[0])
    ax.set_ylabel(ax_names[1])
    ax.set_title("")
    ax.set_xticks([])
    ax.set_yticks([])
    out = os.path.join(FIG_DIR, "fig6_umap_tai.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    ok(out)


# ---------- 图 7：PCA 零假设诊断 ----------
def fig7():
    """两条信息：
       左 —— 真实前3累计 vs 三条零假设基线的分布（说明 61.9% 不是机械效应）
       右 —— 逐 PC 对比（真实 EVR vs 基线 EVR），突出 PC5/PC6 的 z 为**负**
             → 结构集中在前 4 个成分，高阶成分反而比零假设更平
    """
    sum_path = os.path.join(ANALYSIS_DIR, "pca_summary.tsv")
    diag_path = os.path.join(ANALYSIS_DIR, "pca_diagnosis.tsv")
    if not (os.path.exists(sum_path) and os.path.exists(diag_path)):
        info("跳过图7：先跑 08_主成分诊断.py")
        return

    ps = pd.read_csv(sum_path, sep="\t", index_col=0).iloc[0]
    dg = pd.read_csv(diag_path, sep="\t")

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5.4))

    # ---- 左：真实值 vs 三条基线 ----
    real = float(ps["cum3_real"]) * 100
    bases = [
        ("逐列打乱\n（破坏层间共变）", float(ps["cum3_perm_mean"]) * 100,
         float(ps["cum3_perm_sd"]) * 100, "#4c72b0"),
        ("独立正态\n（绝对无结构）", float(ps["cum3_norm_mean"]) * 100,
         float(ps["cum3_norm_mean"]) * 0 + float(ps["cum3_perm_sd"]) * 100, "#55a868"),
        ("Dirichlet 成分\n（同秩对照）", float(ps["cum3_dirichlet_mean"]) * 100,
         float(ps["cum3_dirichlet_sd"]) * 100, "#c44e52"),
    ]
    xpos = 0
    ticks, ticklabels = [], []
    for name, mu, sd, color in bases:
        a1.bar(xpos, mu, width=0.62, color=color, alpha=0.82,
               yerr=sd, capsize=4, ecolor="#333333")
        a1.text(xpos, mu + sd + 1.6, f"{mu:.2f}%", ha="center",
                fontsize=9, fontweight="bold")
        ticks.append(xpos)
        ticklabels.append(name)
        xpos += 1
    a1.bar(xpos, real, width=0.62, color="#111111", alpha=0.9)
    a1.text(xpos, real + 1.6, f"{real:.1f}%", ha="center",
            fontsize=10, fontweight="bold")
    ticks.append(xpos)
    ticklabels.append("真实数据")

    a1.set_xticks(ticks)
    a1.set_xticklabels(ticklabels, fontsize=8)
    a1.set_ylabel("前 3 个主成分累计解释方差（%）")
    a1.set_ylim(0, real * 1.25)
    a1.axhline(3 / 16 * 100, color="#888888", ls=":", lw=1.2)
    a1.text(0.02, 3 / 16 * 100 + 1.2, "理论上限 18.8%（16 维各自独立）",
            fontsize=7.5, color="#666666")
    a1.set_title(f"左：真实值 vs 三条零假设基线\n"
                 f"z = {ps['z_vs_perm']:.1f}（逐列打乱） / "
                 f"z = {ps['z_vs_dirichlet']:.1f}（Dirichlet 同秩）", fontsize=10)
    a1.grid(axis="y", alpha=0.25)

    # ---- 右：逐 PC ----
    n_show = 8
    g = dg.head(n_show)
    idx = np.arange(n_show)
    a2.bar(idx - 0.19, g["explained_variance_ratio"] * 100, width=0.36,
           color="#111111", label="真实数据")
    a2.bar(idx + 0.19, g["baseline_perm_mean"] * 100, width=0.36,
           color="#4c72b0", alpha=0.85, label="逐列打乱基线均值")
    a2.plot(idx + 0.19, g["baseline_perm_p95"] * 100, "_", color="#222222",
            ms=14, mew=1.6, label="基线 95 分位")
    a2.axhline(0, color="k", lw=0.8)
    for i, z in enumerate(g["z_vs_perm"]):
        a2.text(i, max(g["explained_variance_ratio"].iloc[i],
                       g["baseline_perm_mean"].iloc[i]) * 100 + 1.4,
                f"z={z:+.1f}", ha="center", fontsize=7.5,
                color="#c44e52" if z < 0 else "#333333")
    a2.set_xticks(idx)
    a2.set_xticklabels(g["component"], fontsize=8)
    a2.set_xlabel("主成分")
    a2.set_ylabel("解释方差（%）")
    a2.set_ylim(min(-2.5, float(g["z_vs_perm"].min()) * 0.35), 45)
    a2.set_title("右：逐主成分对比（注意 PC5/PC6 的 z 为负 →\n"
                 "结构集中在前 4 个成分，高阶成分反而比零假设更平）", fontsize=10)
    a2.legend(fontsize=8, loc="upper right")
    a2.grid(axis="y", alpha=0.25)

    fig.tight_layout()
    out = os.path.join(FIG_DIR, "fig7_pca_null_baseline.png")
    fig.savefig(out, dpi=180)
    plt.close(fig)
    ok(out)


# ---------- 图 8：三重稳健性逐层对照 ----------
def fig8():
    rob_path = os.path.join(ANALYSIS_DIR, "robustness_comparison.tsv")
    if not os.path.exists(rob_path):
        info("跳过图8：先跑 06_稳健性对比.py")
        return
    rob = pd.read_csv(rob_path, sep="\t").sort_values("stratum")

    runs = [("A", "A 全量", "#4c72b0"),
            ("B", "B 共识度过滤", "#dd8452"),
            ("C", "C 短蛋白剔除", "#55a868")]

    fig, ax = plt.subplots(figsize=(11, 6))
    x = np.arange(len(rob))
    w = 0.26
    # 层号 → 位置索引的映射（x 轴是 0 起始的位置，不是层号！）
    pos_of = {int(s): i for i, s in enumerate(rob["stratum"])}
    handles = []
    for k, (tag, label, color) in enumerate(runs):
        vals = rob[f"fold_{tag}"].values
        # 不显著的层用浅色，显著层用实色
        cols = []
        for i in range(len(rob)):
            d = rob[f"dir_{tag}"].iloc[i]
            cols.append(color if d != "不显著" else "#dddddd")
        bars = ax.bar(x + (k - 1) * w, vals, width=w, color=cols,
                      edgecolor="none", label=label)
        # 图例必须用**纯色**方块，不能用第一个 bar 的 facecolor
        # （第一个 bar 可能刚好是不显著的浅灰，图例就会显示成灰色）
        handles.append(plt.Rectangle((0, 0), 1, 1, facecolor=color,
                                     edgecolor="none"))
        for i, b in enumerate(bars):
            v = vals[i]
            d = rob[f"dir_{tag}"].iloc[i]
            if d == "富集":
                ax.text(b.get_x() + b.get_width() / 2, v + 0.06, "★",
                        ha="center", fontsize=8, color="#c44e52")
            elif d == "耗竭":
                ax.text(b.get_x() + b.get_width() / 2, v - 0.22, "▼",
                        ha="center", fontsize=8, color="#4c72b0")

    ax.axhline(1, color="k", lw=0.9, ls="--")
    ax.set_xticks(x)
    ax.set_xticklabels([f"S{int(s)}" for s in rob["stratum"]], fontsize=8)
    ax.set_xlabel("基因年龄层（1 = 最古老，16 = 智人）")
    ax.set_ylabel("富集倍数")
    ax.set_ylim(0, float(rob[[f"fold_{t}" for t, _, _ in runs]].values.max()) * 1.2)

    # 标出三次共同显著的层（用位置索引，不是层号）
    for s in [8, 10]:
        i = pos_of.get(s)
        if i is not None:
            ax.axvspan(i - 0.45, i - 0.34, color="#c44e52", alpha=0.45, zorder=0)
    for s in [2]:
        i = pos_of.get(s)
        if i is not None:
            ax.axvspan(i - 0.45, i - 0.34, color="#4c72b0", alpha=0.45, zorder=0)

    ax.set_title("")
    ax.legend(handles=handles, fontsize=8)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    out = os.path.join(FIG_DIR, "fig8_robustness_layers.png")
    fig.savefig(out, dpi=180)
    plt.close(fig)
    ok(out)


# ---------- 图 9：短蛋白前提检验 ----------
def fig9():
    spb_path = os.path.join(ANALYSIS_DIR, "short_protein_by_stratum.tsv")
    if not os.path.exists(spb_path):
        info("跳过图9：先跑 06_稳健性对比.py")
        return
    tab = pd.read_csv(spb_path, sep="\t").sort_values("stratum")
    s = tab["stratum"].values
    pct = tab["pct_short"].values
    med = tab["median_length"].values

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 5))

    # ---- 左：短蛋白占比 ----
    cols = []
    for v in s:
        if v <= 4:
            cols.append("#4c72b0")
        elif v >= 13:
            cols.append("#c44e52")
        else:
            cols.append("#aaaaaa")
    a1.bar(s, pct, color=cols, width=0.72)
    old = float(tab[tab["stratum"] <= 4]["pct_short"].mean())
    young = float(tab[tab["stratum"] >= 13]["pct_short"].mean())
    ymax = float(pct.max()) * 1.42
    a1.set_ylim(0, ymax)
    a1.axhline(old, color="#4c72b0", ls="--", lw=1.2)
    a1.axhline(young, color="#c44e52", ls="--", lw=1.2)
    # 标注放在**图内**（靠右），避免跑出画布被裁切
    a1.text(len(s) - 0.4, old + ymax * 0.03,
            f"老层 S1–S4 均值 {old:.2f}%",
            ha="right", va="bottom", fontsize=9, color="#2b5c9a")
    a1.text(len(s) - 0.4, young + ymax * 0.03,
            f"年轻层 S13–S16 均值 {young:.2f}%（{young/old:.1f} 倍）",
            ha="right", va="bottom", fontsize=9, color="#a83232")
    a1.set_xlabel("基因年龄层（1 = 最古老）")
    a1.set_ylabel("蛋白长度 < 100 aa 的基因占比（%）")
    a1.set_xticks(s)
    a1.set_title("左：短蛋白在各年龄层的占比\n"
                 "（蓝 = 老层 S1–S4，灰 = 中间层，红 = 年轻层 S13–S16）", fontsize=10)
    a1.grid(axis="y", alpha=0.25)

    # ---- 右：中位蛋白长度 ----
    a2.plot(s, med, "o-", color="#333333", lw=1.6, ms=5)
    a2.set_xlabel("基因年龄层（1 = 最古老）")
    a2.set_ylabel("中位蛋白长度（aa）")
    a2.set_xticks(s)
    a2.invert_yaxis()
    a2.set_title("右：各层的中位蛋白长度\n"
                 "（越年轻 → 序列越短，与左图同向）", fontsize=10)
    a2.grid(alpha=0.25)

    fig.suptitle("")
    fig.tight_layout()
    out = os.path.join(FIG_DIR, "fig9_shortprot_premise.png")
    fig.savefig(out, dpi=180)
    plt.close(fig)
    ok(out)


# ---------- 图 10：长蛋白版 TAI 对照（检验 D） ----------
def fig10():
    """
    检验 D 的配图：把「短基因换成同层其余长基因」后，TAI 绝对值怎么变。

    左图：主口径 TAI vs 长蛋白版 TAI 的散点 + y=x 参考线。
          全部 177 点都在 y=x 下方 → 绝对年龄被系统性低估（整体下移）。
    右图：ΔTAI（= 长蛋白版 − 主口径）按主口径 TAI 排序的散点。
          若偏倚是「随年龄梯度增强」，右图应呈明显的上升/下降趋势；
          实测梯度 ρ = −0.1326（p = 0.079，不显著）→ 是**平移**不是**梯度**。
    """
    lpc_path = os.path.join(ANALYSIS_DIR, "longprot_comparison.tsv")
    if not os.path.exists(lpc_path):
        info("跳过图10：先跑 03_细胞类型年龄谱.py --min-length 100 --tag longprot "
             "与 94_长蛋白对照.py")
        return
    tab = pd.read_csv(lpc_path, sep="\t")
    tb = tab["TAI_base"].values
    tl = tab["TAI_longprot"].values
    d = tab["delta"].values

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 5))

    # ---- 左：主口径 vs 长蛋白版 ----
    lo = float(min(tb.min(), tl.min())) - 0.06
    hi = float(max(tb.max(), tl.max())) + 0.06
    a1.plot([lo, hi], [lo, hi], "--", color="#999999", lw=1.2,
            label="y = x（无偏倚）")
    a1.scatter(tb, tl, s=16, alpha=0.62, color="#4c72b0",
               edgecolor="none", zorder=3)
    # 回归线（仅作视觉参考）
    k, b = np.polyfit(tb, tl, 1)
    xs = np.linspace(lo, hi, 50)
    a1.plot(xs, k * xs + b, "-", color="#c44e52", lw=1.6,
            label=f"回归线 y = {k:.3f}x + {b:.3f}")
    a1.set_xlim(lo, hi)
    a1.set_ylim(lo, hi)
    a1.set_xlabel("主口径 TAI（全部年龄注释基因）")
    a1.set_ylabel("长蛋白版 TAI（仅 ≥100 aa 基因）")
    a1.set_title("左：177 种细胞类型的 TAI 对照\n"
                 f"（全部落在 y=x 下方 → 绝对年龄系统性低估；"
                 f"Spearman ρ = 0.9845）", fontsize=10)
    a1.legend(fontsize=8.5, loc="upper left")
    a1.grid(alpha=0.25)

    # ---- 右：ΔTAI vs 主口径 TAI ----
    order = np.argsort(tb)
    xs2 = tb[order]
    ds = d[order]
    a2.axhline(0.0, color="#999999", lw=1.0)
    a2.scatter(xs2, ds, s=16, alpha=0.65, color="#55a868",
               edgecolor="none", zorder=3)
    k2, b2 = np.polyfit(xs2, ds, 1)
    xs3 = np.linspace(xs2.min(), xs2.max(), 50)
    a2.plot(xs3, k2 * xs3 + b2, "-", color="#c44e52", lw=1.6,
            label=f"拟合斜率 {k2:+.4f}")
    a2.set_xlabel("主口径 TAI（从古老到年轻）")
    a2.set_ylabel("ΔTAI = 长蛋白版 − 主口径")
    a2.set_title("右：偏倚是否随年龄增强？\n"
                 "（无明显趋势 → 整体平移而非梯度；ρ = −0.1326，p = 0.079）",
                 fontsize=10)
    a2.legend(fontsize=8.5, loc="lower right")
    a2.grid(alpha=0.25)

    fig.suptitle("")
    fig.tight_layout()
    out = os.path.join(FIG_DIR, "fig10_longprot_comparison.png")
    fig.savefig(out, dpi=180)
    plt.close(fig)
    ok(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--all", action="store_true", help="含 UMAP（需读 h5ad）")
    ap.add_argument("--group-key", default="cell_ontology_class")
    args = ap.parse_args()

    require_file(PROFILE_PATH, "请先运行 03_细胞类型年龄谱.py")
    prof = pd.read_csv(PROFILE_PATH, sep="\t", index_col=0)
    info(f"{len(prof)} 种细胞类型")

    fig1(prof)
    fig2()
    fig3(prof)
    fig4()
    fig5(prof)
    fig7()
    fig8()
    fig9()
    fig10()
    if args.all:
        fig6(prof, args.group_key)

    print("\n" + "=" * 50)
    print(f"图片输出目录：{FIG_DIR}")
    print("=" * 50)


if __name__ == "__main__":
    main()
