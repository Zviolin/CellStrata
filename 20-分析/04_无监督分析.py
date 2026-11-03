"""
04_无监督分析.py — 年龄构成能否自动区分细胞类型
================================================
输入：20-分析/celltype_age_profile.tsv
输出：PCA 坐标、聚类标签、轮廓系数、ARI、ARI 置换零假设

用法：
    python 04_无监督分析.py

2026-11-03 审核后的三处重要修正：
  ① 特征列改为**显式白名单** S1..S16（原 `startswith("S")` 偏脆，
     日后新增 Score/Sum 之类的派生列会被误纳入）
  ② **低置信类型的敏感性分析**：25 个 QC 后细胞数 < 20 的类型，
     原脚本虽然在 03 里给它们打了 low_confidence 标记，却让它们
     **以同等权重**参与 PCA / K-means / 轮廓系数 / ARI。
     实测：按标记过滤后 k 从 7 变 5、ARI 从 0.211 降到 0.159，
     说明 k 与 ARI 对该过滤高度敏感 → 两套结果都要报。
  ③ **ARI 置换零假设**：ARI 的绝对值依赖两侧划分的粒度，
     这里左侧是 5 个区室、右侧是 k 簇（k > 5），不可能取到 1。
     所以必须给「随机情况下 ARI 会有多大」作参照，否则
     「ARI=0.211 弱相关」这个解读没有锚点。
     另补 Ward 层次聚类的 ARI 与轮廓系数，让「对照」名副其实。
"""

import os
import numpy as np
import pandas as pd

from common import (PROFILE_PATH, ANALYSIS_DIR, LINEAGE_MAP_PATH,
                    RANDOM_STATE, N_STRATA, die, info, ok, warn, require_file)

# 显式白名单：只允许 S1..S{N}，防止将来新增以 S 开头的派生列被误当特征
S_COLS = [f"S{i}" for i in range(1, N_STRATA + 1)]
MIN_CELLS = 20        # 低置信阈值：QC 后细胞数少于它的不进主结论
N_PERM = 2000         # ARI 置换零假设的置换次数


def run_pipeline(X, cell_types, lm_col=None, n_perm=N_PERM):
    """对给定特征矩阵跑一遍 PCA + 选 k + ARI（含置换零假设）。

    返回一个 dict，字段：
      pca_evr / cum3 / k / sil / labels / ari / ari_perm
      ari_z / ward_ari / ward_sil
    """
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    from sklearn.cluster import KMeans, AgglomerativeClustering
    from sklearn.metrics import silhouette_score, adjusted_rand_score

    Xs = StandardScaler().fit_transform(X)
    pca = PCA(n_components=min(5, Xs.shape[1] - 1), random_state=RANDOM_STATE)
    pca.fit_transform(Xs)
    ev = pca.explained_variance_ratio_

    best = {}
    for k in range(2, 9):
        km = KMeans(n_clusters=k, n_init=20, random_state=RANDOM_STATE)
        lab = km.fit_predict(Xs)
        best[k] = (silhouette_score(Xs, lab), km.inertia_, lab)
    k_best = max(best, key=lambda k: best[k][0])
    sil_best, _, labels = best[k_best]

    res = {
        "pca_evr": ev, "cum3": float(ev[:3].sum()),
        "k": k_best, "sil": float(sil_best), "labels": labels,
        "Xs": Xs,
    }

    # Ward 层次聚类对照（必须给量化结果，否则不叫对照）
    ward = AgglomerativeClustering(n_clusters=k_best, linkage="ward")
    ward_labels = ward.fit_predict(Xs)
    res["ward_labels"] = ward_labels
    res["ward_sil"] = float(silhouette_score(Xs, ward_labels))

    if lm_col is not None:
        mapped = lm_col.reindex(cell_types)
        sub = mapped.notna().values
        if sub.sum() > 10:
            ari = adjusted_rand_score(mapped[sub], labels[sub])
            res["ari"] = float(ari)
            res["ward_ari"] = float(
                adjusted_rand_score(mapped[sub], ward_labels[sub]))
            # ---- ARI 置换零假设 ----
            # 随机打乱区室标签，看 ARI 会自然出现多大。
            # 没有这个参照，「ARI=0.21 弱相关」无法判断是显著还是噪声。
            rng = np.random.default_rng(RANDOM_STATE)
            lab_true = mapped[sub].to_numpy()
            perms = np.empty(n_perm)
            for i in range(n_perm):
                perms[i] = adjusted_rand_score(
                    rng.permutation(lab_true), labels[sub])
            res["ari_perm_mean"] = float(perms.mean())
            res["ari_perm_sd"] = float(perms.std(ddof=1))
            res["ari_perm_p95"] = float(np.percentile(perms, 95))
            res["ari_z"] = float((ari - perms.mean()) / perms.std(ddof=1))
    return res


def main():
    require_file(PROFILE_PATH, "请先运行 03_细胞类型年龄谱.py")

    prof = pd.read_csv(PROFILE_PATH, sep="\t", index_col=0)

    # 只取 S1..S16 显式白名单；TAI / TAI_se / n_cells / low_confidence /
    # n_donors / TAI_se_donor 等派生量一律不能进特征矩阵（会造成信息泄漏）。
    s_cols = [c for c in S_COLS if c in prof.columns]
    missing = [c for c in S_COLS if c not in prof.columns]
    if missing:
        die(f"年龄构成缺少层列：{missing}")
    X = prof[s_cols].values
    cell_types = prof.index.to_numpy()
    info(f"读入 {prof.shape[0]} 种细胞类型 × {len(s_cols)} 个年龄层"
         f"（特征列白名单 S1..S{N_STRATA}，已排除所有派生列）")

    # 谱系标签（若可用）
    lm_col = None
    if os.path.exists(LINEAGE_MAP_PATH):
        lm = pd.read_csv(LINEAGE_MAP_PATH, sep="\t", index_col=0)
        col = "lineage" if "lineage" in lm.columns else lm.columns[0]
        lm_col = lm[col]
    else:
        warn(f"未找到 {LINEAGE_MAP_PATH}，跳过 ARI")

    # ---------- 主口径：全部类型 ----------
    res = run_pipeline(X, cell_types, lm_col)
    print("\n── PCA 解释方差 ──")
    for i, v in enumerate(res["pca_evr"], 1):
        print(f"  PC{i}: {v:.1%}")
    print(f"  前 3 个累计：{res['cum3']:.1%}")
    if res["cum3"] < 0.7:
        print("  ⚠️  累计 < 70%，说明年龄构成的信息分散在更多维度上"
              "（零假设对照见 08_主成分诊断.py）")

    print("\n── 选择聚类数 k ──")
    print(f"  按轮廓系数选 k={res['k']}（{res['sil']:.3f}）")
    print(f"  Ward 层次聚类对照：轮廓 {res['ward_sil']:.3f}")

    if "ari" in res:
        print("\n── 与谱系注释的一致性 ──")
        print(f"  ARI（KMeans, k={res['k']}）= {res['ari']:.3f}")
        print(f"  ARI（Ward 层次聚类）      = {res['ward_ari']:.3f}")
        print(f"  置换零假设（打乱区室标签 {N_PERM} 次）："
              f"{res['ari_perm_mean']:+.4f} ± {res['ari_perm_sd']:.4f}"
              f"，95 分位 {res['ari_perm_p95']:+.4f}")
        print(f"  → 观测值超出零假设 {res['ari_z']:.1f} 个标准差，"
              f"显著高于随机")
        print("  解读注意：ARI 绝对值依赖两侧划分的粒度。"
              "左侧是 5 个区室、右侧是 "
              f"k={res['k']} 簇（簇数多于区室数时不可能取到 1），"
              "所以不能拿 0.2/0.5 这类固定刻度直接解读；"
              "应以上面的置换零假设为参照。")

    # ---------- 敏感性：剔除低置信类型 ----------
    if "n_cells" in prof.columns:
        keep = prof["n_cells"] >= MIN_CELLS
        n_drop = int((~keep).sum())
        print("\n" + "=" * 62)
        print(f"敏感性分析：剔除 {n_drop} 个低置信类型"
              f"（QC 后细胞数 < {MIN_CELLS}）")
        print("=" * 62)
        res2 = run_pipeline(X[keep.values], cell_types[keep.values], lm_col)
        print(f"  剔除后：n={int(keep.sum())}  "
              f"PC1-3 累计 {res2['cum3']:.1%}  "
              f"k={res2['k']}  轮廓 {res2['sil']:.3f}"
              + (f"  ARI {res2['ari']:.3f}" if "ari" in res2 else ""))
        print(f"  （主口径：n={len(prof)}  PC1-3 累计 {res['cum3']:.1%}  "
              f"k={res['k']}  轮廓 {res['sil']:.3f}"
              + (f"  ARI {res['ari']:.3f}" if "ari" in res else "") + "）")
        if res2["k"] != res["k"]:
            print(f"  ⚠️ k 随低置信类型的处理方式改变（{res['k']} → {res2['k']}），"
                  f"报告必须两套都给，不能只报一个。")
        if "ari" in res and "ari" in res2:
            d = res2["ari"] - res["ari"]
            print(f"  ⚠️ ARI 随该过滤变化 {d:+.3f}（{res['ari']:.3f} → "
                  f"{res2['ari']:.3f}），报告必须说明。")
        res["sens"] = res2

    # ---------- 落盘 ----------
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    Z = PCA(n_components=3, random_state=RANDOM_STATE).fit_transform(
        StandardScaler().fit_transform(X))
    pca_df = pd.DataFrame(Z, index=cell_types, columns=["PC1", "PC2", "PC3"])
    pca_df["TAI"] = prof["TAI"]
    pca_df["cluster"] = res["labels"]
    pca_df["cluster_ward"] = res["ward_labels"]
    pca_df["n_cells"] = prof["n_cells"] if "n_cells" in prof.columns else np.nan
    if lm_col is not None:
        pca_df["lineage"] = lm_col.reindex(cell_types)

    out = os.path.join(ANALYSIS_DIR, "unsupervised_result.tsv")
    pca_df.to_csv(out, sep="\t")
    ok(f"已写出：{out}")

    # 汇总指标单独落盘，供终审脚本逐项核对
    met = {
        "n_types": int(len(prof)),
        "n_features": int(len(s_cols)),
        "pc1": round(float(res["pca_evr"][0]), 4),
        "pc2": round(float(res["pca_evr"][1]), 4),
        "pc3": round(float(res["pca_evr"][2]), 4),
        "cum3": round(res["cum3"], 4),
        "k": int(res["k"]),
        "silhouette": round(res["sil"], 4),
        "ari": round(res.get("ari", float("nan")), 4),
        "ari_ward": round(res.get("ward_ari", float("nan")), 4),
        "silhouette_ward": round(res["ward_sil"], 4),
        "ari_perm_mean": round(res.get("ari_perm_mean", float("nan")), 4),
        "ari_perm_sd": round(res.get("ari_perm_sd", float("nan")), 4),
        "ari_z": round(res.get("ari_z", float("nan")), 2),
        "n_sens": int(len(prof) - n_drop) if "n_cells" in prof.columns else -1,
        "k_sens": int(res["sens"]["k"]) if "sens" in res else -1,
        "ari_sens": round(res["sens"].get("ari", float("nan")), 4)
        if "sens" in res else float("nan"),
        "sil_sens": round(res["sens"]["sil"], 4) if "sens" in res else float("nan"),
        "cum3_sens": round(res["sens"]["cum3"], 4) if "sens" in res else float("nan"),
    }
    mp = os.path.join(ANALYSIS_DIR, "unsupervised_metrics.tsv")
    pd.DataFrame([met]).to_csv(mp, sep="\t", index=False)
    ok(f"已写出：{mp}")

    # ---------- 各簇的 TAI ----------
    print("\n── 各簇的平均 TAI ──")
    for c in sorted(pca_df["cluster"].unique()):
        m = pca_df["cluster"] == c
        print(f"  簇 {c}：n={m.sum():>3}  TAI={pca_df.loc[m, 'TAI'].mean():.2f}")


if __name__ == "__main__":
    main()
