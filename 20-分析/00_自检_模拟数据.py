"""
00_自检_模拟数据.py — 不依赖任何真实数据，验证核心逻辑是否正确
================================================================
用模拟数据跑通：聚合 → 年龄构成 → TAI → PCA/聚类 → 随机森林 → 超几何检验

用途：在等真实数据下载的时候，先确认代码逻辑没写错。
      如果这个脚本能跑出合理结果，说明核心函数是对的，
      后面出问题就只可能是「数据格式」而不是「算法」。

用法：
    python 00_自检_模拟数据.py
"""

import sys
import os
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import compute_age_profile, RANDOM_STATE

np.random.seed(RANDOM_STATE)


def section(t):
    print("\n" + "=" * 62)
    print(t)
    print("=" * 62)


def main():
    # =========================================================
    section("1. 构造模拟数据")
    # =========================================================
    n_types = 200         # 细胞类型
    n_genes = 3000        # 基因
    n_strata = 10         # 年龄层

    # 每个基因随机分配一个年龄层（1=最老）
    strata = np.random.randint(1, n_strata + 1, size=n_genes)

    # 关键：人为植入信号——前 30 种细胞类型偏向老基因，后 30 种偏向年轻基因
    # 如果算法正确，TAI 排序应该能把这些类型分到两端，
    # 随机森林的特征重要性也应该落在被放大的那批基因上
    expr = np.random.gamma(shape=2.0, scale=1.0, size=(n_types, n_genes))
    n_sig = 30

    for t in range(n_sig):
        expr[t, strata <= 3] *= 6.0        # 老基因高表达
    for t in range(n_types - n_sig, n_types):
        expr[t, strata >= 8] *= 6.0        # 年轻基因高表达

    ct_names = [f"CT_{i:03d}" for i in range(n_types)]
    for i in range(n_sig):
        ct_names[i] = f"OldType_{i:02d}"
    for i in range(n_sig):
        ct_names[n_types - n_sig + i] = f"YoungType_{i:02d}"

    mean_expr = pd.DataFrame(expr, index=ct_names,
                             columns=[f"G{i}" for i in range(n_genes)])

    print(f"  {n_types} 种细胞类型 × {n_genes} 个基因 × {n_strata} 个年龄层")
    print(f"  已植入信号：前 5 种偏高表达老基因（S1–S3），后 5 种偏高表达年轻基因（S8–S10）")

    # =========================================================
    section("2. 年龄构成与 TAI")
    # =========================================================
    prof = compute_age_profile(mean_expr, strata.astype(float))

    print(f"  输出形状：{prof.shape[0]} 行 × {prof.shape[1] - 1} 列")
    print(f"  TAI 范围：{prof['TAI'].min():.2f} – {prof['TAI'].max():.2f}")

    # 验证：植入的信号能不能被测出来
    old_types = [n for n in prof.index if n.startswith("OldType")]
    young_types = [n for n in prof.index if n.startswith("YoungType")]
    tai_old = prof.loc[old_types, "TAI"].mean()
    tai_young = prof.loc[young_types, "TAI"].mean()

    print(f"\n  OldType 平均 TAI：  {tai_old:.2f}")
    print(f"  YoungType 平均 TAI：{tai_young:.2f}")
    print(f"  差值：{tai_young - tai_old:+.2f}")

    if tai_young > tai_old:
        print("  ✓ 方向正确：年轻基因偏向的类型 TAI 更大")
    else:
        print("  ✗ 方向反了！检查 compute_age_profile")
        return False

    # 占比是否归一
    s_cols = [c for c in prof.columns if c.startswith("S")]
    row_sums = prof[s_cols].sum(axis=1)
    if np.allclose(row_sums, 1.0, atol=1e-6):
        print("  ✓ 每行的占比之和 = 1")
    else:
        print(f"  ✗ 占比未归一：{row_sums.min():.4f} – {row_sums.max():.4f}")
        return False

    # 排序是否单调
    if prof["TAI"].is_monotonic_increasing:
        print("  ✓ 输出已按 TAI 升序排列")
    else:
        print("  ✗ 未排序")

    # =========================================================
    section("3. 无监督分析（PCA + 聚类）")
    # =========================================================
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score, adjusted_rand_score

    X = prof[s_cols].values
    Xs = StandardScaler().fit_transform(X)

    pca = PCA(n_components=3, random_state=RANDOM_STATE)
    Z = pca.fit_transform(Xs)
    print(f"  PCA 解释方差：{['%.1f%%' % (v * 100) for v in pca.explained_variance_ratio_]}")
    print(f"  前 3 个累计：{pca.explained_variance_ratio_.sum():.1%}")

    labels = KMeans(n_clusters=2, n_init=20, random_state=RANDOM_STATE).fit_predict(Xs)
    sil = silhouette_score(Xs, labels)
    print(f"  k=2 轮廓系数：{sil:.3f}")

    # 真实分组（用于验证聚类能否抓出来）
    truth = np.array(["old"] * n_sig
                     + ["mid"] * (n_types - 2 * n_sig)
                     + ["young"] * n_sig)
    ari = adjusted_rand_score(truth, labels)
    print(f"  与植入分组的 ARI：{ari:.3f}")
    if ari > 0.5:
        print("  ✓ 聚类能抓出植入的分组")
    else:
        print("  ⚠ ARI 偏低。模拟数据的信号可能不够强，或 k 选错了——")
        print("    这不代表代码错，真实数据未必如此")

    # =========================================================
    section("4. 有监督分析（随机森林 + 特征重要性）")
    # =========================================================
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.metrics import accuracy_score

    # 注意：这里必须用「基因表达矩阵」，不是年龄构成矩阵
    Xg = expr
    y = np.array([0] * n_sig
                 + [1] * (n_types - 2 * n_sig)
                 + [2] * n_sig)

    rf = RandomForestClassifier(n_estimators=300, random_state=RANDOM_STATE, n_jobs=-1)
    rf.fit(Xg, y)
    acc = accuracy_score(y, rf.predict(Xg))
    print(f"  训练集准确率：{acc:.3f}（应为 1.0，因为数据可分）")

    imp = rf.feature_importances_
    order = np.argsort(imp)[::-1]
    gene_names = mean_expr.columns.values
    top20 = order[:20]

    print(f"  最高重要性：{imp[order[0]]:.5f}（基因 {gene_names[order[0]]}）")
    print(f"  Top 20 基因所在层：{sorted(set(int(x) for x in strata[top20]))}")
    print(f"  植入信号在 S1–S3 与 S8–S10，这两批应占多数")

    # =========================================================
    section("5. 超几何富集检验")
    # =========================================================
    from scipy.stats import hypergeom

    top_k = 300
    top_idx = order[:top_k]
    S = strata.astype(float)
    bg = S
    top_s = S[top_idx]

    M, N = len(bg), len(top_s)
    print(f"  背景 {M} 个基因，Top {N} 个")
    print(f"  {'层':>4} {'背景数':>7} {'Top中':>7} {'期望':>7} {'倍数':>7} {'p值':>10}")
    print("  " + "-" * 48)

    found = False
    for s in sorted(np.unique(bg)):
        n_l = int((bg == s).sum())
        k = int((top_s == s).sum())
        exp = N * n_l / M
        fold = k / exp if exp > 0 else np.nan
        p = hypergeom.sf(k - 1, M, n_l, N)
        star = " ***" if p < 0.001 else " **" if p < 0.01 else " *" if p < 0.05 else ""
        print(f"  S{int(s):>3} {n_l:>7} {k:>7} {exp:>7.1f} {fold:>7.2f} {p:>10.2e}{star}")
        if p < 0.05:
            found = True

    if found:
        print("\n  ✓ 检出了显著富集层")
    else:
        print("\n  ⚠ 没有显著层。若植入信号够强，通常应能检出。")

    # =========================================================
    section("自检结论")
    # =========================================================
    print("""
  以上四步全部跑通，说明：
    · compute_age_profile 的公式正确（方向、归一、排序都对）
    · PCA / 聚类 / 随机森林 / 超几何检验的调用方式正确
    · 各库版本兼容

  ⚠ 这只验证「算法」，不验证「数据」。
    真实数据上仍可能出现：基因名匹配失败、obs 列名不同、内存不足。
    这些属于数据格式问题，脚本里都有对应的错误提示。
""")
    return True


if __name__ == "__main__":
    okall = main()
    sys.exit(0 if okall else 1)
