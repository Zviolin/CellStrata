"""
05_监督分析.py — 细胞身份由哪个年代的基因定义
==============================================
输入：h5ad + gene_age.tsv
输出：
    20-分析/feature_importance.tsv   每个基因的重要性 + 年龄层
    20-分析/age_enrichment.tsv       各年龄层在 top 基因中的富集（超几何检验）
    20-分析/classification_report.txt 分类性能详情

用法：
    python 05_监督分析.py
    python 05_监督分析.py --n-hvg 3000 --top-k 500 --max-cells 300
    python 05_监督分析.py --min-share 0.8    # 稳健性检验 B：只用高共识度的基因
    python 05_监督分析.py --min-length 100   # 稳健性检验 C：剔除短蛋白

方法要点
    1. 按 donor_id 划分训练/测试，**不能按细胞随机划分** ——
       同一供者的细胞高度相似，随机划分会让测试集信息泄漏到训练集，
       准确率虚高十几个百分点。这是单细胞分析最常见的错误之一。
    2. 高变基因在降采样后的子集上算，用稀疏矩阵的矩估计，
       全程不把 483,152 × 58,604 稠密化。
    3. --min-share 用于响应 Moyers & Zhang 的方法学批评：
       只保留**共识度**（13 个同源推断算法中最优势层占比 mode_share）≥ 阈值的基因，
       重跑一遍，看富集结论是否还成立。
       ⚠️ **不要用 --min-dbs**：原始表的 `NumDBsContributing` 是
       **参与投票的算法总数（分母）**，不是"支持数"。实测
       corr(n_dbs, mode_share) = −0.128（负相关），
       按 n_dbs 过滤筛出来的恰恰是**低共识**基因，方向与"高置信"相反。
"""

import argparse
import os

import numpy as np
import pandas as pd
from scipy import sparse

from common import (H5AD_PATH, ANALYSIS_DIR, DATA_DIR, RANDOM_STATE,
                    die, info, ok, warn, require_file, load_gene_age, map_genes,
                    open_source)
from ts_io import TabulaSapiens


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--n-hvg", type=int, default=3000, help="高变基因数")
    p.add_argument("--top-k", type=int, default=500, help="取重要性前多少基因")
    p.add_argument("--max-cells", type=int, default=300,
                   help="每种细胞类型最多取多少（默认 300，控制内存）")
    p.add_argument("--group-key", default="cell_ontology_class")
    p.add_argument("--donor-key", default="donor_id")
    p.add_argument("--min-genes", type=int, default=200,
                   help="QC：每细胞最少检出基因数（与 03 保持一致）")
    p.add_argument("--max-mito", type=float, default=20.0,
                   help="QC：每细胞线粒体占比上限 %%（与 03 保持一致）")
    p.add_argument("--no-qc", action="store_true",
                   help="⚠ 不做 QC，直接用全部细胞（不推荐：会与 TAI 分析口径不一致）")
    p.add_argument("--min-share", type=float, default=0.0,
                   help="共识度过滤：只保留 mode_share ≥ 此值的基因"
                        "（mode_share = 支持该层号的算法比例，1=全体一致；"
                        "0=不过滤）。需先跑 02_构建年龄表.py")
    p.add_argument("--min-length", type=int, default=0,
                   help="只保留蛋白长度 ≥ 此值 aa 的基因（0=不过滤；"
                        "需先跑 93_获取蛋白长度.py）")
    p.add_argument("--tag", default="", help="输出文件名后缀，便于对比多次运行")
    return p.parse_args()


def bh_fdr(p):
    """Benjamini-Hochberg FDR 校正。

    本项目每次运行要做 16 层 × 2 方向 = 32 次超几何检验，
    跨三次运行共 96 次。不校正的话 p<0.05 会产生大量假阳性
    （实测 A 组原始口径 10 层显著，BH 校正后只剩 7 层）。
    """
    p = np.asarray(p, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    # 从后往前取累积最小值，再还原顺序
    q = np.empty(n)
    q[order] = np.minimum.accumulate(ranked[::-1])[::-1]
    return np.minimum(q, 1.0)


def pick_cells(labels, max_cells, seed):
    rng = np.random.default_rng(seed)
    keep = []
    for ct in pd.unique(labels):
        idx = np.where(labels == ct)[0]
        if len(idx) > max_cells:
            idx = rng.choice(idx, max_cells, replace=False)
        keep.append(idx)
    return np.sort(np.concatenate(keep))


def normalize(X):
    total = np.asarray(X.sum(axis=1)).ravel()
    scale = np.where(total > 0, 1e4 / np.maximum(total, 1e-12), 0.0)
    X = (sparse.diags(scale.astype(np.float32)) @ X).tocsr()
    X.data = np.log1p(X.data)
    return X


def gene_stats(ts, rows, chunk=10000):
    """
    在选中的细胞子集上算每个基因的均值、方差、检出率。

    用 np.bincount 直接在 CSR 的 indices 上累加 —— 比 X.sum(axis=0) /
    X.multiply(X).sum(axis=0) 快一个数量级，内存峰值也低得多
    （后者会额外造出两个同规模的稀疏矩阵）。
    数学上完全等价：CSR 的列和 = 按列号累加非零值。
    """
    n_vars = ts.n_vars
    s1 = np.zeros(n_vars, dtype=np.float64)
    s2 = np.zeros(n_vars, dtype=np.float64)
    det = np.zeros(n_vars, dtype=np.int64)
    n = 0

    for start in range(0, len(rows), chunk):
        idx = rows[start:start + chunk]
        if len(idx) > 1 and idx[-1] - idx[0] + 1 == len(idx):
            X = ts.raw_block(int(idx[0]), int(idx[-1]) + 1)
        else:
            X = ts.raw_rows(idx)
        X = normalize(X)
        X = X.tocsr()
        cols = X.indices
        vals = X.data
        s1 += np.bincount(cols, weights=vals, minlength=n_vars)
        s2 += np.bincount(cols, weights=vals * vals, minlength=n_vars)
        det += np.bincount(cols, minlength=n_vars)
        n += X.shape[0]
        del X, cols, vals

    mean = s1 / n
    var = s2 / n - mean ** 2
    return mean, var, det / n, n


def build_dense(ts, rows, gene_idx, chunk=20000):
    """把 rows × gene_idx 的归一化表达读成稠密 float32 矩阵。"""
    out = np.zeros((len(rows), len(gene_idx)), dtype=np.float32)
    for start in range(0, len(rows), chunk):
        idx = rows[start:start + chunk]
        if len(idx) > 1 and idx[-1] - idx[0] + 1 == len(idx):
            X = ts.raw_block(int(idx[0]), int(idx[-1]) + 1)
        else:
            X = ts.raw_rows(idx)
        X = normalize(X)
        out[start:start + len(idx)] = X[:, gene_idx].toarray()
    return out


def main():
    args = parse_args()
    age_df = load_gene_age()

    import sklearn
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.metrics import accuracy_score, f1_score
    from scipy.stats import hypergeom
    info(f"scikit-learn {sklearn.__version__}")

    with open_source() as ts:
        cols = ts.obs_columns()
        for k in (args.group_key, args.donor_key):
            if k not in cols:
                die(f"obs 缺 '{k}'。现有列：{cols}")

        labels = np.asarray(ts.obs(args.group_key)).astype(str)
        donors = np.asarray(ts.obs(args.donor_key)).astype(str)
        info(f"全量：{ts.n_obs:,} 细胞 × {ts.n_vars:,} 基因")

        # ---------- 1. QC + 降采样 ----------
        # ⚠️ **必须与 03 用同一批细胞**（2026-11-03 审核发现的致命口径问题）：
        #   原脚本不做 QC，直接对全部 28,017 个细胞建模，而 TAI/年龄谱是
        #   QC 后的 22,195 个。两套数字来自不同细胞集却写在同一篇报告里。
        #   这里直接读 03 落盘的 QC 通过行索引，保证口径一致。
        qc_rows_path = os.path.join(ANALYSIS_DIR, "qc_passed_rows.npy")
        if args.no_qc:
            rows_all = np.arange(ts.n_obs)
            info(f"⚠️ 已禁用 QC：使用全部 {len(rows_all):,} 个细胞"
                 f"（与 TAI 分析口径不一致，报告需明确说明）")
        elif os.path.exists(qc_rows_path):
            rows_all = np.load(qc_rows_path)
            info(f"QC 后细胞：{len(rows_all):,} / {ts.n_obs:,} "
                 f"= {len(rows_all)/ts.n_obs:.1%}（复用 03 的 QC 结果）")
        else:
            die(f"找不到 {qc_rows_path}\n"
                f"请先运行 03_细胞类型年龄谱.py（它会落盘 QC 通过的行索引）")

        rows = pick_cells(labels[rows_all], args.max_cells, RANDOM_STATE)
        rows = rows_all[rows]
        y = labels[rows]
        d = donors[rows]
        info(f"建模细胞：{len(rows):,}，{len(np.unique(y))} 种类型，"
             f"{len(np.unique(d))} 个供者")

        # ---------- 2. 高变基因 ----------
        info("统计基因均值/方差（分块，稀疏运算）……")
        mean, var, det, n = gene_stats(ts, rows)
        # 只在有一定检出率的基因里选，避免选到全零基因
        ok_gene = det >= 0.01
        disp = np.where(ok_gene, var, -1)
        hvg = np.argsort(disp)[::-1][:args.n_hvg]
        hvg = np.sort(hvg[disp[hvg] > 0])
        info(f"高变基因：{len(hvg)} 个")

        # ---------- 3. 基因年龄对齐 ----------
        symbols_all = ts.symbols()
        matched_all, stratum_all = map_genes(symbols_all, age_df)
        stratum = stratum_all[hvg]
        gene_names = symbols_all[hvg]
        info(f"HVG 中有年龄注释的：{int((~np.isnan(stratum)).sum()):,}")

        # ---------- 4. 共识度过滤（可选） ----------
        # ⚠️ 口径修正（2026-11-03 审核发现）：
        #   原参数 `--min-dbs` 的语义是「参与投票的算法数 ≥ N」，
        #   而报告把它写成「≥N 个算法支持该年龄调用」——**方向是反的**。
        #   实测 n_dbs 与共识度负相关（corr = −0.128）：
        #     n_dbs ≤ 4  的基因 mode_share 中位 1.000（全体一致）
        #     n_dbs ≥ 10 的基因 mode_share 中位 0.692（一致性更差）
        #   真正的共识度是 mode_share = 支持该层号的算法比例。
        #   现在改用 `--min-share`，语义与报告表述一致。
        if args.min_share > 0:
            if "mode_share" not in age_df.columns:
                die("年龄表缺 mode_share 列，请先运行 02_构建年龄表.py")
            conf = dict(zip(age_df["symbol"].astype(str).str.upper(),
                            age_df["mode_share"]))
            ms = np.array([conf.get(str(g).upper(), np.nan) for g in gene_names],
                          dtype=float)
            keep_g = (~np.isnan(ms)) & (ms >= args.min_share)
            info(f"共识度过滤（mode_share ≥ {args.min_share}）："
                 f"{int(keep_g.sum())} / {len(hvg)} 个基因保留")
            if keep_g.sum() < 500:
                die("过滤后基因太少，降低 --min-share")
            hvg, gene_names, stratum = hvg[keep_g], gene_names[keep_g], stratum[keep_g]

        # ---------- 4b. 短蛋白剔除（可选，Moyers & Zhang 建议） ----------
        # phylostratigraphy 的误差非随机：序列短的基因更容易找不到远缘同源，
        # 从而被判为"年轻"。Moyers & Zhang 建议剔除短蛋白后重跑，
        # 看已发表的结论是否消失或反转。
        if args.min_length > 0:
            plen_path = os.path.join(DATA_DIR, "protein_length.tsv")
            if not os.path.exists(plen_path):
                die(f"缺少蛋白长度表：{plen_path}\n"
                    f"请先运行 93_获取蛋白长度.py")
            pl = pd.read_csv(plen_path, sep="\t")
            length_map = dict(zip(pl["symbol"].astype(str).str.upper(),
                                  pl["length"]))
            ln = np.array([length_map.get(str(g).upper(), np.nan)
                           for g in gene_names], dtype=float)
            keep_g = (~np.isnan(ln)) & (ln >= args.min_length)
            n_short = int(((~np.isnan(ln)) & (ln < args.min_length)).sum())
            n_nolen = int(np.isnan(ln).sum())
            info(f"短蛋白剔除（长度 ≥ {args.min_length} aa）："
                 f"{int(keep_g.sum())} / {len(hvg)} 个基因保留"
                 f"（{int((~keep_g).sum())} 个因过短或无长度信息被剔除）")
            # 明细落盘：报告里「剔除 N 个过短或无长度信息」这句话必须可追溯
            reason = np.where(np.isnan(ln), "无长度信息",
                              np.where(ln < args.min_length, "过短", "保留"))
            drop_df = pd.DataFrame({
                "gene": gene_names,
                "length": ln,
                "reason": reason,
            }).sort_values(["reason", "length"])
            drop_path = os.path.join(ANALYSIS_DIR,
                                     f"filtered_genes_{args.tag}.tsv")
            drop_df.to_csv(drop_path, sep="\t", index=False)
            info(f"短蛋白剔除明细：过短 {n_short} 个 / 无长度信息 {n_nolen} 个 "
                 f"→ {drop_path}")
            if keep_g.sum() < 500:
                die("过滤后基因太少，降低 --min-length")
            hvg, gene_names, stratum = hvg[keep_g], gene_names[keep_g], stratum[keep_g]

        # ---------- 5. 稠密化 ----------
        info(f"构建稠密矩阵：{len(rows):,} × {len(hvg):,} "
             f"≈ {len(rows)*len(hvg)*4/1e9:.2f} GB")
        X = build_dense(ts, rows, hvg)

    # ---------- 6. 按供者划分 ----------
    donors_u = np.unique(d)
    rng = np.random.default_rng(RANDOM_STATE)
    rng.shuffle(donors_u)
    n_test = max(1, int(round(len(donors_u) * 0.25)))
    test_donors = set(donors_u[:n_test])
    is_test = np.isin(d, list(test_donors))

    info(f"供者划分：训练 {len(donors_u)-n_test} 人 / 测试 {n_test} 人")
    info(f"细胞划分：训练 {(~is_test).sum():,} / 测试 {is_test.sum():,}")
    print(f"  测试集供者：{sorted(str(x) for x in test_donors)}")

    if len(np.unique(y[is_test])) < 2:
        die("测试集里类别太少，无法评估。提高 --max-cells。")

    # ---------- 7. 训练 ----------
    info("训练随机森林（500 棵树）……")
    rf = RandomForestClassifier(
        n_estimators=500, min_samples_leaf=2, n_jobs=-1,
        random_state=RANDOM_STATE, class_weight="balanced_subsample",
    )
    rf.fit(X[~is_test], y[~is_test])
    pred = rf.predict(X[is_test])
    acc = accuracy_score(y[is_test], pred)

    # ⚠️ 宏平均 F1 的两种口径（2026-11-03 审核发现）
    # sklearn 的 f1_score(average="macro") 默认 labels=None，
    # 标签集 = y_true ∪ y_pred，可能包含**测试集里根本不存在**的类
    #（模型把测试细胞误判成了训练里有、测试集没有的类）。
    # 这些类 F1=0 却计入分母 → 系统性稀释。
    # 两种口径都报，标明各自含义。
    te_labels = np.unique(y[is_test])
    all_labels = np.unique(y)
    f1m_te = f1_score(y[is_test], pred, average="macro",
                       labels=te_labels, zero_division=0)
    f1m_all = f1_score(y[is_test], pred, average="macro",
                       labels=all_labels, zero_division=0)

    # ⚠️ 正确的随机基线与准确率天花板（2026-11-03 审核发现）
    # 原脚本用 1/177（**全体细胞**的类别数）作随机基线，但测试集里
    # 只有 106 类，且其中 25 类在训练集中从未出现 → 这部分细胞必然全错。
    train_labels = set(np.unique(y[~is_test]))
    unseen = sorted(set(te_labels) - train_labels)
    n_unseen_cells = int(np.isin(y[is_test], unseen).sum())
    acc_ceiling = 1 - n_unseen_cells / len(y[is_test])
    rand_te = 1 / len(te_labels)
    maj = float(pd.Series(y[is_test]).value_counts().iloc[0] / len(y[is_test]))

    lines = []
    def say(s=""):
        print(s)
        lines.append(s)

    say("\n" + "=" * 60)
    say("分类性能（按供者划分，无泄漏）")
    say("=" * 60)
    say(f"  准确率：{acc:.3f}")
    say(f"  宏平均 F1（测试集实际出现的 {len(te_labels)} 类）：{f1m_te:.3f}")
    say(f"  宏平均 F1（全部 {len(all_labels)} 类口径）：{f1m_all:.3f}")
    say(f"  类别数：{len(all_labels)}（训练集 {len(train_labels)} 类 / "
        f"测试集 {len(te_labels)} 类）")
    say(f"  特征数：{len(hvg)}")
    if args.min_share > 0:
        say(f"  共识度过滤：mode_share ≥ {args.min_share}")
    if args.min_length > 0:
        say(f"  短蛋白过滤：蛋白长度 ≥ {args.min_length} aa")

    say("\n  ── 基线与天花板（务必与准确率一起看）──")
    say(f"  随机基线（测试集 {len(te_labels)} 类均匀猜）：{rand_te:.4f}")
    say(f"  多数类基线（恒猜最高频类）：{maj:.4f}   ← 更贴近实际")
    say(f"  准确率是多数类基线的 {acc / maj:.1f} 倍")
    say(f"  ⚠️ 测试集中有 {len(unseen)} 类在训练集中从未出现，"
        f"涉及 {n_unseen_cells:,} 个细胞（占测试集 "
        f"{n_unseen_cells / len(y[is_test]):.1%}），模型不可能预测对")
    say(f"  准确率天花板（剔除这些不可预测细胞后）：{acc_ceiling:.3f}"
        f"，实际达到天花板的 {acc / acc_ceiling:.1%}")

    # ---------- 8. 特征重要性 ----------
    imp = rf.feature_importances_
    order = np.argsort(imp)[::-1]
    top_idx = order[:args.top_k]

    imp_df = pd.DataFrame({
        "gene": gene_names[order],
        "importance": imp[order],
        "stratum": stratum[order],
    })
    tag = f"_{args.tag}" if args.tag else ""
    imp_path = os.path.join(ANALYSIS_DIR, f"feature_importance{tag}.tsv")
    imp_df.to_csv(imp_path, sep="\t", index=False)
    ok(f"已写出：{imp_path}")

    say(f"\n── Top 20 重要基因 ──")
    for _, r in imp_df.head(20).iterrows():
        s = f"S{int(r['stratum'])}" if not np.isnan(r["stratum"]) else "无年龄"
        say(f"  {r['importance']:.5f}  {r['gene']:<18} {s}")

    # ---------- 9. 富集检验 ----------
    bg = stratum[~np.isnan(stratum)]
    top_s = stratum[top_idx]
    top_s = top_s[~np.isnan(top_s)]

    M, N = len(bg), len(top_s)
    if M == 0 or N == 0:
        die("背景或 Top 集为空，无法做富集检验")

    say("\n" + "=" * 60)
    say("各年龄层在 Top 基因中的富集（超几何检验，双单尾 + BH-FDR）")
    say("=" * 60)
    say(f"  背景：{M} 个基因   Top 集合：{N} 个基因")
    say(f"  {'层':>4}  {'背景':>6}  {'Top中':>6}  {'期望':>6}  {'倍数':>7}  "
        f"{'p富集':>10}  {'q富集':>9}  {'p耗竭':>10}  {'q耗竭':>9}")
    say("  " + "-" * 80)

    # 超几何检验用「双单尾」而不是双尾：
    #   p_up   = P(X ≥ k)  检验该层是否被富集（比随机更多）
    #   p_down = P(X ≤ k)  检验该层是否被耗竭（比随机更少）
    # 只用右尾 sf 时，fold<1 的层会得到 p≈1，那只说明"没有证据表明富集"，
    # 不能说成"显著回避"。方向判定必须各自配对检验。
    res_rows = []
    for s in np.sort(np.unique(bg)):
        n_layer = int((bg == s).sum())
        k_hit = int((top_s == s).sum())
        exp = N * n_layer / M
        fold = k_hit / exp if exp > 0 else np.nan
        p_up = hypergeom.sf(k_hit - 1, M, n_layer, N)
        p_down = hypergeom.cdf(k_hit, M, n_layer, N)
        res_rows.append((int(s), n_layer, k_hit, exp, fold, p_up, p_down))

    res = pd.DataFrame(res_rows, columns=["stratum", "n_background",
                                          "n_top", "expected", "fold",
                                          "p_value", "p_depletion"])

    # ---------- BH-FDR 多重检验校正 ----------
    # ⚠️ 16 层 × 2 方向 = 32 次检验（跨三次运行共 96 次），
    #   不做校正的话 p<0.05 会大量假阳性。
    #   做法：富集方向单独做一次 BH，耗竭方向单独做一次 BH
    #   （两个方向是事先设定的不同假设，不该混在一起校正）。
    res["q_value"] = bh_fdr(res["p_value"].to_numpy())
    res["q_depletion"] = bh_fdr(res["p_depletion"].to_numpy())

    for _, r in res.iterrows():
        flag = ""
        if r["q_value"] < 0.05 and r["fold"] > 1:
            flag = " ★富集(BH)"
        elif r["q_depletion"] < 0.05 and r["fold"] < 1:
            flag = " ▼耗竭(BH)"
        say(f"  S{int(r['stratum']):>3}  {int(r['n_background']):>6}  "
            f"{int(r['n_top']):>6}  {r['expected']:>6.1f}  "
            f"{r['fold']:>7.2f}  {r['p_value']:>10.2e}  "
            f"{r['q_value']:>9.2e}  {r['p_depletion']:>10.2e}  "
            f"{r['q_depletion']:>9.2e}{flag}")

    res_path = os.path.join(ANALYSIS_DIR, f"age_enrichment{tag}.tsv")
    res.to_csv(res_path, sep="\t", index=False)
    ok(f"已写出：{res_path}")

    say("\n  解读（★ 用 BH 校正后的 q<0.05 判定，这是本项目的正式口径）：")
    sig_up = res[(res["q_value"] < 0.05) & (res["fold"] > 1)]
    sig_dn = res[(res["q_depletion"] < 0.05) & (res["fold"] < 1)]
    raw_up = res[(res["p_value"] < 0.05) & (res["fold"] > 1)]
    raw_dn = res[(res["p_depletion"] < 0.05) & (res["fold"] < 1)]
    if len(sig_up) == 0 and len(sig_dn) == 0:
        say("    没有任何年龄层经 FDR 校正后显著富集或显著耗竭。")
        say("    → 区分细胞类型并不依赖特定年代的基因。这是有效结论，如实写。")
    else:
        for _, r in sig_up.iterrows():
            say(f"    S{int(r['stratum'])}：显著富集，倍数 {r['fold']:.2f}"
                f"（p={r['p_value']:.1e}, q={r['q_value']:.1e}）")
        for _, r in sig_dn.iterrows():
            say(f"    S{int(r['stratum'])}：显著耗竭，倍数 {r['fold']:.2f}"
                f"（p={r['p_depletion']:.1e}, q={r['q_depletion']:.1e}）")
    if len(raw_up) != len(sig_up) or len(raw_dn) != len(sig_dn):
        say(f"    （未校正口径下富集 {len(raw_up)} 层 / 耗竭 {len(raw_dn)} 层；"
            f"BH 校正后 {len(sig_up)} / {len(sig_dn)} 层。"
            f"被校正剔除的富集层："
            f"{sorted(set(raw_up.stratum) - set(sig_up.stratum))}）")

    say("\n" + "=" * 60)
    say("稳健性检验（应对 Moyers & Zhang 的方法学批评）")
    say("=" * 60)
    if args.min_share == 0 and args.min_length == 0:
        say("  稳健性检验 B（共识度过滤）：")
        say("      python 05_监督分析.py --min-share 0.8 --tag confident")
        say("  稳健性检验 C（剔除短蛋白，Moyers & Zhang 建议）：")
        say("      python 93_获取蛋白长度.py")
        say("      python 05_监督分析.py --min-length 100 --tag shortprot")
        say("  若三次运行的富集方向一致，说明结论既不依赖低质量年龄调用，")
        say("  也不依赖短序列基因——这正是对 Moyers & Zhang 批评的直接回应。")
    elif args.min_share > 0:
        say(f"  本次已启用共识度过滤（mode_share ≥ {args.min_share}）。")
        say("  与不加过滤的结果对比：feature_importance*.tsv 里的 stratum 分布是否同向。")
    else:
        say(f"  本次已启用短蛋白剔除（长度 ≥ {args.min_length} aa）。")
        say("  与全量结果对比：feature_importance*.tsv 里的 stratum 分布是否同向。")
    say("=" * 60)

    rep = os.path.join(ANALYSIS_DIR, f"classification_report{tag}.txt")
    with open(rep, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    ok(f"已写出：{rep}")

    # ---------- 分类指标单独落盘（供终审脚本逐项核对）----------
    met = pd.DataFrame([{
        "tag": args.tag,
        "n_cells": int(len(rows)),
        "n_types_all": int(len(all_labels)),
        "n_types_train": int(len(train_labels)),
        "n_types_test": int(len(te_labels)),
        "n_unseen_classes": int(len(unseen)),
        "n_unseen_cells": int(n_unseen_cells),
        "pct_unseen_cells": round(n_unseen_cells / len(y[is_test]), 4),
        "accuracy": round(float(acc), 4),
        "accuracy_ceiling": round(float(acc_ceiling), 4),
        "pct_of_ceiling": round(float(acc / acc_ceiling), 4),
        "macro_f1_test": round(float(f1m_te), 4),
        "macro_f1_all": round(float(f1m_all), 4),
        "rand_baseline_test": round(float(rand_te), 4),
        "majority_baseline": round(float(maj), 4),
        "x_majority": round(float(acc / maj), 2),
        "n_features": int(len(hvg)),
        "n_donors_train": int(len(donors_u) - n_test),
        "n_donors_test": int(n_test),
        "n_cells_train": int((~is_test).sum()),
        "n_cells_test": int(is_test.sum()),
    }])
    mp = os.path.join(ANALYSIS_DIR, f"classification_metrics{tag}.tsv")
    met.to_csv(mp, sep="\t", index=False)
    ok(f"已写出：{mp}")


if __name__ == "__main__":
    main()
