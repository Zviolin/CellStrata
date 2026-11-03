"""
06_稳健性对比.py — 三重稳健性检验的统一对照
===============================================
用途：
    应对 Moyers & Zhang (2015/2017/2018) 对 phylostratigraphy 的方法学批评。
    他们的核心指控是：**误差非随机**——序列短、快速演化的基因更容易找不到远缘同源，
    从而被判为"年轻"。他们明确建议：剔除易错基因后重跑，
    看已发表的结论是否消失甚至反转。

    本项目做三重检验：
      A. 全量          —— 3,000 个 HVG，不做任何过滤（基线）
      B. 共识度过滤    —— 只保留 Trigos/Liebeskind 年龄表中「多算法共识度高」
                          （mode_share ≥ 0.8，即 ≥80% 的同源推断算法给出同一层号）
                          的基因（905 个）
      C. 短蛋白剔除    —— 只用蛋白长度 ≥ 100 aa 的基因（2,405 个）

    ⚠ B 组用的是 **mode_share（共识度）**，不是 `NumDBsContributing`。
      后者是**参与投票的算法数（分母）**，与共识度**负相关**
      （实测 corr = −0.128）：参与算法多的基因反而常常更分散。
      按 `n_dbs ≥ 10` 过滤会筛掉高共识的基因，方向与"高置信"完全相反。

    另加一项**前提检验**（Moyers & Zhang 指控的直接检验）：
      短蛋白在各年龄层的分布是否真的偏斜？
      若短蛋白确实集中在年轻层，说明剔除它们是必要的修正；
      若均匀分布，说明这个担忧在本数据上不成立。

用法：
    python 06_稳健性对比.py

输出：
    robustness_comparison.tsv      三次运行的逐层方向 / 倍数 / p 值对照
    short_protein_by_stratum.tsv    各年龄层的长度分布（前提检验）
"""

import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import ANALYSIS_DIR, DATA_DIR, die, info, ok, warn  # noqa: E402

# 每个运行 = (显示名, 富集结果文件, 分类报告文件)
# ⚠ 富集文件名与分类报告文件的命名规则不同（age_enrichment_confident.tsv
#   对应的是 classification_report_confident.txt，不是 classification_report.txt），
#   所以必须显式配对，不能用字符串 replace 去猜。
RUNS = [
    ("A 全量", "age_enrichment.tsv", "classification_report.txt"),
    ("B 共识度过滤", "age_enrichment_confident.tsv", "classification_report_confident.txt"),
    ("C 短蛋白剔除", "age_enrichment_shortprot.tsv", "classification_report_shortprot.txt"),
]

# 显著性阈值：BH 校正后的 q < 0.05（本项目的正式口径）
ALPHA = 0.05


def direction(fold, q_up, q_down):
    """方向判定必须各自配对检验：
       右尾 q 只支持"富集"结论，左尾 q 才支持"耗竭"结论。
       两者都不显著 → 不显著（不能说成回避）。

       ⚠️ 用 **q（BH 校正后）** 而不是原始 p：每次运行有 16 层 × 2 方向
       = 32 次检验，跨三次运行共 96 次。不做多重检验校正的话，
       p<0.05 会产生大量假阳性（实测 A 组原始口径 8 层富集，
       BH 校正后只剩 6 层）。"""
    if q_up < ALPHA and fold > 1:
        return "富集"
    if q_down < ALPHA and fold < 1:
        return "耗竭"
    return "不显著"


def parse_metrics(path):
    """从分类报告文本里抽取关键指标。

    ⚠ 正则必须**行首锚定**并写全限定语：
      - 「准确率：」也有「准确率天花板：」「准确率是多数类基线的 …」等兄弟行，
        不锚定会抓到错的数；
      - 「宏平均 F1」现在有两个口径（测试集实际类别 / 全部 177 类口径），
        必须用带括号限定语的那一行，否则会误取。
    """
    m = {}
    if not os.path.exists(path):
        return m
    pats = [
        ("accuracy", r"^\s*准确率：([\d.]+)\s*$"),
        ("macro_f1", r"^\s*宏平均 F1（测试集实际出现的 \d+ 类）：([\d.]+)\s*$"),
        ("macro_f1_all", r"^\s*宏平均 F1（全部 \d+ 类口径）：([\d.]+)\s*$"),
        ("n_features", r"^\s*特征数：(\d+)\s*$"),
        ("n_types", r"^\s*类别数：(\d+)"),
        ("accuracy_ceiling", r"准确率天花板（剔除这些不可预测细胞后）：([\d.]+)"),
        ("majority_baseline", r"多数类基线（恒猜最高频类）：([\d.]+)"),
        ("x_majority", r"准确率是多数类基线的 ([\d.]+) 倍"),
        ("pct_of_ceiling", r"实际达到天花板的 ([\d.]+)%"),
    ]
    for line in open(path, encoding="utf-8"):
        line = line.rstrip()
        for key, pat in pats:
            if key in m:
                continue
            mm = re.search(pat, line)
            if mm:
                m[key] = int(mm.group(1)) if key in ("n_features", "n_types") \
                    else float(mm.group(1))
    return m


def main():
    # ---------- 载入三次运行 ----------
    tables = {}
    for name, enr_f, _ in RUNS:
        p = os.path.join(ANALYSIS_DIR, enr_f)
        if not os.path.exists(p):
            die(f"缺少 {name} 的富集结果：{p}\n请先运行对应的 05_监督分析.py 命令")
        tables[name] = pd.read_csv(p, sep="\t").set_index("stratum")

    print("=" * 76)
    print("① 三次运行的分类性能")
    print("=" * 76)
    print(f"  {'运行':<16}{'准确率':>10}{'宏F1':>10}{'特征数':>10}"
          f"{'达天花板':>10}")
    for name, _, rep_f in RUNS:
        m = parse_metrics(os.path.join(ANALYSIS_DIR, rep_f))
        if not m:
            die(f"缺少 {name} 的分类报告：{rep_f}\n"
                f"请先运行对应的 05_监督分析.py 命令")
        if "accuracy" not in m or "macro_f1" not in m:
            die(f"{name} 的分类报告格式对不上，解析不出准确率/宏F1。"
                f"\n文件：{rep_f}\n解析到的键：{sorted(m)}")
        print(f"  {name:<16}{m['accuracy']:>10.3f}{m['macro_f1']:>10.3f}"
              f"{m.get('n_features', 0):>10}{m.get('pct_of_ceiling', 0):>9.1f}%")

    # ---------- 逐层对照 ----------
    print()
    print("=" * 76)
    print("② 各年龄层方向对照（fold 倍数 + 方向，显著性按 BH 校正后的 q<0.05）")
    print("=" * 76)
    hdr = f"  {'层':<5}" + "".join(f"{name:>22}" for name, _, _ in RUNS)
    print(hdr)

    rows = []
    all_strata = sorted(set.intersection(*[set(t.index) for t in tables.values()]))
    for s in all_strata:
        line = f"  S{s:<4}"
        rec = {"stratum": int(s)}
        dirs = []
        for name, _, _ in RUNS:
            t = tables[name]
            f = float(t.loc[s, "fold"])
            pu = float(t.loc[s, "p_value"])
            pd_ = float(t.loc[s, "p_depletion"]) if "p_depletion" in t.columns else 1.0
            qu = float(t.loc[s, "q_value"]) if "q_value" in t.columns else pu
            qd = (float(t.loc[s, "q_depletion"])
                  if "q_depletion" in t.columns else pd_)
            d = direction(f, qu, qd)
            dirs.append(d)
            line += f"{f:>7.2f} {d:<7}"
            rec[f"fold_{name[0]}"] = round(f, 4)
            rec[f"p_enrich_{name[0]}"] = pu
            rec[f"p_deplete_{name[0]}"] = pd_
            rec[f"q_enrich_{name[0]}"] = qu
            rec[f"q_deplete_{name[0]}"] = qd
            rec[f"dir_{name[0]}"] = d
        conflict = ("富集" in dirs and "耗竭" in dirs)
        rec["conflict"] = "YES" if conflict else "no"
        tag = "✗冲突" if conflict else ("✓一致" if len(set(dirs)) == 1 else "—部分")
        print(line + f"   {tag}")
        rows.append(rec)

    out = pd.DataFrame(rows)
    out_path = os.path.join(ANALYSIS_DIR, "robustness_comparison.tsv")
    out.to_csv(out_path, sep="\t", index=False)
    ok(f"已写出：{out_path}")

    # ---------- 显著层汇总 ----------
    print()
    print("=" * 76)
    print("③ 显著层清单汇总（BH 校正后 q<0.05 为正式口径）")
    print("=" * 76)
    sig_up, sig_dn = {}, {}
    for name, _, _ in RUNS:
        t = tables[name]
        qc = "q_value" if "q_value" in t.columns else "p_value"
        qdc = "q_depletion" if "q_depletion" in t.columns else "p_depletion"
        up = sorted(t[(t[qc] < ALPHA) & (t["fold"] > 1)].index.tolist())
        dn = sorted(t[(t[qdc] < ALPHA) & (t["fold"] < 1)].index.tolist())
        # 同时记下未校正口径，供报告说明差异
        raw_up = sorted(t[(t["p_value"] < ALPHA) & (t["fold"] > 1)].index.tolist())
        raw_dn = sorted(t[(t["p_depletion"] < ALPHA) & (t["fold"] < 1)].index.tolist())
        sig_up[name], sig_dn[name] = up, dn
        print(f"  {name:<16} 富集 {up}")
        print(f"  {'':<16} 耗竭 {dn}")
        if raw_up != up or raw_dn != dn:
            print(f"  {'':<16}（未校正口径：富集 {raw_up} / 耗竭 {raw_dn}）")

    core = set.intersection(*[set(v) for v in sig_up.values()])
    if core:
        union = set.union(*[set(v) for v in sig_up.values()])
        print(f"\n  → 三次运行共同显著富集的层：{sorted(core)}")
        print(f"    Jaccard 相似度 {len(core) / len(union):.3f}")
    dn_common = set.intersection(*[set(v) for v in sig_dn.values()])
    print(f"  → 三次运行共同显著耗竭的层：{sorted(dn_common)}")
    n_conf = int((out["conflict"] == "YES").sum())
    print(f"  → 方向冲突的层数：{n_conf}")

    # ---------- 前提检验：短蛋白是否真的偏年轻 ----------
    print()
    print("=" * 76)
    print("④ 前提检验：短蛋白在各年龄层的分布（Moyers & Zhang 指控的直接检验）")
    print("=" * 76)
    plen_path = os.path.join(DATA_DIR, "protein_length.tsv")
    if not os.path.exists(plen_path):
        warn(f"缺少 {plen_path}，跳过前提检验")
        return
    pl = pd.read_csv(plen_path, sep="\t")
    age = pd.read_csv(os.path.join(DATA_DIR, "gene_age.tsv"), sep="\t")
    m = pl.merge(age[["symbol", "stratum"]], on="symbol", how="inner")
    m["short"] = m["length"] < 100

    rows2 = []
    for s in sorted(m["stratum"].unique()):
        sub = m[m["stratum"] == s]
        rows2.append({
            "stratum": int(s),
            "n_genes": len(sub),
            "median_length": int(sub["length"].median()),
            "mean_length": round(float(sub["length"].mean()), 1),
            "n_short_lt100": int(sub["short"].sum()),
            "pct_short": round(float(sub["short"].mean()) * 100, 2),
        })
    tab = pd.DataFrame(rows2)
    out2 = os.path.join(ANALYSIS_DIR, "short_protein_by_stratum.tsv")
    tab.to_csv(out2, sep="\t", index=False)
    ok(f"已写出：{out2}")

    print(f"  {'层':<6}{'基因数':>8}{'中位长度':>10}{'均值':>10}{'<100aa':>10}{'占比%':>9}")
    for _, r in tab.iterrows():
        print(f"  S{int(r['stratum']):<5}{int(r['n_genes']):>8}"
              f"{int(r['median_length']):>10}{r['mean_length']:>10.1f}"
              f"{int(r['n_short_lt100']):>10}{r['pct_short']:>9.2f}")

    from scipy.stats import spearmanr
    rho, p = spearmanr(tab["stratum"], tab["pct_short"])
    print(f"\n  层号 vs 短蛋白占比：Spearman ρ = {rho:+.3f}（p = {p:.3f}）")
    print(f"  ⚠ 该相关只用 16 个点（16 个年龄层），且各层基因数差异巨大"
          f"（S2 有 {int(tab[tab['stratum']==2]['n_genes'].iloc[0]):,} 个，"
          f"S16 只有 {int(tab[tab['stratum']==16]['n_genes'].iloc[0])} 个），"
          f"秩相关隐含偏向大层，p 值不应被单独引用。")
    young = tab[tab["stratum"] >= 13]["pct_short"].mean()
    old = tab[tab["stratum"] <= 4]["pct_short"].mean()
    print(f"  老层 S1–S4 平均短蛋白占比 {old:.2f}%，"
          f"年轻层 S13–S16 平均 {young:.2f}%（相差 {young/old:.1f} 倍）")
    if rho > 0 and p < 0.05:
        print("  → 短蛋白确实更集中在年轻层：**年龄表本身**带有系统性「年轻偏倚」，")
        print("    Moyers & Zhang 描述的机制在本数据上成立。")
        print("  ⚠ 但方向要说清楚：偏倚是**已经发生**在 Trigos 表里的，")
        print("    检验 C 是在这张**已经偏斜的表**上把短蛋白从背景里剔掉，")
        print("    它**不能修复**既有偏倚，只能让背景集合更干净、")
        print("    同时**降低**受影响层（尤其 S13/S15）的统计功效。")
        print("    所以检验 C 的正确定位是**保守性压力测试**"
              "（结论不依赖短蛋白），")
        print("    不是「已修正偏差」——它的实际效果是削弱而非增强年轻层信号。")
    else:
        print("  → 短蛋白在各层分布无显著偏斜：该担忧在本数据的年龄表上不成立，")
        print("    C 组结果与 A 组一致也就不足为奇。")

    print()
    print("=" * 76)
    print("结论")
    print("=" * 76)
    print(f"  · 三次运行的方向冲突层数：{n_conf}")
    if n_conf == 0:
        print("  → 核心结论不依赖于低置信度调用，也不依赖于短蛋白基因。")
    if core:
        print(f"  · 三次共同显著富集：{sorted(core)}")
    print(f"  · 三次共同显著耗竭：{sorted(dn_common)}")


if __name__ == "__main__":
    main()
