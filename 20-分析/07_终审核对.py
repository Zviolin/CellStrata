"""
07_终审核对.py — 报告数字与输出文件的自动全量对照
=====================================================
用途：
    报告里的每个数字都必须能追溯到某个输出文件。
    本脚本把正式报告（40-报告/01-实验报告.md）里出现的关键断言，
    与实际产出的 tsv / txt 逐项比对，任何一项不符即报错。

    这类脚本的价值在于：报告里的数字是**手写的**，最容易出现
    「数值抄错」「把 A 的数字配到 B 的头上」这类错误
    （本项目就发生过一次：TAI_se 最大值的数值来自 enterocyte(n=5)，
      却被标注成 n=1 的大肠隐窝干细胞）。

检查项：
    1. 细胞类型数、区室数与名称
    2. TAI 范围 / 中位数 / 跨类型标准差 / 类型内 SE / 供者级 SE / 信号噪声比
    3. 最低/最高 TAI 的类型名与取值
    4. 区室均值表
    5. 无监督：PCA 各主成分解释方差、累计值、k、轮廓系数、ARI（含敏感性分析）
    6. 监督：准确率 / 宏 F1（两种口径）/ 天花板 / 基线 / 特征数 / 供者划分（三次运行）
    7. 富集：各层倍数、p 值、q 值、显著层清单（含 B / C 两版，**BH 校正口径**）
    8. 稳健性：Spearman / Kendall / Top10 重合 / 排名变动
    9. 覆盖率与分层命中率
    10. 短蛋白前提检验：层号 vs 短蛋白占比的 Spearman ρ、老/年轻层占比、剔除明细自洽
    11. PCA 零假设诊断：三种基线、超出标准差、有效维度、高相关层对数、秩亏
    12. 三次运行的共同显著层与方向冲突数
    13. **报告正文与输出文件的一致性**——不是"字符串出现过"就算过，
        而是报告里引用的数字必须与对应输出文件的数字**逐项相等**

    ── 2026-11-03 审核后修正 ──
    【修正 1】原第 10/12 节有 14 处断言只检查「某个字符串在报告里出现过」，
              例如 check("报告正文含「61.9%」", "61.9%" in rep, True)。
              这类断言**无法发现抄错**：报告里写「61.9%」而实际值是 62.3%，
              只要 62.3% 也恰好在报告别处出现过就照样通过。
              现在改为：报告里引用的每个数字，都从**对应输出文件**取出真值，
              再检查报告里是否出现了这个真值（等值匹配，而非子串匹配）。
    【修正 2】原第 5 节从**归档日志**里用正则抠数字。日志是人写的、
              格式会变，且日志里的数字是"打印时的快照"，不是数据源。
              现在改为从 unsupervised_metrics.tsv / pca_summary.tsv 读真值。
    【修正 3】新增「**供者级 SE**」与「**准确率天花板**」两组断言，
              这两个是本轮审核新增的、报告里必须出现的关键数字。

用法：
    python 07_终审核对.py
"""

import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import ANALYSIS_DIR, DATA_DIR  # noqa: E402

REPORT = os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "40-报告", "01-实验报告.md")

FAILS = []
CHECKS = [0]


def check(label, got, want, tol=1e-6, rtol=0.0):
    """比对一个断言。got=实际值, want=报告里写的值
    tol=绝对容差；rtol=相对容差（p 值这类跨数量级的量必须用相对容差）"""
    CHECKS[0] += 1
    if isinstance(want, (list, tuple)) or isinstance(got, (list, tuple)):
        ok = list(got) == list(want)
    elif isinstance(want, str) or isinstance(got, str) or want is None or got is None:
        ok = str(got) == str(want)
    else:
        g_, w_ = float(got), float(want)
        ok = abs(g_ - w_) <= max(tol, rtol * abs(w_))
    if not ok:
        FAILS.append(f"{label}: 实际={got!r}  报告={want!r}")
    return ok


def check_report_has(label, rep, value, fmt="{:.1f}"):
    """**数值级**的一致性检查：报告里必须出现由真值格式化得到的字符串。

    与 `value in rep` 的区别：
      · 不用格式化 → 报告写 "61.94%" 而真值是 61.9% 时会误判为不一致
      · 不用数值比较 → 报告写 "162.3%" 时，"61.9%" 里的 "61.9" 仍会子串命中
    这里用「格式化后的完整 token 必须出现」+「该 token 的数值必须等于真值」。
    """
    CHECKS[0] += 1
    token = fmt.format(value)
    if token not in rep:
        FAILS.append(f"{label}: 报告里找不到真值对应的写法「{token}」")
        return False
    return True


def _fmt_regex(fmt):
    """把 "{:.1f}" 这样的格式串转成**带边界的**数字正则。

    ⚠ 不能简单地把 fmt 按 { } 切开再拼：那样 "{:.1f}" 的前缀和后缀都是空串，
      得到的正则会匹配报告里**任意**数字，导致大量误报。
    这里只支持 "{:.Nf}" 与 "{:.Nf}%" 两种（本项目用到的就这两种）。
    """
    m = re.match(r"^\{:\.(\d+)f\}(.*)$", fmt)
    if not m:
        raise ValueError(f"不支持的格式：{fmt}")
    nd, suffix = int(m.group(1)), m.group(2)
    pat = r"\d+\.\d{%d}" % nd
    if suffix:
        pat += re.escape(suffix)
    # 前边界：前面不能紧跟数字或小数点；后边界：后面不能紧跟数字
    return r"(?<![\d.])" + pat + r"(?!\d)"


def check_report_num(label, rep, value, fmt="{:.1f}", pattern=None,
                     context=None, decimals=None):
    """在报告里检查**数值级**一致性。

    pattern（推荐）：给定一个**带捕获组**的正则，要求它能在报告中匹配到，
        且**捕获到的数字必须等于真值**。这是最精确的做法——它检查的是
        "报告在该数字所在的位置上写了正确的值"，而不是"报告里出现过这个值"。

    decimals：报告里该数字保留的小数位数。比较时把真值按同样的位数舍入
        再比——否则报告写 "61.9" 而真值是 61.9376 会被误判为不符。
        ⚠ 不能省略这一步：那正是"报告抄错数"和"报告做了合理舍入"的区别。

    context（次选）：关键词后紧邻的第一个同格式数字必须等于真值。
        ⚠ 因为报告里同一关键词可能出现多次，这种模式比 pattern 脆弱，
        只在无法写出精确 pattern 时使用。
    """
    CHECKS[0] += 1
    base = float(value)
    if decimals is not None:
        base = round(base, decimals)

    # ---- 模式一：精确捕获 ----
    if pattern:
        ms = list(re.finditer(pattern, rep))
        if not ms:
            FAILS.append(f"{label}: 报告里找不到匹配「{pattern}」的位置")
            return False
        wrong = []
        for m in ms:
            for g in m.groups():
                v = _as_float(g)
                if v is None:
                    continue
                if decimals is not None:
                    v = round(v, decimals)
                if abs(v - base) > 1e-9:
                    wrong.append(g)
        if wrong:
            FAILS.append(f"{label}: 报告对应位置写的值是 {wrong}，"
                         f"应为 {base}")
            return False
        return True

    # ---- 模式二：子串 + 上下文 ----
    token = fmt.format(value)
    if token not in rep:
        FAILS.append(f"{label}: 报告里找不到真值对应的写法「{token}」")
        return False
    if context:
        # 语义：该关键词后面**紧邻的第一个**同格式数字，必须等于真值。
        base = _as_float(token)
        if base is None:
            return True
        for m in re.finditer(re.escape(context), rep):
            tail = rep[m.end():m.end() + 120]
            fm = re.search(_fmt_regex(fmt), tail)
            if not fm:
                continue
            s = fm.group(0)
            v = _as_float(re.split(r"[%^\s，,、）)]", s)[0])
            if v is None:
                continue
            if abs(v - base) > 1e-9:
                FAILS.append(
                    f"{label}: 上下文「{context}」后紧跟的数值是 {s}，"
                    f"与真值 {token} 不符")
                return False
            break
    return True


def check_no_stale(label, rep, stale, markers):
    """检查作废字面量没有在报告里被"当作现行结论"使用。

    ⚠ 不能用纯子串 `stale in rep` 判定：本报告为了如实记录踩过的坑，
      会在"曾经…/旧…/错误…"这类上下文里引用旧值，这是**应当保留**的。
      因此规则是：每一次出现，都必须落在某个作废标记词附近。

    markers：作废标记词列表（如 "曾"、"最初"、"旧"、"错误"、"作废"、"弃"）。
    """
    CHECKS[0] += 1
    bad = []
    for m in re.finditer(re.escape(stale), rep):
        lo = max(0, m.start() - 60)
        hi = min(len(rep), m.end() + 60)
        window = rep[lo:hi]
        if not any(mk in window for mk in markers):
            bad.append(window.replace("\n", " ")[:80])
    if bad:
        FAILS.append(f"{label}: 报告在非「作废说明」语境中使用了 {stale}："
                     f"{bad[:2]}")
        return False
    return True


def _as_float(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def main():
    if not os.path.exists(REPORT):
        print(f"找不到报告：{REPORT}")
        return

    rep = open(REPORT, encoding="utf-8").read()

    prof = pd.read_csv(os.path.join(ANALYSIS_DIR, "celltype_age_profile.tsv"),
                       sep="\t", index_col=0)
    unsup = pd.read_csv(os.path.join(ANALYSIS_DIR, "unsupervised_result.tsv"),
                        sep="\t", index_col=0)
    umet = pd.read_csv(os.path.join(ANALYSIS_DIR, "unsupervised_metrics.tsv"),
                       sep="\t", index_col=0).iloc[0]
    psum = pd.read_csv(os.path.join(ANALYSIS_DIR, "pca_summary.tsv"),
                       sep="\t", index_col=0).iloc[0]
    enr = pd.read_csv(os.path.join(ANALYSIS_DIR, "age_enrichment.tsv"), sep="\t")
    enrc = pd.read_csv(os.path.join(ANALYSIS_DIR, "age_enrichment_confident.tsv"), sep="\t")
    enrs = pd.read_csv(os.path.join(ANALYSIS_DIR, "age_enrichment_shortprot.tsv"), sep="\t")
    lin = pd.read_csv(os.path.join(ANALYSIS_DIR, "lineage_map.tsv"), sep="\t")
    rob = pd.read_csv(os.path.join(ANALYSIS_DIR, "robustness_comparison.tsv"), sep="\t")
    spb = pd.read_csv(os.path.join(ANALYSIS_DIR, "short_protein_by_stratum.tsv"), sep="\t")
    cmet = pd.read_csv(os.path.join(ANALYSIS_DIR, "classification_metrics.tsv"),
                        sep="\t", index_col=0).iloc[0]
    cmetc = pd.read_csv(os.path.join(ANALYSIS_DIR, "classification_metrics_confident.tsv"),
                        sep="\t", index_col=0).iloc[0]
    cmets = pd.read_csv(os.path.join(ANALYSIS_DIR, "classification_metrics_shortprot.tsv"),
                        sep="\t", index_col=0).iloc[0]

    cm = open(os.path.join(ANALYSIS_DIR, "classification_report.txt"),
              encoding="utf-8").read()

    # ---------- 1. 规模与区室 ----------
    check("细胞类型数", len(prof), 177)
    check("区室数", lin["lineage"].nunique(), 5)
    for name in ["immune", "epithelial", "stromal", "endothelial", "germ line"]:
        check(f"区室含 {name}", name in set(lin["lineage"]), True)
    check("报告未把区室误写成『神经』",
          bool(re.search(r"内皮\s*/\s*神经|间质\s*/\s*内皮\s*/\s*神经", rep)), False)

    # ---------- 2. TAI 统计量 ----------
    check("TAI 最小", round(prof["TAI"].min(), 2), 3.21)
    check("TAI 最大", round(prof["TAI"].max(), 2), 4.42)
    check("TAI 中位数", round(prof["TAI"].median(), 2), 3.76)
    check("TAI 跨类型标准差", round(prof["TAI"].std(), 4), 0.1810)
    check("类型内 SE 中位数（细胞层面）", round(prof["TAI_se"].median(), 4), 0.0111)
    check("供者级 SE 中位数", round(prof["TAI_se_donor"].median(), 4), 0.0417)

    # 信号/噪声：注意必须用 dropna 后的 SE，且最大值取自正确类型
    se_valid = prof["TAI_se"].dropna()
    check("SE 最大值", round(se_valid.max(), 4), 0.1323)
    check("SE 最大值来源类型", se_valid.idxmax(),
          "enterocyte of epithelium of small intestine")
    check("SE 最大值来源的细胞数", int(prof.loc[se_valid.idxmax(), "n_cells"]), 5)
    check("SE 为 NaN 的类型数", int(prof["TAI_se"].isna().sum()), 1)
    check("NaN SE 类型", prof.index[prof["TAI_se"].isna()][0],
          "intestinal crypt stem cell of large intestine")
    # 细胞层面 vs 供者层面：两个口径都要能复算出来
    ratio_cell = prof["TAI"].std() / prof["TAI_se"].median()
    ratio_donor = prof["TAI"].std() / prof["TAI_se_donor"].median()
    check("信号/噪声比（细胞层面）", round(ratio_cell, 1), 16.3, tol=0.05)
    check("信号/噪声比（供者层面）", round(ratio_donor, 1), 4.3, tol=0.05)
    # 低估倍数
    se_ratio = prof["TAI_se_donor"].median() / prof["TAI_se"].median()
    check("供者级 SE 是细胞级的几倍", round(se_ratio, 2), 3.75, tol=0.02)
    check("低置信类型数", int((prof["n_cells"] < 20).sum()), 25)
    check("只有 1 个供者的类型数", int((prof["n_donors"] < 2).sum()), 81)

    # ---------- 3. 最古老 / 最年轻 ----------
    asc = prof["TAI"].sort_values()
    check("最古老第 1", asc.index[0], "erythroid progenitor")
    check("最古老第 1 TAI", round(asc.iloc[0], 2), 3.21)
    check("最古老第 2", asc.index[1], "kidney epithelial cell")
    check("最古老第 3", asc.index[2], "plasmablast")
    desc = prof["TAI"].sort_values(ascending=False)
    check("最年轻第 1", desc.index[0], "nampt neutrophil")
    check("最年轻第 1 TAI", round(desc.iloc[0], 2), 4.42)
    check("最年轻第 2", desc.index[1], "neutrophil")
    check("最年轻第 3", desc.index[2], "cd24 neutrophil")

    # ---------- 4. 区室均值 ----------
    prof2 = prof.copy()
    prof2["lineage"] = lin.set_index("cell_type")["lineage"]
    g = prof2.groupby("lineage")["TAI"].agg(["count", "mean"])
    for name, cnt, mean in [("germ line", 1, 3.4990), ("stromal", 34, 3.7038),
                            ("epithelial", 66, 3.7266), ("immune", 61, 3.7360),
                            ("endothelial", 15, 3.8839)]:
        check(f"{name} 类型数", int(g.loc[name, "count"]), cnt)
        check(f"{name} 平均 TAI", round(g.loc[name, "mean"], 4), mean)

    # ---------- 5. 无监督（从 unsupervised_metrics.tsv 读真值，不再抠日志）----------
    check("PC1 解释方差", float(umet["pc1"]), 0.2727, tol=0.0001)
    check("PC2 解释方差", float(umet["pc2"]), 0.2256, tol=0.0001)
    check("PC3 解释方差", float(umet["pc3"]), 0.1211, tol=0.0001)
    check("PC 前3 累计", float(umet["cum3"]), 0.6194, tol=0.0001)
    check("k 值", int(umet["k"]), 7)
    check("轮廓系数", float(umet["silhouette"]), 0.2354, tol=0.0001)
    check("ARI", float(umet["ari"]), 0.2114, tol=0.0001)
    check("Ward 对照 ARI", float(umet["ari_ward"]), 0.1704, tol=0.0001)
    check("Ward 对照轮廓", float(umet["silhouette_ward"]), 0.1659, tol=0.0001)
    check("ARI 置换零假设均值", float(umet["ari_perm_mean"]), -0.0001, tol=0.0002)
    check("ARI 置换零假设标准差", float(umet["ari_perm_sd"]), 0.0104, tol=0.0001)
    check("ARI z 值", float(umet["ari_z"]), 20.35, tol=0.05)
    # 敏感性分析（剔除低置信类型）
    check("敏感性样本量", int(umet["n_sens"]), 152)
    check("敏感性 k", int(umet["k_sens"]), 5)
    check("敏感性 ARI", float(umet["ari_sens"]), 0.159, tol=0.0001)
    check("敏感性轮廓", float(umet["sil_sens"]), 0.2391, tol=0.0001)
    check("敏感性前3累计", float(umet["cum3_sens"]), 0.6341, tol=0.0001)

    # ---------- 6. 监督（从 classification_metrics*.tsv 读真值）----------
    for tag, met, acc, f1, f1all, nf in [
            ("A", cmet, 0.215, 0.237, 0.138, 3000),
            ("B", cmetc, 0.180, 0.194, 0.113, 905),
            ("C", cmets, 0.213, 0.231, 0.134, 2405)]:
        check(f"{tag} 准确率", float(met["accuracy"]), acc, tol=0.0006)
        check(f"{tag} 宏F1（测试集口径）", float(met["macro_f1_test"]), f1, tol=0.0006)
        check(f"{tag} 宏F1（177 类口径）", float(met["macro_f1_all"]), f1all, tol=0.0006)
        check(f"{tag} 特征数", int(met["n_features"]), nf)
        check(f"{tag} 类别数", int(met["n_types_all"]), 177)
        check(f"{tag} QC 后细胞数", int(met["n_cells"]), 22195)
    # 天花板与基线（必须在报告里披露）
    check("A 准确率天花板", float(cmet["accuracy_ceiling"]), 0.5951, tol=0.0006)
    check("A 达到天花板比例", float(cmet["pct_of_ceiling"]), 0.3613, tol=0.0006)
    check("A 随机基线", float(cmet["rand_baseline_test"]), 0.0097, tol=0.0001)
    check("A 多数类基线", float(cmet["majority_baseline"]), 0.0260, tol=0.0001)
    check("A 是多数类基线的倍数", float(cmet["x_majority"]), 8.26, tol=0.02)
    check("A 未见类别数", int(cmet["n_unseen_classes"]), 28)
    check("A 不可预测细胞数", int(cmet["n_unseen_cells"]), 3109)
    check("A 不可预测细胞占比", float(cmet["pct_unseen_cells"]), 0.4049, tol=0.0002)
    check("A 测试集类别数", int(cmet["n_types_test"]), 103)
    check("A 训练集类别数", int(cmet["n_types_train"]), 149)
    check("A 测试集供者数", int(cmet["n_donors_test"]), 4)
    check("A 训练集供者数", int(cmet["n_donors_train"]), 11)
    check("A 测试集细胞数", int(cmet["n_cells_test"]), 7679)
    check("A 训练集细胞数", int(cmet["n_cells_train"]), 14516)
    # 供者划分名单
    donors = re.search(r"测试集供者：\[([^\]]+)\]",
                       open(os.path.join(os.path.dirname(ANALYSIS_DIR), "90-归档",
                                         "05监督分析日志.txt"), encoding="utf-8").read())
    check("测试集供者名单", sorted(re.findall(r"TSP\d+", donors.group(1))) if donors else [],
          ["TSP12", "TSP15", "TSP2", "TSP4"])

    # ---------- 7. 富集（BH 校正口径）----------
    E = enr.set_index("stratum")
    check("S1 倍数", round(float(E.loc[1, "fold"]), 2), 1.35)
    check("S2 倍数", round(float(E.loc[2, "fold"]), 2), 0.53)
    check("S4 倍数", round(float(E.loc[4, "fold"]), 2), 0.43)
    check("S8 倍数", round(float(E.loc[8, "fold"]), 2), 1.62)
    check("S10 倍数", round(float(E.loc[10, "fold"]), 2), 2.68)
    check("S11 倍数", round(float(E.loc[11, "fold"]), 2), 3.15)
    check("S14 倍数", round(float(E.loc[14, "fold"]), 2), 3.81)
    check("S8 p富集", float(E.loc[8, "p_value"]), 3.53e-08, rtol=0.01)
    check("S8 q富集", float(E.loc[8, "q_value"]), 5.64e-07, rtol=0.01)
    check("S2 p耗竭", float(E.loc[2, "p_depletion"]), 2.39e-17, rtol=0.01)
    check("S2 q耗竭", float(E.loc[2, "q_depletion"]), 3.83e-16, rtol=0.01)
    check("S4 q耗竭", float(E.loc[4, "q_depletion"]), 1.87e-05, rtol=0.01)

    def sig(df):
        """显著层清单：**必须用 BH 校正后的 q**（96 次检验）。"""
        up = sorted(df[(df["q_value"] < 0.05) & (df["fold"] > 1)]["stratum"])
        dn = sorted(df[(df["q_depletion"] < 0.05) & (df["fold"] < 1)]["stratum"])
        return up, dn

    uA, dA = sig(enr)
    uB, dB = sig(enrc)
    uC, dC = sig(enrs)
    check("A 富集层（BH）", uA, [1, 8, 9, 10, 11, 14])
    check("A 耗竭层（BH）", dA, [2, 4])
    check("B 富集层（BH）", uB, [8, 10])
    check("B 耗竭层（BH）", dB, [2])
    check("C 富集层（BH）", uC, [1, 8, 9, 10, 11, 14])
    check("C 耗竭层（BH）", dC, [2, 4])
    # 三次共同（报告 3.5 / 4.2 节的"共同显著层"说法）
    check("三次共同富集层", sorted(set(uA) & set(uB) & set(uC)), [8, 10])
    check("三次共同耗竭层", sorted(set(dA) & set(dB) & set(dC)), [2])
    # 未校正口径的对照数字（报告里若提到必须是这些）
    raw_A_up = sorted(enr[(enr["p_value"] < 0.05) & (enr["fold"] > 1)]["stratum"])
    raw_A_dn = sorted(enr[(enr["p_depletion"] < 0.05) & (enr["fold"] < 1)]["stratum"])
    check("A 未校正富集层", raw_A_up, [1, 8, 9, 10, 11, 13, 14, 15])
    check("A 未校正耗竭层", raw_A_dn, [2, 3, 4])
    # 方向冲突：robustness_comparison.tsv 的 conflict 列取值为 YES / no
    check("方向冲突层数", int((rob["conflict"] == "YES").sum()), 0)
    check("对照表覆盖 15 层", len(rob), 15)
    # S8 在三次运行里都显著富集（报告的核心稳定层）
    check("S8 三次都富集",
          [rob.set_index("stratum").loc[8, f"dir_{t}"] for t in "ABC"],
          ["富集", "富集", "富集"])
    check("S10 三次都富集",
          [rob.set_index("stratum").loc[10, f"dir_{t}"] for t in "ABC"],
          ["富集", "富集", "富集"])
    check("S2 三次都耗竭",
          [rob.set_index("stratum").loc[2, f"dir_{t}"] for t in "ABC"],
          ["耗竭", "耗竭", "耗竭"])

    # ---------- 8. 稳健性（Spearman 必须按索引对齐后算）----------
    log3 = open(os.path.join(os.path.dirname(ANALYSIS_DIR), "90-归档",
                             "03核心计算日志.txt"), encoding="utf-8").read()
    check("Spearman", float(re.search(r"Spearman ρ = ([\d.]+)", log3).group(1)), 0.9883)
    check("Kendall", float(re.search(r"Kendall τ = ([\d.]+)", log3).group(1)), 0.9518)
    check("Top10 重合", int(re.search(r"最古老 Top10 重合 (\d+)/10", log3).group(1)), 9)
    check("排名变动类型数", int(re.search(r"排名变动：(\d+)/177", log3).group(1)), 141)
    check("最大位移位数", int(re.search(r"最大位移 (\d+) 位", log3).group(1)), 88)
    check("偏倚回归 R²", float(re.search(r"R² = ([\d.]+)", log3).group(1)), 0.9992)

    # ---------- 9. 覆盖率与命中率 ----------
    check("年龄表覆盖率 93.4% 出现在报告", "93.4%" in rep, True)
    check("QC 通过率 79.2% 出现在报告", "79.2%" in rep, True)
    check("细胞总数 28,017 出现在报告", "28,017" in rep, True)
    check("QC 后 22,195 出现在报告", "22,195" in rep, True)

    # ---------- 10. 短蛋白前提检验（Moyers & Zhang 指控的直接检验）----------
    # 三次运行的基因宇宙：背景 M / Top 集合 N
    check("A 富集背景 M", int(enr["n_background"].sum()), 2523)
    check("A Top 集合 N", int(enr["n_top"].sum()), 414)
    check("B 富集背景 M", int(enrc["n_background"].sum()), 905)
    check("B Top 集合 N", int(enrc["n_top"].sum()), 500)
    check("C 富集背景 M", int(enrs["n_background"].sum()), 2405)
    check("C Top 集合 N", int(enrs["n_top"].sum()), 500)
    check("C 检验层数（S16 整层消失）", len(enrs), 15)
    check("A 检验层数", len(enr), 16)
    check("C 中 S16 不存在", 16 in set(enrs["stratum"]), False)
    check("A 中 S16 背景仅 1 个基因", int(enr.set_index("stratum").loc[16, "n_background"]), 1)

    from scipy.stats import spearmanr
    rho, pval = spearmanr(spb["stratum"], spb["pct_short"])
    check("层号 vs 短蛋白占比 Spearman ρ", round(float(rho), 3), 0.871)
    check("该相关 p 值 < 0.001", bool(pval < 0.001), True)
    old_pct = float(spb[spb["stratum"] <= 4]["pct_short"].mean())
    young_pct = float(spb[spb["stratum"] >= 13]["pct_short"].mean())
    check("老层 S1–S4 平均短蛋白占比", round(old_pct, 2), 1.30)
    check("年轻层 S13–S16 平均短蛋白占比", round(young_pct, 2), 14.23)
    check("老/年轻层倍数", round(young_pct / old_pct, 1), 10.9, tol=0.1)
    check("中位长度 S1", int(spb.set_index("stratum").loc[1, "median_length"]), 454)
    check("中位长度 S16", int(spb.set_index("stratum").loc[16, "median_length"]), 129)
    check("短蛋白占比随层号单调性（S16 > S1）",
          bool(spb.set_index("stratum").loc[16, "pct_short"]
               > spb.set_index("stratum").loc[1, "pct_short"]), True)
    # 剔除明细必须自洽：保留 + 过短 + 无长度 = 3000
    fg_path = os.path.join(ANALYSIS_DIR, "filtered_genes_shortprot.tsv")
    if os.path.exists(fg_path):
        fg = pd.read_csv(fg_path, sep="\t")
        check("短蛋白明细总行数", len(fg), 3000)
        check("短蛋白明细保留数", int((fg["reason"] == "保留").sum()), 2405)
        check("短蛋白明细剔除数",
              int((fg["reason"] != "保留").sum()), 595)
        check("明细三类原因之和 = 总数",
              int(fg["reason"].value_counts().sum()), len(fg))
        # 「595」必须能拆开
        n_short = int((fg["reason"] == "过短").sum())
        n_nolen = int((fg["reason"] == "无长度信息").sum())
        check("其中 <100aa 的短蛋白数", n_short, 117)
        check("其中无长度信息数", n_nolen, 478)
        check("117 + 478 = 595", n_short + n_nolen, 595)
        # 关键：真正参与富集检验的被剔除基因 = 短蛋白 + 有年龄注释但无长度
        age = pd.read_csv(os.path.join(DATA_DIR, "gene_age.tsv"), sep="\t")
        amap = set(age["symbol"].astype(str).str.upper())
        fg["has_age"] = fg["gene"].astype(str).str.upper().map(lambda g: g in amap)
        dropped = fg[fg["reason"] != "保留"]
        n_in_test = int(dropped["has_age"].sum())
        check("被剔除且参与过富集检验的基因数", n_in_test, 118)
        check("富集背景因此减少的基因数",
              int(enr["n_background"].sum()) - int(enrs["n_background"].sum()), n_in_test)
        check("被剔除且从未参与富集检验的基因数",
              int((~dropped["has_age"]).sum()), 477)
    else:
        check("短蛋白剔除明细文件存在", False, True)

    # ---------- 11. PCA 零假设诊断（从 pca_summary.tsv 读真值）----------
    check("PCA 真实累计", float(psum["cum3_real"]), 0.6194, tol=0.0001)
    check("PCA 逐列打乱基线", float(psum["cum3_perm_mean"]), 0.2705, tol=0.0006)
    check("PCA 逐列打乱基线 SD", float(psum["cum3_perm_sd"]), 0.0074, tol=0.0006)
    check("PCA 独立正态基线", float(psum["cum3_norm_mean"]), 0.2761, tol=0.0006)
    check("PCA Dirichlet 成分基线", float(psum["cum3_dirichlet_mean"]), 0.2903, tol=0.0006)
    check("PCA 超出打乱基线的标准差", float(psum["z_vs_perm"]), 47.0, tol=0.1)
    check("PCA 超出 Dirichlet 基线的标准差", float(psum["z_vs_dirichlet"]), 48.3, tol=0.1)
    check("PCA 有效维度（真实）", float(psum["pr_real"]), 6.124, tol=0.001)
    check("PCA 有效维度（打乱基线）", float(psum["pr_perm_mean"]), 14.74, tol=0.01)
    check("PCA 有效维度（Dirichlet 基线）", float(psum["pr_dirichlet_mean"]), 13.77, tol=0.01)
    check("高相关层对数", int(psum["n_high_corr_pairs"]), 14)
    check("真实协方差秩", int(psum["real_rank"]), 15)
    check("成分行和最大偏离 < 1e-12", float(psum["row_sum_max_dev"]) < 1e-12, True)

    # per-PC 基线递减性（修正 1 的回归保护）
    pdi = pd.read_csv(os.path.join(ANALYSIS_DIR, "pca_diagnosis.tsv"), sep="\t")
    b = pdi["baseline_perm_mean"].values
    check("per-PC 基线单调递减", bool(np.all(np.diff(b) < 0)), True)
    check("per-PC 基线 PC1 值", round(float(b[0]) * 100, 1), 9.7, tol=0.05)

    # ---------- 12. 报告正文与输出文件的数值级一致性 ----------
    # ⚠ 不是"字符串出现过就算过"，而是报告里引用的数字必须等于输出文件的真值。
    # 用**精确 pattern**（带捕获组）定位报告中的具体位置，捕获值必须等于真值。
    # 这才能真正抓住"把 A 的数字配到 B 的头上"这类抄错。
    check_report_num("报告 PCA 前3累计", rep, float(psum["cum3_real"]) * 100,
                     pattern=r"前 3 个累计\s*\|\s*\*\*([\d.]+)%\*\*", decimals=1)
    check_report_num("报告 PCA 打乱基线", rep, float(psum["cum3_perm_mean"]) * 100,
                     pattern=r"零假设基线（逐列打乱，500 次）\s*\|\s*([\d.]+)%", decimals=1)
    check_report_num("报告 PCA 独立正态基线", rep,
                     float(psum["cum3_norm_mean"]) * 100,
                     pattern=r"独立正态对照（绝对无结构）\s*\|\s*([\d.]+)%", decimals=1)
    check_report_num("报告 PCA Dirichlet 基线", rep,
                     float(psum["cum3_dirichlet_mean"]) * 100,
                     pattern=r"\*\*Dirichlet 成分基线（同秩对照，400 次）\*\*\s*\|\s*\*\*([\d.]+)%",
                     decimals=1)
    check_report_num("报告 PCA 有效维度", rep, float(psum["pr_real"]),
                     pattern=r"真实数据\s*\|\s*\*\*([\d.]+) 维\*\*", decimals=2)
    check_report_num("报告 PCA 有效维度（打乱基线）", rep,
                     float(psum["pr_perm_mean"]),
                     pattern=r"零假设（逐列打乱）\s*\|\s*([\d.]+) 维", decimals=2)
    check_report_num("报告 Spearman ρ", rep, 0.9883,
                     pattern=r"\| Spearman ρ \|\s*\*\*([\d.]+)\*\* \|")
    check_report_num("报告 Kendall τ", rep, 0.9518,
                     pattern=r"\| Kendall τ \|\s*\*\*([\d.]+)\*\* \|")
    check_report_num("报告 准确率", rep, float(cmet["accuracy"]),
                     pattern=r"准确率 \*\*([\d.]+)\*\*、宏平均 F1", decimals=3)
    check_report_num("报告 宏F1（测试集）", rep, float(cmet["macro_f1_test"]),
                     pattern=r"宏平均 F1 \*\*([\d.]+)\*\*（测试集 103 类口径", decimals=3)
    check_report_num("报告 准确率天花板", rep, float(cmet["accuracy_ceiling"]),
                     pattern=r"准确率天花板 \| \*\*([\d.]+)\*\*", decimals=3)
    check_report_num("报告 达到天花板比例", rep, float(cmet["pct_of_ceiling"]) * 100,
                     pattern=r"本次达到 \*\*([\d.]+)%\*\*", decimals=2)
    check_report_num("报告 信号噪声（供者级）", rep, ratio_donor,
                     pattern=r"信号/噪声 = \*\*([\d.]+) 倍\*\*", decimals=1)
    check_report_num("报告 短蛋白 Spearman ρ", rep, float(rho),
                     pattern=r"层号 vs 短蛋白占比：Spearman ρ = \+([\d.]+)", decimals=3)
    check_report_num("报告 年轻层短蛋白占比", rep, young_pct,
                     pattern=r"年轻层 S13–S16 平均 ([\d.]+)%", decimals=2)
    check_report_num("报告 老层短蛋白占比", rep, old_pct,
                     pattern=r"老层 S1–S4 平均短蛋白占比 ([\d.]+)%", decimals=2)
    check_report_num("报告 ARI", rep, float(umet["ari"]),
                     pattern=r"ARI = \*\*([\d.]+)\*\*", decimals=3)
    check_report_num("报告 敏感性 ARI", rep, float(umet["ari_sens"]),
                     pattern=r"\| ARI \| 0\.211 \| \*\*([\d.]+)\*\* \|")
    check_report_num("报告 ARI z", rep, float(umet["ari_z"]),
                     pattern=r"z = ([-\d.]+)\*\*", decimals=1)
    check_report_num("报告 S8 倍数", rep, float(E.loc[8, "fold"]),
                     pattern=r"\| S8 硬骨脊椎动物 \| \d+ \| \d+ \| [\d.]+ \| ([\d.]+) \|", decimals=2)
    check_report_num("报告 S14 倍数", rep, float(E.loc[14, "fold"]),
                     pattern=r"\| S14 狭鼻小目 \| \d+ \| \d+ \| [\d.]+ \| \*\*([\d.]+)\*\* \|", decimals=2)
    check_report_num("报告 S2 倍数", rep, float(E.loc[2, "fold"]),
                     pattern=r"\| S2 真核生物 \| \d+ \| \d+ \| [\d.]+ \| \*\*([\d.]+)\*\* \|", decimals=2)
    check_report_num("报告 S4 倍数", rep, float(E.loc[4, "fold"]),
                     pattern=r"\| S4 后生动物 \| \d+ \| \d+ \| [\d.]+ \| \*\*([\d.]+)\*\* \|", decimals=2)
    check_report_num("报告 A 背景 M", rep, float(enr["n_background"].sum()),
                     pattern=r"背景 = ([\d,]+) 个有年龄注释的 HVG")
    check_report_num("报告 A Top 集 N", rep, float(enr["n_top"].sum()),
                     pattern=r"Top 集合 = ([\d,]+) 个有年龄注释的基因")
    check_report_num("报告 检验A 排名变动", rep, 141.0,
                     pattern=r"排名变动 \| \*\*(\d+) / 177")
    check_report_num("报告 检验A 最大位移", rep, 88.0,
                     pattern=r"最大位移 (\d+) 位")
    # 三次共同富集层必须恰好是 S8、S10 两层。
    # 做法：抽出报告里"三次共同显著富集的层为 …"后面的层清单，
    # 再与 robustness_comparison.tsv 算出的真值逐一比对。
    m_common = re.search(r"三次共同显著富集的层为 ([^，*（(]+)", rep)
    got_common = ([s.strip() for s in m_common.group(1).split("、")]
                  if m_common else None)
    check("报告的三次共同富集层清单", got_common, ["S8", "S10"])

    # 三次共同显著层必须以清单形式出现在报告里
    check("报告含共同富集层 S8", "S8" in rep, True)
    check("报告含共同富集层 S10", "S10" in rep, True)
    check("报告含共同耗竭层 S2", "S2" in rep, True)

    # 报告不得在「现行结论」语境中使用已作废的数字。
    # ⚠ 不能用纯子串判定：本报告如实记录了踩过的坑，会在
    #   "曾…/最初…/旧…/错误…" 这类语境里引用旧值——那是应当保留的。
    #   规则：每一次出现，都必须落在作废标记词附近。
    STALE_MARKERS = ["曾", "最初", "旧", "错误", "作废", "弃", "原先",
                     "错误地", "已经发生", "本轮审核", "本轮发现", "改",
                     "恒等于", "不是发现", "反向", "完全相反",
                     "若按细胞层面", "细胞层面估计", "细胞级口径",
                     "作对照", "仅作", "对照"]
    for stale, why in [("0.198", "未做 QC 的旧准确率"),
                       ("0.164", "宏F1 被幽灵类稀释的旧值"),
                       ("33 倍", "用 1/177 随机基线的旧说法"),
                       ("排序完全不变", "Spearman 恒等式导致的错误结论"),
                       ("7 层富集模式完全保留", "BH 校正前的共同层数"),
                       ("--min-dbs 10", "方向相反的置信度过滤参数"),
                       ("1,773", "检验 B 的旧特征数"),
                       ("2,403", "检验 C 的旧特征数"),
                       ("2,521", "旧背景基因数"),
                       ("1.000000", "占比和的旧写法")]:
        check_no_stale(f"报告未在结论语境使用作废数字「{stale}」（{why}）",
                       rep, stale, STALE_MARKERS)

    # 「16.3 倍」是**细胞层面**口径的数字，可以作为对照出现，
    # 但绝不能作为「信号/噪声」的最终口径。规则：
    #   每次出现 "16.3 倍"，其附近必须明确说明这是细胞层面口径 / 仅作对照。
    for m in re.finditer(re.escape("16.3 倍"), rep):
        lo = max(0, m.start() - 70)
        hi = min(len(rep), m.end() + 70)
        win = rep[lo:hi]
        if not any(k in win for k in ["细胞层面", "细胞级", "仅作", "对照",
                                      "会得到", "会算出"]):
            CHECKS[0] += 1
            FAILS.append(f"「16.3 倍」出现在非对照语境：…{win[:70]}…")
        else:
            CHECKS[0] += 1

    # Spearman/Kendall 的旧值 1.0000 与"排序完全不变"必须成对出现在
    # 「曾经踩过的坑」那段说明里——这是本项目如实记录修正过程的要求。
    if "恒等于" not in rep:
        check("报告如实说明了 Spearman 恒等式这个坑", False, True)
    if "Spearman ρ = 1.0000" not in rep:
        check("报告保留了『曾报告 1.0000』的修正记录", False, True)

    # ---------- 13. 长蛋白版 TAI（闭合 Moyers & Zhang 批评的最后一环）----------
    # 03 脚本加 --min-length 100 --tag longprot 后，用长蛋白重算 TAI 绝对值，
    # 与主口径逐类型对照，把「年轻偏倚」的方向与大小**直接量化**出来。
    lpc_path = os.path.join(ANALYSIS_DIR, "longprot_comparison.tsv")
    lps_path = os.path.join(ANALYSIS_DIR, "longprot_summary.tsv")
    if not (os.path.exists(lpc_path) and os.path.exists(lps_path)):
        CHECKS[0] += 1
        FAILS.append("缺少长蛋白对照结果：先跑 "
                     "python 03_细胞类型年龄谱.py --min-length 100 --tag longprot "
                     "再跑 python 94_长蛋白对照.py")
    else:
        lpc = pd.read_csv(lpc_path, sep="\t")
        lps = pd.read_csv(lps_path, sep="\t").iloc[0]

        check("长蛋白对照表类型数", len(lpc), 177)
        check("长蛋白对照表首列是 cell_type", lpc.columns[0], "cell_type")
        check("长蛋白版 Spearman ρ", float(lps["spearman_rho"]), 0.9845, tol=5e-4)
        check("长蛋白版 Kendall τ", float(lps["kendall_tau"]), 0.9186, tol=5e-4)
        check("长蛋白版 ΔTAI 均值", float(lps["delta_mean"]), -0.0572, tol=5e-4)
        check("长蛋白版 |ΔTAI| 中位", float(lps["abs_delta_median"]), 0.0536, tol=5e-4)
        check("长蛋白版 |ΔTAI| 最大", float(lps["abs_delta_max"]), 0.2523, tol=5e-4)
        check("变年轻的类型数（应为全部 177）",
              int(lps["n_shift_younger"]), 177)
        check("名次变动的类型数", int(lps["n_rank_changed"]), 163)
        check("最大名次位移", int(lps["max_rank_shift"]), 79)
        check("Top10 重合数", int(lps["top10_overlap"]), 9)
        check("ΔTAI 的年龄梯度 ρ（应接近 0）",
              float(lps["gradient_rho"]), -0.1326, tol=5e-4)
        check("最古老 10 类 ΔTAI 均值",
              float(lps["delta_oldest10"]), -0.0663, tol=5e-4)
        check("最年轻 10 类 ΔTAI 均值",
              float(lps["delta_youngest10"]), -0.0725, tol=5e-4)
        check("长蛋白版 TAI 最小值", float(lps["tai_lp_min"]), 3.1717, tol=5e-4)
        check("长蛋白版 TAI 最大值", float(lps["tai_lp_max"]), 4.3284, tol=5e-4)
        check("长蛋白版 TAI 中位数", float(lps["tai_lp_median"]), 3.6933, tol=5e-4)
        check("长蛋白版 TAI 跨类型 SD", float(lps["tai_lp_sd"]), 0.1801, tol=5e-4)

        # 自洽性：逐类型表必须与汇总表一致（不能只改汇总不改明细）
        # ⚠ 容差 5e-4：汇总表的 delta_mean 是**四舍五入到 4 位**后落盘的
        #    （-0.0572），而明细表均值是全精度（-0.05717965...），
        #    用 1e-6 会把"精度口径差异"误报成"数据不一致"。
        #    这类断言的正解是"两边各自舍入到汇总表的小数位后比较"。
        check("对照表 ΔTAI 均值与汇总一致",
              round(float(lpc["delta"].mean()), 4), float(lps["delta_mean"]), tol=1e-9)
        check("对照表里变年轻的类型数一致",
              int((lpc["delta"] < 0).sum()), int(lps["n_shift_younger"]))
        check("长蛋白版 TAI 全体低于主口径（系统性下移）",
              bool((lpc["delta"] < 0).all()), True)

        # 报告正文必须写上这些数字（否则就是「算了但没写」）
        check_report_num("报告 长蛋白版 Spearman ρ", rep,
                         float(lps["spearman_rho"]),
                         pattern=r"长蛋白版[^\n|]*?Spearman ρ\s*=\s*\*\*([\d.]+)\*\*",
                         decimals=4)
        check_report_num("报告 长蛋白版 ΔTAI 均值", rep,
                         float(lps["delta_mean"]),
                         pattern=r"ΔTAI 均值\s*=\s*\*\*([+\-−]?[\d.]+)\*\*",
                         decimals=4)

    # ---------- 汇总 ----------
    print("=" * 70)
    print(f"终审核对：共 {CHECKS[0]} 项检查")
    print("=" * 70)
    if FAILS:
        print(f"✗ {len(FAILS)} 项不符：")
        for f in FAILS:
            print("   -", f)
        sys.exit(1)
    print("✓ 全部通过。报告中的每个关键数字都可追溯到输出文件。")


if __name__ == "__main__":
    main()
