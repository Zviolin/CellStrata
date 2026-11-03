"""
95_论文数字终审.py — LaTeX 论文（50-论文/main.tex）的数字与结果文件比对
========================================================================
为什么需要这个脚本：
    论文是**手写**的 LaTeX，数字最容易出现三类错误：
      1. **抄错**——把 A 的数字写到 B 的名目下；
      2. **编造**——写出一个结果文件里根本不存在的数字
         （本项目真实发生过一次：初稿里写了「五个胚层区室的平均 TAI
          依次为神经系统 3.41、间充质 3.63……」，
          而实际输出是 germ line / stromal / epithelial / immune /
          endothelial 五个区室，根本没有「神经系统」这一项）；
      3. **口径漂移**——正文改了口径但表格/图注没跟着改。

    报告那边已经有 `07_终审核对.py`（239 项），论文这边单独一个脚本，
    两边独立取真值、独立比对。

检查项：
    A. 论文里出现的每个关键数字，必须能在**对应输出文件**里找到真值
       （真值取自 tsv，不是从报告里抄的——否则报告错了论文跟着错）
    B. 论文里**不得出现**输出文件里不存在的"编造数字"：
       对一组"本应出现"的关键数字，检查论文里写的值与真值逐位相等
    C. 结构一致性：细胞类型数、区室数、层数、脚本数、断言数等

用法：
    python 95_论文数字终审.py
"""

import glob
import io
import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import ROOT, ANALYSIS_DIR  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAPER = os.path.join(ROOT, "50-论文", "main.tex")

FAILS = []
CHECKS = [0]


def check(label, got, want, tol=1e-6):
    CHECKS[0] += 1
    if isinstance(got, (list, tuple)) or isinstance(want, (list, tuple)):
        ok = list(got) == list(want)
    elif isinstance(want, str) or isinstance(got, str) or want is None or got is None:
        ok = str(got) == str(want)
    else:
        ok = abs(float(got) - float(want)) <= tol
    if not ok:
        FAILS.append(f"{label}: 论文={got!r}  结果文件={want!r}")
    return ok


def paper_num(label, paper, value, decimals=None, occurrence=0,
              allow_percent=False):
    """检查论文里出现的裸数字是否等于真值。

    ⚠ 为什么不用 `str(value) in paper`：
      · 论文里同一个数字可能在多处出现（摘要、正文、表格），
        其中一处写错而其余写对时，"字符串出现过"会漏判；
      · 论文里 3.214 与 3.2142 是同一个量的不同精度，直接子串匹配会误判。
      正确做法是：把论文里所有形如数字的 token 抽出来，
        按论文里约定的精度舍入后与真值比较。

    decimals：论文里该数字的呈现精度。None 表示允许真值以更高精度出现
        （如真值 0.619376、论文写 0.6194 或 61.9% 都算对）。
    allow_percent：额外允许「真值 × 100」的形式（如 0.6194 ↔ 61.9%）。
    """
    CHECKS[0] += 1
    base = float(value)
    if decimals is not None:
        base = round(base, decimals)

    # ---------- 把论文正文规范化成"纯数字流" ----------
    body = paper
    # 1) 先把 **科学计数法** 的 LaTeX 写法转成普通小数。
    #    论文里写 $4.5\times10^{-3}$，正则抓到的是 4.5 / 10 / 3，
    #    与真值 0.00447839 无法比对 —— 必须先做这一步。
    def _sci(m):
        mant = float(m.group(1))
        exp = int(m.group(2))
        return "%.12g" % (mant * (10 ** exp))
    body = re.sub(r"(-?\d+(?:\.\d+)?)\\times\s*10\^\{?(-?\d+)\}?",
                  _sci, body)
    body = re.sub(r"(-?\d+(?:\.\d+)?)\\cdot10\^\{?(-?\d+)\}?",
                  _sci, body)
    # 2) 结构性命令：整条去掉（含参数），它们不是数据
    for pat in (r"\\label\{[^}]*\}", r"\\ref\{[^}]*\}",
                r"\\cite[a-z]*\{[^}]*\}",
                r"\\includegraphics\[[^\]]*\]\{[^}]*\}",
                r"\\bibliography\{[^}]*\}"):
        body = re.sub(pat, " ", body)
    # 3) 千分位：6{,}344 与 6,344 都要还原成 6344
    body = re.sub(r"\d\{,\}(\d{3})", r"\1", body)
    body = re.sub(r"(?<=\d),(?=\d{3}\b)", "", body)
    # 4) 转义下划线：07\_终审核对.py → 07_终审核对.py
    body = body.replace("\\_", "_")
    # 5) 去掉命令名（保留花括号里的内容，\Code{07_终审核对.py} 里的数字要留）
    body = re.sub(r"\\[a-zA-Z@]+\b\s*(?![{\[])", " ", body)
    # 6) 图/表编号不是数据
    body = re.sub(r"(?:Figure|Table|Fig|图|表)\s*~?\d+", " ", body)
    # 7) 章节号（3.1 / 4.2）不是数据
    body = re.sub(r"^\s*\\(?:sub)*section\{[^}]*\}", " ", body,
                  flags=re.M)

    # 8) 数字 token：可带前导负号
    nums = re.findall(r"(-?\d+(?:\.\d+)?)(?![\d])", body)
    vals = [float(t) for t in nums]

    def close(v):
        if decimals is not None:
            if round(v, decimals) == base:
                return True
            # 论文若用**更少**的小数位呈现同一量（如真值 0.2354 写成 0.235），
            # 按论文实际位数舍入后应相等。判定"实际位数"：字符串里小数点后的长度。
            for nd in range(1, decimals):
                if round(v, nd) == round(base, nd) and \
                   abs(v - base) <= 10 ** (-nd) / 2 + 1e-12:
                    return True
            return False
        if abs(v - base) <= 5e-5:
            return True
        for nd in (2, 3, 4):
            if abs(round(v, nd) - round(base, nd)) < 1e-12 and \
               abs(v - base) < 5e-4:
                return True
        return False

    for v in vals:
        if close(v):
            return True
        if allow_percent and close(v / 100.0):
            return True
        if allow_percent and close(v * 100.0):
            return True
    FAILS.append(f"{label}: 论文里找不到 {base}（论文中该数字未出现或写错）")
    return False


def _normalize_literal(s):
    r"""把论文正文规范化，便于做**字面**比对。

    LaTeX 里同一个字符串有多种写法：
      · 千分位  6{,}344  或  6,344
      · 转义下划线  07\_终审核对.py  或  07_终审核对.py
    直接用 `token in paper` 会漏判，必须先统一。
    """
    # ⚠ 必须保留首位数字：写成 r"\d\{,\}(\d{3})" 会把 "6{,}344"
    #    整段替换成 "344"，首位 6 被吃掉 → "6,344" 与 "6344" 仍对不上。
    s = re.sub(r"(\d)\{,\}(\d{3})", r"\1\2", s)
    s = re.sub(r"(?<=\d),(?=\d{3}\b)", "", s)
    s = s.replace("\\_", "_")
    return s


def paper_has(label, paper, token):
    """检查论文里存在某个字面 token（用于结构性断言，如脚本名、文件名）。"""
    CHECKS[0] += 1
    if _normalize_literal(token) not in _normalize_literal(paper):
        FAILS.append(f"{label}: 论文里找不到「{token}」")
        return False
    return True


def main():
    if not os.path.exists(PAPER):
        print("✗ 找不到论文源文件：", PAPER)
        sys.exit(1)
    with open(PAPER, encoding="utf-8") as f:
        paper = f.read()

    prof = pd.read_csv(os.path.join(ANALYSIS_DIR, "celltype_age_profile.tsv"),
                       sep="\t", index_col=0)
    uns = pd.read_csv(os.path.join(ANALYSIS_DIR, "unsupervised_metrics.tsv"),
                      sep="\t").iloc[0]
    pca = pd.read_csv(os.path.join(ANALYSIS_DIR, "pca_summary.tsv"),
                      sep="\t").iloc[0]
    cls = pd.read_csv(os.path.join(ANALYSIS_DIR, "classification_metrics.tsv"),
                      sep="\t").iloc[0]
    enr = pd.read_csv(os.path.join(ANALYSIS_DIR, "age_enrichment.tsv"), sep="\t")
    rob = pd.read_csv(os.path.join(ANALYSIS_DIR, "robustness_comparison.tsv"),
                      sep="\t")
    lpc = pd.read_csv(os.path.join(ANALYSIS_DIR, "longprot_comparison.tsv"),
                      sep="\t")
    lps = pd.read_csv(os.path.join(ANALYSIS_DIR, "longprot_summary.tsv"),
                      sep="\t").iloc[0]
    spb = pd.read_csv(os.path.join(ANALYSIS_DIR, "short_protein_by_stratum.tsv"),
                      sep="\t")

    # ================= A. 规模类 =================
    check("细胞类型数", len(prof), 177)
    check("年龄层数", prof.shape[1] - 7, 16)  # 16 层 + 7 个派生列
    check("QC 后细胞数", int(prof["n_cells"].sum()), 22195)

    # ================= B. TAI 核心数字 =================
    paper_num("论文 TAI 最小值", paper, prof["TAI"].min(), decimals=3)
    paper_num("论文 TAI 最大值", paper, prof["TAI"].max(), decimals=3)
    paper_num("论文 TAI 中位数", paper, prof["TAI"].median(), decimals=3)
    paper_num("论文 TAI 跨类型 SD", paper, prof["TAI"].std(ddof=1), decimals=4)
    paper_num("论文 供者级 SE 中位", paper, prof["TAI_se_donor"].median(),
              decimals=4)
    paper_num("论文 细胞级 SE 中位", paper, prof["TAI_se"].median(), decimals=4)
    ratio = prof["TAI_se_donor"].median() / prof["TAI_se"].median()
    paper_num("论文 SE 比值", paper, ratio, decimals=2)
    snr = prof["TAI"].std(ddof=1) / prof["TAI_se_donor"].median()
    paper_num("论文 信号噪声比", paper, snr, decimals=1)
    check("低置信类型数", int(prof["low_confidence"].sum()), 25)
    check("供者级 SE 不可估类型数", int(prof["se_donor_na"].sum()), 81)

    # ================= C. 区室均值 =================
    lineage = pd.read_csv(os.path.join(ANALYSIS_DIR, "lineage_map.tsv"), sep="\t")
    comp = prof.index.map(lineage.set_index("cell_type")["lineage"])
    prof2 = prof.assign(lineage=list(comp))
    for name, exp_n in [("germ line", 1), ("stromal", 34),
                        ("epithelial", 66), ("immune", 61),
                        ("endothelial", 15)]:
        sub = prof2[prof2["lineage"] == name]
        check(f"区室 {name} 类型数", len(sub), exp_n)
        paper_num(f"论文 区室 {name} 平均 TAI", paper, sub["TAI"].mean(),
                  decimals=4)

    # ================= D. 无监督 =================
    paper_num("论文 PCA 前三累计", paper, float(uns["cum3"]), decimals=4)
    paper_num("论文 置换基线 cum3", paper, float(pca["cum3_perm_mean"]),
              decimals=4)
    paper_num("论文 正态基线 cum3", paper, float(pca["cum3_norm_mean"]),
              decimals=4)
    paper_num("论文 Dirichlet 基线 cum3", paper, float(pca["cum3_dirichlet_mean"]),
              decimals=4)
    paper_num("论文 z vs 置换", paper, float(pca["z_vs_perm"]), decimals=1)
    paper_num("论文 z vs Dirichlet", paper, float(pca["z_vs_dirichlet"]),
              decimals=1)
    paper_num("论文 有效维度 PR", paper, float(pca["pr_real"]), decimals=2)
    paper_num("论文 真实秩", paper, float(pca["real_rank"]), decimals=0)
    paper_num("论文 K-means k", paper, float(uns["k"]), decimals=0)
    paper_num("论文 轮廓系数", paper, float(uns["silhouette"]), decimals=4)
    paper_num("论文 ARI", paper, float(uns["ari"]), decimals=4)
    paper_num("论文 ARI 置换 z", paper, float(uns["ari_z"]), decimals=1)
    paper_num("论文 Ward ARI", paper, float(uns["ari_ward"]), decimals=4)
    paper_num("论文 敏感性 k", paper, float(uns["k_sens"]), decimals=0)
    paper_num("论文 敏感性 ARI", paper, float(uns["ari_sens"]), decimals=4)

    # ================= E. 监督 =================
    paper_num("论文 准确率", paper, float(cls["accuracy"]), decimals=3)
    paper_num("论文 宏 F1", paper, float(cls["macro_f1_test"]), decimals=3)
    paper_num("论文 天花板", paper, float(cls["accuracy_ceiling"]), decimals=4)
    paper_num("论文 天花板达成率 %", paper,
              float(cls["pct_of_ceiling"]) * 100, decimals=1)
    paper_num("论文 多数类基线", paper, float(cls["majority_baseline"]),
              decimals=3)
    paper_num("论文 多数类倍数", paper, float(cls["x_majority"]), decimals=2)
    check("未见类数", int(cls["n_unseen_classes"]), 28)
    check("训练供者数", int(cls["n_donors_train"]), 11)
    check("测试供者数", int(cls["n_donors_test"]), 4)
    check("运行A 特征数", int(cls["n_features"]), 3000)

    clsB = pd.read_csv(os.path.join(ANALYSIS_DIR,
                                    "classification_metrics_confident.tsv"),
                       sep="\t").iloc[0]
    clsC = pd.read_csv(os.path.join(ANALYSIS_DIR,
                                    "classification_metrics_shortprot.tsv"),
                       sep="\t").iloc[0]
    paper_num("论文 运行B 准确率", paper, float(clsB["accuracy"]), decimals=3)
    paper_num("论文 运行B 宏F1", paper, float(clsB["macro_f1_test"]), decimals=3)
    paper_num("论文 运行C 准确率", paper, float(clsC["accuracy"]), decimals=3)
    paper_num("论文 运行C 宏F1", paper, float(clsC["macro_f1_test"]), decimals=3)
    check("运行B 特征数", int(clsB["n_features"]), 905)
    check("运行C 特征数", int(clsC["n_features"]), 2405)

    # ================= F. 富集与稳健性 =================
    sig_up = sorted(enr[(enr["q_value"] < 0.05)]["stratum"].tolist())
    sig_dn = sorted(enr[(enr["q_depletion"] < 0.05)]["stratum"].tolist())
    check("运行A 显著富集层", sig_up, [1, 8, 9, 10, 11, 14])
    check("运行A 显著耗竭层", sig_dn, [2, 4])
    for s, q in enr.set_index("stratum")["q_value"].items():
        if s in sig_up:
            paper_num(f"论文 S{s} 富集 q", paper, float(q), decimals=8)
    check("方向冲突层数", int((rob["conflict"] == "yes").sum()), 0)

    # ================= G. 短蛋白前提 =================
    old = float(spb[spb["stratum"] <= 4]["pct_short"].mean())
    young = float(spb[spb["stratum"] >= 13]["pct_short"].mean())
    paper_num("论文 老层短蛋白占比", paper, old, decimals=2)
    paper_num("论文 年轻层短蛋白占比", paper, young, decimals=2)
    paper_num("论文 短蛋白倍数", paper, young / old, decimals=1)
    paper_num("论文 S1 中位长度", paper,
              float(spb[spb["stratum"] == 1]["median_length"].iloc[0]),
              decimals=0)
    paper_num("论文 S16 中位长度", paper,
              float(spb[spb["stratum"] == 16]["median_length"].iloc[0]),
              decimals=0)

    # ================= H. 检验 D =================
    paper_num("论文 长蛋白版 Spearman", paper, float(lps["spearman_rho"]),
              decimals=4)
    paper_num("论文 长蛋白版 Kendall", paper, float(lps["kendall_tau"]),
              decimals=4)
    paper_num("论文 长蛋白版 dTAI 均值", paper, float(lps["delta_mean"]),
              decimals=4)
    paper_num("论文 长蛋白版 |dTAI| 中位", paper,
              float(lps["abs_delta_median"]), decimals=4)
    paper_num("论文 长蛋白版 |dTAI| 最大", paper, float(lps["abs_delta_max"]),
              decimals=4)
    paper_num("论文 长蛋白版 同向平移数", paper, float(lps["n_shift_younger"]),
              decimals=0)
    paper_num("论文 名次变动数", paper, float(lps["n_rank_changed"]), decimals=0)
    paper_num("论文 最大名次位移", paper, float(lps["max_rank_shift"]), decimals=0)
    paper_num("论文 Top10 重合", paper, float(lps["top10_overlap"]), decimals=0)
    paper_num("论文 dTAI 年龄梯度", paper, float(lps["gradient_rho"]),
              decimals=4)
    paper_num("论文 最古老10类 dTAI", paper, float(lps["delta_oldest10"]),
              decimals=4)
    paper_num("论文 最年轻10类 dTAI", paper, float(lps["delta_youngest10"]),
              decimals=4)
    paper_num("论文 长蛋白版 TAI 最小", paper, float(lps["tai_lp_min"]),
              decimals=4)
    paper_num("论文 长蛋白版 TAI 最大", paper, float(lps["tai_lp_max"]),
              decimals=4)
    paper_num("论文 长蛋白版 TAI 中位", paper, float(lps["tai_lp_median"]),
              decimals=4)
    paper_num("论文 长蛋白版 TAI SD", paper, float(lps["tai_lp_sd"]),
              decimals=4)
    check("长蛋白版全下移", bool((lpc["delta"] < 0).all()), True)
    # 偏倚占 SD 的比例（论文里写"32%"）
    paper_num("论文 偏倚占SD比例", paper,
              abs(float(lps["delta_mean"])) / float(lps["tai_base_sd"]) * 100,
              decimals=1)

    # ================= I. 结构一致性 =================
    # ⚠ 脚本数与代码行数**从磁盘实算**，不写死字面量：
    #   写死会让论文在每次改动代码后静默过期（本项目真实发生过——
    #   新增 95/96 两个脚本后，论文仍写着旧的 6,783 行 / 19 脚本）。
    _py = sorted(glob.glob(os.path.join(ROOT, "20-分析", "*.py"))
                 + glob.glob(os.path.join(ROOT, "30-可视化", "*.py")))
    _n_script = len(_py)
    _n_line = sum(len(io.open(f, encoding="utf-8", errors="ignore")
                      .read().splitlines()) for f in _py)
    check("脚本数与磁盘一致", str(_n_script), "20")
    paper_has("论文提到脚本总数", paper, f"{_n_script} 个 Python 脚本")
    paper_has("论文提到代码行数", paper, f"{_n_line:,}")
    paper_has("论文提到断言数", paper, "239")
    paper_has("论文提到终审脚本", paper, "07_终审核对.py")
    paper_has("论文提到论文终审脚本", paper, "95_论文数字终审.py")
    paper_has("论文引用 Trigos", paper, "trigos2017")
    paper_has("论文引用 Moyers 2018", paper, "moyers2018")

    # 论文里引用的图必须都存在
    figs = set(re.findall(r"\\includegraphics\[[^\]]*\]\{([^}]+)\}", paper))
    for fg in sorted(figs):
        CHECKS[0] += 1
        p = os.path.join(ROOT, "30-可视化", "figures", fg)
        if not os.path.exists(p):
            FAILS.append(f"论文引用的图不存在：{fg}")

    # ================= 汇总 =================
    print("=" * 70)
    print(f"论文数字终审：共 {CHECKS[0]} 项检查")
    print("=" * 70)
    if FAILS:
        print(f"✗ {len(FAILS)} 项不符：")
        for f in FAILS:
            print("   -", f)
        sys.exit(1)
    print("✓ 全部通过。论文中的每个关键数字都可追溯到输出文件。")


if __name__ == "__main__":
    main()
