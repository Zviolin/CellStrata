"""
96_导出报告PDF.py — 把正式报告 Markdown 导出为排版良好的 PDF
==================================================================
路线：Markdown → HTML（自绘样式）→ Word COM → PDF
    （本机无 pandoc；Word COM 已验证可用，ExportAsFixedFormat(pdf, 17)）

为什么走 HTML → Word 而不是直接生成 docx：
    · Word 对 HTML 的表格/图片/样式支持最好，能保留中文字体与代码块底色
    · docx 转 PDF 时页面尺寸/页边距控制不如 Word 直接排版稳

用法：
    python 96_导出报告PDF.py
输出：
    40-报告/01-实验报告.pdf
"""

import os
import re
import sys
import base64
import hashlib
import shutil
import subprocess

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import ROOT, info, ok, die  # noqa: E402

import markdown  # noqa: E402

REPORT_MD = os.path.join(ROOT, "40-报告", "01-实验报告.md")
FIG_DIR = os.path.join(ROOT, "30-可视化", "figures")
OUT_HTML = os.path.join(ROOT, "40-报告", "_report_preview.html")
OUT_PDF = os.path.join(ROOT, "40-报告", "01-实验报告.pdf")
MATH_DIR = os.path.join(ROOT, "40-报告", "_math")
MATH_DPI = 600

# ---------- 图片：按正文出现的"（**图 N**）"锚点插入 ----------
# 报告正文用「（**图 1**、**图 3**）」这样的形式引用图。
# 导出 PDF 时把对应图片**内嵌**到该处，让 PDF 既是文字报告也是图册。
FIG_ORDER = [
    (1, "fig1_strata.png", "图 1  细胞地层主图：177 种细胞类型按 TAI 排序，每根按 16 个年龄层堆叠"),
    (2, "fig2_pca.png", "图 2  PCA 前两主成分散点，按 K-means 簇着色"),
    (3, "fig3_heatmap.png", "图 3  细胞类型 × 年龄层的表达占比热图"),
    (4, "fig4_importance_age.png", "图 4  各年龄层基因数与 Top 500 重要基因的富集倍数"),
    (5, "fig5_tai_rank.png", "图 5  177 种细胞类型的 TAI 排序与双口径误差棒"),
    (6, "fig6_umap_tai.png", "图 6  细胞的 PC1/PC2 投影，按所属细胞类型的 TAI 着色"
                             "（切片源无 UMAP 坐标，用稀疏 SVD 的 PC1/PC2 代替）"),
    (7, "fig7_pca_null_baseline.png", "图 7  PCA 零假设诊断：三条基线对比"),
    (8, "fig8_robustness_layers.png", "图 8  三重稳健性逐层对照"),
    (9, "fig9_shortprot_premise.png", "图 9  短蛋白前提检验"),
    (10, "fig10_longprot_comparison.png", "图 10  长蛋白版 TAI 绝对值对照（检验 D）"),
]

# =====================================================================
#  版式：学术中文报告体例
#  ⚠ 三处反复踩过的坑，改样式前先读：
#   1) `line-height: 1.65` 这种无单位倍数，Word 的 HTML 导入器会把它解释成
#      **固定 12 pt 行距**（10.5 pt 字只有 1.14 倍），页面会又闷又挤。
#      真正的行距必须在 COM 里用 LineSpacingRule=wdLineSpaceMultiple 设置。
#      这里的 CSS 只作兜底，COM 里的 _polish_doc() 会覆盖它。
#   2) Word 会用 **class 名建段落样式**（实测 figcap/figimg 成了独立样式），
#      所以 class 是有意义的钩子，但 CSS 属性只被支持一部分（无 flex、无
#      max-width、无 box-shadow、无 :last-child）。
#   3) 中文正文用宋体、西文用 Times New Roman 靠 font-family 列表**不可靠**，
#      必须在 COM 里分别设 Font.NameFarEast / NameAscii。
# =====================================================================
CSS = """
@page { size: A4; margin: 2.5cm 2.2cm 2.4cm 2.2cm; }
body {
  font-family: "SimSun", "Times New Roman", serif;
  font-size: 10.5pt; line-height: 1.55; color: #000;
  text-align: justify; margin: 0;
}

/* ---------- 封面 ---------- */
h1.doctitle {
  font-family: "SimHei", "Microsoft YaHei", sans-serif;
  font-size: 21pt; font-weight: bold; text-align: center;
  line-height: 1.45; margin: 0 0 20pt 0; text-indent: 0;
  letter-spacing: 0.5pt;
}
p.meta {
  font-size: 10pt; text-align: center; text-indent: 0;
  line-height: 1.75; margin: 0; color: #222;
}

/* ---------- 标题层级 ---------- */
h2 {
  font-family: "SimHei", "Microsoft YaHei", sans-serif;
  font-size: 15pt; font-weight: bold; text-align: left; text-indent: 0;
  margin: 20pt 0 10pt 0; padding: 0 0 3pt 0;
  border-bottom: 0.9pt solid #000; line-height: 1.35;
  page-break-after: avoid;
}
h3 {
  font-family: "SimHei", "Microsoft YaHei", sans-serif;
  font-size: 12.5pt; font-weight: bold; text-align: left; text-indent: 0;
  margin: 15pt 0 6pt 0; line-height: 1.4; page-break-after: avoid;
}
h4 {
  font-family: "SimHei", "Microsoft YaHei", sans-serif;
  font-size: 11pt; font-weight: bold; text-align: left; text-indent: 0;
  margin: 11pt 0 5pt 0; line-height: 1.4; page-break-after: avoid;
}

/* ---------- 摘要区 ---------- */
h2.absh {
  border-bottom: none; text-align: center; font-size: 13pt;
  letter-spacing: 6pt; margin: 16pt 0 12pt 0;
}
p.abs { font-size: 10pt; line-height: 1.6; text-indent: 2em; margin: 0 0 6pt 0; }
p.kw  { font-size: 10pt; line-height: 1.6; text-indent: 0; margin: 10pt 0 0 0; }
p.pbk { page-break-before: always; margin: 0; font-size: 1pt; line-height: 1pt; }

/* ---------- 正文 ---------- */
p { margin: 0 0 0.28em 0; text-indent: 2em; text-align: justify; }
li p, td p, th p { text-indent: 0; }
ul, ol { margin: 4pt 0 8pt 0; padding-left: 2.6em; }
li { text-indent: 0; line-height: 1.55; margin: 0 0 2pt 0; }

/* ---------- 表格：三线表 ---------- */
table {
  border-collapse: collapse; width: 100%;
  margin: 9pt 0 13pt 0; font-size: 9pt; text-indent: 0;
  line-height: 1.35; page-break-inside: avoid;
}
th, td { padding: 3.5pt 5pt; vertical-align: top; text-align: left; }
th {
  font-family: "SimHei", "Microsoft YaHei", sans-serif;
  font-weight: bold; border-top: 1.2pt solid #000;
  border-bottom: 0.6pt solid #000; background: none;
}
td { border: none; }
/* 底部的收口线：Word 不支持 :last-child，由 COM 的 _polish_doc() 补 */

code {
  font-family: "Consolas", "Courier New", monospace; font-size: 9pt;
  background: #f1f2f4;
}
pre {
  font-family: "Consolas", "Courier New", monospace;
  background: #f7f8fa; border: 0.5pt solid #c8ccd2;
  padding: 7pt 9pt; font-size: 8.5pt; line-height: 1.4;
  text-indent: 0; page-break-inside: avoid; margin: 8pt 0;
}
pre code { background: none; font-size: 8.5pt; }

/* ---------- 注记（引用块）：收敛为浅灰 ---------- */
blockquote {
  margin: 9pt 0; padding: 7pt 10pt; background: #f4f5f7;
  border-left: 3pt solid #7b8494; text-indent: 0;
  font-size: 9.8pt; line-height: 1.5;
}
blockquote p { margin: 0 0 3pt 0; text-indent: 0; text-align: left; }
blockquote p:last-child { margin-bottom: 0; }

hr { border: none; border-top: 0.5pt solid #b8b8b8; margin: 12pt 0; }

/* ---------- 图 ---------- */
p.figimg { text-align: center; text-indent: 0; margin: 10pt 0 0 0;
           page-break-after: avoid; }
p.figimg img { border: 0.5pt solid #d5d5d5; }
p.figcap {
  text-align: center; text-indent: 0; font-size: 9pt;
  color: #222; margin: 4pt 0 14pt 0; line-height: 1.4;
}

/* ---------- 公式 ---------- */
.math { text-align: center; text-indent: 0; margin: 10pt 0; line-height: 1.6; }
.math img { vertical-align: middle; }
.math sub, .math sup { font-size: 72%; }
p.mathimg {
  text-align: center; text-indent: 0; margin: 10pt 0;
  line-height: 1; page-break-inside: avoid;
}
sub, sup { font-size: 72%; }
a { color: #000; text-decoration: none; }
strong { font-weight: bold; }
"""



# =====================================================================
#  极简 LaTeX → HTML 公式转换
#  报告里的公式只有 \frac / \sum / \text / \mathrm / \ge / \le / \cdot 这几类，
#  而 Word 的 HTML 导入器不认 MathML、也装不了 MathJax（离线备份场景），
#  所以在这里就把公式降级成"带括号的分式 + 上下标"，保证任何环境都能读。
#  ⚠ 一旦报告新增复杂公式（矩阵、开方、多层分式嵌套），这里必须同步扩展，
#    否则会退化成原样输出 LaTeX 源码——这正是修复前 PDF 里出现 `$$ ... $$` 的原因。
# =====================================================================
_SYMBOLS = {
    r"\ge": "\u2265", r"\le": "\u2264", r"\cdot": "\u00b7",
    r"\times": "\u00d7", r"\in": "\u2208", r"\pm": "\u00b1",
    r"\approx": "\u2248", r"\to": "\u2192", r"\qquad": "\u2003\u2003",
    r"\quad": "\u2003", r"\%": "%", r"\,": "\u2009", r"\;": "\u2005",
    r"\alpha": "\u03b1", r"\beta": "\u03b2", r"\rho": "\u03c1",
    r"\Delta": "\u0394", r"\sigma": "\u03c3",
}


def _take_group(s, i):
    """从 s[i] 起取一个 {...} 组（或单个字符），返回 (内容, 结束位置)。"""
    while i < len(s) and s[i] == " ":
        i += 1
    if i >= len(s):
        return "", i
    if s[i] != "{":
        return s[i], i + 1
    depth, j = 0, i
    while j < len(s):
        if s[j] == "{":
            depth += 1
        elif s[j] == "}":
            depth -= 1
            if depth == 0:
                return s[i + 1:j], j + 1
        j += 1
    return s[i + 1:], len(s)


def latex_to_html(s):
    """把一段行内/行间 LaTeX 转成 Word 能显示的 HTML 片段。"""
    out = []
    i = 0
    n = len(s)
    while i < n:
        if s.startswith(r"\frac", i):
            num, i = _take_group(s, i + 5)
            den, i = _take_group(s, i)
            out.append("(" + latex_to_html(num) + ")/(" + latex_to_html(den) + ")")
        elif s.startswith(r"\sum", i):
            i += 4
            sub = sup = ""
            for _ in range(2):
                while i < n and s[i] == " ":
                    i += 1
                if i < n and s[i] == "_":
                    sub, i = _take_group(s, i + 1)
                elif i < n and s[i] == "^":
                    sup, i = _take_group(s, i + 1)
            frag = "\u03a3"
            if sub:
                frag += "<sub>" + latex_to_html(sub) + "</sub>"
            if sup:
                frag += "<sup>" + latex_to_html(sup) + "</sup>"
            out.append(frag)
        elif (s.startswith(r"\text", i) or s.startswith(r"\mathrm", i)
              or s.startswith(r"\mathbf", i) or s.startswith(r"\mathit", i)):
            for cmd in (r"\mathbf", r"\mathrm", r"\mathit", r"\text"):
                if s.startswith(cmd, i):
                    body, i = _take_group(s, i + len(cmd))
                    txt = latex_to_html(body)
                    out.append(f"<b>{txt}</b>" if cmd == r"\mathbf" else txt)
                    break
        elif s[i] == "\\":
            for sym in sorted(_SYMBOLS, key=len, reverse=True):
                if s.startswith(sym, i):
                    out.append(_SYMBOLS[sym])
                    i += len(sym)
                    break
            else:                      # 未知命令：去掉反斜杠保留名字
                j = i + 1
                while j < n and (s[j].isalpha() or s[j] == "\\"):
                    j += 1
                out.append(s[i + 1:j])
                i = j
        elif s[i] in "_^":
            tag = "sub" if s[i] == "_" else "sup"
            body, i = _take_group(s, i + 1)
            # ⚠ Word 不支持嵌套上/下标（<sub><sub>…</sub></sub> 会排错），
            #   故先把内层的 sub/sup 标签压平，只保留文字。
            inner = re.sub(r"</?su[bp]>", "", latex_to_html(body))
            out.append(f"<{tag}>{inner}</{tag}>")
        elif s[i] == "{":
            body, i = _take_group(s, i)
            out.append(latex_to_html(body))
        else:
            out.append(s[i])
            i += 1
    return "".join(out)


# ---------------------------------------------------------------------
#  行间公式：改用 XeLaTeX 渲成位图
#  ── 为什么不再用 latex_to_html 降级 ───────────────────────────────
#  降级路径把 \frac{a}{b} 写成 (a)/(b)，并把 \sum_{g \in S_s} 里嵌套的
#  下标压平（Word 不支持 <sub><sub>），实测渲染成
#      p(t,s) = (Σ g ∈ SsE(t,g))/(Σ gE(t,g))
#  ——分子分母混在一行、上下标挤成正文大小，读者无法辨认。报告的 5 个
#  行间公式全是分式，必须真排。行内公式（$t$、$\rho$）仍走 HTML 降级，
#  因为它们是单个符号，降级后与正文同基线，反而比贴图更协调。
# ---------------------------------------------------------------------
_TEX_BIN = None


def _tex_bin():
    """定位 xelatex：先查 PATH，再查本机固定安装位置。"""
    global _TEX_BIN
    if _TEX_BIN is None:
        cand = shutil.which("xelatex")
        if not cand:
            for p in (r"D:/Apps/Texlive/texlive/2026/bin/windows/xelatex.exe",
                      r"C:/texlive/2026/bin/windows/xelatex.exe"):
                if os.path.exists(p):
                    cand = p
                    break
        _TEX_BIN = cand or ""
    return _TEX_BIN


def render_math_png(tex_src: str):
    """把一段行间公式渲成 PNG，返回 (文件名, 宽 pt, 高 pt)；失败返回 None。

    按公式源码的 md5 做缓存：同一公式第二次导出直接复用，不重复编译。
    ⚠ 必须带 ctex —— 第 5 个公式的 \\text{} 里是中文（"TAI 跨类型标准差"），
      纯 amsmath 会把它排成空白。
    """
    exe = _tex_bin()
    if not exe:
        return None
    os.makedirs(MATH_DIR, exist_ok=True)
    key = hashlib.md5(tex_src.encode("utf-8")).hexdigest()[:12]
    png = os.path.join(MATH_DIR, key + ".png")
    pdf = os.path.join(MATH_DIR, key + ".pdf")
    if not os.path.exists(png):
        body = ("\\documentclass[border=3pt]{standalone}\n"
                "\\usepackage{amsmath,amssymb}\n"
                "\\usepackage{ctex}\n"
                "\\pagestyle{empty}\n"
                "\\begin{document}\n"
                "$\\displaystyle " + tex_src + "$\n"
                "\\end{document}\n")
        src = os.path.join(MATH_DIR, key + ".tex")
        with open(src, "w", encoding="utf-8") as f:
            f.write(body)
        try:
            r = subprocess.run([exe, "-interaction=nonstopmode", "-halt-on-error",
                                key + ".tex"],
                               cwd=MATH_DIR, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, timeout=180)
            if r.returncode != 0 or not os.path.exists(pdf):
                return None
            import pymupdf
            d = pymupdf.open(pdf)
            pix = d[0].get_pixmap(dpi=MATH_DPI, alpha=False,
                                  colorspace=pymupdf.csGRAY)
            pix.save(png)
            d.close()
        except Exception:
            return None
        finally:
            for ext in (".aux", ".log", ".out"):
                try:
                    os.remove(os.path.join(MATH_DIR, key + ext))
                except OSError:
                    pass
    try:
        from PIL import Image
        w_px, h_px = Image.open(png).size
    except Exception:
        return None
    return key + ".png", w_px * 72.0 / MATH_DPI, h_px * 72.0 / MATH_DPI


def _display_math(tex_src: str) -> str:
    """行间公式 → HTML 片段：优先贴 LaTeX 位图，TeX 不可用时回退文本降级。"""
    got = render_math_png(tex_src)
    if got is None:
        info("公式无法用 LaTeX 渲染，回退为文本降级（排版会变差）")
        return '<p class="math">%s</p>' % latex_to_html(tex_src)
    name, w_pt, h_pt = got
    # 版心宽 21 - 2×2.2 = 16.6 cm ≈ 470 pt，超宽按比例缩
    if w_pt > 470.0:
        h_pt *= 470.0 / w_pt
        w_pt = 470.0
    # HTML 的 width/height 单位是 px@96dpi，1 px = 0.75 pt
    return ('<p class="mathimg"><img src="_math/%s" width="%.0f" height="%.0f"'
            ' alt="公式"/></p>' % (name, w_pt / 0.75, h_pt / 0.75))


def convert_math(md_text):
    """把 Markdown 里的 $$...$$ 与 $...$ 换成 HTML（必须在 markdown 解析前做）。"""
    md_text = re.sub(
        r"^\$\$\s*\n(.+?)\n\s*\$\$$",
        lambda m: "\n" + _display_math(m.group(1).strip()) + "\n",
        md_text, flags=re.S | re.M)
    md_text = re.sub(
        r"(?<![\\$])\$(?!\$)([^$\n]+?)(?<!\\\$)\$(?!\$)",
        lambda m: latex_to_html(m.group(1)),
        md_text)
    return md_text


def _cover_meta_html(md_text: str) -> str:
    """把封面那三行元信息换成三个 <p class="meta"> 块（必须在 markdown 之前做）。

    ⚠ 为什么要在 markdown 之前处理：封面三行是**软换行**分隔的同一段，
    去掉 nl2br 之后它们会被并成一段。而一段里含强制换行的中文，
    在 Word 的两端对齐下会被**分散对齐**，渲染成「课 程 ： 生 物 信 息 学」
    这种逐字分离的灾难（实测旧版就是这样）。
    所以这里在 markdown 之前把它们拆成三个独立段落，
    并自己把 **粗体** 转成 <strong>（raw HTML 块不会被 markdown 再解析）。
    """
    m = re.search(r"\n\n(\*{2}[^\n]+\n\*{2}[^\n]+\n\*{2}[^\n]+)\n\n",
                  md_text)
    if not m:
        return md_text
    lines = [x.strip() for x in m.group(1).split("\n") if x.strip()]
    if not (2 <= len(lines) <= 6):
        return md_text
    out = []
    for ln in lines:
        ln = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", ln)
        out.append(f'<p class="meta">{ln}</p>')
    return md_text[:m.start(1)] + "\n\n".join(out) + md_text[m.end(1):]


def _structure_front(html: str) -> str:
    """把「主标题 + 元信息 + 摘要 + 关键词」整理成封面区，并在摘要后强制分页。

    ⚠ 为什么在 HTML 层面做、而不改 Markdown：
      报告 `.md` 是 07_终审核对.py 的唯一真源（239 项断言按行号/pattern 定位），
      为了排版去动它会把终审全部打乱。所以版式差异一律在渲染层解决。
    ⚠ markdown 的 toc 扩展给每个标题加了 id，匹配时必须容忍 `id="..."`。
    """
    html = re.sub(r'<h1 id="[^"]*">', '<h1 class="doctitle">', html, count=1)

    # 摘要标题
    html = re.sub(r'<h2 id="[^"]*">摘要</h2>',
                  '<h2 class="absh">摘　要</h2>', html, count=1)

    # 摘要正文：关键词段单独一类，其余按摘要体
    m = re.search(r'(<h2 class="absh">摘　要</h2>)([\s\S]*?)\s*<hr\s*/?>', html)
    if m:
        def _cls(mm):
            inner = mm.group(1)
            cls = "kw" if inner.startswith("<strong>关键词</strong>") else "abs"
            return f'<p class="{cls}">{inner}</p>'

        body = re.sub(r"<p>([\s\S]*?)</p>", _cls, m.group(2))
        # 摘要之后另起一页开始正文。分页由 _polish_doc 对「1 引言」设
        # PageBreakBefore 完成（空段落的 page-break-before 在 Word 里不生效）。
        html = html[:m.start()] + m.group(1) + body + html[m.end():]

    return html


def build_html(md_text: str) -> str:
    """Markdown → HTML，并把图片内嵌进去。"""
    # ⚠ 公式必须在 markdown 解析**之前**降级：否则 $$...$$ 会原样穿过管道，
    #    Word 把它当普通文本排进 PDF（修复前实测第 15 页直接印出了 \frac 源码）。
    md_text = convert_math(md_text)
    # ⚠ **不要开 nl2br**：报告 md 里行尾两空格出现 0 次，即**没有一处有意的
    #   硬换行**，所有换行都只是作者排版时的软换行。nl2br 会把每一个软换行
    #   都变成 <br />，于是每一段都被切成源文件里的行宽——参数极不均匀，
    #   左对齐时像没排版、右对齐时被逐字拉伸。去掉它，段落才会自然回流。
    #   封面那三行元信息改在 markdown **之前**直接生成 <p class="meta">。
    md_text = _cover_meta_html(md_text)
    html = markdown.markdown(
        md_text,
        extensions=["tables", "fenced_code", "codehilite", "toc",
                    "sane_lists"],
        extension_configs={
            "codehilite": {"guess_lang": False, "noclasses": True},
            "toc": {"permalink": False},
        },
    )

    # ---- 1) 内嵌图片，并把图与图注成对插入 ----
    # 位置策略（一次性算完再倒序插入，避免边插边算导致偏移错乱）：
    #   · 锚点 = 提到"图 N"的**最早**块级元素；标题优先（"顺带一提"不算讨论）
    #   · 同一锚点有多张图时：第一张紧跟锚点，其余的挪到**该小节末尾**
    #     （否则两张大图连排，第二张放不下就会在页尾留下大片空白）
    #   · 小节末尾 = 下一个同级或更高级标题之前
    def _plain(x):
        return re.sub(r"<[^>]+>", "", x)

    blocks = [(m.start(), m.end(), m.group(1), _plain(m.group(0)))
              for m in re.finditer(
                  r"<(p|li|h2|h3|h4|blockquote)[^>]*>.*?</\1>",
                  html, flags=re.S)]
    head_idx = [i for i, b in enumerate(blocks) if b[2] in ("h2", "h3", "h4")]
    lvl = {"h2": 2, "h3": 3, "h4": 4}

    def _section_end(bi):
        """小节末尾的位置：下一个同级或更高级标题的块首。"""
        cur = lvl.get(blocks[bi][2], 9)
        for j in head_idx:
            if j > bi and lvl.get(blocks[j][2], 9) <= cur:
                return blocks[j][0]
        return len(html)

    # ---- 先算出每张图的锚点与插入位置 ----
    plan = []                       # (pos, num, block_html)
    anchor_of = {}                  # num -> anchor block index
    for num, fname, cap in FIG_ORDER:
        path = os.path.join(FIG_DIR, fname)
        if not os.path.exists(path):
            info(f"跳过图 {num}：{fname} 不存在")
            continue
        cands = [i for i, b in enumerate(blocks)
                 if re.search(r"图\s*%d\b" % num, b[3])]
        if not cands:
            info(f"图 {num} 未在正文找到锚点，追加到文末")
            anchor_of[num] = None
            continue
        heads = [i for i in cands if blocks[i][2] in ("h2", "h3", "h4")]
        anchor_of[num] = (heads or cands)[0]

    # 同一锚点的图：第一张跟锚点，其余放小节末尾
    seen = {}
    for num, fname, cap in FIG_ORDER:
        if num not in anchor_of:
            continue
        ai = anchor_of[num]
        if ai is None:
            # 没找到锚点：追加到文末（并保持后面的偏移有效）
            pass
        else:
            k = seen.get(ai, 0)
            seen[ai] = k + 1
        path = os.path.join(FIG_DIR, fname)
        rel = os.path.relpath(path, os.path.dirname(OUT_HTML)).replace("\\", "/")
        # ⚠ 必须给**固定宽度**（HTML width 属性，单位 px @96dpi）：
        #   Word 的 HTML 导入器不认 CSS 的 max-width，超宽图会被直接裁掉
        #   （实测图 10 右半幅被切）。正文宽 21 - 2*2.2 = 16.6 cm ≈ 627 px，
        #   取 620 px 留一点余量。
        # ⚠ 光限宽不够：**竖长图必须同时限高**。fig5（1282×1937）在 620 px
        #   宽时高 937 px = 703 pt，正好等于版心高度（29.7-2.5-2.4=24.8 cm），
        #   Word 会把它连同图注一起挤出页面、下半截直接消失（实测第 7 页只剩
        #   半张图）。高度上限取 780 px = 585 pt，给图注和上下间距留余量。
        try:
            from PIL import Image
            iw, ih = Image.open(path).size
        except Exception:
            iw = ih = 0
        if iw > 0:
            w_px = 620.0
            h_px = w_px * ih / iw
            MAX_H_PX = 780.0
            if h_px > MAX_H_PX:
                w_px = MAX_H_PX * iw / ih
                h_px = MAX_H_PX
            size_attr = ' width="%.0f" height="%.0f"' % (w_px, h_px)
        else:
            size_attr = ' width="620"'
        # ⚠ 不要用 <figure>/<figcaption>：Word 的 HTML 导入器不认这两个标签，
        #    会把图注拆成两个段落（实测图 1 的"图"被单独挤到一行、其余文字另起一段）。
        #    改用两个带内联样式的 <p>，Word 对 <p> 的解析一直可靠。
        block = (f'<p class="figimg"><img src="{rel}"{size_attr} alt="图 {num}"/></p>'
                 f'<p class="figcap">{cap}</p>')
        if ai is None:
            pos = len(html)
        elif k == 0:
            pos = blocks[ai][1]
        else:
            pos = _section_end(ai)
        plan.append((pos, num, block))

    # 倒序插入：位置从大到小、同一位置按图号从大到小，
    # 这样最终同一位置上的图是**编号递增**的（3.2 的图 2、图 3 不会颠倒）。
    for pos, num, block in sorted(plan, key=lambda x: (-x[0], -x[1])):
        html = html[:pos] + block + html[pos:]

    # ---- 2) 表格里的竖直分隔符会被 markdown 吞掉吗？不，tables 扩展保留 ----
    # ---- 3) 封面 / 摘要区结构化 ----
    # ⚠ 必须在 _structure_front **之后**再删 <hr>：它靠 <hr /> 定位元信息段与
    #   摘要段的边界，先删就找不到锚点了。
    # ⚠ 为什么要全删：md 里的 `---` 是章节分隔线，Word 会把它排成一个带下边框
    #   的**空段落**。当它正好落在分页处（实测第 14 页，`---` 紧跟 `## 3 结果`），
    #   页面顶端就会出现「页眉线 + 一条多余横线」的双线，非常像排版事故。
    #   h2 本身已有下边框，够用。
    html = _structure_front(html)
    html = re.sub(r"\s*<hr\s*/?>\s*", "\n", html)

    # ---- 4) 页面标题 ----
    return f"""<!DOCTYPE html>
<html lang="zh-CN"><head><meta charset="utf-8"/>
<title>细胞地层 CellStrata 实验报告</title>
<style>{CSS}</style></head><body>
{html}
</body></html>"""


def _text_w(txt: str, size: float = 9.0) -> float:
    """估一段文字的排版宽度（pt）。中文按 1 个字宽，西文按 0.5 个字宽。"""
    return sum(size if ord(ch) > 127 else size * 0.5 for ch in txt)


def _table_col_pt(t: object, avail_pt: float = 470.0, size: float = 9.0):
    """算一张表各列的宽度（**磅值**，不是百分比）。

    为什么不用「按最长单元格宽度开方」那种简单权重：
      简单权重让长列吃满、短列只剩 2~3 个汉字宽。9 列表里
      「p（富集）」这种表头会被拆成「p（富 / 集）」，而「1.1×10⁻³」会在
      「10⁻」与「³」之间断开——两者都是读者一眼能看出的排版事故。
    这里的做法：
      1) 每列先算**不可断单元**的最大宽度（连续西文算一个整体，中文逐字可断），
         作为该列的硬下限——保证「p（富集）」不被拆开、数字不被腰斩；
      2) 再按「整格文字宽度」作为理想值；
      3) 先满足所有下限，剩余宽度按（理想 − 下限）的比例分配。
    内边距同时压到 3 pt（默认 5.4 pt）：9 列表里能省出近 50 pt 可用宽度。
    """
    ncol = t.Columns.Count
    pad = 6.0 if ncol >= 7 else 8.0
    try:
        t.LeftPadding = pad / 2
        t.RightPadding = pad / 2
    except Exception:
        pass

    full = [0.0] * ncol
    tok = [0.0] * ncol
    for row in t.Rows:
        cells = row.Cells
        for ci in range(min(ncol, cells.Count)):
            txt = cells(ci + 1).Range.Text
            txt = txt.replace("\r", "").replace("\x07", "").strip()
            if not txt:
                continue
            full[ci] = max(full[ci], _text_w(txt, size))
            for unit in txt.split():
                for run in re.findall(r"[\x21-\x7e]+|.", unit):
                    tok[ci] = max(tok[ci], _text_w(run, size))

    desired = [min(f + pad, 260.0) for f in full]
    mins = [min(max(k + pad, 26.0), 96.0) for k in tok]

    alloc = list(mins)
    rest = avail_pt - sum(alloc)
    if rest > 0:
        wsum = sum(max(d - m, 0.5) for d, m in zip(desired, mins)) or 1.0
        alloc = [m + rest * max(d - m, 0.5) / wsum
                 for m, d in zip(mins, desired)]
    else:                                  # 连下限都放不下：按比例压缩
        s = sum(alloc) or 1.0
        alloc = [a / s * avail_pt for a in alloc]
    total = sum(alloc) or 1.0
    return [a / total * avail_pt for a in alloc]


def _polish_doc(doc):
    """Word 侧排版精修：处理 HTML/CSS 管不到的部分。

    必须在这里做、而不能只靠 CSS 的四件事：
      1) 行距 —— CSS 的 `line-height: 1.55` 会被 Word 当成「固定 12 pt」，
         只有 LineSpacingRule=wdLineSpaceMultiple 才是真正的多倍行距；
      2) 中西文字体分工 —— 靠 font-family 列表不可靠，要分设
         Font.NameFarEast（宋体）与 Font.NameAscii（Times New Roman）；
      3) 三线表 —— border-collapse + :last-child 在 Word 里不生效，边框得逐个设；
      4) 页眉页脚 —— 属于节属性，CSS 的 @page 表达不了。
    """
    ALIGN_LEFT, ALIGN_CENTER, ALIGN_JUSTIFY = 0, 1, 3
    LS_SINGLE, LS_MULTIPLE = 0, 5
    BORDER_TOP, BORDER_LEFT, BORDER_BOTTOM, BORDER_RIGHT = -1, -2, -3, -4
    BORDER_H, BORDER_V = -5, -6
    LINE_SINGLE = 1
    W_050, W_075, W_100, W_150 = 4, 6, 8, 12   # 磅宽枚举
    # CSS 里 blockquote 与 pre 共用的底色 #f7f8fa；Word 用 BGR 存 → 16250356
    NOTE_SHADE = 16250356

    def _no_indent(fmt):
        # ⚠ 两个属性都必须清，缺一不可。
        #   CharacterUnitFirstLineIndent（字符单位）与 FirstLineIndent（磅）在
        #   Word 里**不是**同一个底层值：把字符单位设成 0 之后，磅值仍可能停
        #   在 18 pt（实测表格单元格 CFirst=0.00 而 First=18.00）。
        #   只清前者 → 单元格首行仍缩进 2 个字符、首行可用宽度只剩 1 个汉字，
        #   于是「组织」被拆成「组」/「织」两行、「层面」竖排——这正是表格
        #   看起来"列宽全乱"的真正原因（列宽其实是对的）。
        try:
            fmt.CharacterUnitFirstLineIndent = 0
        except Exception:
            pass
        try:
            fmt.FirstLineIndent = 0
        except Exception:
            pass

    def _first_line_2(fmt):
        try:
            fmt.CharacterUnitFirstLineIndent = 2
        except Exception:
            fmt.FirstLineIndent = 21.0

    def _scale(fmt, multiple):
        fmt.LineSpacingRule = LS_MULTIPLE
        fmt.LineSpacing = multiple * 12.0

    # ---------- 1. 页面与页眉页脚 ----------
    CM = 28.3465
    ps = doc.PageSetup
    ps.PageWidth, ps.PageHeight = 21.0 * CM, 29.7 * CM
    ps.TopMargin, ps.BottomMargin = 2.5 * CM, 2.4 * CM
    ps.LeftMargin, ps.RightMargin = 2.2 * CM, 2.2 * CM
    ps.HeaderDistance, ps.FooterDistance = 1.4 * CM, 1.5 * CM
    ps.DifferentFirstPageHeaderFooter = True    # 封面（第 1 页）不显示页眉页码

    for sec in doc.Sections:
        hdr = sec.Headers(1)
        hdr.Range.Text = "细胞地层 CellStrata　·　人体细胞类型的基因年龄构成"
        hp = hdr.Range.ParagraphFormat
        hp.Alignment = ALIGN_CENTER
        hp.SpaceAfter = 0
        hp.LineSpacingRule = LS_SINGLE
        hp.Borders(BORDER_BOTTOM).LineStyle = LINE_SINGLE
        hp.Borders(BORDER_BOTTOM).LineWidth = W_050
        hdr.Range.Font.Size = 9
        hdr.Range.Font.NameFarEast = "宋体"
        hdr.Range.Font.NameAscii = "Times New Roman"

        ftr = sec.Footers(1)
        ftr.Range.Text = ""
        fp = ftr.Range.ParagraphFormat
        fp.Alignment = ALIGN_CENTER
        fp.LineSpacingRule = LS_SINGLE
        ftr.Range.Font.Size = 9
        ftr.Range.Font.NameFarEast = "宋体"
        ftr.Range.Font.NameAscii = "Times New Roman"

        def _append(text=None, field=None):
            r = ftr.Range
            r.Collapse(0)                      # 折到末尾
            if text is not None:
                r.InsertAfter(text)
            if field is not None:
                r2 = ftr.Range
                r2.Collapse(0)
                r2.Fields.Add(r2, field)       # 33=页码, 26=总页数

        _append("第 ")
        _append(field=33)
        _append(" 页　共 ")
        _append(field=26)
        _append(" 页")

    # ---------- 2. 全局正文格式 ----------
    C = doc.Content
    C.Font.NameFarEast = "宋体"
    C.Font.NameAscii = "Times New Roman"
    C.Font.NameOther = "Times New Roman"
    C.Font.Size = 10.5
    pf = C.ParagraphFormat
    _scale(pf, 1.55)
    pf.Alignment = ALIGN_JUSTIFY
    pf.SpaceBefore, pf.SpaceAfter = 0, 0
    pf.LeftIndent, pf.RightIndent = 0, 0
    _first_line_2(pf)

    # ---------- 3. 逐段按样式覆盖 ----------
    # 标题字号/间距在这张表里统一维护，改版式只动这里。
    HEAD = {
        "标题 1": dict(size=21, align=ALIGN_CENTER, before=0,   after=22, scale=1.40),
        "标题 2": dict(size=15, align=ALIGN_LEFT,   before=18,  after=9,  scale=1.30),
        "标题 3": dict(size=12.5, align=ALIGN_LEFT, before=15,  after=6,  scale=1.35),
        "标题 4": dict(size=11, align=ALIGN_LEFT,   before=11,  after=5,  scale=1.40),
        "标题 5": dict(size=10.5, align=ALIGN_LEFT, before=9,   after=4,  scale=1.45),
    }

    for p in doc.Paragraphs:
        try:
            st = p.Style.NameLocal
        except Exception:
            continue
        fmt, f = p.Format, p.Range.Font

        if st in HEAD:
            cfg = HEAD[st]
            f.NameFarEast, f.NameAscii = "黑体", "Times New Roman"
            f.NameOther = "Times New Roman"
            f.Bold = True
            f.Size = cfg["size"]
            _no_indent(fmt)
            fmt.LeftIndent = 0
            fmt.Alignment = cfg["align"]
            fmt.SpaceBefore, fmt.SpaceAfter = cfg["before"], cfg["after"]
            _scale(fmt, cfg["scale"])
            fmt.KeepWithNext = True          # 标题不落在页脚处单独成行
            txt = p.Range.Text.strip()
            # 「摘　要」这一级标题居中、加字距，与正文章节标题区分
            if txt == "摘　要":
                fmt.Alignment = ALIGN_CENTER
                f.Size = 13
                fmt.SpaceBefore = 26
                fmt.SpaceAfter = 12
            # ⚠ 这里**不设** PageBreakBefore。本报告的摘要长约 1.3 页，
            #   一旦强制分页，摘要尾部会把第 2 页占掉一半、下半页全空，
            #   紧接着又是一次空页（实测 56 页版第 2 页空白 56%）。
            #   改为让正文紧接摘要排下去，整页不留空。
            continue

        if st == "mathimg":                   # 行间公式（LaTeX 位图）
            _no_indent(fmt)
            fmt.Alignment = ALIGN_CENTER
            fmt.SpaceBefore, fmt.SpaceAfter = 10, 10
            fmt.LineSpacingRule = LS_SINGLE
            fmt.KeepWithNext = False
            continue

        if st == "math":                      # 行间公式（TeX 不可用时的文本降级）
            _no_indent(fmt)
            fmt.Alignment = ALIGN_CENTER
            fmt.SpaceBefore, fmt.SpaceAfter = 8, 8
            _scale(fmt, 1.3)
            continue

        if st == "figcap":                    # 图注
            f.NameFarEast, f.NameAscii = "宋体", "Times New Roman"
            f.Size = 9
            f.Bold = False
            _no_indent(fmt)
            fmt.Alignment = ALIGN_CENTER
            fmt.SpaceBefore, fmt.SpaceAfter = 4, 14
            _scale(fmt, 1.35)
            continue

        if st == "figimg":                    # 图本体
            _no_indent(fmt)
            fmt.Alignment = ALIGN_CENTER
            fmt.SpaceBefore, fmt.SpaceAfter = 10, 0
            fmt.LineSpacingRule = LS_SINGLE
            fmt.KeepWithNext = True           # 图与图注必须同页
            continue

        if st in ("HTML 预设格式", "HTML Preformatted"):   # 代码块
            f.NameFarEast, f.NameAscii = "宋体", "Consolas"
            f.Size = 8.5
            _no_indent(fmt)
            fmt.Alignment = ALIGN_LEFT
            fmt.SpaceBefore, fmt.SpaceAfter = 0, 0
            fmt.LineSpacingRule = LS_SINGLE
            continue

        # 注记块（引用块）：Word **没有**把 CSS 的 border-left 导入成段落边框
        #   （实测边框枚举读回来是 0），但**底纹**导入了（#f7f8fa → 16250356）。
        #   所以认注记块只能靠底纹，不能靠左边框——认错的后果是整块文字被
        #   两端对齐拉伸出满行字距（实测页 11 的引用块字距被撑到 1.5 倍）。
        # ⚠ 列表项必须**先于**注记块判定：引用块里嵌着 <ul>，它的 <li>
        #   同时带底纹和列表标记，两条规则都会命中。
        try:
            in_list = p.Range.ListFormat.ListType != 0
        except Exception:
            in_list = False
        try:
            shade = fmt.Shading.BackgroundPatternColor
        except Exception:
            shade = -16777216
        try:
            is_note = (shade == NOTE_SHADE) or (fmt.Borders(BORDER_LEFT).LineStyle != 0)
        except Exception:
            is_note = False

        if in_list:
            # ⚠ 列表项要的是**悬挂缩进**（负的首行缩进），不是 0：
            #   把首行缩进清成 0 后，项目符号会与正文第一个字落在同一位置，
            #   Word 判定重叠，就把整段正文挤到下一行——实测页 11 的
            #   「• Spearman 相关是同义反复」圆点独占一行、文字掉到第二行。
            #   正确做法：左缩进 2 字符、首行缩进 −2 字符。
            try:
                fmt.CharacterUnitFirstLineIndent = -2
            except Exception:
                pass
            fmt.LeftIndent = 21.0
            fmt.FirstLineIndent = -21.0
            fmt.Alignment = ALIGN_LEFT
            fmt.SpaceBefore, fmt.SpaceAfter = 0, 2
            _scale(fmt, 1.45)
            if is_note:
                f.Size = 9.8
            continue

        if is_note:
            _no_indent(fmt)
            f.Size = 9.8
            fmt.Alignment = ALIGN_LEFT
            fmt.SpaceBefore, fmt.SpaceAfter = 2, 2
            _scale(fmt, 1.45)
            continue

    # ---------- 4. 表格：三线表 ----------
    for t in doc.Tables:
        for i in (BORDER_TOP, BORDER_LEFT, BORDER_BOTTOM, BORDER_RIGHT, BORDER_H, BORDER_V):
            try:
                t.Borders(i).LineStyle = 0
            except Exception:
                pass

        # 宽表（≥8 列）降到 8 pt：10 列挤在 16.6 cm 里，9 pt 时每列只有
        # 约 4 个汉字宽，「1.1×10⁻³」这种数字会被断成两行。8 pt 是表格式
        # 排版里常见的下限，正文仍是 10.5 pt，反差可接受。
        ncol0 = t.Columns.Count
        tsize = 8 if ncol0 >= 8 else 9

        r = t.Range
        r.Font.NameFarEast = "宋体"
        r.Font.NameAscii = "Times New Roman"
        r.Font.Size = tsize
        tp = r.ParagraphFormat
        _no_indent(tp)
        tp.Alignment = ALIGN_LEFT
        tp.SpaceBefore, tp.SpaceAfter = 0, 0
        _scale(tp, 1.25)
        tp.LeftIndent, tp.RightIndent = 0, 0

        # 顶线 / 底线
        t.Borders(BORDER_TOP).LineStyle = LINE_SINGLE
        t.Borders(BORDER_TOP).LineWidth = W_150
        t.Borders(BORDER_BOTTOM).LineStyle = LINE_SINGLE
        t.Borders(BORDER_BOTTOM).LineWidth = W_150

        # 表头行：黑体加粗、下细线、跨页重复
        try:
            head_row = t.Rows(1)
            head_row.HeadingFormat = True
            head_row.Range.Font.Bold = True
            head_row.Range.Font.NameFarEast = "黑体"
            for c in head_row.Cells:
                c.Borders(BORDER_BOTTOM).LineStyle = LINE_SINGLE
                c.Borders(BORDER_BOTTOM).LineWidth = W_075
        except Exception:
            pass

        # 行内不断页：否则一个单元格会被劈成上下两半
        try:
            t.Rows.AllowBreakAcrossPages = False
        except Exception:
            pass

        # ---------- 列宽 ----------
        # ⚠ 必须自己算：Word 的自动列宽会**平均分配**，把「层面」「p（富集）」
        #   这种短表头压成竖直排列。算法与踩坑说明见 _table_col_pt()。
        # ⚠ 只能用 Column.Width（磅值），**不能**用 PreferredWidth：
        #   对 Table 设 PreferredWidthType 的合法值是 3（磅）/2（百分比）/1（自动），
        #   传 0 会直接抛「数值超出范围」；而设成百分比时 Word 又把表宽改回自动
        #   (9999999)、按自己的规则重排——实测算出 8.2% 的列被排成 5.1%。
        #   一次性把每列 Width 都设到目标值，Word 会保持总宽 470 pt 不变。
        try:
            ncol = t.Columns.Count
            alloc = _table_col_pt(t, avail_pt=470.0, size=tsize)
            try:
                t.AllowAutoFit = False
            except Exception:
                pass
            t.AutoFitBehavior(0)          # wdAutoFitFixed
            for ci in range(ncol):
                t.Columns(ci + 1).Width = alloc[ci]
        except Exception as e:
            info(f"某一表格的列宽未能调整（不影响导出）：{e}")

    # ---------- 5. 表后留白 ----------
    # ⚠ Word 的表格没有「段后距」概念，表格底部的收口线会紧贴下一段正文
    #   （实测 3.3 节的表格与"本研究的差异化在于三点"只隔约 4 pt）。
    #   表格后面那一段是**独立段落**，单独给它加 SpaceBefore 即可。
    for t in doc.Tables:
        try:
            tail = doc.Range(t.Range.End, t.Range.End)
            if tail.Tables.Count == 0:
                tail.ParagraphFormat.SpaceBefore = 10
        except Exception:
            continue


def html_to_pdf(html_path: str, pdf_path: str):
    """用 Word COM 把 HTML 转成 PDF。"""
    import win32com.client as win32

    html_path = os.path.abspath(html_path)
    pdf_path = os.path.abspath(pdf_path)
    if os.path.exists(pdf_path):
        os.remove(pdf_path)

    word = win32.DispatchEx("Word.Application")
    word.Visible = False
    word.DisplayAlerts = 0
    doc = None
    try:
        doc = word.Documents.Open(html_path, ConfirmConversions=False,
                                  ReadOnly=False)
        # 排版精修：页面 / 页眉页脚 / 字体分工 / 行距 / 三线表，全在 _polish_doc 里
        try:
            _polish_doc(doc)
        except Exception as e:
            info(f"排版精修部分失败（PDF 仍会导出，但版式可能退化）：{e}")

        doc.ExportAsFixedFormat(pdf_path, 17)  # wdExportFormatPDF
        pages = doc.ComputeStatistics(2)  # wdStatisticPages
        words = doc.ComputeStatistics(0)
        doc.Close(SaveChanges=0)
        doc = None
        return pages, words
    finally:
        try:
            if doc is not None:
                doc.Close(SaveChanges=0)
        except Exception:
            pass
        word.Quit()


def main():
    if not os.path.exists(REPORT_MD):
        die(f"找不到报告源文件：{REPORT_MD}")

    with open(REPORT_MD, encoding="utf-8") as f:
        md_text = f.read()

    info(f"读取报告：{len(md_text):,} 字符")
    html = build_html(md_text)
    with open(OUT_HTML, "w", encoding="utf-8") as f:
        f.write(html)
    info(f"生成 HTML：{os.path.getsize(OUT_HTML)/1024:.0f} KB")

    try:
        pages, words = html_to_pdf(OUT_HTML, OUT_PDF)
        ok(f"{OUT_PDF}（{pages} 页 / 约 {words:,} 字）")
    except ImportError:
        print("\n✗ 缺少 pywin32，无法调用 Word COM。", file=sys.stderr)
        print("  备选：直接用浏览器打印 HTML（Ctrl+P → 另存为 PDF）。", file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()
