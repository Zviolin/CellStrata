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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import ROOT, info, ok, die  # noqa: E402

import markdown  # noqa: E402

REPORT_MD = os.path.join(ROOT, "40-报告", "01-实验报告.md")
FIG_DIR = os.path.join(ROOT, "30-可视化", "figures")
OUT_HTML = os.path.join(ROOT, "40-报告", "_report_preview.html")
OUT_PDF = os.path.join(ROOT, "40-报告", "01-实验报告.pdf")

# ---------- 图片：按正文出现的"（**图 N**）"锚点插入 ----------
# 报告正文用「（**图 1**、**图 3**）」这样的形式引用图。
# 导出 PDF 时把对应图片**内嵌**到该处，让 PDF 既是文字报告也是图册。
FIG_ORDER = [
    (1, "fig1_strata.png", "图 1  细胞地层主图：177 种细胞类型按 TAI 排序，每根按 16 个年龄层堆叠"),
    (2, "fig2_pca.png", "图 2  PCA 前两主成分散点，按 K-means 簇着色"),
    (3, "fig3_heatmap.png", "图 3  细胞类型 × 年龄层的表达占比热图"),
    (4, "fig4_importance_age.png", "图 4  各年龄层基因数与 Top 500 重要基因的富集倍数"),
    (5, "fig5_tai_rank.png", "图 5  177 种细胞类型的 TAI 排序与双口径误差棒"),
    (7, "fig7_pca_null_baseline.png", "图 7  PCA 零假设诊断：三条基线对比"),
    (8, "fig8_robustness_layers.png", "图 8  三重稳健性逐层对照"),
    (9, "fig9_shortprot_premise.png", "图 9  短蛋白前提检验"),
    (10, "fig10_longprot_comparison.png", "图 10  长蛋白版 TAI 绝对值对照（检验 D）"),
]

CSS = """
@page { size: A4; margin: 2.2cm 2.0cm 2.2cm 2.0cm; }
body {
  font-family: "Microsoft YaHei", "SimSun", "PingFang SC", serif;
  font-size: 10.5pt; line-height: 1.65; color: #1a1a1a;
}
h1 { font-size: 20pt; text-align: center; margin: 1.2em 0 0.6em;
     border-bottom: 2.5px solid #333; padding-bottom: 0.35em; }
h2 { font-size: 15pt; margin: 1.6em 0 0.7em; padding-left: 0.42em;
     border-left: 5px solid #2b5c9a; line-height: 1.35; }
h3 { font-size: 12.5pt; margin: 1.2em 0 0.5em; color: #1f3f6b; }
h4 { font-size: 11pt; margin: 1em 0 0.4em; color: #333; }
p  { margin: 0.5em 0; text-align: justify; }
table {
  border-collapse: collapse; width: 100%; margin: 0.9em 0;
  font-size: 9pt; page-break-inside: avoid;
}
th, td { border: 1px solid #999; padding: 4px 7px; text-align: left;
         vertical-align: top; }
th { background: #eef2f8; font-weight: bold; }
tr:nth-child(even) td { background: #fafbfd; }
code {
  font-family: "Consolas", "Cascadia Mono", monospace; font-size: 9pt;
  background: #f2f4f7; padding: 1px 4px; border-radius: 3px;
}
pre {
  background: #f5f6f8; border: 1px solid #d5d9e0; border-left: 4px solid #2b5c9a;
  padding: 9px 12px; overflow-x: auto; font-size: 8.8pt; line-height: 1.45;
  page-break-inside: avoid;
}
pre code { background: none; padding: 0; font-size: 8.8pt; }
blockquote {
  margin: 0.9em 0; padding: 0.7em 1em; background: #fff8e6;
  border-left: 4px solid #d9a441; color: #4a3720;
}
blockquote p { margin: 0.35em 0; }
ul, ol { margin: 0.5em 0 0.7em 1.8em; }
li { margin: 0.22em 0; }
hr { border: none; border-top: 1px solid #ccc; margin: 1.5em 0; }
figure { margin: 1.3em 0; text-align: center; page-break-inside: avoid; }
figure img { max-width: 100%; height: auto; border: 1px solid #ddd; }
figcaption { font-size: 9pt; color: #555; margin-top: 0.45em; text-align: center;
             font-style: normal; }
strong { font-weight: bold; }
a { color: #1a4d8f; text-decoration: none; }
"""


def build_html(md_text: str) -> str:
    """Markdown → HTML，并把图片内嵌进去。"""
    html = markdown.markdown(
        md_text,
        extensions=["tables", "fenced_code", "codehilite", "toc",
                    "nl2br", "sane_lists"],
        extension_configs={
            "codehilite": {"guess_lang": False, "noclasses": True},
            "toc": {"permalink": False},
        },
    )

    # ---- 1) 内嵌图片（base64），并把图片段落包进 <figure> ----
    for num, fname, cap in FIG_ORDER:
        path = os.path.join(FIG_DIR, fname)
        if not os.path.exists(path):
            info(f"跳过图 {num}：{fname} 不存在")
            continue
        # ⚠ 用**本地文件路径**而不是 base64 data URI：
        #   Word 打开本地 HTML 时能正确解析相对路径的图片，
        #   而 data URI 会被 Word 统一挪到文档末尾（实测 9 张图全挤在最后 9 页）。
        #   代价是 HTML 不能脱离本目录独立分发，但它只是中间产物，PDF 才是交付物。
        rel = os.path.relpath(path, os.path.dirname(OUT_HTML)).replace("\\", "/")
        # ⚠ 必须给**固定宽度**（HTML width 属性，单位 px @96dpi）：
        #   Word 的 HTML 导入器不认 CSS 的 max-width，超宽图会被直接裁掉
        #   （实测图 10 右半幅被切）。正文宽 21 - 2*2.0 = 17 cm ≈ 643 px，
        #   取 620 px 留一点余量。
        block = (f'<figure><img src="{rel}" width="620" '
                 f'alt="图 {num}"/><figcaption>{cap}</figcaption></figure>')
        # 定位策略：找**最后一个**提到"图 N"的**块级元素**（p / li / h2 / h3 / h4），
        # 把图插在它之后。
        # ⚠ 不能用简单的 "<p>...图 N...</p>" 正则：markdown 转换后正文里是
        #    「（<strong>图 1</strong>、<strong>图 3</strong>）」，段落内部含标签，
        #    [^<]* 匹配不到 → 所有图都会被追加到文末（实测 9 张图全挤在最后 9 页）。
        def _plain(x):
            return re.sub(r"<[^>]+>", "", x)

        # 找所有块级元素的完整范围
        blocks = [(m.start(), m.end(), _plain(m.group(0)))
                  for m in re.finditer(
                      r"<(p|li|h2|h3|h4|blockquote)[^>]*>.*?</\1>",
                      html, flags=re.S)]
        # ⚠ 必须取**最早**的候选，不能取最后一个：
        #   附录 C 的图表清单与附录 A 的命令块里也提到"图 N"，
        #   取最后 → 所有图都被插到附录之后（实测 9 张图全挤在最后 9 页）。
        cands = [b for b in blocks if re.search(r"图\s*%d\b" % num, b[2])]
        if cands:
            pos = cands[0][1]
            # 同一锚点下多张图（如 3.1 同时讲图 1 与图 5）要保持编号顺序：
            # 把新图插到**已插入的同锚点图之前**。
            # ⚠ 只在**同一个锚点之后**才调整顺序：必须校验 nxt > pos。
            #    否则当锚点位置比已插入的 figure 更靠后时，html[pos:nxt] 是空串，
            #    len < 400 成立 → pos 被错误地前移，图 10 被插到图 1 之前。
            tag = '<figure><img src="'
            nxt = html.find(tag, pos)
            while nxt != -1:
                between = html[pos:nxt]
                if len(between) < 400 and "figure" not in between:
                    pos = nxt
                    break
                nxt = html.find(tag, nxt + 1)
            html = html[:pos] + block + html[pos:]
        else:
            html += block
            info(f"图 {num} 未在正文找到锚点，追加到文末")

    # ---- 2) 表格里的竖直分隔符会被 markdown 吞掉吗？不，tables 扩展保留 ----
    # ---- 3) 页面标题 ----
    return f"""<!DOCTYPE html>
<html lang="zh-CN"><head><meta charset="utf-8"/>
<title>细胞地层 CellStrata 实验报告</title>
<style>{CSS}</style></head><body>
{html}
</body></html>"""


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
        # 页面设置：A4 + 合适页边距
        # ⚠ 不要用 word.CentimetersToPoints()：在 Documents.Open 之后的
        #   某些上下文里它会抛 com_error(-2147467259)。
        #   直接用**磅值**（1 cm = 28.3465 pt）更稳，也少一层 COM 调用。
        CM = 28.3465
        ps = doc.PageSetup
        ps.PageWidth = 21.0 * CM
        ps.PageHeight = 29.7 * CM
        ps.TopMargin = 2.2 * CM
        ps.BottomMargin = 2.2 * CM
        ps.LeftMargin = 2.0 * CM
        ps.RightMargin = 2.0 * CM

        # 页脚加页码
        try:
            for sec in doc.Sections:
                footer = sec.Footers(1)  # wdHeaderFooterPrimary
                footer.Range.Fields.Add(footer.Range, 33)  # wdFieldPage
                footer.Range.ParagraphFormat.Alignment = 1  # 居中
                footer.Range.Font.Size = 9
        except Exception as e:
            info(f"页码添加失败（不影响 PDF）：{e}")

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
