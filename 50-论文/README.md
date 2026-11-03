# 50-论文 · LaTeX 小论文

本目录是《生物信息学》课程项目的**独立小论文**版本：用 LaTeX 把项目的完整方法与结果
写成一篇符合**中文科技期刊**格式的论文，并与实验报告**逐位一致**（由 `95_论文数字终审.py` 自动保证）。

## 版式：中文期刊模板 cjc

论文使用 **[`latexstudio/chinesejournal`](https://github.com/latexstudio/chinesejournal)** 的
`cjc.cls`（LPPL-1.3c，作者 Zeping Lee）。这是以《计算机学报》要求编写的中文期刊模板框架，
多数中文期刊模板都由它派生。模板底层是 `ctexart + twocolumn + zihao=5`，即中文科技期刊的标准版式：

| 要素 | 规格 |
|---|---|
| 题目 | 小二黑体，通栏居中 |
| 作者 / 单位 | 三号 / 小五（带括号） |
| 摘要 / 关键词 | 小五宋体，"摘　要："黑体起头 |
| 分类号 | 中图法分类号 + 文献标识码，通栏 |
| 正文 | 五号宋体，**双栏**，栏距 0.75 cm |
| 一级标题 | 四号黑体，二级五号黑体，三级五号宋体 |
| 表格 | 三线表，表题在上（黑体），图题在下 |
| 页眉 | 刊名 + 短题 + 页码，首页为卷期样式 |
| 参考文献 | `cjc.bst`（GB/T 7714 风格），正文用**右上角标** |

模板原件（**未打补丁**的 `cjc.cls`、`cjc.bst`、`cjc.dtx` 与官方示例）归档在 `journal-template/`；
本目录下的 `cjc.cls` 是打过补丁的工作副本。

> ⚠ **模板必须打一个补丁**：`cjc.cls` 用 `\csname str_if_eq_x:nnTF\endcsname` 做字符串比较，
> 而 `str_if_eq_x:nnTF` 在 TeX Live 2021 之后已被 expl3 移除。缺失时 `\csname` 会把它变成 `\relax`，
> 于是**每一个判断分支都会被执行**（mac 的 `Songti SC`、`Menlo` 也会被加载），
> 在本机直接报 `fontspec Error` 并在 l.52 处 `Missing \begin{document}`。
> **官方示例在 TeX Live 2026 上同样是坏的**，不是本项目的用法问题。
> 补丁：`cjc.cls` 中 `str_if_eq_x:nnTF` → `str_if_eq:eeTF`（未改动的原件留存于 `journal-template/cjc.cls`）。

## 文件

| 文件 | 说明 |
|---|---|
| `main.tex` | 论文源文件（中文正文，cjc 文档类 + XeLaTeX） |
| `cjc.cls` / `cjc.bst` | 中文期刊模板的类文件与参考文献样式（编译必需） |
| `references.bib` | 参考文献库（21 条；正文实际引用 12 条，未引条目不进 PDF） |
| `main.pdf` | 编译产物（13 页，A4，中文期刊双栏版式） |
| `journal-template/` | 模板原始归档（含官方示例 `cjc-example.tex`） |
| `main_oup_1col_backup.tex` | 上一版：OUP 官方模板单栏版 |
| `main_ctexart_backup.tex` | 更早一版：ctexart 自建期刊体例 |

## 编译

```bash
cd 50-论文
xelatex main.tex      # 第 1 遍：生成 aux
bibtex  main           # 生成 main.bbl（参考文献）
xelatex main.tex      # 第 2 遍：插入引用与页码
xelatex main.tex      # 第 3 遍：交叉引用收敛
```

必须用 **XeLaTeX**：`cjc.cls` 里有 `\RequireXeTeX`，且用 `unicode-math` + XITSMath 排数学。
需要 TeX Live 完整版（含 ctex、xeCJK、unicode-math、natbib、booktabs、tabularx）。

## 图片路径

论文用 `\graphicspath{{../30-可视化/figures/}}` 引用十张图，
因此**必须从本目录编译**（保持目录结构），否则找不到图。

## 图表放大的处理

双栏版式下单栏宽只有 8.1 cm，多面板图会缩成看不清的小块。本项目按图的长宽比分类处理：

| 类别 | 处理 | 图 |
|---|---|---|
| 横宽图（宽高比 ≥ 1.3） | `figure*[t]` + `width=\textwidth`（17 cm，约为单栏的 2.1 倍） | 图 1、3、4、5、7、8 |
| 近方图（约 1.0） | `figure[t]` + `width=\columnwidth` | 图 6、10 |
| 竖长图（< 0.8） | `figure*[p]` + `height=0.86\textheight`，独占整页 | 图 2（177 类型排序）、图 9（热图） |
| 表格 | 宽表 `table*[t]` + `tabularx{\textwidth}`；单栏表 `table[t]` + `\columnwidth` | 6 张 |

> ⚠ 双栏下 `figure*` / `table*` **只支持 `[t]` 与 `[p]`**，不能 `[H]`/`[h]`/`[b]`。
> 且模板默认 `dbltopfraction=0.66`，意味着一页顶部最多放 16 cm 的跨栏浮动体——
> 按默认值全部通栏图只能一页一张，会被推成一串半空的浮动页。故在导言区放宽到 0.85。

## 论文结构

```
题目/作者/单位 → 中文摘要 → 关键词 → 中图法分类号（通栏）
英文题目/作者/单位 → Abstract → Key words（通栏）
1 引言 → 2 材料与方法（2.1–2.9）→ 3 结果（3.1–3.9）→ 4 讨论（4.1–4.4）→ 5 结论
→ 声明 → 参考文献 → 附录 A 算法流程 → 附录 B 关键参数
```

其中 **§4.3 方法学教训** 专门写了本项目**自己踩过的七个坑**（双单尾右尾退化、派生列泄漏、
`numpy.argmax` 遇 NaN、PCA 零假设置换对象、Spearman 同义反复、宏 F1 幽灵类、
`NumDBsContributing` 语义误读），这些是本项目最有方法学价值的产出。

## 与报告的关系

- 论文 = **对外发布版本**（结构化、文献引用完整、中文期刊双栏版式、13 页）
- 实验报告 = **对内交付版本**（含全部附录 A–E、逐项算法说明、59 页）

两者数字**必须一致**。项目里有两个独立的终审脚本：

| 脚本 | 检查对象 | 断言数 |
|---|---|---|
| `20-分析/07_终审核对.py` | `40-报告/01-实验报告.md` | 239 |
| `20-分析/95_论文数字终审.py` | `50-论文/main.tex` | 104 |

两个脚本**各自独立从结果文件取真值**（不互相引用），因此报告错了论文不会跟着错。
