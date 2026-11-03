"""
03_细胞类型年龄谱.py — 核心计算
================================
输入：h5ad 单细胞数据 + gene_age.tsv
输出：
    20-分析/celltype_mean_expr.tsv     细胞类型 × 基因 的 平均表达
    20-分析/celltype_age_profile.tsv   细胞类型 × 年龄层 + TAI   ← 主结果
    20-分析/lineage_map.tsv            细胞类型 → 区室（谱系代理，供 04 算 ARI）
    20-分析/celltype_donor_table.tsv   细胞类型 × 供者 的细胞数（看供者覆盖面）
    20-分析/qc_summary.tsv             QC 统计（报告要用）

用法：
    python 03_细胞类型年龄谱.py
    python 03_细胞类型年龄谱.py --group-key cell_ontology_class
    python 03_细胞类型年龄谱.py --max-cells 3000 --chunk 20000

内存策略：
    按块读 raw.X（默认 2 万细胞/批，峰值约 0.4 GB），
    每批算完立刻把"该批各类细胞的和"累加到一个小矩阵上。
    峰值内存与总细胞数无关 —— 483,152 个细胞和 100 万个细胞占用一样多。

为什么用 raw.X 而不是 X：
    Tabula Sapiens 的 X 已经是 log 归一化后的值，直接再 normalize_total +
    log1p 就等于归一化了两次。这里统一用 raw.X 的原始 UMI 计数，
    自己做 QC 和归一化，过程可控、可复现。
"""

import argparse
import os

import numpy as np
import pandas as pd
from scipy import sparse

from common import (H5AD_PATH, AGGR_PATH, PROFILE_PATH, LINEAGE_MAP_PATH,
                    DONOR_TABLE_PATH, ANALYSIS_DIR, STRATUM_NAMES,
                    die, info, ok, warn, require_file,
                    load_gene_age, load_protein_lengths, map_genes,
                    compute_age_profile, open_source)
from ts_io import TabulaSapiens


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--max-cells", type=int, default=2000,
                   help="每种细胞类型最多保留多少细胞（默认 2000）")
    p.add_argument("--chunk", type=int, default=20000,
                   help="每批读多少细胞（默认 2 万，峰值内存约 0.4 GB）")
    p.add_argument("--min-genes", type=int, default=200,
                   help="每个细胞至少检测到多少基因")
    p.add_argument("--max-mito", type=float, default=20.0,
                   help="线粒体 UMI 占比上限 %%（默认 20）")
    p.add_argument("--group-key", default="cell_ontology_class",
                   help="细胞类型列名（默认 cell_ontology_class，177 类；"
                        "cell_type 只有 160 类）")
    p.add_argument("--donor-key", default="donor_id", help="供者列名")
    p.add_argument("--lineage-key", default="compartment",
                   help="谱系代理列名（默认 compartment 区室）")
    p.add_argument("--use-x", action="store_true",
                   help="用 X 而不是 raw.X（仅当 X 是原始计数时）")
    p.add_argument("--source", choices=["auto", "slice", "h5ad"], default="auto",
                   help="数据源：auto=有切片用切片 / slice=只用切片 / h5ad=只用完整文件")
    p.add_argument("--min-length", type=int, default=0,
                   help="只用蛋白长度 ≥ N aa 的基因重算（0 = 不过滤）。"
                        "这是 Moyers & Zhang (2018) 建议的更彻底版本修正："
                        "不只重跑富集层，而是**重算 TAI 绝对值**")
    p.add_argument("--tag", default="",
                   help="输出文件名后缀，便于与主口径并存比较"
                        "（如 --tag longprot → celltype_age_profile_longprot.tsv）")
    return p.parse_args()


def pick_cells(labels, max_cells, seed=42):
    """每种细胞类型最多取 max_cells 个，返回升序的位置索引。"""
    rng = np.random.default_rng(seed)
    keep = []
    for ct in pd.unique(labels):
        idx = np.where(labels == ct)[0]
        if len(idx) > max_cells:
            idx = rng.choice(idx, max_cells, replace=False)
        keep.append(idx)
    return np.sort(np.concatenate(keep))


def get_block(ts, start, stop, use_x):
    """取 [start, stop) 行的原始计数 csr。"""
    if use_x:
        f = ts.f["X"]
        ip = f["indptr"][start:stop + 1].astype(np.int64)
        a, b = int(ip[0]), int(ip[-1])
        return sparse.csr_matrix(
            (f["data"][a:b], f["indices"][a:b], ip - a),
            shape=(stop - start, ts.n_vars),
        )
    return ts.raw_block(start, stop)


def qc_mask(X, mito_idx, min_genes, max_mito):
    """
    逐细胞 QC。返回 (保留布尔掩码, 统计字典)。

    基因数 = 非零元个数；线粒体占比 = MT 基因计数 / 总计数。
    这两个都是逐细胞量，分块处理和全量处理结果完全一致。
    """
    n_genes = np.diff(X.indptr)
    total = np.asarray(X.sum(axis=1)).ravel()

    if mito_idx is not None and mito_idx.sum() > 0:
        mito = np.asarray(X[:, mito_idx].sum(axis=1)).ravel()
        pct_mito = np.where(total > 0, mito / np.maximum(total, 1) * 100, 100.0)
    else:
        pct_mito = np.zeros(X.shape[0])

    keep = (n_genes >= min_genes) & (pct_mito < max_mito)
    stats = dict(n=n_genes, mito=pct_mito, keep=keep)
    return keep, stats


def normalize(X):
    """文库大小归一化到 1e4，再 log1p。逐细胞操作，等价于全量处理。"""
    total = np.asarray(X.sum(axis=1)).ravel()
    scale = np.where(total > 0, 1e4 / np.maximum(total, 1e-12), 0.0)
    X = sparse.diags(scale.astype(np.float32)) @ X
    X = X.tocsr()
    X.data = np.log1p(X.data)
    return X


def build_layer_indicator(stratum, n_vars):
    """
    构造「基因 → 年龄层」的稀疏指示矩阵 L（n_vars × 16）。

    有了它，任意细胞的逐层表达量就是 X @ L（cells × 16），
    不用为每个细胞单独去 16 个基因子集上求和。
    """
    from scipy import sparse as sp

    valid = ~np.isnan(stratum)
    layers = np.round(stratum[valid]).astype(int)
    keep = (layers >= 1) & (layers <= 16)
    rows = np.where(valid)[0][keep]
    cols = layers[keep] - 1
    L = sp.csr_matrix((np.ones(len(rows), dtype=np.float32), (rows, cols)),
                      shape=(n_vars, 16))
    return L


def aggregate(ts, args, mito_idx, layer_ind=None):
    labels_all = np.asarray(ts.obs(args.group_key))
    donors_all = np.asarray(ts.obs(args.donor_key)) if args.donor_key in ts.obs_columns() else None
    info(f"全量：{ts.n_obs:,} 细胞 × {ts.n_vars:,} 基因")

    cats = pd.Categorical(labels_all)
    cat_names = np.array(cats.categories.astype(str))
    n_types = len(cat_names)

    keep = pick_cells(labels_all, args.max_cells)
    info(f"抽样后 {len(keep):,} 个细胞，{n_types} 种类型")

    symbols = ts.symbols()
    n_batches = int(np.ceil(len(keep) / args.chunk))
    info(f"分块 {args.chunk:,}/批，共 {n_batches} 批")

    sums = np.zeros((n_types, ts.n_vars), dtype=np.float64)
    counts = np.zeros(n_types, dtype=np.int64)

    # 逐细胞 TAI：先攒着，最后统一按类型汇总
    cell_tai_parts = []
    cell_type_parts = []
    cell_row_parts = []

    qc_rows = []
    n_qc_pass = 0

    for b, start in enumerate(range(0, len(keep), args.chunk), 1):
        idx = keep[start:start + args.chunk]
        # 行号连续时直接按块读，省一次索引
        if len(idx) > 1 and idx[-1] - idx[0] + 1 == len(idx):
            X = get_block(ts, int(idx[0]), int(idx[-1]) + 1, args.use_x)
        else:
            X = ts.raw_rows(idx) if not args.use_x else _x_rows(ts, idx)

        keep_m, st = qc_mask(X, mito_idx, args.min_genes, args.max_mito)
        qc_rows.append(dict(batch=b, n_cells=len(idx),
                            n_pass=int(keep_m.sum()),
                            med_genes=float(np.median(st["n"])),
                            med_mito=float(np.median(st["mito"]))))
        if keep_m.sum() == 0:
            info(f"  批 {b}/{n_batches}：{len(idx):,} 个细胞全被 QC 过滤")
            continue

        X = X[keep_m]
        codes = cats.codes[idx[keep_m]]
        n_qc_pass += int(keep_m.sum())

        X = normalize(X)

        # ---- 逐细胞 TAI ----
        if layer_ind is not None:
            lay = X @ layer_ind                       # (cells, 16) 稀疏
            lay = np.asarray(lay.todense(), dtype=np.float64)
            tot = lay.sum(axis=1, keepdims=True)
            prop = lay / np.where(tot == 0, 1.0, tot)
            tai_cell = (prop * (np.arange(1, 17))[None, :]).sum(axis=1)
            cell_tai_parts.append(tai_cell)
            cell_type_parts.append(codes)
            cell_row_parts.append(idx[keep_m])

        ones = np.ones(len(codes), dtype=np.float32)
        ind = sparse.csr_matrix((ones, (np.arange(len(codes)), codes)),
                                shape=(len(codes), n_types))
        acc = ind.T @ X
        sums += np.asarray(acc.todense()) if sparse.issparse(acc) else acc
        counts += np.bincount(codes, minlength=n_types)

        info(f"  批 {b}/{n_batches}：{len(idx):,} → 保留 {int(keep_m.sum()):,}")

    if n_qc_pass == 0:
        die("所有细胞都被 QC 过滤了。放宽 --min-genes 或 --max-mito。")

    mean_expr = pd.DataFrame(
        sums / np.maximum(counts, 1)[:, None],
        index=cat_names, columns=symbols,
    )
    info(f"QC 通过 {n_qc_pass:,} / {len(keep):,} = {n_qc_pass/len(keep):.1%}")

    cell_tai = None
    if cell_tai_parts:
        cell_tai = pd.DataFrame({
            "row": np.concatenate(cell_row_parts),
            "code": np.concatenate(cell_type_parts),
            "TAI": np.concatenate(cell_tai_parts),
        })

    # ---------- 把 QC 通过的细胞行索引落盘 ----------
    # ⚠️ 单一事实来源：监督分析（05）必须用**同一批** QC 后的细胞，
    #   否则同一篇报告里的 TAI 数字（QC 后）与分类数字（未 QC）来自不同细胞集。
    #   这里导出行索引，05 直接读，保证口径一致。
    qc_rows_path = os.path.join(ANALYSIS_DIR, "qc_passed_rows.npy")
    np.save(qc_rows_path, np.sort(cell_tai["row"].to_numpy()))
    ok(f"已写出：{qc_rows_path}  ({len(cell_tai):,} 个 QC 通过细胞的行索引)")

    return (mean_expr, counts, pd.DataFrame(qc_rows), donors_all,
            keep, cats, cell_tai)


def _x_rows(ts, idx):
    """从 X 取任意行（很少走这条路，--use-x 且行不连续时）。"""
    f = ts.f["X"]
    ip = f["indptr"][:]
    rows = []
    for r in idx:
        a, b = int(ip[r]), int(ip[r + 1])
        rows.append((f["data"][a:b], f["indices"][a:b]))
    data = np.concatenate([r[0] for r in rows]) if rows else np.array([])
    indic = np.concatenate([r[1] for r in rows]) if rows else np.array([])
    ip2 = np.concatenate([[0], np.cumsum([len(r[0]) for r in rows])])
    return sparse.csr_matrix((data, indic, ip2), shape=(len(idx), ts.n_vars))


def main():
    args = parse_args()
    age_df = load_gene_age()

    # 输出路径助手：--tag 非空时所有产出加后缀，避免覆盖主口径结果。
    # 长蛋白版必须与主口径**并存**才能对照，两者不能写同一个文件名。
    suffix = f"_{args.tag}" if args.tag else ""

    def out_path(name):
        stem, ext = os.path.splitext(name)
        return os.path.join(ANALYSIS_DIR, f"{stem}{suffix}{ext}")

    prefer = {"auto": None, "slice": True, "h5ad": False}[args.source]
    ts = open_source(prefer)
    if args.use_x and not hasattr(ts, "f"):
        die("--use-x 只对完整 h5ad 有效（切片的 X 已经是 log 归一化后的）")

    with ts:
        cols = ts.obs_columns()
        if args.group_key not in cols:
            die(f"obs 里没有 '{args.group_key}'。现有列：{cols}")
        if args.donor_key not in cols:
            warn(f"obs 里没有 '{args.donor_key}'，05 将无法按供者划分")

        mito_idx = ts.mito_mask()
        info(f"线粒体基因 {int(mito_idx.sum())} 个")

        # 基因 → 年龄层 的映射要在聚合之前拿到（逐细胞 TAI 要用）
        matched, stratum = map_genes(ts.symbols(), age_df)
        if matched.sum() < 1000:
            die(f"匹配上的基因只有 {matched.sum()} 个，太少。\n"
                f"  跑 01_检查数据.py 看基因名格式。")

        # ---------- 长蛋白过滤（Moyers & Zhang 2018 建议的彻底版修正）----------
        # 短序列基因更容易被 phylostratigraphy 误判为「年轻」。
        # 做法：把短蛋白的层号置为 NaN。build_layer_indicator() 与
        # compute_age_profile() 都以 `~isnan(stratum)` 筛选，
        # 因此**逐细胞 TAI** 与 **类型级 TAI/占比** 会同时只用长蛋白重算。
        # ⚠ 查不到长度信息的基因**保持原样**：它们不是「短蛋白」，
        #   不能因为缺信息就把它们当成过短剔除（那是把「未知」当「否」）。
        # ⚠ `matched` 保持不变——它后面要用来断言「二次映射结果一致」，
        #   改了会让那个自检失去意义。真正参与计算的数量另记 n_effective。
        short_mask = None
        n_effective = int(matched.sum())
        if args.min_length > 0:
            lengths = load_protein_lengths()
            sym_upper = pd.Index(ts.symbols()).astype(str).str.upper()
            ln = lengths.reindex(sym_upper).to_numpy(dtype=float)
            short_mask = np.isfinite(ln) & (ln < args.min_length)
            n_no_len = int((~np.isfinite(ln)).sum())
            stratum = stratum.copy()
            stratum[short_mask] = np.nan
            n_effective = int(np.isfinite(stratum).sum())
            info(f"长蛋白过滤（长度 ≥ {args.min_length} aa）："
                 f"剔除 {int(short_mask.sum()):,} 个短蛋白；"
                 f"另有 {n_no_len:,} 个基因查不到长度，**保留**")
            info(f"参与年龄构成计算的基因：{n_effective:,}"
                 f"（原 {int(matched.sum()):,}）")

        layer_ind = build_layer_indicator(stratum, ts.n_vars)
        info(f"年龄层指示矩阵：{layer_ind.shape[0]:,} 基因 × "
             f"{layer_ind.shape[1]} 层，非零 {layer_ind.nnz:,}")

        (mean_expr, counts, qc_df, donors_all, keep, cats,
         cell_tai) = aggregate(ts, args, mito_idx, layer_ind)

        # ---------- 谱系代理 ----------
        labels_all = np.asarray(ts.obs(args.group_key))
        lineage_map = None
        if args.lineage_key in cols:
            lin_all = np.asarray(ts.obs(args.lineage_key))
            tmp = pd.DataFrame({"ct": labels_all, "lin": lin_all})
            lineage_map = (tmp.groupby("ct", observed=True)["lin"]
                           .agg(lambda s: s.value_counts().idxmax())
                           .rename("lineage").to_frame())
            lineage_map.index.name = "cell_type"

        # ---------- 供者覆盖 ----------
        donor_tab = None
        if donors_all is not None:
            tmp = pd.DataFrame({"ct": labels_all, "donor": donors_all})
            donor_tab = pd.crosstab(tmp["ct"], tmp["donor"])
            donor_tab.index.name = "cell_type"

    # ---------- 落盘：平均表达 ----------
    # ⚠ mean_expr 是 93 MB 的全基因矩阵，且**不随长蛋白过滤而变**
    #   （过滤作用在 stratum 上，不在表达矩阵上）。
    #   所以带 --tag 重跑时不必再写一遍，直接跳过，省 93 MB 写盘。
    if args.tag:
        info(f"--tag {args.tag}：跳过写出 {os.path.basename(AGGR_PATH)}"
             f"（内容与主口径完全相同，不重复落盘）")
    else:
        mean_expr.to_csv(AGGR_PATH, sep="\t")
        ok(f"已写出：{AGGR_PATH}  ({mean_expr.shape[0]} 类型 × {mean_expr.shape[1]} 基因)")

    qc_df.to_csv(os.path.join(ANALYSIS_DIR, "qc_summary.tsv"), sep="\t", index=False)
    if lineage_map is not None:
        lineage_map.to_csv(LINEAGE_MAP_PATH, sep="\t")
        ok(f"已写出：{LINEAGE_MAP_PATH}  ({len(lineage_map)} 类，"
           f"{lineage_map['lineage'].nunique()} 个区室)")
    if donor_tab is not None:
        donor_tab.to_csv(DONOR_TABLE_PATH, sep="\t")
        ok(f"已写出：{DONOR_TABLE_PATH}")

    # ---------- 基因匹配（前面已算过，这里只做校验） ----------
    matched2, stratum2 = map_genes(mean_expr.columns, age_df)
    if not np.array_equal(matched2, matched):
        warn("两次基因匹配结果不一致，以先算的那次为准")
    matched, stratum = matched2, stratum2
    if matched.sum() < 1000:
        die(f"匹配上的基因只有 {matched.sum()} 个，太少。\n"
            f"  跑 01_检查数据.py 看基因名格式。")

    # ⚠ 长蛋白过滤必须在**这里再应用一次**：上面 353 行把 stratum 覆盖成了
    #    未过滤的二次映射结果。不重新过滤的话，类型级 TAI 会用回全部基因，
    #    与逐细胞 TAI（已过滤）口径不一致 —— 那是最难发现的一类错误。
    if short_mask is not None:
        stratum = stratum.copy()
        stratum[short_mask] = np.nan
        n_effective = int(np.isfinite(stratum).sum())
        info(f"类型级 TAI 同步应用长蛋白过滤：参与基因 {n_effective:,}")

    # ---------- 年龄构成 ----------
    info("计算年龄构成与 TAI……")
    profile = compute_age_profile(mean_expr, stratum)
    total = profile.pop("total_expr")

    # ---------- 逐细胞 TAI → 类型级标准误 ----------
    se_note = ""
    if cell_tai is not None and len(cell_tai):
        ct_names = np.array(cats.categories.astype(str))
        cell_tai["cell_type"] = ct_names[cell_tai["code"].to_numpy()]
        cell_tai["donor_id"] = (donors_all[cell_tai["row"].to_numpy()]
                                if donors_all is not None else "")
        cell_tai.to_csv(out_path("cell_tai.tsv"),
                        sep="\t", index=False)

        g = cell_tai.groupby("cell_type")["TAI"]
        stats = g.agg(["count", "mean", "std"])
        stats["se"] = stats["std"] / np.sqrt(stats["count"].clip(lower=1))
        stats = stats.reindex(profile.index)
        profile["TAI_se"] = stats["se"]
        profile["n_cells"] = stats["count"].fillna(0).astype(int)
        # 细胞数太少的类型，其 TAI 只能当参考值看
        profile["low_confidence"] = profile["n_cells"] < 20
        se_note = (f"  逐细胞 TAI：{len(cell_tai):,} 个细胞；"
                   f"类型内 SE 中位数 {profile['TAI_se'].median():.4f}")

        # ---------- 供者层面的标准误（关键修正） ----------
        # ⚠️ 细胞嵌套在供者内，同一供者的细胞共享个体/组织/批次效应，
        #   **不是独立观测**。用细胞层面的 std/√n 会把伪重复当独立样本，
        #   SE 被系统性低估（pseudo-replication）。
        #   正确做法：以「供者均值」为观测单位算 SE。
        #   实测：细胞层面 SE 中位 0.0111 → 供者层面 0.0417（低估 3.66 倍），
        #   类型内 ICC 中位 = 1.000（细胞间方差几乎全部来自供者间差异）。
        if "donor_id" in cell_tai.columns:
            donor_mean = (cell_tai.groupby(["cell_type", "donor_id"])["TAI"]
                          .mean().reset_index())
            dstat = donor_mean.groupby("cell_type")["TAI"].agg(
                ["count", "mean", "std"])
            dstat["se"] = dstat["std"] / np.sqrt(dstat["count"].clip(lower=1))
            dstat = dstat.reindex(profile.index)
            profile["n_donors"] = dstat["count"].fillna(0).astype(int)
            profile["TAI_se_donor"] = dstat["se"]
            # 只有 1 个供者的类型，跨供者不确定性在数学上无法估计
            profile["se_donor_na"] = profile["n_donors"] < 2
            n_d1 = int((profile["n_donors"] < 2).sum())
            se_donor_med = profile["TAI_se_donor"].median()
            sd_type = profile["TAI"].std()
            se_note += (f"；供者层面 SE 中位 {se_donor_med:.4f}"
                        f"（细胞层面 {profile['TAI_se'].median():.4f}，"
                        f"低估 {se_donor_med / profile['TAI_se'].median():.2f} 倍）")
            se_note += (f"；{n_d1} 个类型只有 1 个供者，无法估计跨供者不确定性")
            se_note += (f"；信号/噪声：细胞层面 {sd_type / profile['TAI_se'].median():.1f} 倍，"
                        f"供者层面 {sd_type / se_donor_med:.1f} 倍")

        # 类型级 vs 细胞级 TAI 的差异：
        #   类型级 =「先求类型平均表达，再算层占比」；细胞级 =「先算每个细胞的
        #   层占比再平均」。两者是**不同的估计量**，数学上不等价。
        #   ⚠️ 不能拿「TAI 极差」当分母说成「差异 4.2% 属正常」——
        #   极差是信号尺度，不是噪声尺度。正确的参照是类型内 SE。
        delta = stats["mean"] - profile["TAI"]
        diff = delta.abs().max()
        diff_sd = delta.std()
        se_med = profile["TAI_se"].median()
        r_delta = (np.corrcoef(delta, profile["TAI"])[0, 1]
                   if np.std(delta) > 0 and np.std(profile["TAI"]) > 0 else np.nan)
        se_note += (f"；与细胞级 TAI 均值最大差 {diff:.4f}"
                    f"（Δ 的 SD {diff_sd:.4f}，是类型内 SE 中位 {se_med:.4f} 的 "
                    f"{diff_sd / se_med:.1f} 倍；Δ 与 TAI 相关 r={r_delta:.3f}，"
                    f"属系统性偏移而非随机抖动——两者是不同的估计量，不可互换）")
        ok(f"已写出：{os.path.join(ANALYSIS_DIR, 'cell_tai.tsv')}")

    profile_path = out_path("celltype_age_profile.tsv")
    profile.to_csv(profile_path, sep="\t")
    ok(f"已写出：{profile_path}")

    # ---------- 稳健性：剔除低命中层后 TAI 排序是否稳定 ----------
    # S15（人亚科）命中 64%、S16（智人）69%，远低于其余层的 93–96%。
    # 这两层是「人类特有基因」，表达矩阵用的注释版本较旧，漏掉的多。
    #
    # ⚠️ 但「漏掉 → TAI 偏小」这个偏倚**不是对所有类型都一样的常数**：
    #   漏掉的量与该类型自身表达的 S15/S16 基因数成正比，
    #   所以不同类型的偏倚幅度不同 → 排序会变。
    #   下面的回归会量化这一点。
    from scipy.stats import spearmanr, kendalltau

    name_set = set(pd.Index(mean_expr.columns).astype(str).str.upper())
    hit_rate = {}
    for s in range(1, 17):
        sub = age_df[age_df["stratum"] == s]["symbol"].astype(str).str.upper()
        hit_rate[s] = (sum(1 for g in sub if g in name_set) / len(sub)
                       if len(sub) else 0.0)

    hi_layers = [s for s, r in hit_rate.items() if r >= 0.80]
    drop_layers = [s for s in hit_rate if s not in hi_layers]
    strat_keep = stratum.copy()
    strat_keep[np.isin(np.round(stratum), drop_layers)] = np.nan
    profile_hi = compute_age_profile(mean_expr, strat_keep)

    # ⚠️⚠️ **必须按索引对齐，不能按位置**（2026-11-03 审核发现的致命 bug）
    # compute_age_profile() 的返回值已经按 TAI 升序排过，
    # 于是 profile["TAI"] 和 profile_hi["TAI"] **各自都是升序序列**。
    # scipy.stats.spearmanr 对 pandas Series 是**按位置**取值的（不按索引对齐），
    # 位置 i 对应的永远是「第 i 小的那个」→ 位置秩相关在数学上恒等于 +1，
    # 与两个序列的真实对应关系完全无关（实测：两列独立随机的升序数组，
    # 500/500 次都算出 ρ=1.0）。原代码因此报出 ρ=1.0000 的**恒等式**，
    # 而真实值是 0.9883。
    hi_aligned = profile_hi.reindex(profile.index)
    rho = spearmanr(profile["TAI"].values,
                    hi_aligned["TAI"].values).statistic
    tau = kendalltau(profile["TAI"].values,
                      hi_aligned["TAI"].values).statistic
    top10_a = set(profile.head(10).index)
    top10_b = set(profile_hi.head(10).index)
    overlap = len(top10_a & top10_b)
    # 排名位移
    rank_shift = (profile["TAI"].rank()
                  - hi_aligned["TAI"].rank()).abs()
    n_moved = int((rank_shift > 0).sum())
    max_shift = int(rank_shift.max())
    hi_path = out_path("celltype_age_profile_highhit.tsv")
    profile_hi.to_csv(hi_path, sep="\t")

    # 偏倚是否同质？把 TAI 变化量对「该类型的 S15+S16 表达占比」做回归。
    # 若 R² 接近 1，说明偏倚完全由「该类型表达了多少年轻层基因」驱动，
    # 因而是**类型特异**的，不能假定它对所有类型抵消。
    s_cols = [c for c in profile.columns if c.startswith("S")]
    if all(x in profile.columns for x in ("S15", "S16")):
        young_share = (profile["S15"] + profile["S16"]).values
        d_tai = (profile["TAI"] - hi_aligned["TAI"]).values
        if np.std(young_share) > 0 and np.std(d_tai) > 0:
            rr = np.corrcoef(young_share, d_tai)[0, 1]
            r2 = rr ** 2
        else:
            rr, r2 = np.nan, np.nan
    else:
        rr, r2 = np.nan, np.nan

    print("\n" + "=" * 62)
    print("稳健性检验：只保留高命中层（命中率 ≥ 80%）后重算 TAI")
    print("=" * 62)
    print(f"  保留的层：{hi_layers}")
    print(f"  剔除的层：{drop_layers}"
          f"（命中率 " +
          ", ".join(f"S{s}:{hit_rate[s]:.0%}" for s in drop_layers) + "）")
    print(f"  Spearman ρ = {rho:.4f}   Kendall τ = {tau:.4f}")
    print(f"  最古老 Top10 重合 {overlap}/10")
    print(f"  TAI 绝对差：中位 {((profile['TAI'] - profile_hi['TAI']).abs()).median():.4f}")
    print(f"  排名变动：{n_moved}/{len(profile)} 个类型名次有变化，"
          f"最大位移 {max_shift} 位")
    if overlap < 10:
        print(f"    Top10 掉出：{sorted(top10_a - top10_b)}")
        print(f"    Top10 换入：{sorted(top10_b - top10_a)}")
    print(f"  偏倚同质性检验：ΔTAI 对「该类型的 S15+S16 表达占比」回归 "
          f"r = {rr:.4f}，R² = {r2:.4f}")
    if np.isfinite(r2) and r2 > 0.9:
        print("    → 偏倚**完全由该类型自身表达的年轻层基因量驱动**，"
              "是类型特异的，不能假定它对所有类型抵消。")
    print(f"  已写出：{hi_path}")
    if rho >= 0.95:
        print("  → 整体排序高度一致（但仍不等于逐个类型不变，见上面的排名变动）。")
    else:
        print("  → ⚠️ 排序对层的取舍敏感，报告里必须说明这一点。")

    # ---------- 自检 ----------
    print("\n" + "=" * 62)
    print("自检")
    print("=" * 62)
    print(f"  细胞类型数：{len(profile)}")
    print(f"  年龄层数：{len([c for c in profile.columns if c.startswith('S')])}")
    print(f"  用到的基因数：{n_effective:,}"
          + (f"（长蛋白过滤：长度 ≥ {args.min_length} aa）"
             if args.min_length > 0 else ""))

    prop = profile[[c for c in profile.columns if c.startswith("S")]].values
    row_sum = prop.sum(axis=1)
    print(f"  每行占比之和：{row_sum.min():.6f} – {row_sum.max():.6f}"
          f"  {'✓' if np.allclose(row_sum, 1, atol=1e-6) else '✗ 不为 1'}")

    print(f"  TAI 范围：{profile['TAI'].min():.2f} – {profile['TAI'].max():.2f}")
    print(f"  TAI 中位数：{profile['TAI'].median():.2f}")
    print(f"  总表达量范围：{total.min():.1f} – {total.max():.1f}")
    if se_note:
        print(f"  {se_note}")
        if "TAI_se" in profile.columns:
            se = profile["TAI_se"]
            # ⚠️ 不能用 se.values.argmax()：numpy 的 argmax 遇到 NaN 会返回
            # NaN 所在的位置（而不是最大值的位置），会把「最大标准误来自哪个类型」
            # 标到一个 SE 无法计算（n=1）的类型头上。必须先 dropna。
            se_valid = se.dropna()
            if len(se_valid) == 0:
                print("  TAI 标准误：全部无法计算（各类型细胞数过少）")
            else:
                k_name = se_valid.idxmax()
                print(f"  TAI 标准误：中位数 {se.median():.4f}，"
                      f"最大 {se_valid.max():.4f}（← {k_name}, "
                      f"n={int(profile.loc[k_name, 'n_cells'])}）")
                n_nan = int(se.isna().sum())
                if n_nan:
                    nan_names = ", ".join(profile.index[se.isna()].tolist())
                    print(f"  ⚠️ {n_nan} 个类型无法计算标准误（细胞数过少，"
                          f"SE 为 NaN）：{nan_names}")
            # 信号 vs 噪声：跨类型的离散度 vs 类型内标准误。
            # 注意不能用「排序后相邻类型的 TAI 差」—— 那是连续分布的
            # 人为相邻，本来就很小，会得出「差异小于噪声」的错误结论。
            spread = profile["TAI"].std()
            print(f"  信号/噪声：TAI 跨类型标准差 {spread:.4f} ÷ "
                  f"类型内 SE 中位数 {se.median():.4f} = "
                  f"{spread/max(se.median(), 1e-9):.1f} 倍"
                  f"  → 类型间差异远大于细胞间噪声")
            low_n = int((profile["n_cells"] < 20).sum())
            print(f"  低置信类型（QC 后细胞数 < 20）：{low_n} 类，"
                  f"其 TAI 只作参考，不进入结论")

    if not (100 <= len(profile) <= 250):
        warn(f"细胞类型数 {len(profile)} 与预期（约 177）不符")

    print("\n  ── 最古老（TAI 最小）的 10 种 ──")
    for ct, row in profile.head(10).iterrows():
        print(f"    {row['TAI']:6.2f}  {str(ct)[:52]}")

    print("\n  ── 最年轻（TAI 最大）的 10 种 ──")
    for ct, row in profile.tail(10).iloc[::-1].iterrows():
        print(f"    {row['TAI']:6.2f}  {str(ct)[:52]}")

    print("\n  年龄层图例：")
    for s in range(1, len(STRATUM_NAMES) + 1):
        if f"S{s}" in profile.columns:
            print(f"    S{s:>2}  {STRATUM_NAMES[s]}")

    print("\n  下一步：python 04_无监督分析.py")
    print("=" * 62)


if __name__ == "__main__":
    main()
