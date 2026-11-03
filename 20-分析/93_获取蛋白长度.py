"""
93_获取蛋白长度.py — 批量获取 UniProt 蛋白序列长度
====================================================
用途：
    Moyers & Zhang (2015/2017/2018) 指出 phylostratigraphy 的误差**非随机**：
    序列短、快速演化的基因更容易找不到远缘同源，从而被判为"年轻"。
    他们明确建议：剔除短蛋白（如 < 100 aa）后重跑，看已发表的结论是否消失或反转。

    本项目的年龄表（Trigos et al. 2017）本身**不含序列长度**，
    但它的原始表带 UniProt accession（Entry 列），可据此批量查询长度。

    本脚本把这些长度缓存到 `10-数据/protein_length.tsv`，
    供 `05_监督分析.py --min-length` 做短蛋白剔除检验。

用法：
    python 93_获取蛋白长度.py              # 增量补齐（已缓存的跳过）
    python 93_获取蛋白长度.py --force      # 忽略缓存全部重查

输出：
    10-数据/protein_length.tsv   symbol / uniprot / length
"""

import os
import sys
import time

import pandas as pd
import requests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import DATA_DIR, die, info, ok, warn  # noqa: E402

TPPG_CSV = os.path.join(DATA_DIR, "gene_age_原始_TPPG17.csv")
OUT_TSV = os.path.join(DATA_DIR, "protein_length.tsv")
API = "https://rest.uniprot.org/uniprotkb/search"
# 批量大小：UniProt 的 query 是 URL 参数，OR 拼接 400 个 accession 时
# URL 长度约 4,100 字符会被服务端以 400 拒绝；100 个时约 1,100 字符，稳定。
BATCH = 100
MAX_RETRY = 3


def load_accessions():
    """从 Trigos 原始表拿 symbol → UniProt accession 映射"""
    if not os.path.exists(TPPG_CSV):
        die(f"缺少年龄表原始文件：{TPPG_CSV}")
    df = pd.read_csv(TPPG_CSV)
    need = {"GeneID", "Entry"}
    if not need.issubset(df.columns):
        die(f"{TPPG_CSV} 缺少列 {need}，实际列：{list(df.columns)}")
    df = df[["GeneID", "Entry"]].copy()
    df.columns = ["symbol", "uniprot"]
    df["symbol"] = df["symbol"].astype(str).str.strip().str.upper()
    df["uniprot"] = df["uniprot"].astype(str).str.strip()
    # 去掉空值与占位符
    df = df[(df["uniprot"].str.len() >= 6) & df["uniprot"].str.match(r"^[A-Z0-9]+$")]
    df = df.drop_duplicates(subset="symbol")
    return df


def query_batch(session, accs):
    """查一批 accession 的长度，返回 {accession: length}

    失败时返回 (结果, 是否应降批重试)：
    - 4xx（除 429）是请求本身有问题（如 URL 过长），降批才有用，重试无益
    - 5xx / 网络异常才值得原样重试
    """
    q = " OR ".join(accs)
    params = {
        "query": f"accession:({q})",
        "fields": "accession,length",
        "format": "tsv",
        "size": str(len(accs) + 20),
    }
    for attempt in range(MAX_RETRY):
        try:
            r = session.get(API, params=params, timeout=90)
            if r.status_code == 200:
                out = {}
                for line in r.text.strip().split("\n")[1:]:
                    parts = line.split("\t")
                    if len(parts) >= 2 and parts[1].strip().isdigit():
                        out[parts[0].strip()] = int(parts[1])
                return out, False
            if r.status_code in (400, 413, 414):
                # URL 太长或参数非法：降批，别再退避重试
                return {}, True
            if r.status_code == 429:
                time.sleep(5 * (attempt + 1))
                continue
            time.sleep(2 * (attempt + 1))
        except Exception:
            time.sleep(2 * (attempt + 1))
    return {}, False


def main():
    df = load_accessions()
    info(f"Trigos 表可用基因 {len(df):,} 个（symbol → UniProt accession）")

    cached = {}
    if os.path.exists(OUT_TSV) and "--force" not in sys.argv:
        old = pd.read_csv(OUT_TSV, sep="\t")
        for _, r in old.iterrows():
            try:
                cached[str(r["uniprot"]).strip().upper()] = int(r["length"])
            except Exception:
                pass
        info(f"已有缓存 {len(cached):,} 条（--force 可忽略）")
    else:
        info("--force：忽略已有缓存，全部重查")

    todo = [a for a in df["uniprot"].unique() if a not in cached]
    info(f"待查 {len(todo):,} 个 accession")

    session = requests.Session()
    session.trust_env = False          # 绕开本机代理，UniProt 直连很快

    t0 = time.time()
    done = 0
    failed = 0
    batch_size = BATCH
    if todo:
        i = 0
        while i < len(todo):
            batch = list(todo[i:i + batch_size])
            got, need_smaller = query_batch(session, batch)
            if need_smaller and batch_size > 10:
                batch_size = max(10, batch_size // 2)
                warn(f"批量过大被拒，降到 {batch_size} 个/批后重试")
                continue
            cached.update(got)
            done += len(batch)
            i += len(batch)
            if need_smaller:
                failed += len(batch)
            if (i // 100) % 10 == 0 or i >= len(todo):
                el = time.time() - t0
                rate = done / el if el > 0 else 0
                eta = (len(todo) - done) / rate if rate > 0 else 0
                info(f"  {done:,}/{len(todo):,}  已取回 {len(cached):,} 条  "
                     f"{rate:.0f} acc/s  预计剩余 {eta:.0f}s")
    if failed:
        warn(f"{failed:,} 个 accession 查询失败（已跳过）")
    # 组装输出
    rows = []
    for _, r in df.iterrows():
        acc = r["uniprot"].upper()
        if acc in cached:
            rows.append({"symbol": r["symbol"], "uniprot": acc,
                         "length": cached[acc]})
    out = pd.DataFrame(rows).sort_values("symbol").reset_index(drop=True)
    out.to_csv(OUT_TSV, sep="\t", index=False)
    ok(f"已写出：{OUT_TSV}  （{len(out):,} 个基因）")

    # 统计
    if len(out):
        L = out["length"]
        info(f"长度分布：中位 {L.median():.0f} aa，均值 {L.mean():.0f} aa，"
             f"范围 {L.min():.0f} – {L.max():.0f} aa")
        for thr in (50, 100, 150):
            n = int((L < thr).sum())
            info(f"  长度 < {thr} aa 的基因：{n:,} 个（{n / len(out):.1%}）")


if __name__ == "__main__":
    main()
