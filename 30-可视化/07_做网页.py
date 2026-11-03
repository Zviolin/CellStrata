"""
07_做网页.py — 生成可交互的「细胞地层图」网页
================================================
输入：20-分析/celltype_age_profile.tsv
      （03_细胞类型年龄谱.py 的产出）
输出：30-可视化/cellstrata.html

为什么是单文件 HTML：
    不依赖任何 CDN、不依赖服务器。双击就能开，发给老师也能开。
    数据以 JSON 内嵌进页面，所以离线可用。

页面上能做什么：
    - 每行一个细胞类型，横向堆叠 16 个年龄层（深红=最老 → 亮蓝=最新）
    - 按 TAI 排序，一眼看出谁古老谁年轻
    - 点任意一行 → 右侧展开该细胞类型的 16 层明细
    - 搜索框过滤细胞类型名
    - 切换「按 TAI / 按名字 / 按区室」排序

用法：
    python 07_做网页.py
    python 07_做网页.py --open        # 生成后自动用浏览器打开
"""

import argparse
import json
import os
import sys
import webbrowser

# 允许从 30-可视化/ 目录直接跑脚本时 import 到 20-分析/common.py
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "20-分析"))

from common import (PROFILE_PATH, LINEAGE_MAP_PATH, FIG_DIR,  # noqa: E402
                    STRATUM_NAMES, info, ok, require_file)

VIS_DIR = os.path.dirname(FIG_DIR)


TEMPLATE = """<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>细胞地层 CellStrata</title>
<style>
  :root{
    --bg:#0f1218; --panel:#171b24; --panel2:#1e2330; --line:#2a3040;
    --fg:#e8ecf4; --dim:#93a0b8; --accent:#5aa9e6;
  }
  *{box-sizing:border-box}
  body{margin:0;background:var(--bg);color:var(--fg);
       font:14px/1.55 -apple-system,"Segoe UI","Microsoft YaHei",sans-serif}
  header{padding:22px 28px 14px;border-bottom:1px solid var(--line)}
  h1{margin:0 0 6px;font-size:20px;font-weight:600;letter-spacing:.3px}
  .sub{color:var(--dim);font-size:13px}
  .stats{display:flex;gap:26px;margin-top:14px;flex-wrap:wrap}
  .stat b{display:block;font-size:19px;font-weight:600;color:var(--accent)}
  .stat span{font-size:12px;color:var(--dim)}
  .toolbar{display:flex;gap:12px;align-items:center;padding:12px 28px;
           border-bottom:1px solid var(--line);flex-wrap:wrap}
  input[type=search]{background:var(--panel2);border:1px solid var(--line);
    color:var(--fg);padding:7px 11px;border-radius:7px;width:230px;outline:none}
  input[type=search]:focus{border-color:var(--accent)}
  .seg{display:flex;border:1px solid var(--line);border-radius:7px;overflow:hidden}
  .seg button{background:transparent;border:0;color:var(--dim);padding:7px 12px;
    cursor:pointer;font-size:13px}
  .seg button.on{background:var(--accent);color:#08121e;font-weight:600}
  main{display:grid;grid-template-columns:1fr 330px;gap:0;height:calc(100vh - 168px)}
  #list{overflow-y:auto;padding:10px 18px 40px}
  .row{display:grid;grid-template-columns:230px 1fr 54px;gap:12px;
       align-items:center;padding:3px 8px;border-radius:6px;cursor:pointer}
  .row:hover{background:var(--panel)}
  .row.sel{background:var(--panel2);box-shadow:inset 0 0 0 1px var(--accent)}
  .name{font-size:12px;color:#c8d3e6;white-space:nowrap;overflow:hidden;
        text-overflow:ellipsis}
  .bar{display:flex;height:15px;border-radius:3px;overflow:hidden}
  .bar i{display:block;height:100%}
  .tai{font-size:11px;color:var(--dim);text-align:right;font-variant-numeric:tabular-nums}
  .lc{display:inline-block;margin-left:6px;padding:1px 5px;border-radius:4px;
      font-size:9px;line-height:1.5;color:#c9a227;background:rgba(201,162,39,.13);
      border:1px solid rgba(201,162,39,.32);vertical-align:1px;cursor:help}
  .warn{margin:0 0 14px;padding:8px 10px;border-radius:7px;font-size:11.5px;
        color:#c9a227;background:rgba(201,162,39,.10);
        border:1px solid rgba(201,162,39,.28);line-height:1.6}
  aside{border-left:1px solid var(--line);background:var(--panel);overflow-y:auto;
        padding:18px}
  aside h2{font-size:14px;margin:0 0 4px;line-height:1.35}
  aside .meta{color:var(--dim);font-size:12px;margin-bottom:16px}
  .det{display:grid;grid-template-columns:12px 1fr 46px;gap:9px;align-items:center;
       margin-bottom:5px}
  .chip{width:12px;height:12px;border-radius:3px}
  .det .lb{font-size:11px;color:#c2ccdd;white-space:nowrap;overflow:hidden;
           text-overflow:ellipsis}
  .det .pv{font-size:11px;color:var(--dim);text-align:right;
           font-variant-numeric:tabular-nums}
  .track{height:4px;background:var(--panel2);border-radius:2px;margin-top:2px;
         grid-column:2/4}
  .track i{display:block;height:100%;border-radius:2px}
  .legend{display:flex;align-items:center;gap:8px;font-size:11px;color:var(--dim);
          margin-top:16px}
  .legend .ramp{flex:1;height:9px;border-radius:4px;
    background:linear-gradient(90deg,#b40426,#dddddd,#3b4cc0)}
  footer{color:var(--dim);font-size:11px;padding:10px 28px;border-top:1px solid var(--line)}
  #empty{color:var(--dim);padding:40px 8px}
</style>
</head>
<body>
<header>
  <h1>细胞地层 CellStrata</h1>
  <div class="sub">人体每一种细胞，是由哪个年代的基因搭起来的</div>
  <div class="stats">
    <div class="stat"><b>__N_TYPES__</b><span>细胞类型</span></div>
    <div class="stat"><b>__N_STRATA__</b><span>基因年龄层</span></div>
    <div class="stat"><b>__TAI_MIN__ – __TAI_MAX__</b><span>TAI 范围（越小越古老）</span></div>
    <div class="stat"><b>__N_GENES__</b><span>有年龄注释的基因</span></div>
    <div class="stat"><b>__N_CELLS__</b><span>QC 后细胞（母体 483,152 → 切片 28,017）</span></div>
    <div class="stat"><b>__SE_DONOR__</b><span>供者级标准误中位数</span></div>
  </div>
</header>
<div class="toolbar">
  <input id="q" type="search" placeholder="搜索细胞类型…">
  <div class="seg" id="sortseg">
    <button data-k="tai" class="on">按 TAI</button>
    <button data-k="name">按名字</button>
    <button data-k="lineage">按区室</button>
  </div>
  <span class="sub" id="count"></span>
</div>
<main>
  <div id="list"></div>
  <aside id="detail">
    <h2>点左侧任意一行</h2>
    <div class="meta">查看该细胞类型的 16 层基因年龄构成</div>
    <div class="legend"><span>最古老</span><div class="ramp"></div><span>最新</span></div>
  </aside>
</main>
<footer>数据来源：Tabula Sapiens 1.0（母体 483,152 细胞 / 177 种细胞类型）
→ 本研究经 HTTP Range 远程切片抽取 <b>28,017</b> 个细胞 → 质量控制后保留 <b>22,195</b> 个细胞（79.2%）
× Trigos et al. 2017 基因年龄表（17,318 个基因 / 16 个年龄层）</footer>

<script>
const DATA = __DATA__;

const RED=[180,4,38], MID=[221,221,221], BLUE=[59,76,192];
function ramp(i,n){
  const t = n<=1 ? 0 : i/(n-1);
  const a = t<0.5 ? RED : MID, b = t<0.5 ? MID : BLUE;
  const u = t<0.5 ? t*2 : (t-0.5)*2;
  const c = a.map((v,k)=>Math.round(v+(b[k]-v)*u));
  return `rgb(${c[0]},${c[1]},${c[2]})`;
}

const N = DATA.strata.length;
const COLORS = Array.from({length:N},(_,i)=>ramp(i,N));
const ROWS = DATA.cell_types.map((ct,idx)=>({
  idx, ct,
  vals: DATA.profile[idx],
  tai: DATA.tai[idx],
  lineage: DATA.lineage[ct] || "—",
  // 标准误可能为 null（细胞数过少时数学上无法估计），不要当成 0
  se: (DATA.tai_se && DATA.tai_se[idx]!=null) ? DATA.tai_se[idx] : null,
  // 供者级标准误：推断单位是供者，这一口径才是保守的
  seDonor: (DATA.tai_se_donor && DATA.tai_se_donor[idx]!=null)
           ? DATA.tai_se_donor[idx] : null,
  nDonors: (DATA.n_donors && DATA.n_donors[idx]!=null) ? DATA.n_donors[idx] : null,
  n: (DATA.n_cells && DATA.n_cells[idx]!=null) ? DATA.n_cells[idx] : null,
  low: (DATA.n_cells && DATA.n_cells[idx]!=null) ? DATA.n_cells[idx] < 20 : false
}));

let sortKey="tai", query="", selected=null;

function visible(){
  let r = ROWS.filter(x=>x.ct.toLowerCase().includes(query.toLowerCase()));
  if(sortKey==="tai") r.sort((a,b)=>a.tai-b.tai);
  else if(sortKey==="name") r.sort((a,b)=>a.ct.localeCompare(b.ct));
  else r.sort((a,b)=> (a.lineage||"").localeCompare(b.lineage||"") || a.tai-b.tai);
  return r;
}

function render(){
  const rows = visible();
  document.getElementById("count").textContent = `${rows.length} / ${ROWS.length} 种`;
  const list = document.getElementById("list");
  list.innerHTML = "";
  if(!rows.length){ list.innerHTML='<div id="empty">没有匹配的细胞类型</div>'; return; }
  const frag = document.createDocumentFragment();
  rows.forEach(r=>{
    const div = document.createElement("div");
    div.className = "row" + (selected===r.idx?" sel":"");
    const pct = r.vals.map(v=>(v*100)).map(v=>v.toFixed(2));
    const bar = r.vals.map((v,i)=>
      `<i style="width:${(v*100).toFixed(3)}%;background:${COLORS[i]}"></i>`).join("");
    div.innerHTML =
      `<div class="name" title="${r.ct}">${r.ct}`+
        (r.low?`<sup class="lc" title="QC 后细胞数不足 20，该类型的 TAI 仅作参考，不进入结论">低置信</sup>`:"")+
      `</div>`+
      `<div class="bar" title="${pct.join(' / ')}">${bar}</div>`+
      `<div class="tai">${r.tai.toFixed(2)}</div>`;
    div.onclick = ()=>{ selected=r.idx; render(); showDetail(r); };
    frag.appendChild(div);
  });
  list.appendChild(frag);
}

function showDetail(r){
  const rows = r.vals.map((v,i)=>{
    const p = (v*100);
    return `<div class="det">
      <div class="chip" style="background:${COLORS[i]}"></div>
      <div class="lb" title="${DATA.strata[i]}">S${i+1} ${DATA.strata[i]}</div>
      <div class="pv">${p.toFixed(1)}%</div>
      <div class="track"><i style="width:${Math.min(100,p*2.2).toFixed(1)}%;background:${COLORS[i]}"></i></div>
    </div>`;
  }).join("");
  document.getElementById("detail").innerHTML = `
    <h2>${r.ct}</h2>
    <div class="meta">
      区室：${r.lineage} ·
      TAI = <b style="color:var(--accent)">${r.tai.toFixed(3)}</b>
      ${r.se!=null
        ? `· 细胞级标准误 ±${r.se.toFixed(3)}`
        : `· <span style="color:#c9a227">标准误无法估计（细胞数过少）</span>`}
      ${r.seDonor!=null
        ? `· 供者级标准误 ±${r.seDonor.toFixed(3)}`
        : `· <span style="color:#c9a227">供者级标准误无法估计（仅 1 名供者）</span>`}
      ${r.nDonors!=null?`· 供者数 ${r.nDonors}`:""}
      ${r.n!=null?`· QC 后细胞数 ${r.n}`:""}
    </div>
    ${r.low?`<div class="warn">⚠ 该类型 QC 后细胞数不足 20，TAI 仅供参考，不进入结论。</div>`:""}
    ${(r.nDonors!=null&&r.nDonors<2)?`<div class="warn">⚠ 该类型只有 1 名供者，无法估计跨供者不确定性；供者级标准误不可估。</div>`:""}
    ${rows}
    <div class="legend"><span>最古老</span><div class="ramp"></div><span>最新</span></div>`;
}

document.getElementById("q").oninput = e=>{ query=e.target.value; render(); };
document.querySelectorAll("#sortseg button").forEach(b=>{
  b.onclick = ()=>{
    document.querySelectorAll("#sortseg button").forEach(x=>x.classList.remove("on"));
    b.classList.add("on"); sortKey=b.dataset.k; render();
  };
});

render();
</script>
</body>
</html>
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--open", action="store_true", help="生成后用浏览器打开")
    ap.add_argument("--profile", default=None,
                    help="指定 profile 表（默认 20-分析/celltype_age_profile.tsv）")
    ap.add_argument("--out", default=None,
                    help="指定输出 html（默认 30-可视化/cellstrata.html）")
    ap.add_argument("--n-genes", default="17,318", help="页面上显示的基因数文案")
    args = ap.parse_args()

    prof_path = args.profile or PROFILE_PATH
    require_file(prof_path, "请先运行 03_细胞类型年龄谱.py")

    import pandas as pd
    prof = pd.read_csv(prof_path, sep="\t", index_col=0)
    s_cols = [c for c in prof.columns if c.startswith("S")]
    s_cols.sort(key=lambda x: int(x[1:]))
    info(f"读入 {len(prof)} 种细胞类型 × {len(s_cols)} 层")

    lineage = {}
    if os.path.exists(LINEAGE_MAP_PATH):
        lm = pd.read_csv(LINEAGE_MAP_PATH, sep="\t", index_col=0)
        if "lineage" in lm.columns:
            lineage = lm["lineage"].to_dict()
        info(f"并入区室注释：{len(lineage)} 类")
    else:
        info("没有 lineage_map.tsv，区室列显示为空")

    data = {
        "strata": [STRATUM_NAMES.get(int(c[1:]), c) for c in s_cols],
        "cell_types": prof.index.astype(str).tolist(),
        "profile": prof[s_cols].round(6).values.tolist(),
        "tai": prof["TAI"].round(6).tolist(),
        "lineage": {k: str(v) for k, v in lineage.items()},
    }
    # 有逐细胞 TAI 时，把标准误与细胞数一并带上（网页里显示误差棒）
    if "TAI_se" in prof.columns:
        # ⚠️ 不能 fillna(0)：细胞数过少的类型（如 n=1）其标准误在数学上无法估计，
        # 填 0 会被误读成「该类型的 TAI 极其精确」。这里保持为 null。
        se = prof["TAI_se"]
        data["tai_se"] = [None if pd.isna(v) else round(float(v), 6) for v in se]
    if "TAI_se_donor" in prof.columns:
        sed = prof["TAI_se_donor"]
        data["tai_se_donor"] = [None if pd.isna(v) else round(float(v), 6)
                                for v in sed]
    if "n_donors" in prof.columns:
        data["n_donors"] = prof["n_donors"].astype(int).tolist()
    if "n_cells" in prof.columns:
        data["n_cells"] = prof["n_cells"].astype(int).tolist()

    n_cells_total = int(prof["n_cells"].sum()) if "n_cells" in prof.columns else 0
    se_donor_med = (float(prof["TAI_se_donor"].median())
                    if "TAI_se_donor" in prof.columns else 0.0)

    html = (TEMPLATE
            .replace("__DATA__", json.dumps(data, ensure_ascii=False))
            .replace("__N_TYPES__", str(len(prof)))
            .replace("__N_STRATA__", str(len(s_cols)))
            .replace("__TAI_MIN__", f"{prof['TAI'].min():.2f}")
            .replace("__TAI_MAX__", f"{prof['TAI'].max():.2f}")
            .replace("__N_CELLS__", f"{n_cells_total:,}")
            .replace("__SE_DONOR__", f"{se_donor_med:.4f}")
            .replace("__N_GENES__", args.n_genes))

    out = args.out or os.path.join(VIS_DIR, "cellstrata.html")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write(html)
    ok(f"已写出：{out}  ({os.path.getsize(out)/1024:.0f} KB)")

    if args.open:
        webbrowser.open("file:///" + out.replace("\\", "/"))


if __name__ == "__main__":
    main()
