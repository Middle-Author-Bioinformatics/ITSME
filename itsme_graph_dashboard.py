#!/usr/bin/env python3
"""Build one ITSME dashboard from multiple graph-summary directories.

Expected layout:
    graphs/
      itsme_SAMPLE_A/graph_index.html
      itsme_SAMPLE_B/graph_index.html
      ...

Usage:
    python3 itsme_graph_dashboard.py graphs
    python3 itsme_graph_dashboard.py graphs -o graphs/dashboard.html
    python3 itsme_graph_dashboard.py graphs --title "Arctic SKQ23 ITSME graphs"

The dashboard references each existing graph_index.html by relative path, so keep
it with the sample directories. Python >=3.9; standard library only.
"""

import argparse
import html
import json
import re
from pathlib import Path
from urllib.parse import quote


def natural_key(value):
    return [int(x) if x.isdigit() else x.casefold()
            for x in re.split(r"(\d+)", str(value))]


def display_name(dirname):
    return dirname[6:] if dirname.lower().startswith("itsme_") else dirname


def find_samples(root, recursive=False):
    root = root.resolve()
    if recursive:
        indexes = list(root.rglob("graph_index.html"))
    else:
        indexes = [p / "graph_index.html" for p in root.iterdir()
                   if p.is_dir() and (p / "graph_index.html").is_file()]

    samples = [{
        "dir_name": p.parent.name,
        "label": display_name(p.parent.name),
        "index": p.resolve(),
    } for p in indexes]
    samples.sort(key=lambda x: natural_key(x["label"]))
    return samples


def relative_url(path, output_html):
    rel = path.relative_to(output_html.parent.resolve())
    return "/".join(quote(part) for part in rel.parts)


def build_dashboard(samples, output_html, title):
    output_html = output_html.resolve()
    browser_samples = [{
        "label": s["label"],
        "directory": s["dir_name"],
        "url": relative_url(s["index"], output_html),
    } for s in samples]

    samples_json = json.dumps(browser_samples, ensure_ascii=False).replace("</", "<\\/")

    page = '''<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>__TITLE__</title>
<style>
:root{
  --bg:#07111d;--bg2:#0a1624;--panel:#0d1a2a;--line:#263b52;
  --line2:#35516e;--text:#e4edf7;--muted:#93a7bc;--blue:#78bdff;
  --blue2:#a8d5ff;--active:#173a59;
}
*{box-sizing:border-box}html,body{height:100%;margin:0}
body{
  background:radial-gradient(circle at 20% -10%,rgba(48,126,190,.12),transparent 34%),
             linear-gradient(180deg,var(--bg2),var(--bg));
  color:var(--text);font:14px/1.4 system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;
  overflow:hidden
}
.shell{height:100vh;display:grid;grid-template-rows:auto 1fr}
.topbar{z-index:10;background:rgba(7,17,29,.97);border-bottom:1px solid var(--line);
  box-shadow:0 6px 20px rgba(0,0,0,.2);backdrop-filter:blur(14px)}
.topline{min-height:70px;padding:12px 22px;display:flex;align-items:center;gap:18px}
.brand{min-width:245px}.eyebrow{color:var(--blue);font-size:10px;font-weight:800;
  letter-spacing:.18em;text-transform:uppercase;margin-bottom:2px}
h1{margin:0;font-size:19px;line-height:1.15}
.sample-controls{min-width:0;flex:1;display:flex;align-items:center;gap:8px}
.nav-button,.sample-select,.open-link{border:1px solid var(--line2);background:#0e1c2c;
  color:var(--text);border-radius:8px;min-height:37px}
.nav-button{width:40px;cursor:pointer;font-size:17px;font-weight:800}
.nav-button:hover,.open-link:hover{border-color:#5a82a8;background:#14273b}
.nav-button:disabled{opacity:.32;cursor:default}
.sample-select{flex:1;min-width:180px;max-width:520px;padding:0 35px 0 12px;font-weight:700;
  font-size:14px;outline:none}.sample-count{color:var(--muted);white-space:nowrap;font-size:12px}
.open-link{display:inline-flex;align-items:center;padding:0 11px;color:var(--blue2);
  text-decoration:none;white-space:nowrap;font-size:12px;font-weight:700}
.sample-tabs{padding:0 22px 9px;display:flex;gap:6px;overflow-x:auto;scrollbar-width:thin}
.sample-tab{flex:0 0 auto;border:1px solid #263d55;background:#0c1928;color:#a8b9ca;
  border-radius:999px;padding:5px 10px;cursor:pointer;font-size:11px;font-weight:700}
.sample-tab:hover{color:#d9eaff;border-color:#486d92}
.sample-tab.active{background:var(--active);border-color:#65a7e3;color:#dff0ff}
.viewer{position:relative;min-height:0;background:#050c14}
iframe{display:block;width:100%;height:100%;border:0;background:#07111d}
.loading{position:absolute;inset:0;display:flex;align-items:center;justify-content:center;
  background:#07111d;color:var(--muted);pointer-events:none;transition:opacity .12s ease}
.loading.hidden{opacity:0}
@media(max-width:760px){.topline{padding:10px 12px;gap:8px;flex-wrap:wrap}.brand{min-width:100%;width:100%}
  .sample-controls{width:100%}.sample-tabs{padding:0 12px 8px}.sample-count,.open-link{display:none}}
</style>
</head>
<body>
<div class="shell">
<header class="topbar">
  <div class="topline">
    <div class="brand"><div class="eyebrow">ITSME graph dashboard</div><h1>__TITLE__</h1></div>
    <div class="sample-controls">
      <button class="nav-button" id="previous" title="Previous sample">‹</button>
      <select class="sample-select" id="sampleSelect" aria-label="Select sample"></select>
      <button class="nav-button" id="next" title="Next sample">›</button>
      <span class="sample-count" id="sampleCount"></span>
      <a class="open-link" id="openStandalone" target="_blank" rel="noopener">Open sample separately ↗</a>
    </div>
  </div>
  <div class="sample-tabs" id="sampleTabs"></div>
</header>
<main class="viewer">
  <div class="loading" id="loading">Loading sample…</div>
  <iframe id="sampleFrame" title="ITSME sample graph report"></iframe>
</main>
</div>
<script>
const SAMPLES = __SAMPLES__;
const select=document.getElementById('sampleSelect');
const frame=document.getElementById('sampleFrame');
const previous=document.getElementById('previous');
const next=document.getElementById('next');
const count=document.getElementById('sampleCount');
const tabs=document.getElementById('sampleTabs');
const loading=document.getElementById('loading');
const openStandalone=document.getElementById('openStandalone');
let current=0;

function makeControls(){
  SAMPLES.forEach((sample,i)=>{
    const option=document.createElement('option');option.value=i;option.textContent=sample.label;select.appendChild(option);
    const tab=document.createElement('button');tab.className='sample-tab';tab.type='button';tab.dataset.index=i;
    tab.textContent=sample.label;tab.addEventListener('click',()=>showSample(i));tabs.appendChild(tab);
  });
}
function showSample(index,updateHash=true){
  index=Math.max(0,Math.min(SAMPLES.length-1,Number(index)));current=index;
  const sample=SAMPLES[index];select.value=String(index);loading.classList.remove('hidden');
  frame.src=sample.url;openStandalone.href=sample.url;count.textContent=`${index+1} / ${SAMPLES.length}`;
  previous.disabled=index===0;next.disabled=index===SAMPLES.length-1;
  document.querySelectorAll('.sample-tab').forEach((tab,i)=>tab.classList.toggle('active',i===index));
  const active=document.querySelector(`.sample-tab[data-index="${index}"]`);
  if(active) active.scrollIntoView({behavior:'smooth',block:'nearest',inline:'center'});
  if(updateHash) history.replaceState(null,'','#'+encodeURIComponent(sample.directory));
}
function initialSample(){
  const wanted=decodeURIComponent(location.hash.replace(/^#/,''));
  if(!wanted)return 0;
  const i=SAMPLES.findIndex(s=>s.directory===wanted||s.label===wanted);return i>=0?i:0;
}
select.addEventListener('change',()=>showSample(Number(select.value)));
previous.addEventListener('click',()=>showSample(current-1));
next.addEventListener('click',()=>showSample(current+1));
frame.addEventListener('load',()=>loading.classList.add('hidden'));
document.addEventListener('keydown',event=>{
  if(event.target&&['INPUT','SELECT','TEXTAREA'].includes(event.target.tagName))return;
  if(event.key==='ArrowLeft'&&current>0)showSample(current-1);
  else if(event.key==='ArrowRight'&&current<SAMPLES.length-1)showSample(current+1);
});
makeControls();showSample(initialSample(),false);
</script>
</body>
</html>'''

    page = page.replace("__TITLE__", html.escape(title)).replace("__SAMPLES__", samples_json)
    output_html.parent.mkdir(parents=True, exist_ok=True)
    output_html.write_text(page, encoding="utf-8")


def main():
    p = argparse.ArgumentParser(description="Build a multi-sample ITSME graph dashboard.")
    p.add_argument("graphs_dir", type=Path,
                   help="Directory containing one graph-summary directory per sample")
    p.add_argument("-o", "--output", type=Path,
                   help="Output HTML [default: GRAPHS_DIR/dashboard.html]")
    p.add_argument("--title", default="ITSME multi-sample graph viewer")
    p.add_argument("--recursive", action="store_true",
                   help="Recursively search for graph_index.html")
    a = p.parse_args()

    root = a.graphs_dir.resolve()
    if not root.is_dir():
        p.error(f"Not a directory: {root}")

    output = (a.output or (root / "dashboard.html")).resolve()
    samples = find_samples(root, recursive=a.recursive)
    if not samples:
        p.error(f"No graph_index.html files found under: {root}")

    for sample in samples:
        try:
            sample["index"].relative_to(output.parent)
        except ValueError:
            p.error("Place the dashboard in the graph directory or one of its parents so relative links remain portable.")

    build_dashboard(samples, output, a.title)
    print(f"Dashboard: {output}")
    print(f"Samples:   {len(samples)}")
    for s in samples:
        print(f"  {s['label']}: {s['index']}")


if __name__ == "__main__":
    main()
