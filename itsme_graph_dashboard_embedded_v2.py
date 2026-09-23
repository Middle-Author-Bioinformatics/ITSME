#!/usr/bin/env python3
import argparse
import base64
import hashlib
import html
import json
import re
from collections import defaultdict
from pathlib import Path


def natural_key(value):
    parts = re.split(r"(\d+)", str(value))
    return [int(x) if x.isdigit() else x.casefold() for x in parts]


def display_name(dirname):
    return dirname[6:] if dirname.lower().startswith("itsme_") else dirname


def find_samples(root):
    samples = []
    for d in root.iterdir():
        if not d.is_dir():
            continue
        graph_html = d / "graph_index.html"
        if graph_html.is_file():
            raw = graph_html.read_bytes()
            samples.append({
                "directory": d.name,
                "label": display_name(d.name),
                "graph_html": graph_html.resolve(),
                "sha256": hashlib.sha256(raw).hexdigest(),
            })
    samples.sort(key=lambda x: natural_key(x["label"]))
    return samples


def encode_report(sample):
    page = sample["graph_html"].read_text(encoding="utf-8")
    return base64.b64encode(page.encode("utf-8")).decode("ascii")


def report_duplicates(samples):
    groups = defaultdict(list)
    for s in samples:
        groups[s["sha256"]].append(s["label"])
    duplicates = {k: v for k, v in groups.items() if len(v) > 1}

    if not duplicates:
        print("Input check: all graph_index.html files are unique.")
        return

    print("WARNING: identical graph_index.html inputs detected:")
    for digest, labels in duplicates.items():
        print(f"  {digest[:12]} ({len(labels)} samples)")
        for label in labels:
            print(f"    {label}")
    print()


def build_dashboard(samples, output, title):
    payload = [{
        "label": s["label"],
        "directory": s["directory"],
        "sha256": s["sha256"],
        "html_b64": encode_report(s),
    } for s in samples]

    data = json.dumps(payload, ensure_ascii=False)

    page = f'''<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{html.escape(title)}</title>
<style>
:root{{--bg:#07111d;--panel:#0a1624;--panel2:#102033;--line:#29445e;
--line2:#3c5d7c;--text:#e7f0fa;--muted:#91a7bc;--blue:#79bdff}}
*{{box-sizing:border-box}}
html,body{{width:100%;height:100%;margin:0}}
body{{overflow:hidden;background:var(--bg);color:var(--text);
font-family:system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}}
.dashboard{{width:100%;height:100%;display:grid;grid-template-rows:58px minmax(0,1fr)}}
.toolbar{{display:flex;align-items:center;gap:12px;padding:9px 18px;background:var(--panel);
border-bottom:1px solid var(--line);box-shadow:0 4px 15px rgba(0,0,0,.25)}}
.toolbar label{{color:var(--blue);font-size:11px;font-weight:800;letter-spacing:.13em;
text-transform:uppercase;white-space:nowrap}}
.toolbar select{{width:min(470px,65vw);height:38px;padding:0 11px;border:1px solid var(--line2);
border-radius:7px;background:var(--panel2);color:var(--text);
font:700 13px system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;outline:none;cursor:pointer}}
.toolbar select:hover,.toolbar select:focus{{border-color:var(--blue)}}
.source-info{{color:#7890a7;font:600 10px ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;white-space:nowrap}}
.counter{{margin-left:auto;color:var(--muted);font-size:12px;white-space:nowrap}}
.viewer-wrap{{position:relative;width:100%;min-height:0;overflow:hidden;background:#050b12}}
.sample-frame{{display:block;width:100%;height:100%;border:0;background:#050b12}}
.loading{{position:absolute;inset:0;display:grid;place-items:center;background:#07111d;color:var(--muted);font-size:13px;z-index:2}}
</style>
</head>
<body>
<div class="dashboard">
  <div class="toolbar">
    <label for="sampleSelect">Sample</label>
    <select id="sampleSelect"></select>
    <span class="source-info" id="sourceInfo"></span>
    <span class="counter" id="counter"></span>
  </div>
  <div class="viewer-wrap" id="viewerWrap"></div>
</div>

<script>
const SAMPLES = {data};
const selector = document.getElementById("sampleSelect");
const viewerWrap = document.getElementById("viewerWrap");
const counter = document.getElementById("counter");
const sourceInfo = document.getElementById("sourceInfo");

function decodeB64(v) {{
  const b = atob(v);
  const bytes = new Uint8Array(b.length);
  for(let i=0;i<b.length;i++) bytes[i] = b.charCodeAt(i);
  return new TextDecoder("utf-8").decode(bytes);
}}

function populate() {{
  SAMPLES.forEach((s,i) => {{
    const o = document.createElement("option");
    o.value = String(i);
    o.textContent = s.label;
    selector.appendChild(o);
  }});
}}

function showSample(index) {{
  index = Number(index);
  if(!Number.isInteger(index) || index < 0 || index >= SAMPLES.length) index = 0;
  const s = SAMPLES[index];

  selector.value = String(index);
  counter.textContent = `${{index + 1}} / ${{SAMPLES.length}}`;
  sourceInfo.textContent = `source ${{s.sha256.slice(0,10)}}`;

  // Completely destroy the prior viewer and create a fresh one.
  viewerWrap.replaceChildren();

  const loading = document.createElement("div");
  loading.className = "loading";
  loading.textContent = `Loading ${{s.label}}…`;
  viewerWrap.appendChild(loading);

  const frame = document.createElement("iframe");
  frame.className = "sample-frame";
  frame.title = "ITSME graph report — " + s.label;

  frame.addEventListener("load", () => {{
    if(loading.isConnected) loading.remove();
  }}, {{once:true}});

  frame.srcdoc = decodeB64(s.html_b64);
  viewerWrap.appendChild(frame);

  history.replaceState(null, "", "#" + encodeURIComponent(s.directory));
}}

function initialIndex() {{
  const wanted = decodeURIComponent(location.hash.replace(/^#/, ""));
  if(!wanted) return 0;
  const i = SAMPLES.findIndex(s => s.directory === wanted || s.label === wanted);
  return i >= 0 ? i : 0;
}}

selector.addEventListener("change", () => showSample(selector.value));
populate();
showSample(initialIndex());
</script>
</body>
</html>
'''
    output.write_text(page, encoding="utf-8")


def main():
    ap = argparse.ArgumentParser(
        description="Build a self-contained ITSME multi-sample dashboard."
    )
    ap.add_argument("graphs_dir", nargs="?", type=Path, default=Path("."))
    ap.add_argument("-o", "--output", type=Path)
    ap.add_argument("--title", default="ITSME multi-sample graph viewer")
    args = ap.parse_args()

    root = args.graphs_dir.resolve()
    if not root.is_dir():
        ap.error(f"Not a directory: {root}")

    samples = find_samples(root)
    if not samples:
        ap.error(f"No child directories containing graph_index.html found in {root}")

    report_duplicates(samples)

    output = (args.output or (root / "dashboard.html")).resolve()
    build_dashboard(samples, output, args.title)

    print(f"Created: {output}")
    print(f"Samples embedded: {len(samples)}")
    print(f"Dashboard size: {output.stat().st_size / (1024*1024):.2f} MB")


if __name__ == "__main__":
    main()
