#!/usr/bin/env python3
"""
Build one self-contained ITSME multi-sample graph dashboard.

Reads each sample's graph_index.html and embeds the entire HTML directly into
one dashboard.html. Sample switching swaps the embedded report in-place.

Usage:
    python3 itsme_graph_dashboard_embedded.py .
    python3 itsme_graph_dashboard_embedded.py graphs
"""

import argparse
import base64
import html
import json
import re
from pathlib import Path
from urllib.parse import unquote, urlsplit


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
            samples.append({
                "directory": d.name,
                "label": display_name(d.name),
                "graph_html": graph_html.resolve(),
            })
    samples.sort(key=lambda x: natural_key(x["label"]))
    return samples


def rewrite_local_links(page, sample_dir):
    pattern = re.compile(
        r"(?P<prefix>\bhref\s*=\s*)(?P<quote>[\"'])(?P<url>.*?)(?P=quote)",
        flags=re.IGNORECASE,
    )

    def replace(match):
        url = match.group("url").strip()
        low = url.lower()

        if (
            not url
            or url.startswith("#")
            or low.startswith(("http://", "https://", "data:", "javascript:",
                               "mailto:", "file:"))
        ):
            return match.group(0)

        parsed = urlsplit(url)
        if parsed.scheme or parsed.netloc:
            return match.group(0)

        local_path = sample_dir / unquote(parsed.path)
        if not local_path.exists():
            return match.group(0)

        target = local_path.resolve().as_uri()
        if parsed.fragment:
            target += "#" + parsed.fragment

        q = match.group("quote")
        return match.group("prefix") + q + html.escape(target, quote=True) + q

    return pattern.sub(replace, page)


def encode_page(path):
    page = path.read_text(encoding="utf-8")
    page = rewrite_local_links(page, path.parent)
    return base64.b64encode(page.encode("utf-8")).decode("ascii")


def build_dashboard(samples, output, title):
    embedded = []
    for sample in samples:
        embedded.append({
            "label": sample["label"],
            "directory": sample["directory"],
            "html_b64": encode_page(sample["graph_html"]),
        })

    sample_json = json.dumps(embedded, ensure_ascii=False)

    page = '''<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>''' + html.escape(title) + '''</title>
<style>
:root{
  --bg:#07111d;
  --panel:#0a1624;
  --panel2:#102033;
  --line:#29445e;
  --line2:#3c5d7c;
  --text:#e7f0fa;
  --muted:#91a7bc;
  --blue:#79bdff;
}
*{box-sizing:border-box}
html,body{width:100%;height:100%;margin:0}
body{
  overflow:hidden;
  background:var(--bg);
  color:var(--text);
  font-family:system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;
}
.dashboard{
  width:100%;
  height:100%;
  display:grid;
  grid-template-rows:58px minmax(0,1fr);
}
.toolbar{
  position:relative;
  z-index:5;
  display:flex;
  align-items:center;
  gap:12px;
  padding:9px 18px;
  background:var(--panel);
  border-bottom:1px solid var(--line);
  box-shadow:0 4px 15px rgba(0,0,0,.25);
}
.toolbar label{
  color:var(--blue);
  font-size:11px;
  font-weight:800;
  letter-spacing:.13em;
  text-transform:uppercase;
  white-space:nowrap;
}
.toolbar select{
  width:min(460px,70vw);
  height:38px;
  padding:0 11px;
  border:1px solid var(--line2);
  border-radius:7px;
  background:var(--panel2);
  color:var(--text);
  font:700 13px system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;
  outline:none;
  cursor:pointer;
}
.toolbar select:hover,
.toolbar select:focus{border-color:var(--blue)}
.counter{
  margin-left:auto;
  color:var(--muted);
  font-size:12px;
  white-space:nowrap;
}
.viewer-wrap{
  position:relative;
  width:100%;
  min-height:0;
  overflow:hidden;
  background:#050b12;
}
#viewer{
  display:block;
  width:100%;
  height:100%;
  border:0;
  background:#050b12;
}
.loading{
  position:absolute;
  inset:0;
  display:grid;
  place-items:center;
  background:#07111d;
  color:var(--muted);
  font-size:13px;
  pointer-events:none;
  opacity:1;
  transition:opacity .12s linear;
}
.loading.hidden{opacity:0}
@media(max-width:650px){
  .toolbar{padding:8px 10px}
  .toolbar select{width:auto;min-width:0;flex:1}
  .counter{display:none}
}
</style>
</head>
<body>
<div class="dashboard">
  <div class="toolbar">
    <label for="sampleSelect">Sample</label>
    <select id="sampleSelect"></select>
    <span class="counter" id="counter"></span>
  </div>
  <div class="viewer-wrap">
    <iframe id="viewer" title="ITSME graph report"></iframe>
    <div class="loading" id="loading">Loading graph report…</div>
  </div>
</div>

<script>
const SAMPLES = ''' + sample_json + ''';

const selector = document.getElementById("sampleSelect");
const viewer = document.getElementById("viewer");
const counter = document.getElementById("counter");
const loading = document.getElementById("loading");

function decodeUtf8Base64(value){
  const binary = atob(value);
  const bytes = new Uint8Array(binary.length);
  for(let i=0;i<binary.length;i++) bytes[i] = binary.charCodeAt(i);
  return new TextDecoder("utf-8").decode(bytes);
}

function populateSelector(){
  SAMPLES.forEach((sample, i) => {
    const option = document.createElement("option");
    option.value = String(i);
    option.textContent = sample.label;
    selector.appendChild(option);
  });
}

function showSample(index){
  index = Number(index);
  if(!Number.isInteger(index) || index < 0 || index >= SAMPLES.length) index = 0;

  selector.value = String(index);
  counter.textContent = `${index + 1} / ${SAMPLES.length}`;
  loading.classList.remove("hidden");

  viewer.srcdoc = decodeUtf8Base64(SAMPLES[index].html_b64);
  history.replaceState(null, "", "#" + encodeURIComponent(SAMPLES[index].directory));
}

function initialIndex(){
  const wanted = decodeURIComponent(location.hash.replace(/^#/, ""));
  if(!wanted) return 0;
  const i = SAMPLES.findIndex(
    sample => sample.directory === wanted || sample.label === wanted
  );
  return i >= 0 ? i : 0;
}

viewer.addEventListener("load", () => loading.classList.add("hidden"));
selector.addEventListener("change", () => showSample(selector.value));

populateSelector();
showSample(initialIndex());
</script>
</body>
</html>
'''

    output.write_text(page, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(
        description="Embed multiple ITSME graph_index.html reports in one dashboard."
    )
    parser.add_argument(
        "graphs_dir",
        nargs="?",
        type=Path,
        default=Path("."),
        help="Directory containing one graph-summary directory per sample [.]"
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        help="Output HTML [GRAPHS_DIR/dashboard.html]"
    )
    parser.add_argument(
        "--title",
        default="ITSME multi-sample graph viewer",
        help="Dashboard title."
    )
    args = parser.parse_args()

    root = args.graphs_dir.resolve()
    if not root.is_dir():
        parser.error(f"Not a directory: {root}")

    samples = find_samples(root)
    if not samples:
        parser.error(
            f"No immediate child directories containing graph_index.html found in: {root}"
        )

    output = (args.output or (root / "dashboard.html")).resolve()
    build_dashboard(samples, output, args.title)

    size_mb = output.stat().st_size / (1024 * 1024)
    print(f"Created: {output}")
    print(f"Samples embedded: {len(samples)}")
    print(f"Dashboard size: {size_mb:.2f} MB")
    print("All graph_index.html content is embedded directly in dashboard.html.")


if __name__ == "__main__":
    main()
