#!/usr/bin/env python3
"""
Create a simple multi-sample ITSME graph dashboard WITHOUT iframes.

For every immediate child directory containing graph_index.html, creates:
    SAMPLE/dashboard_view.html

That file is simply the original graph_index.html plus a small sticky sample
selector at the top. Because it remains in the same directory as the original,
all relative FASTG/TSV links continue to work unchanged.

Also creates:
    ROOT/dashboard.html

Usage:
    python3 itsme_graph_dashboard_simple.py .
    python3 itsme_graph_dashboard_simple.py /path/to/graphs

Python >= 3.9; standard library only.
"""

import argparse
import html
import re
from pathlib import Path
from urllib.parse import quote

VIEW_NAME = "dashboard_view.html"


def natural_key(value):
    parts = re.split(r"(\d+)", str(value))
    return [int(x) if x.isdigit() else x.casefold() for x in parts]


def display_name(dirname):
    return dirname[6:] if dirname.lower().startswith("itsme_") else dirname


def url_path(parts):
    return "/".join(quote(str(p)) for p in parts)


def find_samples(root):
    samples = []
    for d in root.iterdir():
        if not d.is_dir():
            continue
        source = d / "graph_index.html"
        if source.is_file():
            samples.append({
                "directory": d.name,
                "label": display_name(d.name),
                "source": source,
                "view": d / VIEW_NAME,
            })
    samples.sort(key=lambda x: natural_key(x["label"]))
    return samples


def selector_bar(samples, current_dir):
    option_html = []
    current_label = current_dir
    for s in samples:
        target = url_path(["..", s["directory"], VIEW_NAME])
        selected = " selected" if s["directory"] == current_dir else ""
        if selected:
            current_label = s["label"]
        option_html.append(
            '<option value="{}"{}>{}</option>'.format(
                html.escape(target, quote=True), selected, html.escape(s["label"])
            )
        )

    return """
<style id="itsme-dashboard-style">
#itsme-sample-nav {
  position: sticky;
  top: 0;
  z-index: 999999;
  display: flex;
  align-items: center;
  gap: 12px;
  min-height: 54px;
  padding: 8px 18px;
  background: #0a1420;
  border-bottom: 1px solid #294057;
  box-shadow: 0 4px 16px rgba(0,0,0,.28);
  color: #dce9f6;
  font-family: system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;
}
#itsme-sample-nav .itsme-nav-label {
  color: #79bdff;
  font-size: 11px;
  font-weight: 800;
  letter-spacing: .12em;
  text-transform: uppercase;
  white-space: nowrap;
}
#itsme-sample-nav select {
  min-width: 300px;
  max-width: 560px;
  height: 36px;
  padding: 0 34px 0 10px;
  border: 1px solid #3a5773;
  border-radius: 7px;
  background: #101f30;
  color: #e7f1fb;
  font: 700 13px system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;
  outline: none;
}
#itsme-sample-nav select:hover,
#itsme-sample-nav select:focus {
  border-color: #72b9f5;
}
#itsme-sample-nav .itsme-current {
  margin-left: auto;
  color: #91a6bb;
  font-size: 12px;
  white-space: nowrap;
}
@media(max-width:700px) {
  #itsme-sample-nav { padding: 8px 10px; }
  #itsme-sample-nav select { min-width: 0; flex: 1; }
  #itsme-sample-nav .itsme-current { display: none; }
}
</style>
<div id="itsme-sample-nav">
  <span class="itsme-nav-label">Sample</span>
  <select id="itsme-sample-select" aria-label="Select ITSME sample">
    OPTIONS_PLACEHOLDER
  </select>
  <span class="itsme-current">CURRENT_PLACEHOLDER</span>
</div>
<script id="itsme-dashboard-script">
(function() {
  var selector = document.getElementById('itsme-sample-select');
  if (!selector) return;
  selector.addEventListener('change', function() {
    window.location.href = this.value;
  });
})();
</script>
""".replace("OPTIONS_PLACEHOLDER", "".join(option_html)).replace(
        "CURRENT_PLACEHOLDER", html.escape(current_label)
    )


def inject_nav(source_html, nav_html):
    body = re.search(r"<body\b[^>]*>", source_html, flags=re.IGNORECASE)
    if body:
        return source_html[:body.end()] + nav_html + source_html[body.end():]
    return nav_html + source_html


def write_sample_views(samples):
    for sample in samples:
        original = sample["source"].read_text(encoding="utf-8")
        updated = inject_nav(original, selector_bar(samples, sample["directory"]))
        sample["view"].write_text(updated, encoding="utf-8")


def write_landing(root, samples, title):
    option_html = ['<option value="">Choose a sample...</option>']
    for s in samples:
        target = url_path([s["directory"], VIEW_NAME])
        option_html.append(
            '<option value="{}">{}</option>'.format(
                html.escape(target, quote=True), html.escape(s["label"])
            )
        )

    page = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>TITLE_PLACEHOLDER</title>
<style>
*{box-sizing:border-box}
html,body{height:100%;margin:0}
body{
  display:grid;
  place-items:center;
  background:#07111d;
  color:#e5eef8;
  font-family:system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;
}
.card{
  width:min(620px,calc(100vw - 32px));
  padding:28px;
  border:1px solid #294057;
  border-radius:12px;
  background:#0d1a29;
  box-shadow:0 18px 55px rgba(0,0,0,.3);
}
.eyebrow{
  margin-bottom:5px;
  color:#79bdff;
  font-size:11px;
  font-weight:800;
  letter-spacing:.14em;
  text-transform:uppercase;
}
h1{margin:0 0 8px;font-size:24px}
p{margin:0 0 20px;color:#96a9bc}
select{
  width:100%;
  height:44px;
  padding:0 12px;
  border:1px solid #3a5773;
  border-radius:8px;
  background:#101f30;
  color:#e7f1fb;
  font-size:14px;
  font-weight:700;
}
</style>
</head>
<body>
<div class="card">
  <div class="eyebrow">ITSME graphs</div>
  <h1>TITLE_PLACEHOLDER</h1>
  <p>Select a sample. It opens as a normal graph page, with the same sample selector at the top.</p>
  <select id="sample">
    OPTIONS_PLACEHOLDER
  </select>
</div>
<script>
document.getElementById('sample').addEventListener('change', function() {
  if (this.value) window.location.href = this.value;
});
</script>
</body>
</html>
"""
    page = page.replace("TITLE_PLACEHOLDER", html.escape(title))
    page = page.replace("OPTIONS_PLACEHOLDER", "".join(option_html))
    (root / "dashboard.html").write_text(page, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(
        description="Create a simple no-iframe dashboard for ITSME graph_index.html reports."
    )
    parser.add_argument(
        "graphs_dir", nargs="?", default=".", type=Path,
        help="Directory containing one graph-summary directory per sample [.]"
    )
    parser.add_argument(
        "--title", default="ITSME multi-sample graph viewer",
        help="Title shown on dashboard.html"
    )
    args = parser.parse_args()

    root = args.graphs_dir.resolve()
    if not root.is_dir():
        parser.error("Not a directory: {}".format(root))

    samples = find_samples(root)
    if not samples:
        parser.error(
            "No immediate child directories containing graph_index.html found in: {}".format(root)
        )

    write_sample_views(samples)
    write_landing(root, samples, args.title)

    print("Created: {}".format(root / "dashboard.html"))
    print("Samples: {}".format(len(samples)))
    for sample in samples:
        print("  {} -> {}".format(sample["label"], sample["view"]))


if __name__ == "__main__":
    main()
