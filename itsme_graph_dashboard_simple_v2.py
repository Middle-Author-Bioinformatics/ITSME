#!/usr/bin/env python3
'''
Simple multi-sample ITSME graph selector.

Creates:
  ROOT/dashboard.html
  ROOT/<sample>/dashboard_view.html

The dropdown uses absolute file:// URLs, avoiding iframe and relative-path issues.

Usage:
    python3 itsme_graph_dashboard_simple_v2.py .
'''

import argparse
import html
import re
from pathlib import Path

VIEW_NAME = "dashboard_view.html"


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
        source = d / "graph_index.html"
        if source.is_file():
            samples.append({
                "directory": d.name,
                "label": display_name(d.name),
                "source": source.resolve(),
                "view": (d / VIEW_NAME).resolve(),
            })
    samples.sort(key=lambda x: natural_key(x["label"]))
    return samples


def nav_html(samples, current):
    options = []
    for sample in samples:
        selected = " selected" if sample["directory"] == current["directory"] else ""
        url = sample["view"].as_uri()
        options.append(
            f'<option value="{html.escape(url, quote=True)}"{selected}>'
            f'{html.escape(sample["label"])}</option>'
        )

    options_html = "".join(options)

    return f'''
<style id="itsme-multisample-nav-style">
#itsme-multisample-nav {{
  position: sticky;
  top: 0;
  z-index: 1000000;
  display: flex;
  align-items: center;
  gap: 12px;
  height: 58px;
  padding: 9px 18px;
  background: #0a1522;
  border-bottom: 1px solid #29425a;
  box-shadow: 0 4px 14px rgba(0,0,0,.28);
  font-family: system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;
}}
#itsme-multisample-nav label {{
  color: #79bdff;
  font-size: 11px;
  font-weight: 800;
  letter-spacing: .13em;
  text-transform: uppercase;
}}
#itsme-multisample-nav select {{
  width: min(430px,70vw);
  height: 38px;
  padding: 0 11px;
  border: 1px solid #3c5a77;
  border-radius: 7px;
  background: #102033;
  color: #e8f2fc;
  font: 700 13px system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;
  outline: none;
  cursor: pointer;
}}
#itsme-multisample-nav select:hover,
#itsme-multisample-nav select:focus {{
  border-color: #79bdff;
}}
</style>

<div id="itsme-multisample-nav">
  <label for="itsme-sample-selector">Sample</label>
  <select id="itsme-sample-selector"
          onchange="if(this.value){{ window.location.href = this.value; }}">
    {options_html}
  </select>
</div>
'''


def inject_after_body(source_html, insertion):
    match = re.search(r"<body\b[^>]*>", source_html, flags=re.IGNORECASE)
    if match:
        return source_html[:match.end()] + "\n" + insertion + "\n" + source_html[match.end():]
    return insertion + "\n" + source_html


def write_views(samples):
    for current in samples:
        original = current["source"].read_text(encoding="utf-8")
        output = inject_after_body(original, nav_html(samples, current))
        current["view"].write_text(output, encoding="utf-8")


def write_landing(root, samples, title):
    options = ['<option value="">Choose a sample...</option>']
    for sample in samples:
        options.append(
            f'<option value="{html.escape(sample["view"].as_uri(), quote=True)}">'
            f'{html.escape(sample["label"])}</option>'
        )
    options_html = "".join(options)

    page = f'''<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{html.escape(title)}</title>
<style>
*{{box-sizing:border-box}}
html,body{{height:100%;margin:0}}
body{{
  display:grid;
  place-items:center;
  background:#07111d;
  color:#e7f0fa;
  font-family:system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;
}}
main{{
  width:min(600px,calc(100vw - 30px));
  padding:28px;
  background:#0d1a29;
  border:1px solid #29425a;
  border-radius:12px;
}}
small{{
  display:block;
  margin-bottom:5px;
  color:#79bdff;
  font-weight:800;
  letter-spacing:.13em;
  text-transform:uppercase;
}}
h1{{margin:0 0 18px;font-size:24px}}
select{{
  width:100%;
  height:44px;
  padding:0 12px;
  border:1px solid #3c5a77;
  border-radius:8px;
  background:#102033;
  color:#e7f0fa;
  font-size:14px;
  font-weight:700;
}}
</style>
</head>
<body>
<main>
  <small>ITSME graphs</small>
  <h1>{html.escape(title)}</h1>
  <select onchange="if(this.value){{ window.location.href = this.value; }}">
    {options_html}
  </select>
</main>
</body>
</html>
'''
    (root / "dashboard.html").write_text(page, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(
        description="Build a simple no-iframe ITSME multi-sample graph selector."
    )
    parser.add_argument(
        "graphs_dir",
        nargs="?",
        default=".",
        type=Path,
        help="Directory containing sample graph-summary directories [.]"
    )
    parser.add_argument(
        "--title",
        default="ITSME multi-sample graph viewer",
        help="Landing-page title."
    )
    args = parser.parse_args()

    root = args.graphs_dir.resolve()
    if not root.is_dir():
        parser.error(f"Not a directory: {root}")

    samples = find_samples(root)
    if not samples:
        parser.error(
            f"No immediate child directories containing graph_index.html found in {root}"
        )

    write_views(samples)
    write_landing(root, samples, args.title)

    print(f"Created landing page: {root / 'dashboard.html'}")
    print(f"Updated {len(samples)} sample dashboard_view.html files.")
    for sample in samples:
        print(f"  {sample['label']} -> {sample['view'].as_uri()}")
    print()
    print("NOTE: links use absolute file:// URLs. Re-run this script if the graphs folder moves.")


if __name__ == "__main__":
    main()
