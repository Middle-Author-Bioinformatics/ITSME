#!/usr/bin/env python3
"""
Build a single standalone dashboard for browsing many ITSME graph_index.html files.

By default, sample folders whose report.json indicates zero mapped sequences or zero
graph reports are omitted from the dashboard.

Expected layout:

ROOT/
  itsme_SAMPLE_1/
    graph_index.html
    report.json
    *.summary.fastg
    *.sequences.tsv
    ...
  itsme_SAMPLE_2/
    graph_index.html
    report.json
    ...

The generated dashboard embeds every selected graph_index.html directly into one
master HTML and swaps the entire viewer inside a single iframe when the sample
dropdown changes.

Python >= 3.9; standard library only.
"""

from __future__ import annotations

import argparse
import base64
import hashlib
import html
import json
import os
import re
import sys
from pathlib import Path
from urllib.parse import quote


def natural_key(text: str):
    """Natural sort: sample2 before sample10."""
    return [
        int(piece) if piece.isdigit() else piece.lower()
        for piece in re.split(r"(\d+)", text)
    ]


def add_base_href(page: str, base_href: str) -> str:
    """
    Inject a <base> tag so links inside iframe srcdoc still resolve back to the
    original sample directory (e.g. FASTG and TSV links).
    """
    tag = f'<base href="{html.escape(base_href, quote=True)}">'
    m = re.search(r"<head\b[^>]*>", page, flags=re.IGNORECASE)
    if m:
        return page[: m.end()] + tag + page[m.end() :]
    return tag + page


def load_report(report_path: Path):
    try:
        return json.loads(report_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError, UnicodeError):
        return None


def is_empty_report(report: dict | None) -> bool:
    """
    Treat a sample as empty when report.json clearly indicates nothing was mapped
    or nothing was rendered.
    """
    if not isinstance(report, dict):
        return False
    mapped = report.get("mapped_sequences")
    graphs = report.get("graphs")
    if isinstance(mapped, int) and mapped == 0:
        return True
    if isinstance(graphs, list) and len(graphs) == 0:
        return True
    return False


def discover_samples(root: Path, viewer_name: str, pattern: str | None, include_empty: bool):
    included = []
    skipped_empty = []
    skipped_missing_viewer = []

    for child in root.iterdir():
        if not child.is_dir():
            continue
        if pattern and not child.match(pattern):
            continue

        viewer = child / viewer_name
        if not viewer.is_file():
            skipped_missing_viewer.append(child.name)
            continue

        report_path = child / "report.json"
        report = load_report(report_path) if report_path.is_file() else None

        if not include_empty and is_empty_report(report):
            skipped_empty.append({
                "name": child.name,
                "mapped_sequences": report.get("mapped_sequences"),
                "requested_sequences": report.get("requested_sequences"),
                "graph_count": len(report.get("graphs", [])) if isinstance(report.get("graphs"), list) else None,
            })
            continue

        included.append((child.name, child, viewer, report))

    included.sort(key=lambda x: natural_key(x[0]))
    skipped_empty.sort(key=lambda x: natural_key(x["name"]))
    return included, skipped_empty, skipped_missing_viewer


def build_dashboard(root: Path, output: Path, viewer_name: str, pattern: str | None, include_empty: bool):
    samples, skipped_empty, skipped_missing_viewer = discover_samples(root, viewer_name, pattern, include_empty)
    if not samples:
        raise ValueError(
            "No eligible sample folders were found. "
            "If you want to include empty viewers, use --include-empty."
        )

    payload = []
    hash_groups: dict[str, list[str]] = {}

    for sample_name, sample_dir, viewer_path, report in samples:
        page = viewer_path.read_text(encoding="utf-8")
        sha = hashlib.sha256(page.encode("utf-8")).hexdigest()
        hash_groups.setdefault(sha, []).append(sample_name)

        # Make links inside the embedded child HTML resolve to that sample folder.
        rel_dir = os.path.relpath(sample_dir, output.parent)
        base_href = quote(Path(rel_dir).as_posix().rstrip("/") + "/", safe="/._-")
        page = add_base_href(page, base_href)

        payload.append(
            {
                "name": sample_name,
                "sha": sha[:12],
                "bytes": len(page.encode("utf-8")),
                "mapped_sequences": report.get("mapped_sequences") if isinstance(report, dict) else None,
                "graph_count": len(report.get("graphs", [])) if isinstance(report, dict) and isinstance(report.get("graphs"), list) else None,
                "b64": base64.b64encode(page.encode("utf-8")).decode("ascii"),
            }
        )

    options = "\n".join(
        f'<option value="{i}">{html.escape(item["name"])}</option>'
        for i, item in enumerate(payload)
    )

    data_json = json.dumps(payload, separators=(",", ":"))

    master = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>ITSME multi-sample graph dashboard</title>
<style>
  :root {{
    color-scheme: dark;
    --bg:#070b12;
    --panel:#0d1420;
    --line:#26354a;
    --ink:#e8eef7;
    --muted:#8fa0b5;
    --accent:#6fb5ff;
  }}
  * {{ box-sizing:border-box; }}
  html,body {{ margin:0; height:100%; background:var(--bg); color:var(--ink);
               font:14px/1.4 system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; }}
  body {{ display:flex; flex-direction:column; overflow:hidden; }}
  .bar {{
    flex:0 0 auto;
    display:flex;
    align-items:center;
    gap:10px;
    padding:10px 14px;
    background:var(--panel);
    border-bottom:1px solid var(--line);
  }}
  .title {{
    font-weight:700;
    white-space:nowrap;
    margin-right:4px;
  }}
  select,button {{
    background:#121e2d;
    color:var(--ink);
    border:1px solid #33445b;
    border-radius:7px;
    padding:7px 10px;
    font:inherit;
  }}
  select {{
    min-width:320px;
    max-width:55vw;
  }}
  button {{ cursor:pointer; }}
  button:hover {{ background:#1b2b40; }}
  button:disabled {{ opacity:.4; cursor:default; }}
  .meta {{
    color:var(--muted);
    font-size:12px;
    margin-left:auto;
    white-space:nowrap;
  }}
  .hash {{
    font-family:ui-monospace,SFMono-Regular,Menlo,monospace;
    color:#b9c9dc;
  }}
  iframe {{
    flex:1 1 auto;
    width:100%;
    min-height:0;
    border:0;
    background:var(--bg);
  }}
  @media(max-width:900px) {{
    .title {{ display:none; }}
    select {{ min-width:0; max-width:none; flex:1; }}
    .meta {{ display:none; }}
  }}
</style>
</head>
<body>

<div class="bar">
  <div class="title">ITSME samples</div>
  <button id="prev" type="button" title="Previous sample">←</button>
  <select id="sampleSelect" aria-label="Choose sample">
    {options}
  </select>
  <button id="next" type="button" title="Next sample">→</button>
  <div class="meta">
    <span id="counter"></span>
    &nbsp;·&nbsp;
    graphs <span id="graphCount"></span>
    &nbsp;·&nbsp;
    mapped <span id="mappedCount"></span>
    &nbsp;·&nbsp;
    HTML <span class="hash" id="hash"></span>
  </div>
</div>

<iframe id="viewer" title="ITSME graph viewer"></iframe>

<script>
const SAMPLES = {data_json};

const select = document.getElementById("sampleSelect");
const viewer = document.getElementById("viewer");
const prev = document.getElementById("prev");
const next = document.getElementById("next");
const counter = document.getElementById("counter");
const hash = document.getElementById("hash");
const graphCount = document.getElementById("graphCount");
const mappedCount = document.getElementById("mappedCount");

function decodeUtf8Base64(b64) {{
  const binary = atob(b64);
  const bytes = new Uint8Array(binary.length);
  for (let i = 0; i < binary.length; i++) bytes[i] = binary.charCodeAt(i);
  return new TextDecoder("utf-8").decode(bytes);
}}

function loadSample(index) {{
  index = Math.max(0, Math.min(SAMPLES.length - 1, Number(index) || 0));
  select.value = String(index);

  const item = SAMPLES[index];

  // Replace the ENTIRE child document. Nothing is shared between samples.
  viewer.srcdoc = decodeUtf8Base64(item.b64);

  counter.textContent = `${{index + 1}} / ${{SAMPLES.length}}`;
  hash.textContent = item.sha;
  graphCount.textContent = item.graph_count ?? "NA";
  mappedCount.textContent = item.mapped_sequences ?? "NA";
  prev.disabled = index === 0;
  next.disabled = index === SAMPLES.length - 1;

  document.title = `${{item.name}} — ITSME dashboard`;
}}

select.addEventListener("change", () => loadSample(select.value));

prev.addEventListener("click", () => {{
  loadSample(Number(select.value) - 1);
}});

next.addEventListener("click", () => {{
  loadSample(Number(select.value) + 1);
}});

document.addEventListener("keydown", (event) => {{
  if (event.target === select) return;
  if (event.key === "ArrowLeft") loadSample(Number(select.value) - 1);
  if (event.key === "ArrowRight") loadSample(Number(select.value) + 1);
}});

loadSample(0);
</script>
</body>
</html>
"""

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(master, encoding="utf-8")

    print(f"Wrote: {output}")
    print(f"Included {len(payload)} sample viewer(s):")
    for item in payload:
        extra = []
        if item["graph_count"] is not None:
            extra.append(f"graphs={item['graph_count']}")
        if item["mapped_sequences"] is not None:
            extra.append(f"mapped={item['mapped_sequences']}")
        extra_txt = "  " + "  ".join(extra) if extra else ""
        print(f"  {item['name']}  sha256={item['sha']}  embedded={item['bytes']:,} bytes{extra_txt}")

    if skipped_empty:
        print(f"\nSkipped {len(skipped_empty)} sample(s) with no reconstructed sequences / no graph reports:")
        for row in skipped_empty:
            detail = []
            if row["requested_sequences"] is not None:
                detail.append(f"requested={row['requested_sequences']}")
            if row["mapped_sequences"] is not None:
                detail.append(f"mapped={row['mapped_sequences']}")
            if row["graph_count"] is not None:
                detail.append(f"graphs={row['graph_count']}")
            print(f"  {row['name']}" + (f"  ({', '.join(detail)})" if detail else ""))

    duplicates = [names for names in hash_groups.values() if len(names) > 1]
    if duplicates:
        print("\nWARNING: some INCLUDED graph_index.html files are byte-for-byte identical:")
        for names in duplicates:
            print("  " + ", ".join(names))
        print(
            "If those samples are expected to differ, inspect the per-sample "
            "graph_index.html generation before blaming the dashboard."
        )

    if skipped_missing_viewer:
        print(f"\nIgnored {len(skipped_missing_viewer)} subfolder(s) with no {viewer_name}:")
        for name in sorted(skipped_missing_viewer, key=natural_key):
            print(f"  {name}")


def main():
    p = argparse.ArgumentParser(
        description="Build a single embedded dashboard from many ITSME graph_index.html outputs."
    )
    p.add_argument(
        "root",
        type=Path,
        help="Folder whose immediate subfolders are ITSME summary-graph output folders.",
    )
    p.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Output HTML. Default: ROOT/itsme_multi_sample_dashboard.html",
    )
    p.add_argument(
        "--viewer",
        default="graph_index.html",
        help="Viewer filename inside each sample folder (default: graph_index.html).",
    )
    p.add_argument(
        "--pattern",
        default=None,
        help="Optional subfolder glob, e.g. 'itsme_*'.",
    )
    p.add_argument(
        "--include-empty",
        action="store_true",
        help="Include samples even when report.json shows zero mapped sequences or zero graphs.",
    )
    args = p.parse_args()

    root = args.root.expanduser().resolve()
    if not root.is_dir():
        p.error(f"Not a directory: {root}")

    output = (
        args.output.expanduser().resolve()
        if args.output
        else root / "itsme_multi_sample_dashboard.html"
    )

    if output.exists() and output.is_dir():
        p.error(f"Output is a directory: {output}")

    try:
        build_dashboard(root, output, args.viewer, args.pattern, args.include_empty)
    except (OSError, ValueError, UnicodeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
