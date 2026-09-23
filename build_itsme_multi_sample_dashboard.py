#!/usr/bin/env python3
"""
Build a single standalone dashboard for browsing many ITSME graph_index.html files.

Expected layout:

ROOT/
  itsme_SAMPLE_1/
    graph_index.html
    *.summary.fastg
    *.sequences.tsv
    ...
  itsme_SAMPLE_2/
    graph_index.html
    ...

The generated dashboard embeds every graph_index.html directly into one master HTML
and swaps the entire viewer inside a single iframe when the sample dropdown changes.

Python >= 3.9; standard library only.
"""

from __future__ import annotations

import argparse
import base64
import hashlib
import html
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


def discover_samples(root: Path, viewer_name: str, pattern: str | None):
    samples = []
    for child in root.iterdir():
        if not child.is_dir():
            continue
        if pattern and not child.match(pattern):
            continue
        viewer = child / viewer_name
        if viewer.is_file():
            samples.append((child.name, child, viewer))
    samples.sort(key=lambda x: natural_key(x[0]))
    return samples


def build_dashboard(root: Path, output: Path, viewer_name: str, pattern: str | None):
    samples = discover_samples(root, viewer_name, pattern)
    if not samples:
        raise ValueError(
            f"No immediate subfolders under {root} contained {viewer_name!r}"
            + (f" matching pattern {pattern!r}" if pattern else "")
        )

    payload = []
    hash_groups: dict[str, list[str]] = {}

    for sample_name, sample_dir, viewer_path in samples:
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
                "b64": base64.b64encode(page.encode("utf-8")).decode("ascii"),
            }
        )

    options = "\n".join(
        f'<option value="{i}">{html.escape(item["name"])}</option>'
        for i, item in enumerate(payload)
    )

    # Base64 avoids all nested </script>, quoting, and multiline-HTML problems.
    import json

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
  @media(max-width:800px) {{
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
    print(f"Embedded {len(payload)} sample viewer(s):")
    for item in payload:
        print(f"  {item['name']}  sha256={item['sha']}  embedded={item['bytes']:,} bytes")

    duplicates = [names for names in hash_groups.values() if len(names) > 1]
    if duplicates:
        print("\nWARNING: some graph_index.html files are byte-for-byte identical:")
        for names in duplicates:
            print("  " + ", ".join(names))
        print(
            "If those samples are expected to differ, inspect the per-sample "
            "graph_index.html generation before blaming the dashboard."
        )


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
        build_dashboard(root, output, args.viewer, args.pattern)
    except (OSError, ValueError, UnicodeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
