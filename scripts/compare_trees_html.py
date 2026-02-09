#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import tskit

try:
    import tszip
except Exception:  # pragma: no cover - optional dependency
    tszip = None


def load_ts(path: Path) -> tskit.TreeSequence:
    if path.suffix == ".tsz":
        if tszip is None:
            raise RuntimeError("tszip is required to load .tsz files")
        return tszip.load(str(path))
    return tskit.load(str(path))


def parse_args():
    p = argparse.ArgumentParser(
        description="Render two tree indices from two tree sequences into an HTML report",
    )
    p.add_argument("ts_a", help="First tree sequence file (.ts, .trees, or .tsz)")
    p.add_argument("index_a", type=int, help="Tree index for the first file (0-based)")
    p.add_argument("ts_b", help="Second tree sequence file (.ts, .trees, or .tsz)")
    p.add_argument("index_b", type=int, help="Tree index for the second file (0-based)")
    p.add_argument("--out", default="tree_compare.html", help="Output HTML file")
    return p.parse_args()


def main():
    args = parse_args()
    path_a = Path(args.ts_a)
    path_b = Path(args.ts_b)

    ts_a = load_ts(path_a)
    ts_b = load_ts(path_b)

    tree_a = ts_a.at_index(args.index_a)
    tree_b = ts_b.at_index(args.index_b)

    ascii_a = tree_a.draw_text()
    ascii_b = tree_b.draw_text()

    html = f"""<!doctype html>
<html lang="en">
<meta charset="utf-8">
<title>Tree comparison</title>
<style>
body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; margin: 24px; }}
pre {{ background: #f5f5f5; padding: 12px; overflow-x: auto; }}
.panel {{ margin-bottom: 24px; }}
.meta {{ color: #444; font-size: 13px; margin-bottom: 8px; }}
</style>
<h1>Tree comparison</h1>
<div class="panel">
  <div class="meta">File A: {path_a.name} | Tree index: {args.index_a}</div>
  <pre>{ascii_a}</pre>
</div>
<div class="panel">
  <div class="meta">File B: {path_b.name} | Tree index: {args.index_b}</div>
  <pre>{ascii_b}</pre>
</div>
</html>
"""
    Path(args.out).write_text(html)


if __name__ == "__main__":
    main()
