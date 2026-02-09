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
        description="Render all trees from two tree sequences into an HTML gallery",
    )
    p.add_argument("ts_top", help="Top tree sequence file (.ts, .trees, or .tsz)")
    p.add_argument("ts_bottom", help="Bottom tree sequence file (.ts, .trees, or .tsz)")
    p.add_argument("--out", default="trees_gallery.html", help="Output HTML file")
    return p.parse_args()


def tree_card(tree, index: int, label: str) -> str:
    interval = tree.interval
    ascii_tree = tree.draw_text()
    return f"""<div class="tree-card">
  <div class="tree-meta">{label} | Tree {index} | [{interval.left:.3f}, {interval.right:.3f})</div>
  <pre>{ascii_tree}</pre>
</div>"""


def main():
    args = parse_args()
    path_top = Path(args.ts_top)
    path_bottom = Path(args.ts_bottom)

    ts_top = load_ts(path_top)
    ts_bottom = load_ts(path_bottom)

    n = min(ts_top.num_trees, ts_bottom.num_trees)
    top_cards = []
    bottom_cards = []
    for i in range(n):
        top_cards.append(tree_card(ts_top.at_index(i), i, f"Top: {path_top.name}"))
        bottom_cards.append(tree_card(ts_bottom.at_index(i), i, f"Bottom: {path_bottom.name}"))

    html = f"""<!doctype html>
<html lang="en">
<meta charset="utf-8">
<title>Tree gallery</title>
<style>
body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; margin: 24px; }}
h1 {{ font-size: 20px; }}
h2 {{ font-size: 16px; margin-top: 28px; }}
.tree-row {{ display: flex; flex-direction: row; gap: 8px; overflow-x: auto; padding-bottom: 8px; }}
.tree-card {{ border: 1px solid #ddd; border-radius: 6px; padding: 8px; background: #fafafa; min-width: 320px; }}
.tree-meta {{ color: #444; font-size: 12px; margin-bottom: 6px; }}
pre {{ margin: 0; font-size: 10px; line-height: 1.2; overflow-x: auto; }}
</style>
<h1>Tree gallery</h1>
<div class="meta">Top: {path_top.name} | Trees: {ts_top.num_trees}</div>
<div class="meta">Bottom: {path_bottom.name} | Trees: {ts_bottom.num_trees}</div>
<h2>Top: {path_top.name}</h2>
<div class="tree-row">
{''.join(top_cards)}
</div>
<h2>Bottom: {path_bottom.name}</h2>
<div class="tree-row">
{''.join(bottom_cards)}
</div>
</html>
"""

    Path(args.out).write_text(html)


if __name__ == "__main__":
    main()
