#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from argtest_common import (
    dump_ts,
    load_remove_intervals,
    load_ts,
    name_to_nodes_map,
    validate_trimmed_ts,
)


def parse_remove_list(values):
    # Accept comma-separated lists or repeated flags.
    if not values:
        return []
    paths = []
    for value in values:
        parts = [v.strip() for v in value.split(",")]
        paths.extend([p for p in parts if p])
    return paths


def parse_individuals(values):
    # Parse comma-separated IDs into a list.
    if not values:
        return []
    return [v.strip() for v in values.split(",") if v.strip()]


def merge_intervals(base, extra):
    # Merge per-individual intervals, keeping them sorted.
    merged = {k: {"starts": list(v["starts"]), "ends": list(v["ends"])} for k, v in base.items()}
    for name, spans in extra.items():
        entry = merged.setdefault(name, {"starts": [], "ends": []})
        entry["starts"].extend(spans["starts"])
        entry["ends"].extend(spans["ends"])
    for name, spans in merged.items():
        paired = sorted(zip(spans["starts"], spans["ends"]))
        spans["starts"] = [s for s, _ in paired]
        spans["ends"] = [e for _, e in paired]
    return merged


def parse_args():
    p = argparse.ArgumentParser(
        description="Remove individuals over BED intervals and write a trimmed tree sequence",
    )
    p.add_argument("ts", help="Tree sequence file (.ts, .trees, or .tsz)")
    p.add_argument(
        "--individuals",
        help="Comma-separated individual IDs to remove across the entire sequence",
    )
    p.add_argument(
        "--remove",
        action="append",
        help="BED file(s) of regions to remove per individual (comma-separated or repeated)",
    )
    p.add_argument("--out", help="Output tree sequence path (.ts, .trees, or .tsz)")
    p.add_argument("--suffix-to-strip", default="_anchorwave")
    return p.parse_args()


def remove_ancestry(ts, samples, left, right):
    """
    Removes the ancestry for `samples` over `[left, right)`, by:
    1. Split all edges intersecting "left" and "right"
    2. Remove singleton edges above nodes for which we don't want ancestry over [left, right]
    3. Throw into simplify, which will remove this ancestry
    4. Squash edges to "join" previously split edges
    """

    def split_edges_at(tables, position):
        for i, edge in enumerate(tables.edges):
            if edge.left < position < edge.right:
                tables.edges[i] = edge.replace(right=position)
                tables.edges.append(edge.replace(left=position))

    tables = ts.dump_tables()
    # Split edges so we can drop the target interval exactly.
    split_edges_at(tables, left)
    split_edges_at(tables, right)
    drop_edges = np.logical_and.reduce(
        [
            np.isin(tables.edges.child, samples),
            tables.edges.left >= left,
            tables.edges.right <= right,
        ]
    )
    tables.edges.keep_rows(~drop_edges)
    tables.sort()
    # Simplify drops the removed ancestry and may renumber nodes.
    tables.simplify()
    tables.edges.squash()
    return tables.tree_sequence()


def main():
    args = parse_args()
    ts_path = Path(args.ts)
    ts = load_ts(ts_path)

    remove_intervals = {}
    if args.remove:
        remove_paths = parse_remove_list(args.remove)
        remove_intervals = load_remove_intervals(remove_paths)

    individuals = parse_individuals(args.individuals)
    if individuals:
        # Expand full-length removals to [0, sequence_length).
        full = {
            name: {"starts": [0.0], "ends": [float(ts.sequence_length)]}
            for name in individuals
        }
        remove_intervals = merge_intervals(remove_intervals, full)

    if not remove_intervals:
        raise SystemExit("ERROR: provide --individuals and/or --remove")

    # Natefun-style removal can change sample identities/order via simplify.
    trimmed_ts = ts
    for name, intervals in remove_intervals.items():
        # Rebuild the name->nodes map after each simplify.
        name_to_nodes = name_to_nodes_map(trimmed_ts, suffix_to_strip=args.suffix_to_strip)
        samples = name_to_nodes.get(name, [])
        if not samples:
            continue
        for left, right in zip(intervals["starts"], intervals["ends"]):
            trimmed_ts = remove_ancestry(trimmed_ts, samples, left, right)
    validate_trimmed_ts(trimmed_ts)

    if args.out:
        out_path = Path(args.out)
    else:
        # Default output to results/ with a trimmed suffix.
        out_dir = Path(__file__).resolve().parent.parent / "results"
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{ts_path.stem}_trimmed.tsz"
    dump_ts(trimmed_ts, out_path)


if __name__ == "__main__":
    main()
