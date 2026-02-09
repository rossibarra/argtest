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


def dump_ts(ts: tskit.TreeSequence, out_path: Path) -> None:
    if out_path.suffix == ".tsz":
        if tszip is None:
            raise RuntimeError("tszip is required to write .tsz files")
        tszip.compress(ts, out_path)
        return
    ts.dump(str(out_path))


def parse_regions_bed(paths):
    intervals = []
    for path in paths:
        p = Path(path)
        if not p.exists():
            raise FileNotFoundError(f"Remove BED not found: {p}")
        with open(p, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) < 3:
                    raise ValueError(f"Invalid BED line in {p}: {line}")
                start = float(parts[1])
                end = float(parts[2])
                if end <= start:
                    continue
                intervals.append((start, end))
    intervals.sort()
    return intervals


def merge_intervals(intervals):
    if not intervals:
        return []
    merged = [list(intervals[0])]
    for start, end in intervals[1:]:
        last = merged[-1]
        if start <= last[1]:
            last[1] = max(last[1], end)
        else:
            merged.append([start, end])
    return [(s, e) for s, e in merged]


def clamp_intervals(intervals, sequence_length):
    clamped = []
    for start, end in intervals:
        s = max(0.0, float(start))
        e = min(float(sequence_length), float(end))
        if e > s:
            clamped.append((s, e))
    return clamped


def complement_intervals(remove, sequence_length):
    keep = []
    cursor = 0.0
    for start, end in remove:
        if start > cursor:
            keep.append((cursor, start))
        cursor = max(cursor, end)
    if cursor < sequence_length:
        keep.append((cursor, sequence_length))
    return keep


def parse_args():
    p = argparse.ArgumentParser(
        description="Remove regions from a tree sequence based on BED intervals",
    )
    p.add_argument("ts", help="Tree sequence file (.ts, .trees, or .tsz)")
    p.add_argument(
        "--remove",
        action="append",
        required=True,
        help="BED file(s) of regions to remove (comma-separated or repeated)",
    )
    p.add_argument(
        "--simplify",
        action="store_true",
        help="Simplify the output tree sequence after trimming",
    )
    p.add_argument("--out", help="Output tree sequence path (.ts, .trees, or .tsz)")
    return p.parse_args()


def parse_remove_list(values):
    if not values:
        return []
    paths = []
    for value in values:
        parts = [v.strip() for v in value.split(",")]
        paths.extend([p for p in parts if p])
    return paths


def main():
    args = parse_args()
    ts_path = Path(args.ts)
    ts = load_ts(ts_path)

    remove_paths = parse_remove_list(args.remove)
    remove_intervals = merge_intervals(parse_regions_bed(remove_paths))
    remove_intervals = clamp_intervals(remove_intervals, ts.sequence_length)
    keep = complement_intervals(remove_intervals, ts.sequence_length)
    if not keep:
        tables = ts.dump_tables()
        tables.edges.clear()
        tables.sites.clear()
        tables.mutations.clear()
        tables.sequence_length = ts.sequence_length
        trimmed = tables.tree_sequence()
    else:
        trimmed = ts.keep_intervals(keep, simplify=False)

    if args.simplify:
        trimmed = trimmed.simplify(keep_unary=True)
        tables = trimmed.dump_tables()
        # Preserve original coordinate system after simplify.
        tables.sequence_length = ts.sequence_length
        trimmed = tables.tree_sequence()

    if args.out:
        out_path = Path(args.out)
    else:
        out_dir = Path(__file__).resolve().parent.parent / "results"
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{ts_path.stem}_regions_trimmed.tsz"
    dump_ts(trimmed, out_path)


if __name__ == "__main__":
    main()
