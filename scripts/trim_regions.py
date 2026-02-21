#!/usr/bin/env python
from __future__ import annotations

import argparse
import pickle
import re
from pathlib import Path

import numpy as np

from argtest_common import build_shared_mask, collapse_masked_intervals, dump_ts, load_ts


def infer_mu_base(ts_stem: str) -> list[str]:
    bases = [ts_stem]
    m = re.match(r"^(.+)\.(\d+)$", ts_stem)
    if m:
        bases.append(m.group(1))
    m = re.match(r"^(.+)[_-](\d+)$", ts_stem)
    if m:
        bases.append(m.group(1))
    dedup = []
    seen = set()
    for b in bases:
        if b not in seen:
            dedup.append(b)
            seen.add(b)
    return dedup


def infer_mu_path(ts_path: Path) -> Path:
    bases = infer_mu_base(ts_path.stem)
    search_dirs = [ts_path.parent, ts_path.parent.parent]
    for d in search_dirs:
        for b in bases:
            p = d / f"{b}.mut_rate.p"
            if p.exists():
                return p
    candidates = []
    for d in search_dirs:
        if not d.exists():
            continue
        for p in d.glob("*.mut_rate.p"):
            nm = p.name
            if any(nm.startswith(f"{b}.") or nm == f"{b}.mut_rate.p" for b in bases):
                candidates.append(p)
    candidates = sorted(set(candidates))
    if len(candidates) == 1:
        return candidates[0]
    if len(candidates) > 1:
        raise RuntimeError(
            f"Ambiguous mutation maps for {ts_path.name}: "
            + ", ".join(str(x) for x in candidates)
        )
    raise FileNotFoundError(
        f"Could not infer mutation map for {ts_path}. Tried bases={bases} in {search_dirs}"
    )


def format_num(x: float) -> str:
    if float(x).is_integer():
        return str(int(x))
    return str(x).replace(".", "p")


def output_name(ts_path: Path, window_size: float, cutoff_bp: float) -> str:
    ext = ts_path.suffix
    return (
        f"{ts_path.stem}.collapsed.ws{format_num(window_size)}."
        f"accbp{format_num(cutoff_bp)}{ext}"
    )


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Trim regions by collapsing masked and low-accessibility windows in a directory of tree sequences."
        )
    )
    p.add_argument(
        "--ts-dir",
        required=True,
        type=Path,
        help="Directory containing tree sequence files (.tsz, .ts, .trees).",
    )
    p.add_argument(
        "--window-size",
        required=True,
        type=float,
        help="Window size in bp for accessibility checks.",
    )
    p.add_argument(
        "--cutoff-bp",
        required=True,
        type=float,
        help="Minimum accessible bp in a window to keep it.",
    )
    p.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory for collapsed tree files (default: <ts-dir>/collapsed).",
    )
    p.add_argument(
        "--pattern",
        default="*",
        help="Optional glob pattern to filter tree sequence filenames (default: '*').",
    )
    p.add_argument(
        "--log",
        type=Path,
        default=None,
        help="Optional log file path (default: <out-dir>/collapse_log.txt).",
    )
    return p.parse_args()


def find_tree_files(ts_dir: Path, pattern: str) -> list[Path]:
    if not ts_dir.exists():
        raise FileNotFoundError(f"Tree directory does not exist: {ts_dir}")
    files = sorted(
        [
            p
            for p in ts_dir.glob(pattern)
            if p.is_file() and p.suffix in {".tsz", ".ts", ".trees"}
        ]
    )
    if not files:
        raise RuntimeError(f"No tree files found in {ts_dir} matching pattern '{pattern}'.")
    return files


def main():
    args = parse_args()
    ts_files = find_tree_files(args.ts_dir, args.pattern)
    out_dir = args.out_dir or (args.ts_dir / "collapsed")
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = args.log or (out_dir / "collapse_log.txt")

    first_ts = load_ts(ts_files[0])
    mu_path = infer_mu_path(ts_files[0])
    with open(mu_path, "rb") as f:
        mu = pickle.load(f)
    mask, good_intervals, windows, acc_bp = build_shared_mask(
        ts=first_ts,
        mu=mu,
        window_size=args.window_size,
        cutoff_bp=args.cutoff_bp,
    )
    kept_bp_windows = float(sum(r - l for l, r in good_intervals))
    kept_bp_mask = float(mask.get_cumulative_mass(first_ts.sequence_length))

    with open(log_path, "w") as log:
        log.write("# trim_regions\n")
        log.write(f"# ts_dir={args.ts_dir}\n")
        log.write(f"# out_dir={out_dir}\n")
        log.write(f"# window_size={args.window_size}\n")
        log.write(f"# cutoff_bp={args.cutoff_bp}\n\n")
        log.write("# shared_window_filter\n")
        log.write(f"reference_ts={ts_files[0]}\n")
        log.write(f"mu_file={mu_path}\n")
        log.write(f"sequence_length={first_ts.sequence_length}\n")
        log.write(f"n_windows={len(windows) - 1}\n")
        log.write(f"n_good_windows={int(np.sum(acc_bp >= args.cutoff_bp))}\n")
        log.write(f"kept_bp_from_windows={kept_bp_windows}\n")
        log.write(f"kept_bp_from_mask={kept_bp_mask}\n\n")

        for ts_file in ts_files:
            ts = load_ts(ts_file)
            if float(ts.sequence_length) != float(first_ts.sequence_length):
                raise RuntimeError(
                    f"Sequence length mismatch for {ts_file}: "
                    f"{ts.sequence_length} != {first_ts.sequence_length}"
                )
            inferred_mu = infer_mu_path(ts_file)
            if inferred_mu != mu_path:
                raise RuntimeError(
                    f"Inferred mutation map differs for {ts_file}: {inferred_mu} != {mu_path}. "
                    "This script currently uses one shared window filter for all files."
                )
            ts2 = collapse_masked_intervals(ts, mask)

            out_path = out_dir / output_name(ts_file, args.window_size, args.cutoff_bp)
            dump_ts(ts2, out_path)

            dropped_bp = float(ts.sequence_length - ts2.sequence_length)

            log.write("=" * 80 + "\n")
            log.write(f"ts_file={ts_file}\n")
            log.write(f"out_file={out_path}\n")
            log.write(f"old_L={ts.sequence_length}\n")
            log.write(f"new_L={ts2.sequence_length}\n")
            log.write(f"dropped_bp={dropped_bp}\n")
            log.write("\n")

            print(f"Wrote: {out_path.name}")

    print(f"Done. Log: {log_path}")


if __name__ == "__main__":
    main()
