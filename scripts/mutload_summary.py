#!/usr/bin/env python
from __future__ import annotations

import argparse
import base64
from io import BytesIO
import os
import sys
from pathlib import Path

import numpy as np

from argtest_common import (
    aggregate_by_individual,
    load_ts,
    mutational_load,
    sample_names,
)


def fig_to_data_url(fig) -> str:
    # Encode a matplotlib figure as an inline PNG data URL.
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=160, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode("ascii")
    return f"data:image/png;base64,{b64}"


def plot_single(load, names, title):
    # Simple per-individual bar chart.
    fig, ax = plt.subplots(figsize=(max(6, 0.35 * len(names)), 4))
    ax.bar(names, load, color="#777777")
    ax.set_ylabel("Derived mutational load")
    ax.set_title(title)
    ax.tick_params(axis="x", labelrotation=90, labelsize=8)
    fig.tight_layout()
    return fig


def plot_windows(load, names, windows):
    # Grid of per-window bar charts with shared y-scale.
    nwin = load.shape[0]
    ncols = 4
    nrows = int(np.ceil(nwin / ncols))
    # Use a shared y-axis scale across all window plots
    ymax = float(load.max()) if load.size else 0.0
    ytop = ymax * 1.05 if ymax > 0 else 1.0
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4, nrows * 3), squeeze=False)
    for i in range(nwin):
        r, c = divmod(i, ncols)
        ax = axes[r][c]
        ax.bar(names, load[i], color="#777777")
        ax.set_ylim(0, ytop)
        left = int(windows[i])
        right = int(windows[i + 1])
        ax.set_title(f"{left:,}-{right:,} bp")
        ax.tick_params(axis="x", labelrotation=90, labelsize=6)
    for j in range(nwin, nrows * ncols):
        r, c = divmod(j, ncols)
        axes[r][c].axis("off")
    fig.tight_layout()
    return fig


def plot_outlier_hist(counts):
    # Histogram of how many outlier windows each individual has.
    counts = [int(c) for c in counts]
    max_count = max(counts) if counts else 0
    bins = list(range(0, max_count + 2))
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(counts, bins=bins, color="#777777", edgecolor="#222222")
    ax.set_xlabel("Outlier windows per individual")
    ax.set_ylabel("Count of individuals")
    ax.set_title("Outlier window counts")
    ax.set_xticks(bins)
    fig.tight_layout()
    return fig


def parse_args():
    p = argparse.ArgumentParser(
        description="Summarize mutational load per individual from a tree sequence",
    )
    p.add_argument("ts", help="Tree sequence file (.ts, .trees, or .tsz)")
    p.add_argument("--window-size", type=float, help="Window size in bp")
    p.add_argument(
        "--cutoff",
        type=float,
        default=0.25,
        help="Outlier cutoff as a fraction of the window mean (default: 0.25)",
    )
    p.add_argument("--out", default="mutational_load_summary.html")
    p.add_argument("--suffix-to-strip", default="_anchorwave")
    return p.parse_args()


class Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for s in self.streams:
            s.write(data)
            s.flush()

    def flush(self):
        for s in self.streams:
            s.flush()


def main():
    args = parse_args()
    ts_path = Path(args.ts)
    repo_root = Path(__file__).resolve().parent.parent
    results_dir = repo_root / "results"
    logs_dir = repo_root / "logs"
    results_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    out = results_dir / Path(args.out).name

    log_path = logs_dir / f"{out.stem}.log"
    old_stdout = sys.stdout
    old_stderr = sys.stderr
    with open(log_path, "w") as log_fh:
        try:
            sys.stdout = Tee(old_stdout, log_fh)
            sys.stderr = Tee(old_stderr, log_fh)
            # Ensure matplotlib cache is writable and avoid stderr warnings
            os.environ.setdefault("MPLCONFIGDIR", str(logs_dir / ".matplotlib"))
            global plt  # pylint: disable=global-statement
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt  # noqa: E402

            ts = load_ts(ts_path)

            windows = None
            if args.window_size:
                # Build explicit window edges covering the full sequence.
                L = float(ts.sequence_length)
                windows = np.arange(0, L + args.window_size, args.window_size, dtype=float)
                if windows[-1] > L:
                    windows[-1] = L

            load = mutational_load(ts, windows=windows)
            names = sample_names(ts, suffix_to_strip=args.suffix_to_strip)
            load, unique_names = aggregate_by_individual(load, names)

            outlier_counts = None
            if windows is None:
                fig = plot_single(load, unique_names, "Mutational load")
            else:
                fig = plot_windows(load, unique_names, windows)
                window_means = load.mean(axis=1)
                valid = window_means > 0
                high = (1 + args.cutoff) * window_means
                low = (1 - args.cutoff) * window_means
                # Identify per-individual outliers in each window.
                mask = (load > high[:, None]) | (load < low[:, None])
                mask &= valid[:, None]
                outlier_counts = mask.sum(axis=0).astype(int).tolist()

            img_url = fig_to_data_url(fig)
            hist_html = ""
            if outlier_counts is not None:
                hist_fig = plot_outlier_hist(outlier_counts)
                hist_url = fig_to_data_url(hist_fig)
                hist_html = (
                    "<h2>Outlier window counts</h2>\n"
                    f"<img src=\"{hist_url}\" alt=\"Outlier window counts histogram\">\n"
                )
            html = f"""<!doctype html>
<html lang="en">
<meta charset="utf-8">
<title>Mutational load summary</title>
<style>
body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; margin: 24px; }}
h1 {{ font-size: 20px; }}
.meta {{ color: #444; font-size: 13px; margin-bottom: 16px; }}
img {{ max-width: 100%; height: auto; }}
</style>
<h1>Mutational load summary</h1>
{hist_html}
<div class="meta">Input: {ts_path.name} | Samples: {ts.num_samples} | Individuals: {len(unique_names)}</div>
<img src="{img_url}" alt="Mutational load plot">
"""

            if windows is not None:
                html += f"<div class=\"meta\">Window size: {int(args.window_size)} bp</div>"
                html += f"<div class=\"meta\">Outlier cutoff: {args.cutoff:.3f} of window mean</div>"

            html += "</html>\n"
            out.write_text(html)

            if windows is not None:
                # Write one BED listing outliers per window.
                out_path = results_dir / f"{ts_path.stem}_outliers.bed"
                lines = []
                for w in range(load.shape[0]):
                    if not valid[w]:
                        continue
                    row_mask = mask[w]
                    if not row_mask.any():
                        continue
                    outlier_names = [unique_names[i] for i in np.where(row_mask)[0]]
                    outlier_vals = [f"{load[w, i]:.3f}" for i in np.where(row_mask)[0]]
                    start = int(windows[w])
                    end = int(windows[w + 1])
                    lines.append(
                        f"{ts_path.stem}\t{start}\t{end}\t{','.join(outlier_names)}\t{','.join(outlier_vals)}\t{window_means[w]:.3f}"
                    )
                out_path.write_text("\n".join(lines) + ("\n" if lines else ""))
        finally:
            sys.stdout = old_stdout
            sys.stderr = old_stderr


if __name__ == "__main__":
    main()
