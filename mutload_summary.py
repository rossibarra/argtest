#!/usr/bin/env python
import argparse
import base64
from io import BytesIO
from pathlib import Path
import os
import sys

import numpy as np

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


def get_individual_name(ind, suffix_to_strip="_anchorwave") -> str:
    nm = None
    try:
        if isinstance(ind.metadata, dict):
            nm = ind.metadata.get("id")
    except Exception:
        nm = None
    if nm is None:
        nm = f"ind{ind.id}"
    if isinstance(nm, bytes):
        nm = nm.decode()
    return nm.replace(suffix_to_strip, "")


def drop_individuals_by_name(ts: tskit.TreeSequence, targets, suffix_to_strip="_anchorwave"):
    if not targets:
        return ts

    targets = set(targets)
    drop_nodes = set()
    found = set()
    for ind_id, ind in enumerate(ts.individuals()):
        nm = get_individual_name(ind, suffix_to_strip=suffix_to_strip)
        if nm in targets:
            found.add(nm)
            drop_nodes.update(ts.individual(ind_id).nodes)

    missing = sorted(targets - found)
    if missing:
        raise ValueError(f"Individuals not found: {', '.join(missing)}")

    keep = [u for u in ts.samples() if u not in drop_nodes]
    return ts.simplify(samples=keep, keep_unary=True)


def mutational_load(ts: tskit.TreeSequence, windows: np.ndarray = None) -> np.ndarray:
    genome_windows = np.array([0, ts.sequence_length]) if windows is None else windows
    assert genome_windows[0] == 0 and genome_windows[-1] == ts.sequence_length
    mutations_window = np.digitize(ts.sites_position[ts.mutations_site], genome_windows) - 1
    assert mutations_window.min() >= 0 and mutations_window.max() < genome_windows.size - 1
    load = np.zeros((genome_windows.size - 1, ts.num_samples))
    tree = ts.first(sample_lists=True)
    for s in ts.sites():
        tree.seek(s.position)
        for m in s.mutations:
            if m.edge != tskit.NULL:
                window = mutations_window[m.id]
                samples = list(tree.samples(m.node))
                load[window, samples] += 1.0
    return load.squeeze(0) if windows is None else load


def sample_names(ts: tskit.TreeSequence, suffix_to_strip="_anchorwave"):
    names = []
    for u in ts.samples():
        ind_id = ts.node(u).individual
        if ind_id != tskit.NULL:
            nm = get_individual_name(ts.individual(ind_id), suffix_to_strip=suffix_to_strip)
        else:
            nm = f"node{u}"
        names.append(nm)
    return names


def aggregate_by_individual(load, names):
    unique = []
    idx_map = {}
    for i, nm in enumerate(names):
        if nm not in idx_map:
            idx_map[nm] = len(unique)
            unique.append(nm)

    if load.ndim == 1:
        agg = np.zeros(len(unique), dtype=float)
        for i, nm in enumerate(names):
            agg[idx_map[nm]] += load[i]
        return agg, unique

    agg = np.zeros((load.shape[0], len(unique)), dtype=float)
    for i, nm in enumerate(names):
        agg[:, idx_map[nm]] += load[:, i]
    return agg, unique


def fig_to_data_url(fig) -> str:
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=160, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode("ascii")
    return f"data:image/png;base64,{b64}"


def plot_single(load, names, title):
    fig, ax = plt.subplots(figsize=(max(6, 0.35 * len(names)), 4))
    ax.bar(names, load, color="#777777")
    ax.set_ylabel("Derived mutational load")
    ax.set_title(title)
    ax.tick_params(axis="x", labelrotation=90, labelsize=8)
    fig.tight_layout()
    return fig


def plot_windows(load, names, windows):
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
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(counts, bins=20, color="#777777", edgecolor="#222222")
    ax.set_xlabel("Outlier windows per individual")
    ax.set_ylabel("Count of individuals")
    ax.set_title("Outlier window counts")
    fig.tight_layout()
    return fig


def parse_remove_list(value):
    if value is None:
        return []
    parts = [v.strip() for v in value.split(",")]
    return [p for p in parts if p]


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
    repo_root = Path(__file__).resolve().parent
    results_dir = repo_root / "results"
    logs_dir = repo_root / "logs"
    beds_dir = results_dir / "beds"
    results_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)
    beds_dir.mkdir(parents=True, exist_ok=True)

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
                outlier_counts = []
                for i in range(load.shape[1]):
                    count = 0
                    for w in range(load.shape[0]):
                        if window_means[w] <= 0:
                            continue
                    if load[w, i] > (1 + args.cutoff) * window_means[w] or load[w, i] < (1 - args.cutoff) * window_means[w]:
                        count += 1
                    outlier_counts.append(count)

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
<html lang=\"en\">
<meta charset=\"utf-8\">
<title>Mutational load summary</title>
<style>
body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; margin: 24px; }}
h1 {{ font-size: 20px; }}
.meta {{ color: #444; font-size: 13px; margin-bottom: 16px; }}
img {{ max-width: 100%; height: auto; }}
</style>
<h1>Mutational load summary</h1>
{hist_html}
<div class=\"meta\">Input: {ts_path.name} | Samples: {ts.num_samples} | Individuals: {len(unique_names)}</div>
<img src=\"{img_url}\" alt=\"Mutational load plot\">
"""

            if windows is not None:
                html += f"<div class=\"meta\">Window size: {int(args.window_size)} bp</div>"
                html += f"<div class=\"meta\">Outlier cutoff: {args.cutoff:.3f} of window mean</div>"

            html += "</html>\n"
            out.write_text(html)

            if windows is not None:
                # Write per-individual BED files for windows outside cutoff range of window mean
                for i, name in enumerate(unique_names):
                    safe = "".join(c if c.isalnum() or c in ("-", "_", ".") else "_" for c in name)
                    bed_path = beds_dir / f"{safe}.bed"
                    lines = []
                    for w in range(load.shape[0]):
                        if window_means[w] <= 0:
                            continue
                        if load[w, i] > (1 + args.cutoff) * window_means[w] or load[w, i] < (1 - args.cutoff) * window_means[w]:
                            start = int(windows[w])
                            end = int(windows[w + 1])
                            lines.append(
                                f"{ts_path.stem}\t{start}\t{end}\t{name}\t{window_means[w]:.3f}\t{load[w, i]:.3f}"
                            )
                    bed_path.write_text("\n".join(lines) + ("\n" if lines else ""))
        finally:
            sys.stdout = old_stdout
            sys.stderr = old_stderr


if __name__ == "__main__":
    main()
