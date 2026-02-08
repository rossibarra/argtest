#!/usr/bin/env python
import argparse
import base64
from io import BytesIO
from pathlib import Path
import os
import sys
from bisect import bisect_right

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


def mutational_load(
    ts: tskit.TreeSequence,
    windows: np.ndarray = None,
    remove_intervals=None,
    name_to_nodes=None,
) -> np.ndarray:
    genome_windows = np.array([0, ts.sequence_length]) if windows is None else windows
    assert genome_windows[0] == 0 and genome_windows[-1] == ts.sequence_length
    load = np.zeros((genome_windows.size - 1, ts.num_samples))

    segments = [(0.0, ts.sequence_length, frozenset())]
    if remove_intervals and name_to_nodes:
        segments = build_removal_segments(remove_intervals, ts.sequence_length)

    for left, right, drop_names in segments:
        if right <= left:
            continue
        drop_nodes = set()
        for nm in drop_names:
            drop_nodes.update(name_to_nodes.get(nm, []))

        sub = ts.keep_intervals([(left, right)], simplify=False)
        if drop_nodes:
            keep = [u for u in ts.samples() if u not in drop_nodes]
            sub, node_map = sub.simplify(samples=keep, keep_unary=True, map_nodes=True)
            rev_map = np.full(sub.num_nodes, tskit.NULL, dtype=int)
            for u in keep:
                new_id = node_map[u]
                if new_id != tskit.NULL:
                    rev_map[new_id] = u
        else:
            rev_map = None

        site_windows = None
        for tree in sub.trees(sample_lists=True):
            for s in tree.sites():
                if site_windows is None:
                    site_windows = np.digitize(sub.sites_position, genome_windows) - 1
                window = int(site_windows[s.id])
                if window < 0 or window >= genome_windows.size - 1:
                    continue
                for m in s.mutations:
                    if m.edge != tskit.NULL:
                        samples = list(tree.samples(m.node))
                        if rev_map is not None:
                            samples = [rev_map[u] for u in samples if rev_map[u] != tskit.NULL]
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


def name_to_nodes_map(ts: tskit.TreeSequence, suffix_to_strip="_anchorwave"):
    mapping = {}
    for ind_id, ind in enumerate(ts.individuals()):
        nm = get_individual_name(ind, suffix_to_strip=suffix_to_strip)
        nodes = list(ts.individual(ind_id).nodes)
        mapping[nm] = nodes
    return mapping


def _pos_in_intervals(pos, intervals):
    starts = intervals["starts"]
    ends = intervals["ends"]
    idx = bisect_right(starts, pos) - 1
    return idx >= 0 and pos < ends[idx]


def _interval_overlaps(left, right, intervals):
    starts = intervals["starts"]
    ends = intervals["ends"]
    idx = bisect_right(starts, right) - 1
    if idx < 0:
        return False
    return ends[idx] > left


def build_removal_segments(remove_intervals, sequence_length):
    events = {}
    for name, intervals in remove_intervals.items():
        for start, end in zip(intervals["starts"], intervals["ends"]):
            events.setdefault(start, []).append((name, 1))
            events.setdefault(end, []).append((name, -1))

    positions = sorted(set([0.0, float(sequence_length)] + list(events.keys())))
    active = {}
    segments = []
    for i, pos in enumerate(positions[:-1]):
        for name, delta in events.get(pos, []):
            active[name] = active.get(name, 0) + delta
            if active[name] <= 0:
                active.pop(name, None)
        next_pos = positions[i + 1]
        if next_pos > pos:
            segments.append((pos, next_pos, frozenset(active.keys())))
    return segments


def build_segments_with_drop_nodes(remove_intervals, name_to_nodes, sequence_length):
    segments = build_removal_segments(remove_intervals, sequence_length)
    out = []
    for left, right, drop_names in segments:
        drop_nodes = set()
        for nm in drop_names:
            drop_nodes.update(name_to_nodes.get(nm, []))
        out.append((left, right, drop_nodes))
    return out


def trim_ts_by_intervals(ts: tskit.TreeSequence, remove_intervals, name_to_nodes):
    segments = build_segments_with_drop_nodes(remove_intervals, name_to_nodes, ts.sequence_length)
    seg_lefts = [s[0] for s in segments]
    seg_rights = [s[1] for s in segments]

    tables = ts.dump_tables()
    tables.edges.clear()
    tables.sites.clear()
    tables.mutations.clear()

    # Edges: split by segments, remove if parent/child dropped in that segment
    for edge in ts.edges():
        i = bisect_right(seg_rights, edge.left)
        while i < len(segments) and seg_lefts[i] < edge.right:
            left, right, drop_nodes = segments[i]
            seg_l = max(edge.left, left)
            seg_r = min(edge.right, right)
            if seg_r > seg_l:
                if edge.parent not in drop_nodes and edge.child not in drop_nodes:
                    tables.edges.add_row(
                        left=seg_l,
                        right=seg_r,
                        parent=edge.parent,
                        child=edge.child,
                        metadata=edge.metadata,
                    )
            i += 1

    # Sites and mutations: drop mutations whose node is removed in the segment
    for site in ts.sites():
        i = bisect_right(seg_rights, site.position)
        if i >= len(segments) or seg_lefts[i] > site.position:
            i = bisect_right(seg_lefts, site.position) - 1
        if i < 0 or i >= len(segments):
            continue
        _, _, drop_nodes = segments[i]
        kept = []
        for mut in site.mutations:
            if mut.node in drop_nodes:
                continue
            kept.append(mut)
        if not kept:
            continue
        new_site_id = tables.sites.add_row(
            position=site.position,
            ancestral_state=site.ancestral_state,
            metadata=site.metadata,
        )
        for mut in kept:
            tables.mutations.add_row(
                site=new_site_id,
                node=mut.node,
                derived_state=mut.derived_state,
                parent=tskit.NULL,
                time=mut.time,
                metadata=mut.metadata,
            )

    tables.sort()
    tables.build_index()
    tables.compute_mutation_parents()
    return tables.tree_sequence()


def load_remove_intervals(paths):
    remove = {}
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
                if len(parts) >= 4:
                    raw = parts[3]
                    names = [n.strip() for n in raw.split(",") if n.strip()]
                else:
                    names = [p.stem]
                for name in names:
                    remove.setdefault(name, []).append((start, end))

    intervals = {}
    for name, spans in remove.items():
        spans.sort()
        starts = [s for s, _ in spans]
        ends = [e for _, e in spans]
        intervals[name] = {"starts": starts, "ends": ends}
    return intervals


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


def parse_remove_list(values):
    if not values:
        return []
    paths = []
    for value in values:
        parts = [v.strip() for v in value.split(",")]
        paths.extend([p for p in parts if p])
    return paths


def parse_args():
    p = argparse.ArgumentParser(
        description="Summarize mutational load per individual from a tree sequence",
    )
    p.add_argument("ts", help="Tree sequence file (.ts, .trees, or .tsz)")
    p.add_argument("--window-size", type=float, help="Window size in bp")
    p.add_argument(
        "--remove",
        action="append",
        help="BED file(s) of regions to remove per individual (comma-separated or repeated)",
    )
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
                L = float(ts.sequence_length)
                windows = np.arange(0, L + args.window_size, args.window_size, dtype=float)
                if windows[-1] > L:
                    windows[-1] = L

            remove_paths = parse_remove_list(args.remove)
            remove_intervals = load_remove_intervals(remove_paths) if remove_paths else None
            name_to_nodes = name_to_nodes_map(ts, suffix_to_strip=args.suffix_to_strip)
            if remove_intervals:
                ts = trim_ts_by_intervals(ts, remove_intervals, name_to_nodes)
                if tszip is None:
                    raise RuntimeError("tszip is required to write trimmed .tsz files")
                trimmed_path = results_dir / f"{ts_path.stem}_trimmed.tsz"
                tszip.compress(ts, trimmed_path)
                remove_intervals = None
            load = mutational_load(
                ts,
                windows=windows,
                remove_intervals=remove_intervals,
                name_to_nodes=name_to_nodes,
            )
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
                # Write one BED listing outliers per window
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
