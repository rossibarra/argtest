#!/usr/bin/env python
from __future__ import annotations

import argparse
import os
from pathlib import Path
import warnings

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from argtest_common import load_ts, mutational_load

matplotlib.rcParams["figure.dpi"] = 300


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Generate SINGER-style validation/diagnostic plots from a set of tree sequences."
        )
    )
    p.add_argument(
        "--ts-dir",
        required=True,
        type=Path,
        help="Directory containing tree sequence files (.tsz, .ts, .trees).",
    )
    p.add_argument(
        "--pattern",
        default="*.tsz",
        help="Glob pattern for input trees (default: *.tsz).",
    )
    p.add_argument(
        "--out-dir",
        type=Path,
        default=Path("results/validation_plots"),
        help="Output directory for plots.",
    )
    p.add_argument(
        "--time-bins-file",
        required=True,
        type=Path,
        help=(
            "Text file with explicit time-bin edges (whitespace/comma separated). "
            "If 0 or inf are missing they are added automatically."
        ),
    )
    p.add_argument(
        "--burnin-frac",
        type=float,
        default=0.5,
        help="Fraction of initial trees to ignore when computing posterior means.",
    )
    p.add_argument(
        "--tail-cutoff",
        type=float,
        default=1e-12,
        help="Set coalescence rates to NaN when survival drops below this threshold.",
    )
    p.add_argument(
        "--time-adjust",
        type=float,
        default=1.0,
        help="Divide time-bin x-values by this factor for plotting (default: 1.0).",
    )
    p.add_argument(
        "--log-rates",
        action="store_true",
        help="Use log scale on y-axis for pair coalescence rates.",
    )
    p.add_argument(
        "--year",
        type=float,
        default=None,
        help="Optional x-position for a red dashed vertical marker on the Ne plot.",
    )
    p.add_argument(
        "--prefix",
        default="",
        help="Optional prefix for output plot filenames.",
    )
    p.add_argument(
        "--folded",
        action="store_true",
        help="Plot folded SFS (minor-allele frequency) instead of polarised derived-frequency SFS.",
    )
    p.add_argument(
        "--window-size",
        type=float,
        default=5.0e4,
        help="Window size (bp) for diversity/Tajima's D validation plots (default: 50000).",
    )
    p.add_argument(
        "--mutation-rate",
        type=float,
        required=True,
        help=(
            "User-defined mutation rate used to scale branch diversity for comparison to site diversity; "
            "Ne plot is computed as 1/(2*coal_rate)."
        ),
    )
    return p.parse_args()


def find_tree_files(ts_dir: Path, pattern: str) -> list[Path]:
    files = sorted(
        [p for p in ts_dir.glob(pattern) if p.is_file() and p.suffix in {".tsz", ".ts", ".trees"}]
    )
    if not files:
        raise RuntimeError(f"No tree files found in {ts_dir} matching pattern '{pattern}'.")
    return files


def load_time_windows(path: Path) -> np.ndarray:
    if not path.exists():
        raise FileNotFoundError(f"Time bins file not found: {path}")
    text = path.read_text().replace(",", " ")
    vals = [float(x) for x in text.split() if x.strip()]
    if len(vals) < 2:
        raise ValueError("Time bins file must contain at least two bin edges.")
    windows = np.array(vals, dtype=float)
    if not np.all(np.diff(windows) > 0):
        raise ValueError("Time bins must be strictly increasing.")
    if windows[0] > 0:
        windows = np.append([0.0], windows)
    if not np.isinf(windows[-1]):
        windows = np.append(windows, np.inf)
    return windows


def plottable_interval_mask(time_windows: np.ndarray) -> np.ndarray:
    left = time_windows[:-1]
    right = time_windows[1:]
    return np.logical_and(np.isfinite(right), left > 0)


def genome_windows(sequence_length: float, window_size: float) -> np.ndarray:
    windows = np.arange(0, float(sequence_length) + float(window_size), float(window_size), dtype=float)
    if windows[-1] > float(sequence_length):
        windows[-1] = float(sequence_length)
    return windows


def compute_pair_coal(ts, time_windows: np.ndarray, tail_cutoff: float, keep_mask: np.ndarray):
    pdf = ts.pair_coalescence_counts(time_windows=time_windows, pair_normalise=True)
    rates = ts.pair_coalescence_rates(time_windows=time_windows)
    survival = np.append(1.0, 1.0 - np.cumsum(pdf))
    rates[survival[:-1] <= tail_cutoff] = np.nan
    return pdf[keep_mask], rates[keep_mask]


def safe_nanmean(a: np.ndarray, axis=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(a, axis=axis)


def safe_nanquantile(a: np.ndarray, q, axis=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanquantile(a, q, axis=axis)


def main():
    args = parse_args()
    ts_files = find_tree_files(args.ts_dir, args.pattern)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    if args.time_adjust <= 0:
        raise ValueError("--time-adjust must be > 0")

    time_windows = load_time_windows(args.time_bins_file)
    keep_mask = plottable_interval_mask(time_windows)
    breaks = time_windows[:-1][keep_mask]
    plot_breaks = breaks / args.time_adjust
    burnin = int(np.floor(len(ts_files) * args.burnin_frac))
    keep_idx = np.arange(len(ts_files))
    keep_post = keep_idx[burnin:] if burnin < len(ts_files) else keep_idx[-1:]

    pdf_vals = []
    rate_vals = []
    load_vals = []
    seq_lengths = []
    site_div_vals = []
    branch_div_vals = []
    site_td_vals = []
    branch_td_vals = []
    site_afs_vals = []
    branch_afs_vals = []

    n_samples = None
    stat_windows = None
    for ts_path in ts_files:
        ts = load_ts(ts_path)
        if n_samples is None:
            n_samples = ts.num_samples
            stat_windows = genome_windows(ts.sequence_length, args.window_size)
        elif ts.num_samples != n_samples:
            raise RuntimeError(
                f"Sample count mismatch: {ts_path} has {ts.num_samples}, expected {n_samples}."
            )
        if float(ts.sequence_length) != float(stat_windows[-1]):
            raise RuntimeError(
                f"Sequence length mismatch: {ts_path} has {ts.sequence_length}, "
                f"expected {stat_windows[-1]}."
            )

        pdf, rates = compute_pair_coal(ts, time_windows, args.tail_cutoff, keep_mask)
        pdf_vals.append(pdf)
        rate_vals.append(rates)

        load = mutational_load(ts)
        load = load / float(ts.sequence_length)
        load_vals.append(load)
        seq_lengths.append(float(ts.sequence_length))
        site_div_vals.append(ts.diversity(mode="site", windows=stat_windows))
        branch_div_vals.append(ts.diversity(mode="branch", windows=stat_windows))
        site_td_vals.append(ts.Tajimas_D(mode="site", windows=stat_windows))
        branch_td_vals.append(ts.Tajimas_D(mode="branch", windows=stat_windows))
        site_afs_vals.append(
            ts.allele_frequency_spectrum(mode="site", polarised=not args.folded)
        )
        branch_afs_vals.append(
            ts.allele_frequency_spectrum(mode="branch", polarised=not args.folded)
        )

    pdf_vals = np.stack(pdf_vals, axis=0)
    rate_vals = np.stack(rate_vals, axis=0)
    load_vals = np.stack(load_vals, axis=-1)  # [sample, replicate]
    site_div_vals = np.stack(site_div_vals, axis=-1)  # [window, replicate]
    branch_div_vals = np.stack(branch_div_vals, axis=-1)
    site_td_vals = np.stack(site_td_vals, axis=-1)
    branch_td_vals = np.stack(branch_td_vals, axis=-1)
    site_afs_vals = np.stack(site_afs_vals, axis=-1)  # [freq_bin, replicate]
    branch_afs_vals = np.stack(branch_afs_vals, axis=-1)
    seq_lengths = np.array(seq_lengths)
    coord = stat_windows[:-1] / 2.0 + stat_windows[1:] / 2.0

    mean_pdf = safe_nanmean(pdf_vals[keep_post], axis=0)
    mean_rates = safe_nanmean(rate_vals[keep_post], axis=0)
    mean_load = safe_nanmean(load_vals[:, keep_post], axis=-1)
    q_load = safe_nanquantile(load_vals[:, keep_post], [0.025, 0.975], axis=-1)
    mean_site_div = safe_nanmean(site_div_vals[:, keep_post], axis=-1)
    mean_branch_div = safe_nanmean(branch_div_vals[:, keep_post], axis=-1) * args.mutation_rate
    q_branch_div = safe_nanquantile(branch_div_vals[:, keep_post], [0.025, 0.975], axis=-1) * args.mutation_rate
    mean_site_td = safe_nanmean(site_td_vals[:, keep_post], axis=-1)
    mean_branch_td = safe_nanmean(branch_td_vals[:, keep_post], axis=-1)
    q_branch_td = safe_nanquantile(branch_td_vals[:, keep_post], [0.025, 0.975], axis=-1)
    trace_branch_div = safe_nanmean(branch_div_vals, axis=0) * args.mutation_rate
    trace_branch_td = safe_nanmean(branch_td_vals, axis=0)
    mean_site_afs = safe_nanmean(site_afs_vals[:, keep_post], axis=-1)
    mean_branch_afs = safe_nanmean(branch_afs_vals[:, keep_post], axis=-1) * args.mutation_rate
    q_branch_afs = safe_nanquantile(branch_afs_vals[:, keep_post], [0.025, 0.975], axis=-1) * args.mutation_rate

    reps_kwargs = {"color": "gray", "alpha": 0.15}
    mean_kwargs = {"color": "black", "linewidth": 1.5}

    # Pair coalescence PDF
    fig, ax = plt.subplots(1, 1, figsize=(5, 4), constrained_layout=True)
    for val in pdf_vals:
        ax.step(plot_breaks, val, **reps_kwargs)
    ax.step(plot_breaks, mean_pdf, **mean_kwargs)
    ax.set_xlabel("Adjusted generations in past")
    ax.set_ylabel("Proportion coalescing pairs")
    ax.set_xscale("log")
    pdf_path = args.out_dir / f"{args.prefix}pair-coalescence-pdf.png"
    plt.savefig(pdf_path)
    plt.clf()

    # Pair coalescence rates
    fig, ax = plt.subplots(1, 1, figsize=(5, 4), constrained_layout=True)
    for val in rate_vals:
        ax.step(plot_breaks, val, **reps_kwargs)
    ax.step(plot_breaks, mean_rates, **mean_kwargs)
    ax.set_xlabel("Adjusted generations in past")
    ax.set_ylabel("Pair coalescence rate")
    ax.set_xscale("log")
    if args.log_rates:
        ax.set_yscale("log")
    rate_path = args.out_dir / f"{args.prefix}pair-coalescence-rates.png"
    plt.savefig(rate_path)
    plt.clf()

    # Effective population size through time from pair coalescence rates.
    ne_vals = np.full_like(rate_vals, np.nan, dtype=float)
    valid_rates = np.isfinite(rate_vals) & (rate_vals > 0)
    ne_vals[valid_rates] = 1.0 / (2.0 * rate_vals[valid_rates])
    mean_ne = safe_nanmean(ne_vals[keep_post], axis=0)
    fig, ax = plt.subplots(1, 1, figsize=(5, 4), constrained_layout=True)
    for val in ne_vals:
        ax.step(plot_breaks, val, **reps_kwargs)
    ax.step(plot_breaks, mean_ne, **mean_kwargs)
    if args.year is not None:
        ax.axvline(args.year, color="red", linestyle="--", linewidth=1.2)
    ax.set_xlabel("Adjusted generations in past")
    ax.set_ylabel("Estimated Ne = 1 / (2 * coal. rate)")
    ax.set_xscale("log")
    if args.log_rates:
        ax.set_yscale("log")
    ne_path = args.out_dir / f"{args.prefix}effective-pop-size.png"
    plt.savefig(ne_path)
    plt.clf()

    # Mutational load by sample (posterior summary)
    samples = np.arange(mean_load.size)
    fig, ax = plt.subplots(1, 1, figsize=(max(5, 0.1 * samples.size), 4), constrained_layout=True)
    ax.axhline(y=float(np.nanmean(mean_load)), color="firebrick", linestyle="dashed", label="mean")
    ax.plot(samples, mean_load, "o", color="black", markersize=2)
    ax.vlines(samples, q_load[0], q_load[1], color="black", label="95% interval")
    ax.set_xlabel("Sample ID")
    ax.set_ylabel("Derived mutations / base")
    ax.legend()
    load_path = args.out_dir / f"{args.prefix}mutational-load.png"
    plt.savefig(load_path)
    plt.clf()

    # Mutational load trace across replicates
    fig, ax = plt.subplots(1, 1, figsize=(5, 4), constrained_layout=True)
    for i in range(load_vals.shape[0]):
        ax.plot(keep_idx, load_vals[i], "-", color="black", linewidth=0.5, alpha=0.2)
    ax.axvline(burnin, color="firebrick", linestyle="--", linewidth=1, label="burnin cutoff")
    ax.set_xlabel("Replicate index")
    ax.set_ylabel("Derived mutations / base in each sample")
    ax.legend()
    trace_path = args.out_dir / f"{args.prefix}mutational-load-trace.png"
    plt.savefig(trace_path)
    plt.clf()

    # Diversity: observed(site) vs expected(branch)
    fig, ax = plt.subplots(1, 1, figsize=(5, 4), constrained_layout=True)
    ax.scatter(mean_site_div, mean_branch_div, c="firebrick", s=8)
    x = float(np.nanmean(mean_site_div))
    ax.axline((x, x), slope=1.0, color="black", linestyle="dashed")
    ax.set_xlabel("Site diversity per window")
    ax.set_ylabel("Expected site diversity per window")
    div_scatter_path = args.out_dir / f"{args.prefix}diversity-scatter.png"
    plt.savefig(div_scatter_path)
    plt.clf()

    fig, ax = plt.subplots(1, 1, figsize=(8, 4), constrained_layout=True)
    ax.fill_between(coord, q_branch_div[0], q_branch_div[1], color="firebrick", alpha=0.1)
    ax.plot(coord, mean_branch_div, "-o", c="firebrick", label="branch", markersize=3)
    ax.plot(coord, mean_site_div, "-o", c="black", label="site", markersize=3)
    ax.set_xlabel("Position on chromosome")
    ax.set_ylabel("Diversity")
    ax.legend()
    div_skyline_path = args.out_dir / f"{args.prefix}diversity-skyline.png"
    plt.savefig(div_skyline_path)
    plt.clf()

    fig, ax = plt.subplots(1, 1, figsize=(5, 4), constrained_layout=True)
    ax.plot(keep_idx, trace_branch_div, "-", c="firebrick", label="branch")
    ax.set_xlabel("Replicate index")
    ax.set_ylabel("Expected genome-wide diversity")
    ax.legend()
    div_trace_path = args.out_dir / f"{args.prefix}diversity-trace.png"
    plt.savefig(div_trace_path)
    plt.clf()

    # Tajima's D: observed(site) vs expected(branch)
    fig, ax = plt.subplots(1, 1, figsize=(5, 4), constrained_layout=True)
    ax.scatter(mean_site_td, mean_branch_td, c="firebrick", s=8)
    x = float(np.nanmean(mean_site_td))
    ax.axline((x, x), slope=1.0, color="black", linestyle="dashed")
    ax.set_xlabel("Site Tajima's D per window")
    ax.set_ylabel("Branch Tajima's D per window")
    td_scatter_path = args.out_dir / f"{args.prefix}tajima-d-scatter.png"
    plt.savefig(td_scatter_path)
    plt.clf()

    fig, ax = plt.subplots(1, 1, figsize=(8, 4), constrained_layout=True)
    ax.fill_between(coord, q_branch_td[0], q_branch_td[1], color="firebrick", alpha=0.1)
    ax.plot(coord, mean_branch_td, "-o", c="firebrick", label="branch", markersize=3)
    ax.plot(coord, mean_site_td, "-o", c="black", label="site", markersize=3)
    ax.set_xlabel("Position on chromosome")
    ax.set_ylabel("Tajima's D")
    ax.legend()
    td_skyline_path = args.out_dir / f"{args.prefix}tajima-d-skyline.png"
    plt.savefig(td_skyline_path)
    plt.clf()

    fig, ax = plt.subplots(1, 1, figsize=(5, 4), constrained_layout=True)
    ax.plot(keep_idx, trace_branch_td, "-", c="firebrick", label="branch")
    ax.set_xlabel("Replicate index")
    ax.set_ylabel("Expected genome-wide Tajima's D")
    ax.legend()
    td_trace_path = args.out_dir / f"{args.prefix}tajima-d-trace.png"
    plt.savefig(td_trace_path)
    plt.clf()

    # Site vs branch SFS
    freq = np.arange(1, mean_site_afs.size)
    fig, ax = plt.subplots(1, 1, figsize=(8, 4), constrained_layout=True)
    ax.fill_between(freq, q_branch_afs[0, 1:], q_branch_afs[1, 1:], color="firebrick", alpha=0.1)
    ax.scatter(freq, mean_branch_afs[1:], c="firebrick", label="branch (scaled by mu)", s=8)
    ax.scatter(freq, mean_site_afs[1:], c="black", label="site observed", s=8)
    ax.set_xlabel("Minor allele frequency" if args.folded else "Derived allele frequency")
    ax.set_ylabel("# of variants / base")
    ax.set_yscale("log")
    ax.legend()
    sfs_path = args.out_dir / f"{args.prefix}frequency-spectrum.png"
    plt.savefig(sfs_path)
    plt.clf()

    summary_path = args.out_dir / f"{args.prefix}summary.txt"
    summary_path.write_text(
        "\n".join(
            [
                f"ts_dir={args.ts_dir}",
                f"pattern={args.pattern}",
                f"n_files={len(ts_files)}",
                f"burnin_frac={args.burnin_frac}",
                f"burnin_index={burnin}",
                f"time_bins_file={args.time_bins_file}",
                f"time_windows={time_windows.tolist()}",
                f"time_adjust={args.time_adjust}",
                f"year_marker={args.year}",
                f"tail_cutoff={args.tail_cutoff}",
                f"window_size={args.window_size}",
                f"mutation_rate={args.mutation_rate}",
                f"folded_sfs={args.folded}",
                f"n_samples={n_samples}",
                f"sequence_length_min={float(np.min(seq_lengths))}",
                f"sequence_length_max={float(np.max(seq_lengths))}",
                f"pair_coalescence_pdf_plot={pdf_path}",
                f"pair_coalescence_rates_plot={rate_path}",
                f"effective_pop_size_plot={ne_path}",
                f"mutational_load_plot={load_path}",
                f"mutational_load_trace_plot={trace_path}",
                f"diversity_scatter_plot={div_scatter_path}",
                f"diversity_skyline_plot={div_skyline_path}",
                f"diversity_trace_plot={div_trace_path}",
                f"tajima_d_scatter_plot={td_scatter_path}",
                f"tajima_d_skyline_plot={td_skyline_path}",
                f"tajima_d_trace_plot={td_trace_path}",
                f"frequency_spectrum_plot={sfs_path}",
            ]
        )
        + "\n"
    )

    print(f"Wrote: {pdf_path}")
    print(f"Wrote: {rate_path}")
    print(f"Wrote: {ne_path}")
    print(f"Wrote: {load_path}")
    print(f"Wrote: {trace_path}")
    print(f"Wrote: {div_scatter_path}")
    print(f"Wrote: {div_skyline_path}")
    print(f"Wrote: {div_trace_path}")
    print(f"Wrote: {td_scatter_path}")
    print(f"Wrote: {td_skyline_path}")
    print(f"Wrote: {td_trace_path}")
    print(f"Wrote: {sfs_path}")
    print(f"Wrote: {summary_path}")


if __name__ == "__main__":
    main()
