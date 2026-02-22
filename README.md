# ARG Tree Sequence Utilities and Validation Plotting

Standalone scripts for post-processing, QC, and visualization of ARG tree sequences (`.ts`, `.trees`, `.tsz`).

## Install

```bash
conda env create -f environment.yml
conda activate argtest
```

Core dependencies are in `environment.yml` (`numpy`, `matplotlib`, `tskit`, `tszip`, `msprime`).

## Scripts

### `scripts/mutload_summary.py`
Builds an HTML summary of per-individual mutational load and writes outlier windows as BED when windowing is enabled.

Inputs:
- one tree sequence file

Key outputs:
- `results/<name>.html`
- `results/<ts_stem>_outliers.bed` (if `--window-size` is used)
- `logs/<name>.log`

Example:
```bash
python scripts/mutload_summary.py example_data/maize.tsz --window-size 50000 --cutoff 0.25 --out mutload.html
```

Main options:
- `--window-size`
- `--cutoff`
- `--out`
- `--suffix-to-strip`

### `scripts/trim_samples.py`
Removes selected individuals either genome-wide (`--individuals`) or over BED intervals (`--remove`).

Inputs:
- one tree sequence file
- optional BED(s) with per-individual intervals

Key output:
- trimmed tree sequence (`--out`, or `results/<ts_stem>_trimmed.tsz`)

Example:
```bash
python scripts/trim_samples.py example_data/maize.tsz --individuals B73,Mo17 --out results/maize_trimmed.tsz
```

Main options:
- `--individuals`
- `--remove`
- `--out`
- `--suffix-to-strip`

### `scripts/trim_regions.py`
Applies a shared accessibility-based mask across a directory of tree sequences by:
1. inferring mutation map (`*.mut_rate.p`) from TS names,
2. keeping windows with accessible bp >= `--cutoff-bp`,
3. collapsing masked intervals and writing renamed output TS files.

Inputs:
- directory of tree sequences
- inferred mutation-rate map file(s)

Key outputs:
- collapsed tree sequences in output directory
- one summary log (`collapse_log.txt` by default)

Example:
```bash
python scripts/trim_regions.py \
  --ts-dir /path/to/trees \
  --window-size 50000 \
  --cutoff-bp 2500 \
  --pattern "*.tsz"
```

Main options:
- `--ts-dir`
- `--window-size`
- `--cutoff-bp`
- `--out-dir`
- `--pattern`
- `--log`

### `scripts/validation_plots_from_ts.py`
Generates SINGER-style validation/diagnostic plots directly from a set of TS replicates.

Plots produced:
- `pair-coalescence-pdf.png`
- `pair-coalescence-rates.png`
- `effective-pop-size.png` (`Ne = 1 / (2 * coal_rate)`)
- `mutational-load.png`
- `mutational-load-trace.png`
- `diversity-scatter.png`
- `diversity-skyline.png`
- `diversity-trace.png`
- `tajima-d-scatter.png`
- `tajima-d-skyline.png`
- `tajima-d-trace.png`
- `frequency-spectrum.png`
- `summary.txt`

Notes:
- branch diversity is scaled by `--mutation-rate` for site-vs-branch comparison
- trace plots are branch-only MCMC outcomes
- time-axis bins are read from `--time-bins-file`
- `--time-adjust` divides plotted time-axis values by a factor

Example:
```bash
python scripts/validation_plots_from_ts.py \
  --ts-dir ~/crud/collapsed \
  --pattern "*.tsz" \
  --time-bins-file /path/to/time_bins.txt \
  --window-size 100000 \
  --mutation-rate 3.3e-8 \
  --burnin-frac 0.5 \
  --time-adjust 6.19476 \
  --year 534 \
  --out-dir results/validation_plots \
  --log-rates
```

Main options:
- `--ts-dir`
- `--pattern`
- `--time-bins-file`
- `--mutation-rate`
- `--burnin-frac`
- `--tail-cutoff`
- `--time-adjust`
- `--year`
- `--window-size`
- `--folded`
- `--log-rates`
- `--out-dir`
- `--prefix`

### `scripts/compare_trees_html.py`
Renders one tree index from each of two tree sequences into a single HTML comparison.

Example:
```bash
python scripts/compare_trees_html.py a.tsz 9 b.tsz 9 --out tree_compare.html
```

### `scripts/trees_gallery_html.py`
Renders all trees from two tree sequences as a top/bottom-row scrollable HTML gallery.

Example:
```bash
python scripts/trees_gallery_html.py a.tsz b.tsz --out trees_gallery.html
```

## Shared module

`scripts/argtest_common.py` contains shared tree-sequence helpers used by multiple scripts:
- TS I/O (`load_ts`, `dump_ts`)
- mutational load/stat helpers
- trimming and masking helpers

Use this module for internal script imports.

## Repository notes

- Generated `logs/` and `results/` are git-ignored.
- `.DS_Store` is git-ignored.
