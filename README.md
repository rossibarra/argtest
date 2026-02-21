# Mutational Load Summary Script

`scripts/mutload_summary.py` generates an HTML report of derived mutational load per individual from a tree sequence and writes a BED of outlier windows when a window size is provided.
`scripts/trim_samples.py` removes individuals over BED intervals or by name and writes a trimmed tree sequence.
`scripts/trim_regions.py` collapses low-accessibility windows (from mutation-map accessibility) across a directory of tree sequences and writes renamed output tree files.
`scripts/compare_trees_html.py` renders two specific trees from two tree sequences into a single HTML report.
`scripts/trees_gallery_html.py` renders all trees from two tree sequences into a horizontally scrollable HTML gallery.

## Requirements

- Python 3.8+
- `tskit`
- `numpy`
- `matplotlib`
- `tszip` (only required for `.tsz` inputs)
- `msprime` (required for accessibility-mask collapse workflows)

Install via conda:

```bash
conda env create -f environment.yml
conda activate argtest
```

## Usage

Basic:

```bash
python scripts/mutload_summary.py example_data/maize.tsz
```


Compute per-window load (tiled barplots across the contig):

```bash
python scripts/mutload_summary.py example_data/maize.tsz --window-size 50000 --out load_windows.html
```

Remove individuals from the treesequence only within specific regions (BEDs):

```bash
python scripts/mutload_summary.py example_data/maize.tsz --window-size 1000000 --cutoff 0.5
python scripts/trim_samples.py example_data/maize.tsz --remove results/maize_outliers.bed
```

Remove specific individuals everywhere:

```bash
python scripts/trim_samples.py example_data/maize.tsz --individuals B73,Mo17
```

Collapse low-accessibility windows across a directory of trees (shared window filter, one summary log):

```bash
python scripts/trim_regions.py \
  --ts-dir /path/to/trees \
  --window-size 50000 \
  --cutoff-bp 2500 \
  --pattern "*.tsz"
```

Compare two trees by index:

```bash
python scripts/compare_trees_html.py natefun.tsz 9 vanilla.tsz 9 --out tree_compare.html
```

Render full galleries (top and bottom rows):

```bash
python scripts/trees_gallery_html.py natefun.tsz vanilla.tsz --out trees_gallery.html
```

## Inputs

- Tree sequence file: `.ts`, `.trees`, or `.tsz`.
- Individual IDs are matched against `individual.metadata["id"]` with the suffix `_anchorwave` stripped by default.

## Outputs

- An HTML file with an embedded PNG plot.
- Default output file is `results/mutational_load_summary.html`.

## Notes

- If `--window-size` is omitted, the report shows a single barplot of per-individual load.
- If `--window-size` is provided, the report shows one barplot per window, tiled across the contig.
- When `--window-size` is provided, a single BED file is emitted listing outliers per window
  where an individual's load is greater than (1 + `cutoff`) × the window mean or less than
  (1 - `cutoff`) × the window mean. The cutoff is a fraction of the mean for each window.
  The BED includes columns: `chrom`, `start`, `end`, `outlier_ids`, `outlier_values`, `window_mean`.
  Output is written to `results/` and the run log to `logs/`.
- `trim_samples.py` accepts one or more BED files listing regions where individuals are removed from the tree sequence.
  If the BED has a 4th column, it is used as the individual ID (comma-separated IDs supported); otherwise the filename stem is used.
  Individuals are removed only within the listed regions.
- Shared helpers now live in `scripts/argtest_common.py` (tree I/O, mutational-load helpers, trim helpers, and collapse helpers).
- Internal imports should use `argtest_common` (the previous `mutload_common.py` module has been removed).
- `scripts/trim_regions.py` now infers mutation-map files from TS filenames, computes one shared window mask, applies it to all files, and writes a single summary log.
- Generated `logs/` and `results/` outputs are ignored by git.

## Options

`mutload_summary.py`:

```text
positional arguments:
  ts                    Tree sequence file (.ts, .trees, or .tsz)

options:
  --window-size         Window size in bp
  --cutoff              Outlier cutoff as a fraction of the window mean (default: 0.25)
  --out                 Output HTML file (default: mutational_load_summary.html; written to results/)
  --suffix-to-strip     Suffix removed from individual IDs (default: _anchorwave)
```

`trim_samples.py`:

```text
positional arguments:
  ts                    Tree sequence file (.ts, .trees, or .tsz)

options:
  --individuals         Comma-separated individual IDs to remove across the entire sequence
  --remove              BED file(s) of regions to remove per individual (comma-separated or repeated)
  --out                 Output tree sequence path (.ts, .trees, or .tsz)
  --suffix-to-strip     Suffix removed from individual IDs (default: _anchorwave)
```

`trim_regions.py`:

```text
options:
  --ts-dir              Directory containing tree sequence files (.tsz, .ts, .trees)
  --window-size         Window size in bp for accessibility checks
  --cutoff-bp           Minimum accessible bp per window to keep
  --out-dir             Output directory (default: <ts-dir>/collapsed)
  --pattern             Glob pattern to filter inputs (default: *)
  --log                 Log path (default: <out-dir>/collapse_log.txt)
```

`compare_trees_html.py`:

```text
positional arguments:
  ts_a                  First tree sequence file (.ts, .trees, or .tsz)
  index_a               Tree index for the first file (0-based)
  ts_b                  Second tree sequence file (.ts, .trees, or .tsz)
  index_b               Tree index for the second file (0-based)

options:
  --out                 Output HTML file (default: tree_compare.html)
```

`trees_gallery_html.py`:

```text
positional arguments:
  ts_top                Top tree sequence file (.ts, .trees, or .tsz)
  ts_bottom             Bottom tree sequence file (.ts, .trees, or .tsz)

options:
  --out                 Output HTML file (default: trees_gallery.html)
```
