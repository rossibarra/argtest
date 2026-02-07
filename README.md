# Mutational Load Summary Script

`mutload_summary.py` generates an HTML report of derived mutational load per individual from a tree sequence.

## Requirements

- Python 3.8+
- `tskit`
- `numpy`
- `matplotlib`
- `tszip` (only required for `.tsz` inputs)

Install via conda:

```bash
conda env create -f environment.yml
conda activate argtest
```

## Usage

Basic:

```bash
python mutload_summary.py example_data/maize.tsz
```

Drop individuals by name:

```bash
python mutload_summary.py example_data/maize.tsz --out load.html
```

Drop individuals from a file (one ID per line):

```bash
python mutload_summary.py example_data/maize.tsz
```

Compute per-window load (tiled barplots across the contig):

```bash
python mutload_summary.py example_data/maize.tsz --window-size 50000 --out load_windows.html
```

Mask individuals only within specific regions (BEDs):

```bash
python mutload_summary.py example_data/maize.tsz --window-size 1000000 --remove results/beds/Ki11.bed,results/beds/Ki3.bed
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
- When `--window-size` is provided, per-individual BED files are emitted for windows
  where the individual's load is greater than (1 + `cutoff`) × the window mean or less
  than (1 - `cutoff`) × the window mean. The cutoff is a fraction of the mean for each window.
  BED files are written to `results/beds/` and the run log to `logs/`.
- `--remove` accepts one or more BED files listing regions to remove individuals from.
  If the BED has a 4th column, it is used as the individual ID; otherwise the filename stem is used.
  Individuals are removed only within the listed regions.

## Options

```text
positional arguments:
  ts                    Tree sequence file (.ts, .trees, or .tsz)

options:
  --window-size         Window size in bp
  --remove              BED file(s) of regions to remove per individual (comma-separated or repeated)
  --cutoff              Outlier cutoff as a fraction of the window mean (default: 0.25)
  --out                 Output HTML file (default: mutational_load_summary.html; written to results/)
  --suffix-to-strip     Suffix removed from individual IDs (default: _anchorwave)
```
