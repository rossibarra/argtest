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
conda activate mutload
```

## Usage

Basic:

```bash
python mutload_summary.py path/to/trees.tsz
```

Drop individuals by name:

```bash
python mutload_summary.py path/to/trees.tsz --remove Ky21,Ky22 --out load.html
```

Drop individuals from a file (one ID per line):

```bash
python mutload_summary.py path/to/trees.tsz --remove-file remove.txt
```

Compute per-window load (tiled barplots across the contig):

```bash
python mutload_summary.py path/to/trees.tsz --window-size 50000 --out load_windows.html
```

## Inputs

- Tree sequence file: `.ts`, `.trees`, or `.tsz`.
- Individual IDs are matched against `individual.metadata["id"]` with the suffix `_anchorwave` stripped by default.

## Outputs

- An HTML file with an embedded PNG plot.
- Default output file is `mutational_load_summary.html`.

## Notes

- If `--window-size` is omitted, the report shows a single barplot of per-individual load.
- If `--window-size` is provided, the report shows one barplot per window, tiled across the contig.
- If any requested individuals are not found, the script exits with an error.

## Options

```text
positional arguments:
  ts                    Tree sequence file (.ts, .trees, or .tsz)

options:
  --remove              Comma-separated individual IDs to drop
  --remove-file         File with one individual ID per line
  --window-size         Window size in bp
  --out                 Output HTML file (default: mutational_load_summary.html)
  --suffix-to-strip     Suffix removed from individual IDs (default: _anchorwave)
```
