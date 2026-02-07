# Tree Sequence Simulations

This directory contains 10 tree sequences generated with `msprime`.

Simulation details:
- Method: `msprime.sim_ancestry`
- Replicates: 10
- Sample size: 24 haploid samples
- Population size (Ne): 100,000
- Sequence length: 5,000,000 bp (5 Mb)
- Recombination rate: 1e-8 per bp
- Random seed: replicate index (1..10)

Files are named `sim_###.tsz` (e.g., `sim_001.tsz`).

The directory also contains one tree sequence estimated from 26 maize individuals named `maize.tsz`

Mutational load output:
- `maize_mutload.html` generated from `maize.tsz` using 1,000,000 bp windows.
- Per-individual BED files are emitted alongside the HTML, listing windows where
  the individual's load is >1.25x the mean load for that window.
