# AlphaFold3 protein–protein interface geometry

This repository contains analysis code and selected results for comparing protein–protein interfaces in AlphaFold 3 (AF3) models with experimental structures. The project examines interface contact size, local residue-pair distance geometry, close-contact counts, broad residue-contact classes, interface-size effects, and selected matched structural cases.

## Datasets and scope

Analyses cover Pedro AF3 predictions and the paired FoldBench and SAFES experimental/AF3 benchmarks, with a selected Mycoplasma genitalium experimental-reference analysis. Raw structure archives and benchmark source files are not redistributed here. See [`data/README.md`](data/README.md) for input expectations and the selected reports under [`results/`](results/README.md) for the analyses included in this package.

## Contact definition

Unless a report states otherwise, an inter-chain residue pair is a contact when its minimum distance between non-hydrogen atoms is at most 8 Å. Each residue pair is counted once, regardless of the number of atom pairs within the cutoff. Close contacts are defined as `3 < d_min ≤ 4 Å`. Zero-contact interfaces are retained where the analysis population includes them; reports state when a selected analysis is restricted to nonzero reference interfaces.

## Repository layout

- `src/`: parsing, mapping, contact calculation, validation, and reporting scripts
- `data/`: input-data instructions; raw structures and archives are excluded
- `results/`: selected validated tables, summaries, and figures
- `requirements.txt`: Python package dependencies

## Setup

Use Python 3.9 or newer and install the listed dependencies:

```bash
python -m pip install -r requirements.txt
```

The full Pedro archive parse is computationally intensive and requires the original input archives and HPC resources. Do not run the full parsing scripts as a quick setup check. Benchmark-specific scripts may expect the data layout described in their source comments and accompanying reports. Slurm files under `src/` reflect the original UCPH HPC workflow and may require adaptation for another cluster.

## Reproducibility and interpretation

The distributed outputs are selected summaries, not a replacement for the original raw datasets. Experimental chain pairs and model interfaces are not always independent biological replicates; see individual reports for pairing, mapping, QC, and sampling limitations. Descriptive pooled-contact distributions should not be interpreted as independent-contact statistical tests or as direct measures of binding affinity.
