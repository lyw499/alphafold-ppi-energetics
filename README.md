# alphafold-ppi-energetics
Analysis of AlphaFold predictions and biophysical preferences in protein-protein interactions.

## Project goal

This project investigates whether AlphaFold predictions capture known
biophysical energetic preferences in protein–protein interactions.

## Current workflow

1. Collect candidate protein–protein interactions.
2. Obtain AlphaFold 3 prediction results.
3. Extract AlphaFold confidence metrics.
4. Obtain STRING interaction scores.
5. Select high-confidence interactions using both sources.
6. Analyze structural and energetic properties.

## Repository structure

- `src/`: reusable Python scripts
- `notebooks/`: exploratory analyses
- `data/raw/`: original downloaded data
- `data/processed/`: processed datasets
- `results/`: analysis results
- `figures/`: generated figures

## Status

Project setup and data-selection workflow in progress.