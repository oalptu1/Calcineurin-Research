# ProDy Analysis Workflow

This directory contains scripts used for the structural dynamics comparison of the 1TCO and 6TZ6 calcineurin–FK506 complexes.

## 1. Structure Preparation

The preparation scripts extract FK506 (residue name `FK5`) from the input structures and generate ligand parameter files using AmberTools.

- `prepare_1TCO.sh`
- `prepare_6TZ6.sh`

These scripts perform preparation only. They do not calculate protein fluctuations.

## 2. ProDy Fluctuation Analysis

`prody_compare_1TCO_6TZ6.sh`

This script performs the comparative ProDy GNM analysis of 1TCO and 6TZ6.

It:

- builds GNM models using protein C-alpha atoms,
- calculates residue-level square fluctuations,
- matches homologous chains and residues between the two structures,
- normalizes the fluctuation profiles,
- calculates the fluctuation difference for each matched residue,
- identifies the 20 residues with the largest fluctuation differences,
- identifies the 20 residues with the smallest fluctuation differences.

### Main Outputs

- `fluctuations_1TCO.csv`
- `fluctuations_6TZ6.csv`
- `matched_residue_fluctuation_comparison.csv`
- `top20_most_different_fluctuations.csv`
- `top20_least_different_fluctuations.csv`
- `TOP_FLUCTUATION_DIFFERENCES.txt`
- `RUN_REPORT.txt`

## Usage

```bash
bash prody_compare_1TCO_6TZ6.sh
```

Default input files:

```text
1TCO_FK506.pdb
6TZ6_FK506.pdb
```

## Requirements

- Python 3
- ProDy
- NumPy
- AmberTools for the preparation scripts
