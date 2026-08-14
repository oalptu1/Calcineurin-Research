# Calcineurin-Research
# FK506–Calcineurin Computational Analysis

This repository contains the computational structures, analysis scripts, and AMBER-based workflows used to investigate FK506 and FK506-derived compounds in complex with calcineurin.

## Contents

The repository includes:

- **PDB structures** used for structural and computational analyses
- **ProDy scripts** for structural dynamics and residue-level analysis
- **AMBER scripts** for molecular dynamics system preparation and analysis
- **Antechamber and tleap workflows** for ligand parameterization and system preparation
- **MM/GBSA calculations** for binding-energy analysis
- **Per-residue energy decomposition** for identifying residue-level energetic contributions
- Selected computational results and intermediate files generated during the analyses

## Computational Workflow

The general workflow includes:

1. Preparation and inspection of PDB structures
2. Protein and ligand preparation
3. Ligand parameterization using Antechamber
4. AMBER topology and coordinate generation using tleap
5. Molecular dynamics trajectory preparation and analysis
6. MM/GBSA binding-energy calculations
7. Per-residue energy decomposition
8. Structural and residue-level analysis using ProDy and related tools

## Scientific Objective

The computational analyses are used to investigate how structural modifications of FK506 affect its interactions with calcineurin.

Particular attention is given to residue-level interactions within the calcineurin B subunit and to changes in energetic contributions associated with FK506 derivatives.

## Software

The workflows make use of:

- AMBER / AmberTools
- Antechamber
- tleap
- cpptraj
- ProDy
- Python
- PyMOL

## Repository Status

This repository is under active development. Scripts and computational files may be updated as the analysis progresses.

The current version primarily documents the computational workflow and the structures used in the study.

## License

No license is currently specified for this repository.
