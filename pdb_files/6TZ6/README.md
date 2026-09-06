# 6TZ6–FK506 Prepared Complex

This directory contains the prepared *Candida albicans* calcineurin–FKBP12–FK506 complex derived from the crystal structure with PDB ID **6TZ6**.

## Structure preparation

Inspection of the original crystal structure revealed incomplete side chains in several protein residues. In some cases, only the first side-chain carbon atom was resolved. These missing side-chain heavy atoms were rebuilt during structure preparation with **AmberTools**, using standard amino-acid residue templates.

The experimentally resolved protein and ligand heavy atoms were retained as the structural reference while the rebuilt side-chain atoms and added hydrogen atoms were allowed to relax. The prepared complex was then subjected to restrained energy minimization to remove unfavorable contacts and optimize the newly generated coordinates without substantially altering the experimental binding geometry.

The genuine chain discontinuities present in 6TZ6 were preserved with `TER` records during topology generation. This prevents artificial peptide bonds from being created across unresolved sequence gaps.

## Quality control

The prepared and minimized structure was checked for:

- preservation of the three genuine sequence gaps in the 6TZ6 protein chains;
- absence of artificial peptide bonds across those gaps;
- absence of inter-residue heavy-atom contacts shorter than 1.0 Å;
- preservation of the experimentally resolved protein and FK506 heavy-atom geometry; and
- physically reasonable van der Waals energies after minimization.

The final minimized structure was subsequently processed with the **Protein Preparation Wizard in Maestro** before PLIP interaction analysis. This step assigned chemically appropriate protonation and hydrogen-bonding states and restored the expected protein–ligand hydrogen-bond network in the PLIP results.

## MM/GBSA analysis

The accompanying shell script starts directly from the prepared PDB file and automatically:

1. separates FK506 (`FK5`) from the protein;
2. assigns AM1-BCC charges and GAFF2 parameters to FK506;
3. builds the protein–ligand topology with ff14SB and GAFF2;
4. verifies the expected 6TZ6 chain breaks and performs structural quality-control checks;
5. generates the complex, receptor, and ligand topologies; and
6. performs single-snapshot MM/GBSA and per-residue decomposition calculations.

Example:

```bash
conda activate ambertools

bash run_analysis_6TZ6_from_PDB.sh \
  6TZ6_FK506_MASTER_MAESTRO_PREP.pdb \
  6TZ6_FK506
```

The calculation uses the GB model with `igb=5` and a salt concentration of 0.150 M. No configurational entropy correction is applied. Consequently, the reported values should be interpreted as relative single-structure MM/GBSA estimates rather than absolute binding free energies.

## Important note

The prepared PDB file, rather than the incomplete crystallographic coordinates, should be used as the starting structure for subsequent derivative construction, PLIP analysis, and MM/GBSA calculations. When FK506 is replaced by another derivative, the ligand must be parameterized separately and the resulting complex must pass the same topology, steric-clash, and van der Waals quality-control checks.
