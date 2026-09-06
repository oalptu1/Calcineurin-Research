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
