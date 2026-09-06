# 1TCO MM/GBSA Workflow

This directory contains the files and shell script used for the single-snapshot MM/GBSA analysis of the *Bos taurus* calcineurin–FKBP12–FK506 complex derived from the crystal structure with PDB ID **1TCO**.

## Starting structure

The analysis starts directly from the prepared protein–FK506 complex in PDB format. Hydrogen atoms and PDB connectivity records are removed during preprocessing and are subsequently regenerated according to the selected Amber force fields. The additional `MYR` residue is excluded from the MM/GBSA system.

Unlike the 6TZ6 structure, the 1TCO complex did not require reconstruction of multiple incomplete protein side chains or the extended restrained minimization procedure used for 6TZ6.

## MM/GBSA analysis

The accompanying shell script automatically:

1. cleans the input PDB file and separates FK506 (`FK5`) from the protein;
2. assigns AM1-BCC partial charges and GAFF2 parameters to FK506;
3. constructs the protein–ligand topology using ff14SB for the protein and GAFF2 for FK506;
4. generates a single-frame NetCDF coordinate file;
5. creates the complex, receptor, and ligand topology files; and
6. performs single-snapshot MM/GBSA and per-residue energy-decomposition calculations.

## Ligand residue assignment

In the TLEaP-generated 1TCO topology, FK506 is residue **629**:

```text
Ligand residue name: FK5
Ligand topology residue: 629
```

This assignment is essential. Residue 509 is a protein glycine (`GLY`) and must not be used as the ligand mask. Receptor and ligand topologies are therefore generated using `:629` as the FK506 residue selection.

## Example

Place the input PDB file and the script in the same directory, activate the AmberTools environment, and run:

```bash
conda activate ambertools

bash run_analysis_1TCO.sh
```

The default input filename defined in the script is:

```text
1TCO_FK506_complex.pdb
```

If the PDB or script has been renamed, the corresponding filename in the script must be updated before execution.

## Calculation settings

- Protein force field: ff14SB
- FK506 force field: GAFF2
- FK506 charge method: AM1-BCC
- Generalized Born model: `igb=5`
- Salt concentration: 0.150 M
- Number of analyzed structures: 1
- Per-residue decomposition: `idecomp=1`
- Configurational entropy correction: not applied

Because the calculation is based on a single prepared structure and does not include a configurational entropy term, the resulting values should be interpreted as relative MM/GBSA estimates rather than absolute experimental binding free energies.

## Output files

The principal output files are:

```text
FINAL_RESULTS_1TCO_FK506_MMPBSA.dat
FINAL_RESULTS_1TCO_FK506_DECOMP.dat
```

The first file contains the total MM/GBSA energy components and calculated binding energy. The second contains the per-residue decomposition results used to identify the principal protein residues contributing to FK506 binding.

## Application to FK506 derivatives

For derivative complexes, FK506 must first be replaced with the relevant derivative while preserving the prepared 1TCO protein structure and binding pose as the structural reference. Each derivative must then be parameterized independently with its correct atom composition, protonation state, and total charge. The ligand residue index must be verified in every newly generated topology before MM/GBSA analysis.
