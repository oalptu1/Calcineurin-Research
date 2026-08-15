#!/bin/bash

set -e

# ============================================================
# USER-CONFIGURABLE SETTINGS
# ============================================================

INPUT_PDB="1TCO_FK506_complex.pdb"
LIGAND="FK5"
LIGAND_RESIDUE="629"
PREFIX="1TCO_FK506"

# ============================================================
# FILE NAMES
# ============================================================

CLEAN_PDB="${PREFIX}_clean.pdb"
PROTEIN_PDB="${PREFIX}_protein.pdb"
LIGAND_PDB="${PREFIX}_${LIGAND}.pdb"
LIGAND_MOL2="${PREFIX}_${LIGAND}.mol2"
LIGAND_FRCMOD="${PREFIX}_${LIGAND}.frcmod"

PRMTOP="${PREFIX}.prmtop"
INPCRD="${PREFIX}.inpcrd"
NC="${PREFIX}.nc"

RECEPTOR_PRMTOP="${PREFIX}_receptor.prmtop"
LIGAND_PRMTOP="${PREFIX}_ligand.prmtop"

MMPBSA_IN="mmpbsa_${PREFIX}.in"
DECOMP_IN="mmpbsa_${PREFIX}_decomp.in"

MMPBSA_OUT="FINAL_RESULTS_${PREFIX}_MMPBSA.dat"
DECOMP_OUT="FINAL_RESULTS_${PREFIX}_DECOMP.dat"

# ============================================================
echo ""
echo "============================================================"
echo "             MM/GBSA PIPELINE"
echo "============================================================"
echo ""
echo "Input PDB      : $INPUT_PDB"
echo "Ligand         : $LIGAND"
echo "Ligand residue : $LIGAND_RESIDUE"
echo "Prefix         : $PREFIX"
echo ""

# ============================================================
# 0. INPUT VALIDATION
# ============================================================

if [ ! -f "$INPUT_PDB" ]; then
    echo "ERROR: $INPUT_PDB was not found."
    exit 1
fi

echo "=== 0. Input PDB found ==="

# ============================================================
# 1. REMOVE MYR, HYDROGENS, AND CONECT RECORDS
# ============================================================

echo "=== 1. Cleaning PDB ==="

awk '
BEGIN { OFS="" }
{
    # Remove CONECT records
    if (substr($0,1,6) == "CONECT") next

    # Remove MYR residue
    if (substr($0,18,3) == "MYR") next

    # Remove hydrogen atoms
    atom = substr($0,13,4)
    gsub(/ /,"",atom)

    element = substr($0,77,2)
    gsub(/ /,"",element)

    if (atom ~ /^H/) next
    if (element == "H") next

    print
}
' "$INPUT_PDB" > "$CLEAN_PDB"

if [ ! -s "$CLEAN_PDB" ]; then
    echo "ERROR: Clean PDB could not be created."
    exit 1
fi

echo "OK: $CLEAN_PDB created."

# ============================================================
# 2. SEPARATE FK5
# ============================================================

echo "=== 2. Separating FK5 ==="

grep " ${LIGAND} " "$INPUT_PDB" | grep -v "CONECT" > "$LIGAND_PDB"

if [ ! -s "$LIGAND_PDB" ]; then
    echo "ERROR: ${LIGAND} residue was not found."
    exit 1
fi

grep -v " ${LIGAND} " "$CLEAN_PDB" > "$PROTEIN_PDB"

echo "OK: $LIGAND_PDB"
echo "OK: $PROTEIN_PDB"

echo "FK5 atom count:"
grep " ${LIGAND} " "$CLEAN_PDB" | wc -l

# ============================================================
# 3. ANTECHAMBER
# ============================================================

echo "=== 3. Running Antechamber ==="

antechamber \
-fi pdb \
-fo mol2 \
-i "$LIGAND_PDB" \
-o "$LIGAND_MOL2" \
-c bcc \
-at gaff2 \
-rn "$LIGAND" \
-s 2

if [ ! -s "$LIGAND_MOL2" ]; then
    echo "ERROR: MOL2 file could not be created."
    exit 1
fi

echo "OK: $LIGAND_MOL2"

# ============================================================
# 4. PARMCHK2
# ============================================================

echo "=== 4. Running Parmchk2 ==="

parmchk2 \
-i "$LIGAND_MOL2" \
-f mol2 \
-o "$LIGAND_FRCMOD" \
-s 2

if [ ! -s "$LIGAND_FRCMOD" ]; then
    echo "ERROR: FRCMOD file could not be created."
    exit 1
fi

echo "OK: $LIGAND_FRCMOD"

# ============================================================
# 5. TLEAP
# ============================================================

echo "=== 5. Running TLEaP ==="

cat > tleap_${PREFIX}.in <<EOF

source leaprc.protein.ff14SB
source leaprc.gaff2

loadamberparams ${LIGAND_FRCMOD}

protein = loadpdb ${PROTEIN_PDB}
ligand = loadmol2 ${LIGAND_MOL2}

complex = combine {protein ligand}

saveamberparm complex ${PRMTOP} ${INPCRD}
savepdb complex ${PREFIX}_tleap.pdb

quit
EOF

tleap -f tleap_${PREFIX}.in

if [ ! -s "$PRMTOP" ]; then
    echo "ERROR: PRMTOP file could not be created."
    exit 1
fi

if [ ! -s "$INPCRD" ]; then
    echo "ERROR: INPCRD file could not be created."
    exit 1
fi

echo "OK: $PRMTOP"
echo "OK: $INPCRD"

# ============================================================
# 6. SINGLE-FRAME TRAJECTORY
# ============================================================

echo "=== 6. Creating single-frame trajectory ==="

cat > single_frame_${PREFIX}.in <<EOF
parm ${PRMTOP}
trajin ${INPCRD}
trajout ${NC} netcdf
run
EOF

cpptraj -i single_frame_${PREFIX}.in

if [ ! -s "$NC" ]; then
    echo "ERROR: NetCDF trajectory could not be created."
    exit 1
fi

echo "OK: $NC"

# ============================================================
# 7. ANTE-MMPBSA
# ============================================================

echo "=== 7. Creating receptor/ligand topologies ==="

ante-MMPBSA.py \
-p "$PRMTOP" \
-c "$PRMTOP" \
-r "$RECEPTOR_PRMTOP" \
-l "$LIGAND_PRMTOP" \
-n ":${LIGAND_RESIDUE}"

if [ ! -s "$RECEPTOR_PRMTOP" ]; then
    echo "ERROR: Receptor topology could not be created."
    exit 1
fi

if [ ! -s "$LIGAND_PRMTOP" ]; then
    echo "ERROR: Ligand topology could not be created."
    exit 1
fi

echo "OK: Receptor and ligand topologies created."

# ============================================================
# 8. MM/GBSA INPUT
# ============================================================

echo "=== 8. Preparing MM/GBSA input ==="

cat > "$MMPBSA_IN" <<EOF
&general
  startframe=1,
  endframe=1,
  interval=1,
  verbose=1,
/

&gb
  igb=5,
  saltcon=0.150,
/
EOF

# ============================================================
# 9. MM/GBSA
# ============================================================

echo "=== 9. Running MM/GBSA calculation ==="

MMPBSA.py -O \
-i "$MMPBSA_IN" \
-o "$MMPBSA_OUT" \
-sp "$PRMTOP" \
-cp "$PRMTOP" \
-rp "$RECEPTOR_PRMTOP" \
-lp "$LIGAND_PRMTOP" \
-y "$NC"

if [ ! -s "$MMPBSA_OUT" ]; then
    echo "ERROR: MM/GBSA result file was not created."
    exit 1
fi

echo "OK: $MMPBSA_OUT"

# ============================================================
# 10. DECOMPOSITION INPUT
# ============================================================

echo "=== 10. Preparing per-residue decomposition input ==="

cat > "$DECOMP_IN" <<EOF
&general
  startframe=1,
  endframe=1,
  interval=1,
  verbose=1,
/

&gb
  igb=5,
  saltcon=0.150,
/

&decomp
  idecomp=1,
  dec_verbose=1,
/
EOF

# ============================================================
# 11. PER-RESIDUE DECOMPOSITION
# ============================================================

echo "=== 11. Running per-residue decomposition ==="

MMPBSA.py -O \
-i "$DECOMP_IN" \
-o "$MMPBSA_OUT" \
-do "$DECOMP_OUT" \
-sp "$PRMTOP" \
-cp "$PRMTOP" \
-rp "$RECEPTOR_PRMTOP" \
-lp "$LIGAND_PRMTOP" \
-y "$NC"

if [ ! -s "$DECOMP_OUT" ]; then
    echo "ERROR: Decomposition result file was not created."
    exit 1
fi

echo "OK: $DECOMP_OUT"

# ============================================================
# 12. RESULTS
# ============================================================

echo ""
echo "============================================================"
echo "                 PIPELINE COMPLETED"
echo "============================================================"
echo ""

echo "Generated files:"
echo ""

ls -lh \
"$CLEAN_PDB" \
"$PROTEIN_PDB" \
"$LIGAND_PDB" \
"$LIGAND_MOL2" \
"$LIGAND_FRCMOD" \
"$PRMTOP" \
"$INPCRD" \
"$NC" \
"$RECEPTOR_PRMTOP" \
"$LIGAND_PRMTOP" \
"$MMPBSA_OUT" \
"$DECOMP_OUT"

echo ""
echo "MM/GBSA Differences:"
grep -n "Differences" "$MMPBSA_OUT" || true

echo ""
echo "============================================================"
echo "                     DONE"
echo "============================================================"
