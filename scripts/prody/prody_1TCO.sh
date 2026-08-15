#!/bin/bash

# ============================================================
# FK506 / Calcineurin - AmberTools Preparation
# Input file:
#   1TCO_FK506.pdb
#
# Output files:
#   FK5_1TCO.pdb
#   FK5_1TCO.mol2
#   FK5_1TCO.frcmod
#
# This script does NOT run tleap.
# ============================================================

set -e

INPUT="1TCO_FK506.pdb"
LIGAND_PDB="FK5_1TCO.pdb"
LIGAND_MOL2="FK5_1TCO.mol2"
LIGAND_FRCMOD="FK5_1TCO.frcmod"

echo "=========================================="
echo " FK506 AmberTools Preparation"
echo "=========================================="
echo

# ------------------------------------------------------------
# 1. Validate input file
# ------------------------------------------------------------

if [ ! -f "$INPUT" ]; then
    echo "ERROR: Input file not found:"
    echo "  $INPUT"
    exit 1
fi

echo "[1/7] Input PDB found"
echo "      $INPUT"
echo

# ------------------------------------------------------------
# 2. Inspect protein chains
# ------------------------------------------------------------

echo "[2/7] Protein chain inspection"
echo "------------------------------------------"

grep "^ATOM" "$INPUT" | cut -c22 | sort | uniq -c

echo

# ------------------------------------------------------------
# 3. Inspect hetero residues
# ------------------------------------------------------------

echo "[3/7] HETATM residue inspection"
echo "------------------------------------------"

grep "^HETATM" "$INPUT" | awk '{print $4}' | sort | uniq -c

echo

# ------------------------------------------------------------
# 4. Extract FK506 (residue name: FK5)
# ------------------------------------------------------------

echo "[4/7] Extracting FK506 (FK5)"
echo "------------------------------------------"

grep "^HETATM" "$INPUT" | awk '$4=="FK5" {print}' > "$LIGAND_PDB"

echo "Created:"
ls -lh "$LIGAND_PDB"

echo

# Verify atom count
FK5_COUNT=$(grep "^HETATM" "$LIGAND_PDB" | wc -l)

echo "FK5 atom count: $FK5_COUNT"

if [ "$FK5_COUNT" -ne 126 ]; then
    echo
    echo "WARNING: Expected 126 FK5 atoms, found $FK5_COUNT"
    echo "Check the input PDB before continuing."
    exit 1
fi

echo "FK5 atom count is correct."
echo

# ------------------------------------------------------------
# 5. Generate FK506 MOL2 file with Antechamber
# ------------------------------------------------------------

echo "[5/7] Running Antechamber"
echo "------------------------------------------"

antechamber \
    -i "$LIGAND_PDB" \
    -fi pdb \
    -o "$LIGAND_MOL2" \
    -fo mol2 \
    -c bcc \
    -nc 0 \
    -s 2

echo

echo "Created:"
ls -lh "$LIGAND_MOL2"

echo

# ------------------------------------------------------------
# 6. Generate FRCMOD file with Parmchk2
# ------------------------------------------------------------

echo "[6/7] Running Parmchk2"
echo "------------------------------------------"

parmchk2 \
    -i "$LIGAND_MOL2" \
    -f mol2 \
    -o "$LIGAND_FRCMOD"

echo

echo "Created:"
ls -lh "$LIGAND_FRCMOD"

echo

# ------------------------------------------------------------
# 7. Final verification
# ------------------------------------------------------------

echo "[7/7] Final verification"
echo "=========================================="

echo
echo "Input:"
ls -lh "$INPUT"

echo
echo "FK506 files:"
ls -lh "$LIGAND_PDB" "$LIGAND_MOL2" "$LIGAND_FRCMOD"

echo
echo "FK5 residue in MOL2:"
grep -n "FK5" "$LIGAND_MOL2" | head

echo
echo "=========================================="
echo " Preparation completed successfully."
echo "=========================================="
echo
echo "Next step: TLEaP"
