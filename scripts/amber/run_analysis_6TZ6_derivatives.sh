#!/bin/bash

set -e

# ============================================================
# USER-CONFIGURABLE SETTINGS
# ============================================================

INPUT_PDB="6TZ6_derivative_2_complex.pdb"
LIGAND="FK5"
LIGAND_RESIDUE="619"
PREFIX="6TZ6_derivative_2"

# ============================================================
# FILE NAMES
# ============================================================

CLEAN_PDB="${PREFIX}_clean.pdb"
PROTEIN_PDB="${PREFIX}_protein.pdb"
LIGAND_PDB="${PREFIX}_${LIGAND}.pdb"
LIGAND_MOL2="${PREFIX}_${LIGAND}.mol2"
LIGAND_FRCMOD="${PREFIX}_${LIGAND}.frcmod"

PRMTOP="${PREFIX}_REBUILT.prmtop"
INPCRD="${PREFIX}_REBUILT.inpcrd"
TLEAP_PDB="${PREFIX}_REBUILT_tleap.pdb"

RECEPTOR_PRMTOP="${PREFIX}_RECEPTOR.prmtop"
RECEPTOR_INPCRD="${PREFIX}_RECEPTOR.inpcrd"
LIGAND_PRMTOP="${PREFIX}_LIGAND.prmtop"
LIGAND_INPCRD="${PREFIX}_LIGAND.inpcrd"

NC="${PREFIX}_oneframe.nc"

MMPBSA_IN="mmpbsa_${PREFIX}.in"
DECOMP_IN="mmpbsa_${PREFIX}_decomp.in"
MMPBSA_OUT="FINAL_RESULTS_${PREFIX}_MMPBSA.dat"
DECOMP_OUT="FINAL_RESULTS_${PREFIX}_DECOMP.dat"

# ============================================================
# FIXED ANALYSIS SETTINGS
# ============================================================

VERSION="1.3.0"
PROFILE="6TZ6"

LIGAND_CHARGE="0"

# Residues of interest from the 6TZ6 interaction table.
# Format: CHAIN:RESIDUE:EXPECTED_AMINO_ACID
# Repeated contacts in the source table are intentionally represented once here.
DISTANCE_RESIDUES="C:30:TYR,C:40:PHE,C:50:PHE,C:59:VAL,C:63:TRP,C:97:TYR,C:106:ILE,C:114:PHE,B:121:MET,B:122:VAL,A:401:TRP,A:405:PHE"
HYDROGEN_BOND_RESIDUES="C:30:TYR,C:41:ASP,C:58:GLN,C:60:ILE,C:97:TYR,A:401:TRP,A:404:PRO,A:408:GLU"

# Unique PDB positions used only for mapping/QC.
# Per-residue decomposition itself is performed for ALL residues.
PDB_TARGETS="C:30,C:40,C:50,C:59,C:63,C:97,C:106,C:114,B:121,B:122,A:401,A:405,C:41,C:58,C:60,A:404,A:408"

PROTEIN_LEAPRC="leaprc.protein.ff14SB"
GAFF_TYPE="gaff2"
LIGAND_LEAPRC="leaprc.${GAFF_TYPE}"

IGB="5"
SALTCON="0.15"
SURFTEN="0.0072"
SURFOFF="0.0"

EXPECTED_PROTEIN_RESIDUES="618"
REFERENCE_DELTA_TOTAL="-92.9522"
REFERENCE_TOLERANCE="0.10"
REFERENCE_CHECK="0"

# ============================================================
# HELPER FUNCTIONS
# ============================================================

die() {
    echo "ERROR: $*" >&2
    exit 1
}

warn() {
    echo "WARNING: $*" >&2
}

info() {
    echo "[$(date +'%H:%M:%S')] $*"
}

command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# ----------------------------- PROGRAMS ------------------------------------
REQUIRED_PROGRAMS=(
    python3
    antechamber
    parmchk2
    tleap
    cpptraj
    MMPBSA.py
)

for prog in "${REQUIRED_PROGRAMS[@]}"; do
    command_exists "$prog" || die \
        "Required program not found: $prog. Activate AmberTools first."
done


# ============================================================
# 0. INPUT VALIDATION AND RUN DIRECTORIES
# ============================================================

if [ ! -f "$INPUT_PDB" ]; then
    echo "ERROR: $INPUT_PDB was not found."
    exit 1
fi

INPUT_PDB="$(cd "$(dirname "$INPUT_PDB")" && pwd)/$(basename "$INPUT_PDB")"

BASE_DIR="$(pwd)"
RUN_ID="$(date +%Y%m%d_%H%M%S)"
OUTPUT_ROOT="${BASE_DIR}/${PREFIX}_RESULTS"

RUN_DIR="${OUTPUT_ROOT}/run_${RUN_ID}"
INPUT_DIR="${RUN_DIR}/input"
PREP_DIR="${RUN_DIR}/prepared"
WORK_DIR="${RUN_DIR}/work"
TABLE_DIR="${RUN_DIR}/tables"
LOG_DIR="${RUN_DIR}/logs"

mkdir -p "$INPUT_DIR" "$PREP_DIR" "$WORK_DIR" "$TABLE_DIR" "$LOG_DIR"
cp "$INPUT_PDB" "${INPUT_DIR}/original_input.pdb"

# Bind file names to the current run directories.
CLEAN_PDB="${PREP_DIR}/${CLEAN_PDB}"
PROTEIN_PDB="${PREP_DIR}/${PROTEIN_PDB}"
LIGAND_PDB="${PREP_DIR}/${LIGAND_PDB}"
LIGAND_MOL2="${PREP_DIR}/${LIGAND_MOL2}"
LIGAND_FRCMOD="${PREP_DIR}/${LIGAND_FRCMOD}"

COMPLEX_TOP="${PREP_DIR}/${PRMTOP}"
COMPLEX_INPCRD="${PREP_DIR}/${INPCRD}"
TLEAP_PDB="${PREP_DIR}/${TLEAP_PDB}"
RECEPTOR_TOP="${PREP_DIR}/${RECEPTOR_PRMTOP}"
RECEPTOR_INPCRD="${PREP_DIR}/${RECEPTOR_INPCRD}"
LIGAND_TOP="${PREP_DIR}/${LIGAND_PRMTOP}"
LIGAND_INPCRD="${PREP_DIR}/${LIGAND_INPCRD}"

MAPPING_TSV="${TABLE_DIR}/pdb_mapping_by_protein_order.tsv"
REMOVED_HETERO_TSV="${TABLE_DIR}/removed_hetero.tsv"
TARGET_MAP_TSV="${TABLE_DIR}/requested_pdb_targets.tsv"

ONEFRAME_TRAJ="${WORK_DIR}/${NC}"
MMPBSA_INPUT="${WORK_DIR}/${MMPBSA_IN}"
DECOMP_INPUT="${WORK_DIR}/${DECOMP_IN}"
FINAL_RESULTS="${RUN_DIR}/${MMPBSA_OUT}"
FINAL_DECOMP="${RUN_DIR}/${DECOMP_OUT}"

SUMMARY_TSV="${TABLE_DIR}/MMGBSA_summary.tsv"
DECOMP_TSV="${TABLE_DIR}/decomposition_mapped.tsv"
INTEREST_TSV="${TABLE_DIR}/residues_of_interest.tsv"
RAW_FALLBACK_TSV="${TABLE_DIR}/decomposition_raw_fallback.tsv"
REPORT_TXT="${RUN_DIR}/RUN_REPORT.txt"
METADATA_TXT="${RUN_DIR}/RUN_METADATA.txt"

echo "=================================================================="
echo "AmberTools PDB -> MM/GBSA + decomposition pipeline"
echo "Version       : $VERSION"
echo "Profile       : $PROFILE"
echo "Input         : $INPUT_PDB"
echo "Run directory : $RUN_DIR"
echo "=================================================================="

# ------------------------ STRUCTURE PREPARATION -----------------------------
info "Preparing structure and building original-PDB mapping..."

STRUCTURE_INFO="$(
python3 - \
    "$INPUT_PDB" \
    "$LIGAND" \
    "$CLEAN_PDB" \
    "$PROTEIN_PDB" \
    "$LIGAND_PDB" \
    "$MAPPING_TSV" \
    "$REMOVED_HETERO_TSV" \
    "$TARGET_MAP_TSV" \
    "$PDB_TARGETS" \
    "$EXPECTED_PROTEIN_RESIDUES" <<'PY'
import csv
import re
import sys
from pathlib import Path

(
    input_s,
    ligand_name,
    clean_s,
    protein_s,
    ligand_s,
    mapping_s,
    removed_s,
    target_map_s,
    pdb_targets_s,
    expected_count_s,
) = sys.argv[1:]

input_path = Path(input_s)
clean_path = Path(clean_s)
protein_path = Path(protein_s)
ligand_path = Path(ligand_s)
mapping_path = Path(mapping_s)
removed_path = Path(removed_s)
target_map_path = Path(target_map_s)

ligand_name = ligand_name.upper()
expected_count = int(expected_count_s) if expected_count_s else None

lines = input_path.read_text(errors="replace").splitlines(True)

def atom_fields(line):
    if len(line) < 27:
        return None
    record = line[0:6].strip()
    if record not in {"ATOM", "HETATM"}:
        return None
    atom_name = line[12:16].strip()
    altloc = line[16:17].strip()
    resname = line[17:20].strip().upper()
    chain = line[21:22].strip() or "_"
    resseq = line[22:26].strip()
    icode = line[26:27].strip()
    serial = line[6:11].strip()
    return record, serial, atom_name, altloc, resname, chain, resseq, icode

# Keep blank/A alternate locations only.
protein_lines = []
ligand_lines = []
clean_lines = []
removed = []
ligand_residue_keys = set()

# For mapping: unique protein residues in exact input order.
protein_residues = []
protein_seen = set()

for line in lines:
    f = atom_fields(line)
    if f is None:
        continue

    record, serial, atom_name, altloc, resname, chain, resseq, icode = f
    if altloc not in {"", "A"}:
        continue

    is_ligand = resname == ligand_name

    if is_ligand:
        ligand_residue_keys.add((chain, resseq, icode, resname))

        # Standardize ligand residue to number 1 for antechamber.
        chars = list(line.rstrip("\n").ljust(80))
        chars[17:20] = list(ligand_name[:3].rjust(3))
        chars[21:22] = list("L")
        chars[22:26] = list(f"{1:4d}")
        chars[26:27] = list(" ")
        ligand_line = "".join(chars).rstrip() + "\n"

        ligand_lines.append(ligand_line)
        clean_lines.append(line)
        continue

    if record == "ATOM":
        # Protein hydrogens from prepared PDB files can use atom names that
        # do not match Amber residue templates (e.g. HB1/HB2/HB3, HXT, H1/H2).
        # Match the validated 1TCO workflow: remove protein hydrogens before
        # tleap and let the Amber protein force field rebuild them.
        element = line[76:78].strip().upper() if len(line) >= 78 else ""
        atom_upper = atom_name.upper()
        is_hydrogen = (
            element == "H"
            or atom_upper.startswith("H")
            or (atom_upper[:1].isdigit() and len(atom_upper) > 1 and atom_upper[1] == "H")
        )
        if is_hydrogen:
            continue

        protein_lines.append(line)
        clean_lines.append(line)
        key = (chain, resseq, icode, resname)
        if key not in protein_seen:
            protein_seen.add(key)
            protein_residues.append(key)
    else:
        removed.append((resname, chain, resseq, icode, atom_name, serial))

if not ligand_lines:
    # Give a useful inventory to the user.
    hetero = []
    seen = set()
    for line in lines:
        f = atom_fields(line)
        if f and f[0] == "HETATM":
            key = (f[4], f[5], f[6], f[7])
            if key not in seen:
                seen.add(key)
                hetero.append(key)
    print(
        f"ERROR: ligand residue name {ligand_name!r} was not found. "
        f"HETATM residues seen: {hetero}",
        file=sys.stderr,
    )
    sys.exit(10)

if len(ligand_residue_keys) != 1:
    print(
        f"ERROR: expected exactly one ligand residue named {ligand_name}; "
        f"found {len(ligand_residue_keys)}: {sorted(ligand_residue_keys)}",
        file=sys.stderr,
    )
    sys.exit(11)

if not protein_lines:
    print("ERROR: no protein ATOM records were found.", file=sys.stderr)
    sys.exit(12)

if expected_count is not None and len(protein_residues) != expected_count:
    print(
        f"ERROR: profile expects {expected_count} protein residues, "
        f"but input contains {len(protein_residues)} unique ATOM residues.",
        file=sys.stderr,
    )
    sys.exit(13)

for path, selected in (
    (clean_path, clean_lines),
    (protein_path, protein_lines),
    (ligand_path, ligand_lines),
):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        fh.writelines(selected)
        fh.write("END\n")

with mapping_path.open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(
        ["Topology_Residue", "PDB_Chain", "PDB_Residue", "PDB_ICode", "PDB_Name"]
    )
    for idx, (chain, resseq, icode, resname) in enumerate(protein_residues, 1):
        w.writerow([idx, chain, resseq, icode, resname])

with removed_path.open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["Residue_Name", "Chain", "Residue_Number", "ICode", "Atom", "Serial"])
    w.writerows(removed)

mapping_by_pdb = {}
for idx, (chain, resseq, icode, resname) in enumerate(protein_residues, 1):
    mapping_by_pdb[(chain, resseq, icode)] = (idx, resname)

target_rows = []
target_indices = []
for raw in [x.strip() for x in pdb_targets_s.split(",") if x.strip()]:
    m = re.fullmatch(r"([^:]+):(-?\d+)([A-Za-z]?)", raw)
    if not m:
        print(
            f"ERROR: invalid --pdb-targets item {raw!r}; expected CHAIN:RESNUM "
            f"(example B:115).",
            file=sys.stderr,
        )
        sys.exit(14)
    chain, resseq, icode = m.group(1), m.group(2), m.group(3)
    key = (chain, resseq, icode)
    if key not in mapping_by_pdb:
        print(f"ERROR: requested PDB target not found: {raw}", file=sys.stderr)
        sys.exit(15)
    idx, resname = mapping_by_pdb[key]
    target_indices.append(idx)
    target_rows.append((raw, idx, resname))

with target_map_path.open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["Requested_PDB_Target", "Topology_Residue", "PDB_Name"])
    w.writerows(target_rows)

ligand_key = next(iter(ligand_residue_keys))
print(f"PROTEIN_RES_COUNT={len(protein_residues)}")
print(f"LIGAND_ORIGINAL_CHAIN={ligand_key[0]}")
print(f"LIGAND_ORIGINAL_RESID={ligand_key[1]}")
print(f"TARGET_TOPOLOGY_RESIDUES={','.join(map(str, target_indices))}")
print(f"REMOVED_HETERO_ATOMS={len(removed)}")
PY
)" || die "Structure preparation failed."

# Import only known key=value lines from our own Python output.
PROTEIN_RES_COUNT=""
LIGAND_ORIGINAL_CHAIN=""
LIGAND_ORIGINAL_RESID=""
TARGET_TOPOLOGY_RESIDUES=""
REMOVED_HETERO_ATOMS=""

while IFS='=' read -r key value; do
    case "$key" in
        PROTEIN_RES_COUNT) PROTEIN_RES_COUNT="$value" ;;
        LIGAND_ORIGINAL_CHAIN) LIGAND_ORIGINAL_CHAIN="$value" ;;
        LIGAND_ORIGINAL_RESID) LIGAND_ORIGINAL_RESID="$value" ;;
        TARGET_TOPOLOGY_RESIDUES) TARGET_TOPOLOGY_RESIDUES="$value" ;;
        REMOVED_HETERO_ATOMS) REMOVED_HETERO_ATOMS="$value" ;;
    esac
done <<< "$STRUCTURE_INFO"

[[ -n "$PROTEIN_RES_COUNT" ]] || die "Could not determine protein residue count."
[[ -n "$TARGET_TOPOLOGY_RESIDUES" ]] || die "Could not map requested PDB targets."

CALCULATED_LIGAND_RESIDUE=$((PROTEIN_RES_COUNT + 1))
if [[ "$CALCULATED_LIGAND_RESIDUE" -ne "$LIGAND_RESIDUE" ]]; then
    die "Ligand residue QC failed: expected ${LIGAND_RESIDUE}, calculated ${CALCULATED_LIGAND_RESIDUE}."
fi

RECEPTOR_MASK=":1-${PROTEIN_RES_COUNT}"
LIGAND_MASK=":${LIGAND_RESIDUE}"

# Validate the amino-acid names supplied in the residues-of-interest list
# against the actual input PDB BEFORE the expensive ligand parameterization.
python3 - "$TARGET_MAP_TSV" "$DISTANCE_RESIDUES" "$HYDROGEN_BOND_RESIDUES" <<'PY'
import csv
import sys

target_map_s, distance_s, hbond_s = sys.argv[1:]

actual = {}
with open(target_map_s) as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        actual[row["Requested_PDB_Target"].upper()] = row["PDB_Name"].upper()

expected = {}
for spec_string in (distance_s, hbond_s):
    for token in [x.strip() for x in spec_string.split(",") if x.strip()]:
        chain, residue, name = token.split(":")
        key = f"{chain.upper()}:{residue}"
        if key in expected and expected[key] != name.upper():
            raise SystemExit(
                f"ERROR: conflicting expected names for {key}: "
                f"{expected[key]} vs {name.upper()}"
            )
        expected[key] = name.upper()

errors = []
for key, expected_name in expected.items():
    actual_name = actual.get(key)
    if actual_name is None:
        errors.append(f"{key}: not found in PDB mapping")
    elif actual_name != expected_name:
        errors.append(
            f"{key}: expected {expected_name}, actual PDB residue is {actual_name}"
        )

if errors:
    print("ERROR: residues-of-interest QC failed:", file=sys.stderr)
    for item in errors:
        print("  " + item, file=sys.stderr)
    sys.exit(1)

print(f"Residues-of-interest QC PASS: {len(expected)} unique residues.")
PY

info "Protein residues: $PROTEIN_RES_COUNT"
info "Ligand topology residue: $LIGAND_RESIDUE (QC PASS)"
info "Per-residue decomposition scope: ALL residues"
info "Residues of interest: 17 unique PDB residues"
info "Removed non-ligand HETATM atoms: $REMOVED_HETERO_ATOMS"

# ------------------------ LIGAND PARAMETERIZATION ---------------------------
info "Parameterizing ligand with antechamber (${GAFF_TYPE}, AM1-BCC)..."

# Antechamber/SQM uses fixed temporary names such as sqm.in and sqm.out.
# Run it inside this run's own work directory to prevent collisions between runs.
ANTE_WORK="${WORK_DIR}/antechamber"
mkdir -p "$ANTE_WORK"
cp "$LIGAND_PDB" "${ANTE_WORK}/ligand_input.pdb"

(
    cd "$ANTE_WORK"
    set +e
    antechamber \
        -i ligand_input.pdb -fi pdb \
        -o ligand_output.mol2 -fo mol2 \
        -at "$GAFF_TYPE" \
        -c bcc \
        -nc "$LIGAND_CHARGE" \
        -rn "$LIGAND" \
        -s 2 \
        > "${LOG_DIR}/antechamber.log" 2>&1
    status=$?
    set -e
    exit "$status"
)
ANTE_STATUS=$?

if [[ $ANTE_STATUS -ne 0 || ! -s "${ANTE_WORK}/ligand_output.mol2" ]]; then
    tail -n 50 "${LOG_DIR}/antechamber.log" >&2 || true
    die "antechamber failed."
fi

cp "${ANTE_WORK}/ligand_output.mol2" "$LIGAND_MOL2"

(
    cd "$ANTE_WORK"
    set +e
    parmchk2 \
        -i ligand_output.mol2 -f mol2 \
        -o ligand_output.frcmod \
        -s "$GAFF_TYPE" \
        > "${LOG_DIR}/parmchk2.log" 2>&1
    status=$?
    set -e
    exit "$status"
)
PARM_STATUS=$?

if [[ $PARM_STATUS -ne 0 || ! -s "${ANTE_WORK}/ligand_output.frcmod" ]]; then
    tail -n 50 "${LOG_DIR}/parmchk2.log" >&2 || true
    die "parmchk2 failed."
fi

cp "${ANTE_WORK}/ligand_output.frcmod" "$LIGAND_FRCMOD"

# ------------------------------- TLEAP -------------------------------------
info "Building complex/receptor/ligand topologies with tleap..."

# Work in PREP_DIR so tleap input uses simple filenames and avoids path parsing issues.
cp "$PROTEIN_PDB" "${PREP_DIR}/protein_for_tleap.pdb"
cp "$LIGAND_MOL2" "${PREP_DIR}/ligand_for_tleap.mol2"
cp "$LIGAND_FRCMOD" "${PREP_DIR}/ligand_for_tleap.frcmod"

cat > "${PREP_DIR}/tleap_build.in" <<EOF
source ${PROTEIN_LEAPRC}
source ${LIGAND_LEAPRC}

loadamberparams ligand_for_tleap.frcmod

PROT = loadpdb protein_for_tleap.pdb
LIG = loadmol2 ligand_for_tleap.mol2
COMPLEX = combine { PROT LIG }

check COMPLEX

saveamberparm COMPLEX $(basename "$COMPLEX_TOP") $(basename "$COMPLEX_INPCRD")
savepdb COMPLEX $(basename "$TLEAP_PDB")

saveamberparm PROT $(basename "$RECEPTOR_TOP") $(basename "$RECEPTOR_INPCRD")
saveamberparm LIG $(basename "$LIGAND_TOP") $(basename "$LIGAND_INPCRD")

quit
EOF

(
    cd "$PREP_DIR"
    set +e
    tleap -f tleap_build.in > "${LOG_DIR}/tleap.log" 2>&1
    status=$?
    set -e
    exit "$status"
) || {
    tail -n 80 "${LOG_DIR}/tleap.log" >&2 || true
    die "tleap failed."
}

for f in \
    "$COMPLEX_TOP" "$COMPLEX_INPCRD" "$TLEAP_PDB" \
    "$RECEPTOR_TOP" "$LIGAND_TOP"; do
    [[ -s "$f" ]] || die "tleap did not create expected file: $f"
done

# Validate that combine {PROT LIG} produced protein first and ligand last.
python3 - "$TLEAP_PDB" "$LIGAND" "$PROTEIN_RES_COUNT" <<'PY'
import sys
from pathlib import Path

pdb = Path(sys.argv[1])
ligand_name = sys.argv[2].upper()
protein_count = int(sys.argv[3])

seen = []
keys = set()
for line in pdb.read_text(errors="replace").splitlines():
    if not (line.startswith("ATOM  ") or line.startswith("HETATM")):
        continue
    resname = line[17:20].strip().upper()
    chain = line[21:22].strip() or "_"
    resseq = line[22:26].strip()
    icode = line[26:27].strip()
    key = (chain, resseq, icode, resname)
    if key not in keys:
        keys.add(key)
        seen.append(key)

if len(seen) != protein_count + 1:
    print(
        f"ERROR: tleap PDB has {len(seen)} residues; "
        f"expected {protein_count + 1}.",
        file=sys.stderr,
    )
    sys.exit(20)

last = seen[-1]
if last[3] != ligand_name:
    print(
        f"ERROR: last tleap residue is {last}, expected ligand {ligand_name}.",
        file=sys.stderr,
    )
    sys.exit(21)

try:
    tleap_resid = int(last[1])
except ValueError:
    print(f"ERROR: tleap ligand residue number is not numeric: {last[1]!r}", file=sys.stderr)
    sys.exit(22)

if tleap_resid != protein_count + 1:
    print(
        f"ERROR: ligand tleap residue is {tleap_resid}; "
        f"expected {protein_count + 1}.",
        file=sys.stderr,
    )
    sys.exit(23)

print(
    f"Topology QC PASS: ligand {ligand_name} is residue {tleap_resid} "
    f"after {protein_count} protein residues."
)
PY

info "Topology QC passed."

# ---------------------- ONE-FRAME TRAJECTORY -------------------------------
info "Creating one-frame NetCDF trajectory..."

cat > "${WORK_DIR}/cpptraj_oneframe.in" <<EOF
parm $COMPLEX_TOP
trajin $COMPLEX_INPCRD 1 1
trajout $ONEFRAME_TRAJ netcdf
run
quit
EOF

cpptraj -i "${WORK_DIR}/cpptraj_oneframe.in" \
    > "${LOG_DIR}/cpptraj_oneframe.log" 2>&1

[[ -s "$ONEFRAME_TRAJ" ]] || {
    tail -n 50 "${LOG_DIR}/cpptraj_oneframe.log" >&2 || true
    die "cpptraj could not create the one-frame trajectory."
}

# ------------------------------ MM/GBSA ------------------------------------
# Run the ordinary MM/GBSA calculation separately from decomposition.
# This mirrors the working 1TCO workflow and preserves the overall energy
# summary even if Amber's decomposition writer later sees ******** overflow.
cat > "$MMPBSA_INPUT" <<EOF
${PREFIX} one-frame MM/GBSA
&general
  startframe=1,
  endframe=1,
  interval=1,
  verbose=1,
  keep_files=0,
  receptor_mask="$RECEPTOR_MASK",
  ligand_mask="$LIGAND_MASK",
/
&gb
  igb=$IGB,
  saltcon=$SALTCON,
  surften=$SURFTEN,
  surfoff=$SURFOFF,
/
EOF

MMGBSA_WORK="${WORK_DIR}/mmgbsa"
mkdir -p "$MMGBSA_WORK"

info "Running standard MM/GBSA..."
(
    cd "$MMGBSA_WORK"
    MMPBSA.py -O \
        -i "$MMPBSA_INPUT" \
        -o "$FINAL_RESULTS" \
        -cp "$COMPLEX_TOP" \
        -rp "$RECEPTOR_TOP" \
        -lp "$LIGAND_TOP" \
        -y "$ONEFRAME_TRAJ" \
        > "${LOG_DIR}/MMPBSA.log" 2>&1
) || {
    tail -n 80 "${LOG_DIR}/MMPBSA.log" >&2 || true
    die "Standard MM/GBSA calculation failed."
}

[[ -s "$FINAL_RESULTS" ]] || die "MM/GBSA result file was not created."

# ----------------------- ALL-RESIDUE DECOMPOSITION --------------------------
cat > "$DECOMP_INPUT" <<EOF
${PREFIX} one-frame all-residue decomposition
&general
  startframe=1,
  endframe=1,
  interval=1,
  verbose=1,
  keep_files=2,
  receptor_mask="$RECEPTOR_MASK",
  ligand_mask="$LIGAND_MASK",
/
&gb
  igb=$IGB,
  saltcon=$SALTCON,
  surften=$SURFTEN,
  surfoff=$SURFOFF,
/
&decomp
  idecomp=1,
  dec_verbose=1,
/
EOF

DECOMP_WORK="${WORK_DIR}/decomposition"
mkdir -p "$DECOMP_WORK"

info "Running all-residue per-residue decomposition..."
set +e
(
    cd "$DECOMP_WORK"
    MMPBSA.py -O \
        -i "$DECOMP_INPUT" \
        -o "${RUN_DIR}/DECOMP_MMPBSA_SUMMARY.tmp" \
        -do "$FINAL_DECOMP" \
        -cp "$COMPLEX_TOP" \
        -rp "$RECEPTOR_TOP" \
        -lp "$LIGAND_TOP" \
        -y "$ONEFRAME_TRAJ" \
        > "${LOG_DIR}/MMPBSA_decomp.log" 2>&1
)
DECOMP_STATUS=$?
set -e

RAW_COMPLEX="${DECOMP_WORK}/_MMPBSA_complex_gb.mdout.0"
RAW_RECEPTOR="${DECOMP_WORK}/_MMPBSA_receptor_gb.mdout.0"
RAW_LIGAND="${DECOMP_WORK}/_MMPBSA_ligand_gb.mdout.0"

[[ -s "$RAW_COMPLEX" ]] || die "Raw complex decomposition output was not created."
[[ -s "$RAW_RECEPTOR" ]] || die "Raw receptor decomposition output was not created."

grep -q "PRINT DECOMP - TOTAL ENERGIES" "$RAW_COMPLEX" || \
    die "Raw complex decomposition table was not found."
grep -q "PRINT DECOMP - TOTAL ENERGIES" "$RAW_RECEPTOR" || \
    die "Raw receptor decomposition table was not found."

if [[ $DECOMP_STATUS -ne 0 ]]; then
    if grep -q '\*\*\*\*' "$RAW_COMPLEX" "$RAW_RECEPTOR" 2>/dev/null; then
        warn "Amber decomposition calculations finished, but its final parser encountered ******** overflow."
        warn "Using the raw TDC fallback parser; no expensive calculation will be repeated."
    else
        warn "MMPBSA.py decomposition writer returned exit code $DECOMP_STATUS."
        warn "Raw TDC outputs exist, so the fallback parser will be used."
    fi
fi

# ------------------------- RAW TDC PARSER + MAPPING -------------------------
info "Parsing raw TDC decomposition and mapping back to original PDB residues..."

python3 - \
    "$FINAL_RESULTS" \
    "$MAPPING_TSV" \
    "$SUMMARY_TSV" \
    "$DECOMP_TSV" \
    "$RAW_FALLBACK_TSV" \
    "$INTEREST_TSV" \
    "$REPORT_TXT" \
    "$DISTANCE_RESIDUES" \
    "$HYDROGEN_BOND_RESIDUES" \
    "$PROTEIN_RES_COUNT" \
    "$LIGAND_RESIDUE" \
    "${REFERENCE_DELTA_TOTAL:-}" \
    "$REFERENCE_TOLERANCE" \
    "$REFERENCE_CHECK" \
    "$SURFTEN" \
    "$SURFOFF" \
    "$RAW_COMPLEX" \
    "$RAW_RECEPTOR" \
    "$RAW_LIGAND" \
    "$DECOMP_STATUS" <<'RAWPARSEPY'
import csv
import re
import sys
from collections import OrderedDict
from pathlib import Path

(
    final_results_s,
    mapping_s,
    summary_s,
    decomp_out_s,
    raw_fallback_s,
    interest_out_s,
    report_s,
    distance_s,
    hbond_s,
    protein_count_s,
    ligand_res_s,
    reference_total_s,
    tolerance_s,
    reference_check_s,
    surften_s,
    surfoff_s,
    raw_complex_s,
    raw_receptor_s,
    raw_ligand_s,
    decomp_status_s,
) = sys.argv[1:]

final_results = Path(final_results_s)
mapping_path = Path(mapping_s)
summary_path = Path(summary_s)
decomp_out = Path(decomp_out_s)
raw_fallback_out = Path(raw_fallback_s)
interest_out = Path(interest_out_s)
report_path = Path(report_s)
raw_complex = Path(raw_complex_s)
raw_receptor = Path(raw_receptor_s)
raw_ligand = Path(raw_ligand_s)

protein_count = int(protein_count_s)
ligand_res = int(ligand_res_s)
reference_total = float(reference_total_s) if reference_total_s else None
tolerance = float(tolerance_s)
reference_check = bool(int(reference_check_s))
surften = float(surften_s)
surfoff = float(surfoff_s)
decomp_status = int(decomp_status_s)

interest = OrderedDict()

def add_specs(spec_string, category):
    for token in [x.strip() for x in spec_string.split(",") if x.strip()]:
        chain, residue, expected_name = token.split(":")
        key = (chain.upper(), residue)
        if key not in interest:
            interest[key] = {"Expected_Name": expected_name.upper(), "Categories": []}
        elif interest[key]["Expected_Name"] != expected_name.upper():
            raise SystemExit(f"ERROR: conflicting expected residue names for {chain}:{residue}")
        if category not in interest[key]["Categories"]:
            interest[key]["Categories"].append(category)

add_specs(distance_s, "Distance")
add_specs(hbond_s, "Hydrogen Bond")

wanted = ["VDWAALS", "EEL", "EGB", "ESURF", "DELTA G gas", "DELTA G solv", "DELTA TOTAL"]
summary = {}
if final_results.exists():
    results_text = final_results.read_text(errors="replace")

    # Parse only the binding-energy DIFFERENCES block.
    # Searching the whole file can accidentally capture absolute complex
    # energies (very large values) instead of Complex-Receptor-Ligand.
    marker = "Differences (Complex - Receptor - Ligand):"
    starts = [m.start() for m in re.finditer(re.escape(marker), results_text)]
    candidate_blocks = []
    for i, start in enumerate(starts):
        end = starts[i + 1] if i + 1 < len(starts) else len(results_text)
        block = results_text[start:end]
        parsed = {}
        for label in wanted:
            m = re.search(
                rf"^\s*{re.escape(label)}\s+"
                rf"([-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?)",
                block,
                flags=re.MULTILINE,
            )
            if m:
                parsed[label] = float(m.group(1))
        if "DELTA TOTAL" in parsed:
            candidate_blocks.append(parsed)

    if not candidate_blocks:
        raise SystemExit(
            "ERROR: Could not locate a complete "
            "'Differences (Complex - Receptor - Ligand)' block "
            "in FINAL_RESULTS."
        )

    # Use the last complete differences block; MMPBSA.py may repeat it.
    summary = candidate_blocks[-1]

    # Internal arithmetic QC for the parsed binding-energy components.
    required = {"VDWAALS", "EEL", "EGB", "ESURF",
                "DELTA G gas", "DELTA G solv", "DELTA TOTAL"}
    if required.issubset(summary):
        gas_calc = summary["VDWAALS"] + summary["EEL"]
        solv_calc = summary["EGB"] + summary["ESURF"]
        total_calc = summary["DELTA G gas"] + summary["DELTA G solv"]
        tol = 0.05
        if abs(gas_calc - summary["DELTA G gas"]) > tol:
            raise SystemExit(
                "ERROR: MM/GBSA summary QC failed: "
                "VDWAALS + EEL does not match DELTA G gas."
            )
        if abs(solv_calc - summary["DELTA G solv"]) > tol:
            raise SystemExit(
                "ERROR: MM/GBSA summary QC failed: "
                "EGB + ESURF does not match DELTA G solv."
            )
        if abs(total_calc - summary["DELTA TOTAL"]) > tol:
            raise SystemExit(
                "ERROR: MM/GBSA summary QC failed: "
                "DELTA G gas + DELTA G solv does not match DELTA TOTAL."
            )

summary_path.parent.mkdir(parents=True, exist_ok=True)
with summary_path.open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["Energy_Component", "Average_kcal_mol"])
    for label in wanted:
        if label in summary:
            w.writerow([label, f"{summary[label]:.6f}"])

mapping = {}
with mapping_path.open() as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        mapping[int(row["Topology_Residue"])] = row

def parse_value(token):
    token = token.strip()
    if "*" in token:
        return None
    try:
        return float(token)
    except ValueError:
        return None


def parse_tdc(path):
    data = {}
    if not path.exists():
        return data
    for line in path.read_text(errors="replace").splitlines():
        parts = line.split()
        if len(parts) < 7 or parts[0] != "TDC":
            continue
        try:
            resid = int(parts[1])
        except ValueError:
            continue
        if resid in data:
            continue
        data[resid] = {
            "Internal": parse_value(parts[2]),
            "VDW": parse_value(parts[3]),
            "Electrostatic": parse_value(parts[4]),
            "Polar_Solvation": parse_value(parts[5]),
            "SASA": parse_value(parts[6]),
        }
    return data

complex_tdc = parse_tdc(raw_complex)
receptor_tdc = parse_tdc(raw_receptor)
ligand_tdc = parse_tdc(raw_ligand)


def subtract_components(a, b):
    keys = ["Internal", "VDW", "Electrostatic", "Polar_Solvation", "SASA"]
    if a is None or b is None:
        return None
    if any(a.get(k) is None or b.get(k) is None for k in keys):
        return None
    d = {k: a[k] - b[k] for k in keys}
    d["NonPolar_Solvation"] = surften * d["SASA"]
    d["TOTAL"] = d["Internal"] + d["VDW"] + d["Electrostatic"] + d["Polar_Solvation"] + d["NonPolar_Solvation"]
    return d

rows = []
overflow_residues = []
for resid in range(1, protein_count + 1):
    pdb = mapping.get(resid)
    a = complex_tdc.get(resid)
    b = receptor_tdc.get(resid)
    diff = subtract_components(a, b)
    if a is None or b is None:
        status = "MISSING_RAW_TDC"
    elif diff is None:
        status = "OVERFLOW"
        overflow_residues.append(resid)
    else:
        status = "OK"

    item = {
        "Topology_Residue": resid,
        "Topology_Name": pdb["PDB_Name"].upper() if pdb else "?",
        "PDB_Chain": pdb["PDB_Chain"].upper() if pdb else "?",
        "PDB_Residue": pdb["PDB_Residue"] if pdb else "?",
        "PDB_ICode": pdb["PDB_ICode"] if pdb else "",
        "PDB_Name": pdb["PDB_Name"].upper() if pdb else "?",
        "Mapping_Status": "OK" if pdb else "NO_MAPPING",
        "Raw_Status": status,
        "Residue_of_Interest": "NO",
        "Interest_Category": "",
        "Expected_Name": "",
        "Interest_Name_QC": "",
        "Internal": None,
        "VDW": None,
        "Electrostatic": None,
        "Polar_Solvation": None,
        "SASA_Delta": None,
        "NonPolar_Solvation": None,
        "TOTAL": None,
    }

    key = (item["PDB_Chain"], str(item["PDB_Residue"]))
    if key in interest:
        spec = interest[key]
        item["Residue_of_Interest"] = "YES"
        item["Interest_Category"] = ";".join(spec["Categories"])
        item["Expected_Name"] = spec["Expected_Name"]
        item["Interest_Name_QC"] = "PASS" if item["PDB_Name"] == spec["Expected_Name"] else "MISMATCH"

    if diff is not None:
        for k in ("Internal", "VDW", "Electrostatic", "Polar_Solvation", "NonPolar_Solvation", "TOTAL"):
            item[k] = diff[k]
        item["SASA_Delta"] = diff["SASA"]
    rows.append(item)

ligand_complex = complex_tdc.get(ligand_res)
ligand_isolated = None
if ligand_tdc:
    ligand_isolated = ligand_tdc.get(1)
    if ligand_isolated is None:
        ligand_isolated = ligand_tdc[sorted(ligand_tdc)[0]]
ligand_diff = subtract_components(ligand_complex, ligand_isolated)
ligand_item = {
    "Topology_Residue": ligand_res, "Topology_Name": "FK5", "PDB_Chain": "LIGAND",
    "PDB_Residue": str(ligand_res), "PDB_ICode": "", "PDB_Name": "FK5",
    "Mapping_Status": "LIGAND", "Raw_Status": "OK" if ligand_diff is not None else "UNAVAILABLE",
    "Residue_of_Interest": "NO", "Interest_Category": "", "Expected_Name": "", "Interest_Name_QC": "",
    "Internal": None, "VDW": None, "Electrostatic": None, "Polar_Solvation": None,
    "SASA_Delta": None, "NonPolar_Solvation": None, "TOTAL": None,
}
if ligand_diff is not None:
    for k in ("Internal", "VDW", "Electrostatic", "Polar_Solvation", "NonPolar_Solvation", "TOTAL"):
        ligand_item[k] = ligand_diff[k]
    ligand_item["SASA_Delta"] = ligand_diff["SASA"]
rows.append(ligand_item)

interest_rows = [r for r in rows if r["Residue_of_Interest"] == "YES"]
interest_lookup = {(r["PDB_Chain"], str(r["PDB_Residue"])): r for r in interest_rows}
missing_interest = [key for key in interest if key not in interest_lookup]
interest_name_mismatches = [r for r in interest_rows if r["Interest_Name_QC"] == "MISMATCH"]
interest_raw_failures = [r for r in interest_rows if r["Raw_Status"] != "OK"]
interest_qc_status = "PASS" if (not missing_interest and not interest_name_mismatches and not interest_raw_failures and len(interest_rows) == len(interest)) else "WARNING"
mapping_mismatches = [r for r in rows[:protein_count] if r["Mapping_Status"] != "OK"]
mapping_qc_status = "PASS" if not mapping_mismatches else "WARNING"

reference_status = "NOT_CHECKED"
delta_total = summary.get("DELTA TOTAL")
if reference_check and reference_total is not None and delta_total is not None:
    dref = abs(delta_total - reference_total)
    reference_status = "PASS" if dref <= tolerance else f"WARNING (difference={dref:.4f} kcal/mol)"

fields = ["Topology_Residue", "Topology_Name", "PDB_Chain", "PDB_Residue", "PDB_ICode", "PDB_Name", "Mapping_Status", "Raw_Status", "Residue_of_Interest", "Interest_Category", "Expected_Name", "Interest_Name_QC", "Internal", "VDW", "Electrostatic", "Polar_Solvation", "SASA_Delta", "NonPolar_Solvation", "TOTAL"]

def printable(v):
    if v is None:
        return "NA"
    if isinstance(v, float):
        return f"{v:.6f}"
    return v

for path in (decomp_out, raw_fallback_out):
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for row in rows:
            w.writerow({k: printable(row.get(k)) for k in fields})

interest_fields = ["PDB_Chain", "PDB_Residue", "Expected_Name", "PDB_Name", "Topology_Residue", "Interest_Category", "Interest_Name_QC", "Raw_Status", "VDW", "Electrostatic", "Polar_Solvation", "SASA_Delta", "NonPolar_Solvation", "TOTAL"]
with interest_out.open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=interest_fields, delimiter="\t")
    w.writeheader()
    for key, spec in interest.items():
        r = interest_lookup.get(key)
        if r is None:
            w.writerow({"PDB_Chain": key[0], "PDB_Residue": key[1], "Expected_Name": spec["Expected_Name"], "PDB_Name": "NOT_FOUND", "Topology_Residue": "", "Interest_Category": ";".join(spec["Categories"]), "Interest_Name_QC": "NOT_FOUND", "Raw_Status": "NOT_FOUND"})
            continue
        w.writerow({
            "PDB_Chain": r["PDB_Chain"], "PDB_Residue": r["PDB_Residue"], "Expected_Name": spec["Expected_Name"],
            "PDB_Name": r["PDB_Name"], "Topology_Residue": r["Topology_Residue"], "Interest_Category": ";".join(spec["Categories"]),
            "Interest_Name_QC": r["Interest_Name_QC"], "Raw_Status": r["Raw_Status"], "VDW": printable(r["VDW"]),
            "Electrostatic": printable(r["Electrostatic"]), "Polar_Solvation": printable(r["Polar_Solvation"]),
            "SASA_Delta": printable(r["SASA_Delta"]), "NonPolar_Solvation": printable(r["NonPolar_Solvation"]), "TOTAL": printable(r["TOTAL"]),
        })

with report_path.open("w") as fh:
    fh.write("FINAL ANALYSIS REPORT\n")
    fh.write("=" * 78 + "\n\n")
    fh.write(f"Protein topology residues : 1-{protein_count}\n")
    fh.write(f"Ligand topology residue   : {ligand_res}\n")
    fh.write("Decomposition scope       : ALL protein residues + ligand\n")
    fh.write("Decomposition parser      : RAW TDC fallback parser\n")
    fh.write(f"Amber decomp exit code    : {decomp_status}\n")
    fh.write(f"Raw overflow residues     : {len(overflow_residues)}\n")
    if overflow_residues:
        fh.write("Overflow topology residues: " + ",".join(map(str, overflow_residues)) + "\n")
    fh.write(f"Residues of interest      : {len(interest)} unique residues\n")
    fh.write(f"Residues-of-interest QC   : {interest_qc_status}\n")
    fh.write(f"Mapping QC                : {mapping_qc_status}\n\n")
    fh.write("MM/GBSA SUMMARY\n" + "-" * 78 + "\n")
    for label in wanted:
        if label in summary:
            fh.write(f"{label:18s} {summary[label]:12.6f} kcal/mol\n")
    fh.write("\n")
    if reference_total is not None and reference_check:
        fh.write(f"Regression check against DELTA TOTAL {reference_total:.4f} +/- {tolerance:.2f}: {reference_status}\n\n")
    fh.write("RESIDUES OF INTEREST\n" + "-" * 78 + "\n")
    for key, spec in interest.items():
        r = interest_lookup.get(key)
        category = ";".join(spec["Categories"])
        if r is None:
            fh.write(f"{key[0]}:{key[1]} {spec['Expected_Name']} [{category}] NOT FOUND\n")
            continue
        if r["Raw_Status"] != "OK":
            fh.write(f"{r['PDB_Chain']}:{r['PDB_Residue']} {r['PDB_Name']} / topology {r['Topology_Residue']} [{category}] RAW STATUS={r['Raw_Status']}\n")
            continue
        fh.write(f"{r['PDB_Chain']}:{r['PDB_Residue']} {r['PDB_Name']} / topology {r['Topology_Residue']} [{category}] TOTAL={r['TOTAL']:10.6f}  vdW={r['VDW']:10.6f}  EEL={r['Electrostatic']:10.6f}  Polar={r['Polar_Solvation']:10.6f}  NonPolar={r['NonPolar_Solvation']:10.6f}\n")
    fh.write("\nRAW FALLBACK NOTE\n" + "-" * 78 + "\n")
    fh.write("For protein residues, binding decomposition is calculated as Complex TDC - Receptor TDC.\n")
    fh.write(f"Non-polar contribution uses surften={surften} and surfoff={surfoff}; NonPolar = surften * delta(SASA).\n")
    fh.write("Rows containing Amber fixed-width overflow (********) are retained as OVERFLOW/NA instead of aborting the entire analysis.\n")

print(f"Parsed protein decomposition rows: {protein_count}")
print(f"Overflow rows retained as NA: {len(overflow_residues)}")
print("Mapping QC:", mapping_qc_status)
print("Residues-of-interest QC:", interest_qc_status)
print("Regression check:", reference_status)

if interest_qc_status != "PASS":
    raise SystemExit("ERROR: one or more residues of interest could not be recovered safely.")
if mapping_qc_status != "PASS":
    raise SystemExit("ERROR: PDB-topology mapping QC failed.")
if reference_check and reference_total is not None and reference_status != "PASS":
    raise SystemExit("ERROR: MM/GBSA regression check failed.")
RAWPARSEPY

# ----------------------------- METADATA ------------------------------------
{
    echo "Pipeline version: $VERSION"
    echo "Run date: $(date -Is)"
    echo "Profile: $PROFILE"
    echo "Input PDB: $INPUT_PDB"
    echo "Ligand residue name: $LIGAND"
    echo "Ligand charge: $LIGAND_CHARGE"
    echo "Protein leaprc: $PROTEIN_LEAPRC"
    echo "Ligand force field: $GAFF_TYPE"
    echo "IGB: $IGB"
    echo "Salt concentration: $SALTCON"
    echo "Surface tension: $SURFTEN"
    echo "Surface offset: $SURFOFF"
    echo "Decomposition scope: ALL residues"
    echo "Residues of interest: $PDB_TARGETS"
    echo "Distance residues: $DISTANCE_RESIDUES"
    echo "Hydrogen-bond residues: $HYDROGEN_BOND_RESIDUES"
    echo "Protein residue count: $PROTEIN_RES_COUNT"
    echo "Ligand topology residue: $LIGAND_RESIDUE"
    echo "Python: $(python3 --version 2>&1)"
    echo "MMPBSA.py: $(command -v MMPBSA.py)"
    echo "antechamber: $(command -v antechamber)"
    echo "tleap: $(command -v tleap)"
} > "$METADATA_TXT"

echo
echo "=================================================================="
echo "DONE"
echo "=================================================================="
echo "Run directory:"
echo "  $RUN_DIR"
echo
echo "Main outputs:"
echo "  $FINAL_RESULTS"
echo "  $FINAL_DECOMP"
echo "  $REPORT_TXT"
echo "  $SUMMARY_TSV"
echo "  $DECOMP_TSV"
echo "  $INTEREST_TSV"
echo "  $RAW_FALLBACK_TSV"
echo "  $TARGET_MAP_TSV"
echo "  $MAPPING_TSV"
echo "  $REMOVED_HETERO_TSV"
echo "  $METADATA_TXT"
echo
echo "Preparation outputs:"
echo "  $COMPLEX_TOP"
echo "  $COMPLEX_INPCRD"
echo "  $RECEPTOR_TOP"
echo "  $LIGAND_TOP"
echo "  $TLEAP_PDB"
echo
echo "For 6TZ6 derivative analyses, inspect RUN_REPORT.txt and residues_of_interest.tsv."
echo "The FK506 reference regression check is disabled."
echo "Verify that Mapping QC and Residues-of-interest QC are PASS."
echo "=================================================================="
