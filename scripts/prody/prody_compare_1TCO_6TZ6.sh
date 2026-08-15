#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# ProDy GNM Fluctuation Comparison: 1TCO vs 6TZ6
#
# Purpose
# -------
# Compare residue-level GNM fluctuations between 1TCO and 6TZ6 and generate:
#   1. Per-structure fluctuation tables
#   2. Sequence-based chain/residue mapping
#   3. A full matched-residue comparison table
#   4. The 20 residues with the largest fluctuation differences
#   5. The 20 residues with the smallest fluctuation differences
#
# Method
# ------
# - Protein C-alpha atoms are used to build one GNM for each complete complex.
# - All non-zero GNM modes are included in calcSqFlucts().
# - Because GNM square fluctuations are in arbitrary/relative units, each
#   structure's square-fluctuation profile is normalized to mean = 1.0 before
#   cross-structure comparison.
# - Homologous chains/residues are mapped using ProDy matchChains() with
#   sequence alignment enabled.
# - Ranking uses:
#       Signed_Difference = Normalized_6TZ6 - Normalized_1TCO
#       Absolute_Difference = abs(Signed_Difference)
#
# Default input files
# -------------------
#   1TCO_FK506.pdb
#   6TZ6_FK506.pdb
#
# Usage
# -----
#   bash prody_compare_1TCO_6TZ6.sh
#
# or:
#   bash prody_compare_1TCO_6TZ6.sh path/to/1TCO.pdb path/to/6TZ6.pdb
#
# Requirements
# ------------
#   Python 3
#   ProDy
#   NumPy
#
# The script does not modify either input PDB.
# ==============================================================================

PDB_1TCO="${1:-1TCO_FK506.pdb}"
PDB_6TZ6="${2:-6TZ6_FK506.pdb}"

TOP_N="${TOP_N:-20}"
SEQID="${SEQID:-50}"
OVERLAP="${OVERLAP:-50}"

RESULT_ROOT="${RESULT_ROOT:-ProDy_1TCO_6TZ6_RESULTS}"
RUN_ID="$(date +%Y%m%d_%H%M%S)"
RUN_DIR="${RESULT_ROOT}/run_${RUN_ID}"

mkdir -p "$RUN_DIR"

die() {
    echo "ERROR: $*" >&2
    exit 1
}

[[ -f "$PDB_1TCO" ]] || die "Input PDB not found: $PDB_1TCO"
[[ -f "$PDB_6TZ6" ]] || die "Input PDB not found: $PDB_6TZ6"

command -v python3 >/dev/null 2>&1 || die "python3 was not found."

python3 - <<'PY' >/dev/null 2>&1 || {
import prody
import numpy
PY
    die "ProDy and/or NumPy could not be imported in the active Python environment."
}

echo "=================================================================="
echo "ProDy GNM Fluctuation Comparison: 1TCO vs 6TZ6"
echo "=================================================================="
echo "1TCO input : $PDB_1TCO"
echo "6TZ6 input : $PDB_6TZ6"
echo "Top N      : $TOP_N"
echo "Seq. ID    : $SEQID%"
echo "Overlap    : $OVERLAP%"
echo "Run dir    : $RUN_DIR"
echo "=================================================================="

python3 - \
    "$PDB_1TCO" \
    "$PDB_6TZ6" \
    "$RUN_DIR" \
    "$TOP_N" \
    "$SEQID" \
    "$OVERLAP" <<'PY'
import csv
import math
import sys
from pathlib import Path

import numpy as np
from prody import GNM, calcSqFlucts, matchChains, parsePDB


pdb1_s, pdb2_s, outdir_s, topn_s, seqid_s, overlap_s = sys.argv[1:]

pdb1 = Path(pdb1_s)
pdb2 = Path(pdb2_s)
outdir = Path(outdir_s)
top_n = int(topn_s)
seqid = float(seqid_s)
overlap = float(overlap_s)

outdir.mkdir(parents=True, exist_ok=True)


def fail(message):
    raise SystemExit(f"ERROR: {message}")


def clean_text(value):
    if value is None:
        return ""
    return str(value).strip()


def atom_metadata(atoms):
    """Return residue metadata arrays for a C-alpha AtomGroup/Selection."""
    chids = [clean_text(x) or "_" for x in atoms.getChids()]
    resnums = [int(x) for x in atoms.getResnums()]
    resnames = [clean_text(x).upper() for x in atoms.getResnames()]

    try:
        raw_icodes = atoms.getIcodes()
    except Exception:
        raw_icodes = None

    if raw_icodes is None:
        icodes = [""] * len(resnums)
    else:
        icodes = [clean_text(x) for x in raw_icodes]

    return chids, resnums, icodes, resnames


def residue_key(chain, resnum, icode):
    return (clean_text(chain) or "_", int(resnum), clean_text(icode))


def build_gnm_profile(pdb_path, label):
    print(f"[{label}] Parsing PDB: {pdb_path}")
    atoms = parsePDB(str(pdb_path))
    if atoms is None:
        fail(f"ProDy could not parse {pdb_path}")

    ca = atoms.select("protein and calpha")
    if ca is None or ca.numAtoms() < 3:
        fail(f"No usable protein C-alpha selection was found in {pdb_path}")

    print(f"[{label}] Protein C-alpha atoms: {ca.numAtoms()}")

    gnm = GNM(f"{label}_GNM")
    gnm.buildKirchhoff(ca)
    gnm.calcModes(n_modes=None)

    sq = np.asarray(calcSqFlucts(gnm[:]), dtype=float)
    if len(sq) != ca.numAtoms():
        fail(
            f"{label}: fluctuation vector length ({len(sq)}) does not match "
            f"C-alpha count ({ca.numAtoms()})."
        )

    if not np.all(np.isfinite(sq)):
        fail(f"{label}: non-finite square fluctuations were produced.")

    mean_sq = float(np.mean(sq))
    if mean_sq <= 0:
        fail(f"{label}: mean square fluctuation is not positive.")

    norm_sq = sq / mean_sq
    rms = np.sqrt(np.clip(sq, 0.0, None))

    chids, resnums, icodes, resnames = atom_metadata(ca)

    rows = []
    lookup = {}

    for i, (chain, resnum, icode, resname) in enumerate(
        zip(chids, resnums, icodes, resnames)
    ):
        key = residue_key(chain, resnum, icode)
        row = {
            "Structure": label,
            "CA_Index": i + 1,
            "Chain": chain,
            "Residue": resnum,
            "ICode": icode,
            "Residue_Name": resname,
            "Sq_Fluctuation_Raw": float(sq[i]),
            "RMS_Fluctuation_Raw": float(rms[i]),
            "Sq_Fluctuation_Normalized": float(norm_sq[i]),
        }
        rows.append(row)

        if key in lookup:
            fail(
                f"{label}: duplicate C-alpha residue key detected: "
                f"{chain}:{resnum}{icode}"
            )
        lookup[key] = row

    return atoms, ca, gnm, rows, lookup, mean_sq


def write_csv(path, rows, fields):
    with Path(path).open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def chain_id_from_map(atom_map):
    chids = []
    for chid in atom_map.getChids():
        chid = clean_text(chid) or "_"
        if chid not in chids:
            chids.append(chid)
    if len(chids) != 1:
        fail(f"Expected one chain per AtomMap, found: {chids}")
    return chids[0]


def choose_unique_chain_matches(matches):
    """
    ProDy returns all matching chain pairs sorted by decreasing sequence
    identity. Keep the best one-to-one set so that a chain is not reused.
    """
    selected = []
    used_1 = set()
    used_2 = set()

    for amap1, amap2, identity, ovlp in matches:
        c1 = chain_id_from_map(amap1)
        c2 = chain_id_from_map(amap2)

        if c1 in used_1 or c2 in used_2:
            continue

        selected.append((amap1, amap2, float(identity), float(ovlp), c1, c2))
        used_1.add(c1)
        used_2.add(c2)

    return selected


atoms1, ca1, gnm1, fluct1, lookup1, mean1 = build_gnm_profile(pdb1, "1TCO")
atoms2, ca2, gnm2, fluct2, lookup2, mean2 = build_gnm_profile(pdb2, "6TZ6")

fluct_fields = [
    "Structure",
    "CA_Index",
    "Chain",
    "Residue",
    "ICode",
    "Residue_Name",
    "Sq_Fluctuation_Raw",
    "RMS_Fluctuation_Raw",
    "Sq_Fluctuation_Normalized",
]

write_csv(outdir / "fluctuations_1TCO.csv", fluct1, fluct_fields)
write_csv(outdir / "fluctuations_6TZ6.csv", fluct2, fluct_fields)

print("[Mapping] Matching protein chains by sequence similarity...")

try:
    matches = matchChains(
        atoms1.protein,
        atoms2.protein,
        subset="calpha",
        seqid=seqid,
        overlap=overlap,
        pwalign=True,
    )
except Exception as exc:
    fail(
        "ProDy matchChains() failed during sequence alignment. "
        f"Original error: {exc}"
    )

if not matches:
    fail(
        f"No matching protein chains were found at seqid={seqid}% "
        f"and overlap={overlap}%."
    )

selected_matches = choose_unique_chain_matches(matches)
if not selected_matches:
    fail("No unique one-to-one chain matches could be selected.")

chain_rows = []
comparison_rows = []
seen_pairs = set()

for rank, (amap1, amap2, identity, ovlp, chain1, chain2) in enumerate(
    selected_matches, 1
):
    n1 = amap1.numAtoms()
    n2 = amap2.numAtoms()
    if n1 != n2:
        fail(
            f"Matched AtomMaps have different lengths for "
            f"1TCO chain {chain1} and 6TZ6 chain {chain2}: {n1} vs {n2}"
        )

    chain_rows.append(
        {
            "Match_Rank": rank,
            "Chain_1TCO": chain1,
            "Chain_6TZ6": chain2,
            "Sequence_Identity_pct": f"{identity:.3f}",
            "Sequence_Overlap_pct": f"{ovlp:.3f}",
            "Matched_CA_Atoms": n1,
        }
    )

    c1, r1, i1, n1_names = atom_metadata(amap1)
    c2, r2, i2, n2_names = atom_metadata(amap2)

    for j in range(len(r1)):
        key1 = residue_key(c1[j], r1[j], i1[j])
        key2 = residue_key(c2[j], r2[j], i2[j])

        pair_key = (key1, key2)
        if pair_key in seen_pairs:
            continue
        seen_pairs.add(pair_key)

        row1 = lookup1.get(key1)
        row2 = lookup2.get(key2)

        if row1 is None:
            fail(f"Could not recover 1TCO fluctuation for {key1}")
        if row2 is None:
            fail(f"Could not recover 6TZ6 fluctuation for {key2}")

        signed = (
            row2["Sq_Fluctuation_Normalized"]
            - row1["Sq_Fluctuation_Normalized"]
        )
        absolute = abs(signed)

        if signed > 0:
            direction = "Higher_in_6TZ6"
        elif signed < 0:
            direction = "Higher_in_1TCO"
        else:
            direction = "Equal"

        comparison_rows.append(
            {
                "Chain_1TCO": key1[0],
                "Residue_1TCO": key1[1],
                "ICode_1TCO": key1[2],
                "Residue_Name_1TCO": row1["Residue_Name"],
                "Chain_6TZ6": key2[0],
                "Residue_6TZ6": key2[1],
                "ICode_6TZ6": key2[2],
                "Residue_Name_6TZ6": row2["Residue_Name"],
                "Sequence_Identity_pct": f"{identity:.3f}",
                "Sequence_Overlap_pct": f"{ovlp:.3f}",
                "Raw_SqFluct_1TCO": f"{row1['Sq_Fluctuation_Raw']:.10g}",
                "Raw_SqFluct_6TZ6": f"{row2['Sq_Fluctuation_Raw']:.10g}",
                "Normalized_SqFluct_1TCO": f"{row1['Sq_Fluctuation_Normalized']:.10g}",
                "Normalized_SqFluct_6TZ6": f"{row2['Sq_Fluctuation_Normalized']:.10g}",
                "Signed_Difference_6TZ6_minus_1TCO": f"{signed:.10g}",
                "Absolute_Difference": f"{absolute:.10g}",
                "Direction": direction,
            }
        )

if not comparison_rows:
    fail("No matched residue pairs were recovered.")

chain_fields = [
    "Match_Rank",
    "Chain_1TCO",
    "Chain_6TZ6",
    "Sequence_Identity_pct",
    "Sequence_Overlap_pct",
    "Matched_CA_Atoms",
]
write_csv(outdir / "chain_matches.csv", chain_rows, chain_fields)

comparison_fields = [
    "Chain_1TCO",
    "Residue_1TCO",
    "ICode_1TCO",
    "Residue_Name_1TCO",
    "Chain_6TZ6",
    "Residue_6TZ6",
    "ICode_6TZ6",
    "Residue_Name_6TZ6",
    "Sequence_Identity_pct",
    "Sequence_Overlap_pct",
    "Raw_SqFluct_1TCO",
    "Raw_SqFluct_6TZ6",
    "Normalized_SqFluct_1TCO",
    "Normalized_SqFluct_6TZ6",
    "Signed_Difference_6TZ6_minus_1TCO",
    "Absolute_Difference",
    "Direction",
]

write_csv(
    outdir / "matched_residue_fluctuation_comparison.csv",
    comparison_rows,
    comparison_fields,
)

def absdiff(row):
    return float(row["Absolute_Difference"])

most_different = sorted(
    comparison_rows,
    key=absdiff,
    reverse=True,
)[: min(top_n, len(comparison_rows))]

least_different = sorted(
    comparison_rows,
    key=absdiff,
)[: min(top_n, len(comparison_rows))]

write_csv(
    outdir / f"top{top_n}_most_different_fluctuations.csv",
    most_different,
    comparison_fields,
)
write_csv(
    outdir / f"top{top_n}_least_different_fluctuations.csv",
    least_different,
    comparison_fields,
)

# Human-readable compact lists.
def residue_label(row, which):
    return (
        f"{row[f'Residue_Name_{which}']} "
        f"{row[f'Chain_{which}']}:{row[f'Residue_{which}']}"
        f"{row[f'ICode_{which}']}"
    )

with (outdir / "TOP_FLUCTUATION_DIFFERENCES.txt").open("w") as fh:
    fh.write("ProDy GNM FLUCTUATION COMPARISON: 1TCO vs 6TZ6\n")
    fh.write("=" * 78 + "\n\n")
    fh.write(
        "Ranking is based on absolute difference between mean-normalized "
        "GNM square fluctuations.\n"
    )
    fh.write(
        "Signed difference = Normalized(6TZ6) - Normalized(1TCO).\n\n"
    )

    fh.write(f"TOP {len(most_different)} MOST DIFFERENT RESIDUES\n")
    fh.write("-" * 78 + "\n")
    for rank, row in enumerate(most_different, 1):
        fh.write(
            f"{rank:2d}. {residue_label(row, '1TCO'):18s} <-> "
            f"{residue_label(row, '6TZ6'):18s} "
            f"abs_diff={float(row['Absolute_Difference']):.6f} "
            f"signed={float(row['Signed_Difference_6TZ6_minus_1TCO']):+.6f} "
            f"{row['Direction']}\n"
        )

    fh.write("\n")
    fh.write(f"TOP {len(least_different)} LEAST DIFFERENT RESIDUES\n")
    fh.write("-" * 78 + "\n")
    for rank, row in enumerate(least_different, 1):
        fh.write(
            f"{rank:2d}. {residue_label(row, '1TCO'):18s} <-> "
            f"{residue_label(row, '6TZ6'):18s} "
            f"abs_diff={float(row['Absolute_Difference']):.6f} "
            f"signed={float(row['Signed_Difference_6TZ6_minus_1TCO']):+.6f} "
            f"{row['Direction']}\n"
        )

with (outdir / "RUN_REPORT.txt").open("w") as fh:
    fh.write("ProDy GNM FLUCTUATION COMPARISON REPORT\n")
    fh.write("=" * 78 + "\n\n")
    fh.write(f"1TCO input: {pdb1}\n")
    fh.write(f"6TZ6 input: {pdb2}\n")
    fh.write(f"1TCO protein C-alpha atoms: {ca1.numAtoms()}\n")
    fh.write(f"6TZ6 protein C-alpha atoms: {ca2.numAtoms()}\n")
    fh.write(f"1TCO raw mean square fluctuation: {mean1:.10g}\n")
    fh.write(f"6TZ6 raw mean square fluctuation: {mean2:.10g}\n")
    fh.write("Normalization: each structure divided by its own mean square fluctuation\n")
    fh.write("Normalized mean target: 1.0\n")
    fh.write(f"Sequence identity threshold: {seqid:.1f}%\n")
    fh.write(f"Sequence overlap threshold: {overlap:.1f}%\n")
    fh.write(f"Unique chain matches: {len(selected_matches)}\n")
    fh.write(f"Matched residue pairs: {len(comparison_rows)}\n")
    fh.write(f"Ranking size requested: {top_n}\n\n")

    fh.write("CHAIN MATCHES\n")
    fh.write("-" * 78 + "\n")
    for row in chain_rows:
        fh.write(
            f"1TCO {row['Chain_1TCO']} <-> 6TZ6 {row['Chain_6TZ6']} | "
            f"identity={row['Sequence_Identity_pct']}% | "
            f"overlap={row['Sequence_Overlap_pct']}% | "
            f"matched_CAs={row['Matched_CA_Atoms']}\n"
        )

    fh.write("\nOUTPUTS\n")
    fh.write("-" * 78 + "\n")
    for name in [
        "fluctuations_1TCO.csv",
        "fluctuations_6TZ6.csv",
        "chain_matches.csv",
        "matched_residue_fluctuation_comparison.csv",
        f"top{top_n}_most_different_fluctuations.csv",
        f"top{top_n}_least_different_fluctuations.csv",
        "TOP_FLUCTUATION_DIFFERENCES.txt",
        "RUN_REPORT.txt",
    ]:
        fh.write(name + "\n")

print()
print("==================================================================")
print("ANALYSIS COMPLETED")
print("==================================================================")
print(f"Unique chain matches : {len(selected_matches)}")
print(f"Matched residue pairs: {len(comparison_rows)}")
print(f"Results directory    : {outdir}")
print()
print("Main outputs:")
print(f"  {outdir / f'top{top_n}_most_different_fluctuations.csv'}")
print(f"  {outdir / f'top{top_n}_least_different_fluctuations.csv'}")
print(f"  {outdir / 'matched_residue_fluctuation_comparison.csv'}")
print(f"  {outdir / 'TOP_FLUCTUATION_DIFFERENCES.txt'}")
print(f"  {outdir / 'RUN_REPORT.txt'}")
print("==================================================================")
PY

echo
echo "DONE"
echo "Results:"
echo "  $RUN_DIR"
