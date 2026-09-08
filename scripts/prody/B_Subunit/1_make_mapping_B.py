#!/usr/bin/env python3
"""
Create residue mapping for Calcineurin B subunit (1TCO B vs 6TZ6 B).

Default inputs:
    1TCO_B_Subunit.pdb
    6TZ6_B_Subunit.pdb

Outputs:
    mapping_B_new.csv
    B_sequence_alignment.txt

The mapping CSV format is:
    Bos residue,Residue (B),Candida residue,Residue (B)

Important:
    ProDy selection "chain B and protein and name CA" is used to ensure
    Ca2+ ions or other HETATM records are NOT mistaken for protein C-alpha atoms.
"""

import argparse
import csv
import os
import sys
from collections import Counter

from Bio import Align
from Bio.Align import substitution_matrices
from Bio.SeqUtils import seq1
from prody import parsePDB


CUSTOM_AA = {
    "MSE": "M",
    "HIE": "H",
    "HID": "H",
    "HIP": "H",
    "CYX": "C",
    "CYM": "C",
    "ASH": "D",
    "GLH": "E",
    "LYN": "K",
}


def pretty_three_letter(resname3):
    resname3 = resname3.strip().upper()
    if resname3 in {"HIE", "HID", "HIP"}:
        return "His"
    if resname3 in {"CYX", "CYM"}:
        return "Cys"
    if resname3 == "ASH":
        return "Asp"
    if resname3 == "GLH":
        return "Glu"
    if resname3 == "LYN":
        return "Lys"
    if resname3 == "MSE":
        return "Met"
    return resname3.capitalize()


def read_ca_residues(pdb_file, chain="B"):
    """Return ordered protein C-alpha residues from one specified chain."""
    structure = parsePDB(pdb_file)
    if structure is None:
        raise RuntimeError(f"ProDy could not parse: {pdb_file}")

    ca = structure.select(f"chain {chain} and protein and name CA")
    if ca is None or ca.numAtoms() == 0:
        raise RuntimeError(
            f"No protein C-alpha atoms found for chain {chain} in {pdb_file}"
        )

    resnums = ca.getResnums()
    resnames = ca.getResnames()
    chids = ca.getChids()
    icodes = ca.getIcodes()

    residues = []
    seen_numeric_ids = set()

    for chid, resnum, resname, icode in zip(chids, resnums, resnames, icodes):
        numeric_key = int(resnum)

        if numeric_key in seen_numeric_ids:
            raise ValueError(
                f"Duplicate numeric residue number {numeric_key} in {pdb_file}. "
                "Insertion codes/duplicate numbering must be resolved before comparison."
            )
        seen_numeric_ids.add(numeric_key)

        resname = str(resname).strip().upper()
        one = seq1(resname, custom_map=CUSTOM_AA, undef_code="X")

        residues.append(
            {
                "chain": str(chid).strip(),
                "resnum": numeric_key,
                "icode": str(icode).strip(),
                "resname3": resname,
                "resname1": one,
            }
        )

    return residues


def build_mapping(bos_residues, candida_residues):
    bos_seq = "".join(r["resname1"] for r in bos_residues)
    candida_seq = "".join(r["resname1"] for r in candida_residues)

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(bos_seq, candida_seq)
    if len(alignments) == 0:
        raise RuntimeError("No sequence alignment could be generated.")

    aln = alignments[0]
    mapping = []
    matches = 0

    for bos_block, candida_block in zip(aln.aligned[0], aln.aligned[1]):
        b0, b1 = map(int, bos_block)
        c0, c1 = map(int, candida_block)

        if (b1 - b0) != (c1 - c0):
            raise RuntimeError("Unexpected unequal aligned block lengths.")

        for bi, ci in zip(range(b0, b1), range(c0, c1)):
            b = bos_residues[bi]
            c = candida_residues[ci]

            if b["resname1"] == c["resname1"]:
                matches += 1

            mapping.append(
                {
                    "Bos residue": b["resnum"],
                    "Residue (B)": pretty_three_letter(b["resname3"]),
                    "Candida residue": c["resnum"],
                    "Residue (B)_Candida": pretty_three_letter(c["resname3"]),
                }
            )

    identity = 100.0 * matches / len(mapping) if mapping else 0.0
    return aln, mapping, identity


def write_mapping_csv(mapping, outfile):
    with open(outfile, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow(
            ["Bos residue", "Residue (B)", "Candida residue", "Residue (B)"]
        )
        for row in mapping:
            writer.writerow(
                [
                    row["Bos residue"],
                    row["Residue (B)"],
                    row["Candida residue"],
                    row["Residue (B)_Candida"],
                ]
            )


def compare_with_old(new_mapping, old_file):
    if not old_file or not os.path.exists(old_file):
        return

    with open(old_file, "r", encoding="utf-8-sig", newline="") as f:
        rows = list(csv.reader(f))

    old_rows = []
    for row in rows[1:]:
        if not row or len(row) < 3:
            continue
        old_rows.append((int(float(row[0])), int(float(row[2]))))

    new_rows = [
        (int(r["Bos residue"]), int(r["Candida residue"])) for r in new_mapping
    ]

    print("\nQC against existing mapping:")
    if old_rows == new_rows:
        print(f"  PASS: new residue-number mapping is IDENTICAL to {old_file}")
    else:
        print(f"  WARNING: new mapping differs from {old_file}")
        old_set = set(old_rows)
        new_set = set(new_rows)
        only_old = sorted(old_set - new_set)
        only_new = sorted(new_set - old_set)
        print(f"  Pairs only in old mapping: {len(only_old)}")
        print(f"  Pairs only in new mapping: {len(only_new)}")
        if only_old[:10]:
            print("  First old-only pairs:", only_old[:10])
        if only_new[:10]:
            print("  First new-only pairs:", only_new[:10])


def main():
    ap = argparse.ArgumentParser(
        description="Generate B-subunit residue mapping for 1TCO vs 6TZ6."
    )
    ap.add_argument(
        "--bos",
        default="1TCO_B_Subunit.pdb",
        help="1TCO B-subunit PDB (default: 1TCO_B_Subunit.pdb)",
    )
    ap.add_argument(
        "--candida",
        default="6TZ6_B_Subunit.pdb",
        help="6TZ6 B-subunit PDB (default: 6TZ6_B_Subunit.pdb)",
    )
    ap.add_argument(
        "--chain",
        default="B",
        help="Protein chain ID to map (default: B)",
    )
    ap.add_argument(
        "--out",
        default="mapping_B_new.csv",
        help="Output mapping CSV (default: mapping_B_new.csv)",
    )
    ap.add_argument(
        "--alignment-out",
        default="B_sequence_alignment.txt",
        help="Text alignment output (default: B_sequence_alignment.txt)",
    )
    ap.add_argument(
        "--old-mapping",
        default="mapping_B.csv",
        help="Existing mapping used only for QC (default: mapping_B.csv)",
    )
    args = ap.parse_args()

    for f in (args.bos, args.candida):
        if not os.path.isfile(f):
            print(f"ERROR: file not found: {f}", file=sys.stderr)
            sys.exit(1)

    print("Reading PDB files...")
    bos = read_ca_residues(args.bos, args.chain)
    candida = read_ca_residues(args.candida, args.chain)

    print(f"1TCO B C-alpha residues : {len(bos)}")
    print(f"6TZ6 B C-alpha residues : {len(candida)}")
    print(f"1TCO numbering          : {bos[0]['resnum']} .. {bos[-1]['resnum']}")
    print(
        f"6TZ6 numbering          : "
        f"{candida[0]['resnum']} .. {candida[-1]['resnum']}"
    )

    aln, mapping, identity = build_mapping(bos, candida)
    write_mapping_csv(mapping, args.out)

    with open(args.alignment_out, "w", encoding="utf-8") as f:
        f.write("1TCO B vs 6TZ6 B global sequence alignment\n")
        f.write("Scoring: BLOSUM62, gap open -10, gap extend -0.5\n\n")
        f.write(str(aln))
        f.write("\n")

    print("\nAlignment completed.")
    print(f"Mapped residue pairs    : {len(mapping)}")
    print(f"Sequence identity       : {identity:.2f}%")

    if mapping:
        first = mapping[0]
        last = mapping[-1]
        print(
            "First mapping pair     : "
            f"{first['Bos residue']} {first['Residue (B)']} -> "
            f"{first['Candida residue']} {first['Residue (B)_Candida']}"
        )
        print(
            "Last mapping pair      : "
            f"{last['Bos residue']} {last['Residue (B)']} -> "
            f"{last['Candida residue']} {last['Residue (B)_Candida']}"
        )

        offsets = Counter(
            int(r["Candida residue"]) - int(r["Bos residue"]) for r in mapping
        )
        print(
            "Most common numbering offsets (6TZ6 - 1TCO):",
            offsets.most_common(5),
        )

    print(f"\nSaved mapping           : {args.out}")
    print(f"Saved alignment         : {args.alignment_out}")

    compare_with_old(mapping, args.old_mapping)

    print("\nIMPORTANT:")
    print("  Review mapping_B_new.csv before using it downstream.")
    print("  This script performs residue sequence mapping only; it does NOT run ANM.")


if __name__ == "__main__":
    main()
