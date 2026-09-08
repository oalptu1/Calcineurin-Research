#!/usr/bin/env python3
"""
Create residue mapping for Calcineurin A subunit (1TCO A vs 6TZ6 A).

Default inputs:
    1TCO_A_Subunit.pdb
    6TZ6_A_Subunit.pdb

Outputs:
    mapping_A_new.csv
    A_sequence_alignment.txt

The mapping CSV is compatible with the existing compare_A_crosscorr.py workflow:
    Bos residue,Residue (A),Candida residue,Residue (A)
"""

import argparse
import csv
import os
import sys
from collections import Counter

from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1


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


def read_ca_residues(pdb_file):
    """Return ordered C-alpha-containing protein residues from first model."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(os.path.basename(pdb_file), pdb_file)
    model = next(structure.get_models())

    residues = []
    seen_numeric_ids = set()

    for chain in model:
        for residue in chain:
            if "CA" not in residue:
                continue

            hetflag, resnum, icode = residue.id
            resname = residue.resname.strip().upper()
            one = seq1(resname, custom_map=CUSTOM_AA, undef_code="X")

            # The downstream comparison script uses residue number only.
            # Stop if insertion codes would make that ambiguous.
            numeric_key = int(resnum)
            if numeric_key in seen_numeric_ids:
                raise ValueError(
                    f"Duplicate numeric residue number {numeric_key} in {pdb_file}. "
                    "Insertion codes/duplicate numbering must be resolved before using the old comparison workflow."
                )
            seen_numeric_ids.add(numeric_key)

            residues.append(
                {
                    "chain": chain.id,
                    "resnum": numeric_key,
                    "icode": icode.strip(),
                    "resname3": resname,
                    "resname1": one,
                }
            )

    if not residues:
        raise ValueError(f"No C-alpha residues found in {pdb_file}")

    return residues


def pretty_three_letter(resname3):
    """Convert VAL -> Val, HIE -> His for human-readable mapping."""
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

    # aligned contains ungapped aligned blocks in sequence-index coordinates
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
                    "Residue (A)": pretty_three_letter(b["resname3"]),
                    "Candida residue": c["resnum"],
                    "Residue (A)_Candida": pretty_three_letter(c["resname3"]),
                }
            )

    identity = 100.0 * matches / len(mapping) if mapping else 0.0
    return aln, mapping, identity, bos_seq, candida_seq


def write_mapping_csv(mapping, outfile):
    with open(outfile, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow(["Bos residue", "Residue (A)", "Candida residue", "Residue (A)"])
        for row in mapping:
            writer.writerow(
                [
                    row["Bos residue"],
                    row["Residue (A)"],
                    row["Candida residue"],
                    row["Residue (A)_Candida"],
                ]
            )


def compare_with_old(new_mapping, old_file):
    if not old_file or not os.path.exists(old_file):
        return

    old_rows = []
    with open(old_file, "r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            old_rows.append(
                (int(row["Bos residue"]), int(row["Candida residue"]))
            )

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
        description="Generate A-subunit residue mapping for 1TCO vs corrected 6TZ6."
    )
    ap.add_argument(
        "--bos",
        default="1TCO_A_Subunit.pdb",
        help="1TCO A-subunit PDB (default: 1TCO_A_Subunit.pdb)",
    )
    ap.add_argument(
        "--candida",
        default="6TZ6_A_Subunit.pdb",
        help="Corrected 6TZ6 A-subunit PDB (default: 6TZ6_A_Subunit.pdb)",
    )
    ap.add_argument(
        "--out",
        default="mapping_A_new.csv",
        help="Output mapping CSV (default: mapping_A_new.csv)",
    )
    ap.add_argument(
        "--alignment-out",
        default="A_sequence_alignment.txt",
        help="Text alignment output (default: A_sequence_alignment.txt)",
    )
    ap.add_argument(
        "--old-mapping",
        default="mapping_A.csv",
        help="Existing mapping used only for QC comparison (default: mapping_A.csv)",
    )
    args = ap.parse_args()

    for f in (args.bos, args.candida):
        if not os.path.isfile(f):
            print(f"ERROR: file not found: {f}", file=sys.stderr)
            sys.exit(1)

    print("Reading PDB files...")
    bos = read_ca_residues(args.bos)
    candida = read_ca_residues(args.candida)

    print(f"1TCO A C-alpha residues : {len(bos)}")
    print(f"6TZ6 A C-alpha residues : {len(candida)}")
    print(
        f"1TCO numbering          : {bos[0]['resnum']} .. {bos[-1]['resnum']}"
    )
    print(
        f"6TZ6 numbering          : {candida[0]['resnum']} .. {candida[-1]['resnum']}"
    )

    aln, mapping, identity, bos_seq, candida_seq = build_mapping(bos, candida)

    write_mapping_csv(mapping, args.out)

    with open(args.alignment_out, "w", encoding="utf-8") as f:
        f.write("1TCO A vs 6TZ6 A global sequence alignment\n")
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
            f"{first['Bos residue']} {first['Residue (A)']} -> "
            f"{first['Candida residue']} {first['Residue (A)_Candida']}"
        )
        print(
            "Last mapping pair      : "
            f"{last['Bos residue']} {last['Residue (A)']} -> "
            f"{last['Candida residue']} {last['Residue (A)_Candida']}"
        )

        offsets = Counter(
            int(r["Candida residue"]) - int(r["Bos residue"]) for r in mapping
        )
        common = offsets.most_common(5)
        print("Most common numbering offsets (6TZ6 - 1TCO):", common)

    print(f"\nSaved mapping           : {args.out}")
    print(f"Saved alignment         : {args.alignment_out}")

    compare_with_old(mapping, args.old_mapping)

    print("\nIMPORTANT:")
    print("  Review mapping_A_new.csv before replacing the old mapping_A.csv.")
    print("  This script performs residue sequence mapping only; it does NOT run ANM.")


if __name__ == "__main__":
    main()
