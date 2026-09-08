#!/usr/bin/env python3
"""
Compare 1TCO vs 6TZ6 A-subunit ANM cross-correlation matrices.

Required inputs (default filenames):
    mapping_A_new.csv
    1TCO_A_crosscorr.csv
    6TZ6_A_crosscorr.csv
    1TCO_A_residues.csv
    6TZ6_A_residues.csv

Outputs:
    A_crosscorr_absolute_difference.csv
    A_crosscorr_residue_difference.csv
    A_top20_most_different.csv
    A_top20_least_different.csv

Definition:
    DeltaCC_i = mean_j | CC_6TZ6(i,j) - CC_1TCO(i,j) |

Only homologous residues listed in mapping_A_new.csv are compared.
Both matrices are first reduced/reordered to the same mapped-residue order.
"""

import argparse
import csv
import os
import sys

import numpy as np


def norm_header(s):
    return str(s).strip().lower().replace("_", " ")


def find_column(fieldnames, candidates):
    normalized = {norm_header(x): x for x in fieldnames}
    for candidate in candidates:
        c = norm_header(candidate)
        if c in normalized:
            return normalized[c]
    return None


def read_mapping(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Mapping file not found: {path}")

    with open(path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            raise RuntimeError(f"No header found in mapping file: {path}")

        bos_num_col = find_column(
            reader.fieldnames,
            ["Bos residue", "Bos_residue", "1TCO residue", "1TCO_residue"]
        )
        bos_name_col = find_column(
            reader.fieldnames,
            ["Residue (A)", "Bos residue name", "Bos_residue_name",
             "1TCO residue name", "1TCO_residue_name"]
        )

        # Because the historical mapping has two columns both visually named
        # "Residue (A)", csv.DictReader may collapse duplicate names.
        # Therefore, if duplicate headers exist, re-read positionally below.
        fieldnames = reader.fieldnames

    # Re-read positionally to robustly support the historical four-column file:
    # Bos residue | Residue (A) | Candida residue | Residue (A)
    with open(path, newline="", encoding="utf-8-sig") as f:
        rows = list(csv.reader(f))

    if len(rows) < 2:
        raise RuntimeError("Mapping file contains no residue mappings.")

    header = [x.strip() for x in rows[0]]

    # Preferred historical layout: first four columns.
    if len(header) >= 4:
        h0, h1, h2, h3 = [norm_header(x) for x in header[:4]]
        looks_historical = (
            ("bos" in h0 or "1tco" in h0)
            and ("candida" in h2 or "6tz6" in h2)
        )
    else:
        looks_historical = False

    mappings = []

    if looks_historical:
        for line_no, row in enumerate(rows[1:], start=2):
            if not row or all(not x.strip() for x in row):
                continue
            if len(row) < 4:
                raise RuntimeError(
                    f"Mapping row {line_no} has fewer than 4 columns: {row}"
                )
            try:
                bos_resnum = int(float(row[0].strip()))
                candida_resnum = int(float(row[2].strip()))
            except ValueError as exc:
                raise RuntimeError(
                    f"Invalid residue number in mapping row {line_no}: {row}"
                ) from exc

            mappings.append({
                "bos_resnum": bos_resnum,
                "bos_resname": row[1].strip(),
                "candida_resnum": candida_resnum,
                "candida_resname": row[3].strip(),
            })
    else:
        # Flexible named-column fallback for cleaner future mapping files.
        with open(path, newline="", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f)
            fns = reader.fieldnames or []

            bos_num_col = find_column(
                fns, ["Bos residue", "Bos_residue", "1TCO residue", "1TCO_residue"]
            )
            cand_num_col = find_column(
                fns, ["Candida residue", "Candida_residue",
                      "6TZ6 residue", "6TZ6_residue"]
            )
            bos_name_col = find_column(
                fns, ["Bos residue name", "Bos_residue_name",
                      "1TCO residue name", "1TCO_residue_name"]
            )
            cand_name_col = find_column(
                fns, ["Candida residue name", "Candida_residue_name",
                      "6TZ6 residue name", "6TZ6_residue_name"]
            )

            if not bos_num_col or not cand_num_col:
                raise RuntimeError(
                    "Could not identify residue-number columns in mapping file.\n"
                    f"Headers found: {fns}"
                )

            for line_no, row in enumerate(reader, start=2):
                if not any(str(v).strip() for v in row.values() if v is not None):
                    continue
                mappings.append({
                    "bos_resnum": int(float(row[bos_num_col])),
                    "bos_resname": row.get(bos_name_col, "").strip()
                        if bos_name_col else "",
                    "candida_resnum": int(float(row[cand_num_col])),
                    "candida_resname": row.get(cand_name_col, "").strip()
                        if cand_name_col else "",
                })

    if not mappings:
        raise RuntimeError("No residue mappings were read.")

    return mappings


def read_residue_index(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Residue list not found: {path}")

    with open(path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        required = {"Matrix_index", "Chain", "Residue_number", "Residue_name"}
        if not reader.fieldnames or not required.issubset(set(reader.fieldnames)):
            raise RuntimeError(
                f"Unexpected residue-list columns in {path}.\n"
                f"Expected: {sorted(required)}\n"
                f"Found: {reader.fieldnames}"
            )

        by_resnum = {}
        rows = []
        for row in reader:
            idx0 = int(row["Matrix_index"]) - 1
            resnum = int(row["Residue_number"])
            resname = row["Residue_name"].strip()
            chain = row["Chain"].strip()

            if resnum in by_resnum:
                raise RuntimeError(
                    f"Duplicate residue number {resnum} in {path}. "
                    "This script expects one C-alpha per residue in one chain."
                )

            info = {
                "index": idx0,
                "resnum": resnum,
                "resname": resname,
                "chain": chain,
            }
            by_resnum[resnum] = info
            rows.append(info)

    return by_resnum, rows


def load_matrix(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Cross-correlation matrix not found: {path}")
    matrix = np.loadtxt(path, delimiter=",")
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise RuntimeError(f"Matrix is not square: {path} -> {matrix.shape}")
    return matrix


def aa_equal(a, b):
    if not a or not b:
        return True
    return a.strip().upper() == b.strip().upper()


def write_difference_matrix(path, matrix):
    np.savetxt(path, matrix, delimiter=",", fmt="%.10f")


def write_table(path, rows):
    columns = [
        "Rank",
        "Bos_residue",
        "Bos_residue_name",
        "Candida_residue",
        "Candida_residue_name",
        "DeltaCC",
    ]
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        for rank, row in enumerate(rows, start=1):
            writer.writerow({
                "Rank": rank,
                "Bos_residue": row["bos_resnum"],
                "Bos_residue_name": row["bos_resname"],
                "Candida_residue": row["candida_resnum"],
                "Candida_residue_name": row["candida_resname"],
                "DeltaCC": f'{row["delta_cc"]:.10f}',
            })


def main():
    ap = argparse.ArgumentParser(
        description="Calculate residue-wise DeltaCC and top/bottom 20 for A subunit."
    )
    ap.add_argument("--mapping", default="mapping_A_new.csv")
    ap.add_argument("--bos-matrix", default="1TCO_A_crosscorr.csv")
    ap.add_argument("--candida-matrix", default="6TZ6_A_crosscorr.csv")
    ap.add_argument("--bos-residues", default="1TCO_A_residues.csv")
    ap.add_argument("--candida-residues", default="6TZ6_A_residues.csv")
    ap.add_argument("--top-n", type=int, default=20)
    args = ap.parse_args()

    mappings = read_mapping(args.mapping)

    bos_index, bos_rows = read_residue_index(args.bos_residues)
    cand_index, cand_rows = read_residue_index(args.candida_residues)

    bos_cc = load_matrix(args.bos_matrix)
    cand_cc = load_matrix(args.candida_matrix)

    if bos_cc.shape[0] != len(bos_rows):
        raise RuntimeError(
            f"1TCO matrix size {bos_cc.shape[0]} does not match residue list "
            f"length {len(bos_rows)}."
        )
    if cand_cc.shape[0] != len(cand_rows):
        raise RuntimeError(
            f"6TZ6 matrix size {cand_cc.shape[0]} does not match residue list "
            f"length {len(cand_rows)}."
        )

    bos_indices = []
    cand_indices = []
    enriched = []
    name_warnings = []

    seen_bos = set()
    seen_cand = set()

    for m in mappings:
        br = m["bos_resnum"]
        cr = m["candida_resnum"]

        if br not in bos_index:
            raise RuntimeError(f"Mapped 1TCO residue {br} not found in residue list.")
        if cr not in cand_index:
            raise RuntimeError(f"Mapped 6TZ6 residue {cr} not found in residue list.")

        if br in seen_bos:
            raise RuntimeError(f"Duplicate mapped 1TCO residue: {br}")
        if cr in seen_cand:
            raise RuntimeError(f"Duplicate mapped 6TZ6 residue: {cr}")
        seen_bos.add(br)
        seen_cand.add(cr)

        binfo = bos_index[br]
        cinfo = cand_index[cr]

        # Use actual PDB residue names in outputs; compare to mapping as QC.
        if m["bos_resname"] and not aa_equal(m["bos_resname"], binfo["resname"]):
            name_warnings.append(
                f"1TCO {br}: mapping={m['bos_resname']} PDB={binfo['resname']}"
            )
        if m["candida_resname"] and not aa_equal(
            m["candida_resname"], cinfo["resname"]
        ):
            name_warnings.append(
                f"6TZ6 {cr}: mapping={m['candida_resname']} PDB={cinfo['resname']}"
            )

        bos_indices.append(binfo["index"])
        cand_indices.append(cinfo["index"])

        enriched.append({
            "bos_resnum": br,
            "bos_resname": binfo["resname"],
            "candida_resnum": cr,
            "candida_resname": cinfo["resname"],
        })

    # Reduce and reorder both matrices to exactly the same homologous-residue order.
    bos_aligned = bos_cc[np.ix_(bos_indices, bos_indices)]
    cand_aligned = cand_cc[np.ix_(cand_indices, cand_indices)]

    if bos_aligned.shape != cand_aligned.shape:
        raise RuntimeError(
            f"Aligned matrices have different shapes: "
            f"{bos_aligned.shape} vs {cand_aligned.shape}"
        )

    # Historical DeltaCC definition:
    # mean absolute cross-correlation difference across all mapped partners.
    absolute_difference = np.abs(cand_aligned - bos_aligned)
    residue_difference = absolute_difference.mean(axis=1)

    for row, value in zip(enriched, residue_difference):
        row["delta_cc"] = float(value)

    # Save full mapped absolute-difference matrix.
    diff_matrix_file = "A_crosscorr_absolute_difference.csv"
    write_difference_matrix(diff_matrix_file, absolute_difference)

    # Save all mapped residues in mapping order.
    all_file = "A_crosscorr_residue_difference.csv"
    write_table(all_file, enriched)

    top_n = min(args.top_n, len(enriched))

    most = sorted(
        enriched, key=lambda x: x["delta_cc"], reverse=True
    )[:top_n]

    least = sorted(
        enriched, key=lambda x: x["delta_cc"]
    )[:top_n]

    most_file = "A_top20_most_different.csv"
    least_file = "A_top20_least_different.csv"
    write_table(most_file, most)
    write_table(least_file, least)

    # QC summary.
    max_asym = float(np.max(np.abs(absolute_difference - absolute_difference.T)))
    diag_mean = float(np.mean(np.diag(absolute_difference)))

    print("\n" + "=" * 72)
    print("A-SUBUNIT DELTACC COMPARISON")
    print("=" * 72)
    print(f"1TCO matrix shape          : {bos_cc.shape}")
    print(f"6TZ6 matrix shape          : {cand_cc.shape}")
    print(f"Mapped homologous residues : {len(enriched)}")
    print(f"Aligned matrix shape       : {bos_aligned.shape}")
    print(f"Difference max asymmetry   : {max_asym:.3e}")
    print(f"Mean diagonal difference   : {diag_mean:.3e}")
    print(f"Residue-name QC warnings   : {len(name_warnings)}")

    if name_warnings:
        print("\nResidue-name warnings (first 10):")
        for w in name_warnings[:10]:
            print("  WARNING:", w)
        if len(name_warnings) > 10:
            print(f"  ... plus {len(name_warnings) - 10} more")

    print("\nTop 5 MOST different:")
    for rank, row in enumerate(most[:5], start=1):
        print(
            f"{rank:2d}. 1TCO {row['bos_resname']}{row['bos_resnum']}  "
            f"<-> 6TZ6 {row['candida_resname']}{row['candida_resnum']}  "
            f"DeltaCC={row['delta_cc']:.6f}"
        )

    print("\nTop 5 LEAST different:")
    for rank, row in enumerate(least[:5], start=1):
        print(
            f"{rank:2d}. 1TCO {row['bos_resname']}{row['bos_resnum']}  "
            f"<-> 6TZ6 {row['candida_resname']}{row['candida_resnum']}  "
            f"DeltaCC={row['delta_cc']:.6f}"
        )

    print("\nCreated:")
    print(f"  {diff_matrix_file}")
    print(f"  {all_file}")
    print(f"  {most_file}")
    print(f"  {least_file}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        sys.exit(1)
