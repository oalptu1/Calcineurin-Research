#!/usr/bin/env python3
"""
ProDy ANM cross-correlation calculation for FKBP12 / C subunit.

Inputs (default):
    1TCO_C_Subunit.pdb
    6TZ6_C_Subunit.pdb

Outputs:
    1TCO_C_crosscorr.csv
    6TZ6_C_crosscorr.csv
    1TCO_C_residues.csv
    6TZ6_C_residues.csv

Method:
    - protein C-alpha atoms only
    - chain C
    - ANM cutoff = 15.0 Angstrom
    - gamma = 1.0
    - 20 non-zero ANM modes
    - ProDy calcCrossCorr()
"""

import argparse
import csv
import os
import sys

import numpy as np
from prody import ANM, calcCrossCorr, parsePDB


def analyze_structure(pdbfile, label, chain, cutoff, gamma, n_modes, outfile):
    print("\n" + "=" * 72)
    print(f"Processing: {label}")
    print(f"Input     : {pdbfile}")
    print("=" * 72)

    if not os.path.isfile(pdbfile):
        raise FileNotFoundError(f"Input PDB not found: {pdbfile}")

    structure = parsePDB(pdbfile)
    if structure is None:
        raise RuntimeError(f"ProDy could not parse: {pdbfile}")

    ca = structure.select(f"chain {chain} and protein and name CA")
    if ca is None or ca.numAtoms() == 0:
        raise RuntimeError(
            f"No protein C-alpha atoms found for chain {chain} in {pdbfile}"
        )

    n_ca = ca.numAtoms()
    print(f"Chain              : {chain}")
    print(f"C-alpha atoms      : {n_ca}")
    print(f"Residues           : {ca.numResidues()}")
    print(f"ANM cutoff (A)     : {cutoff}")
    print(f"ANM gamma          : {gamma}")
    print(f"Requested modes    : {n_modes}")

    anm = ANM(f"{label}_ANM")
    anm.buildHessian(ca, cutoff=cutoff, gamma=gamma)
    anm.calcModes(n_modes=n_modes, zeros=False)

    actual_modes = anm.numModes()
    print(f"Calculated modes   : {actual_modes}")

    if actual_modes != n_modes:
        raise RuntimeError(
            f"Expected {n_modes} ANM modes but ProDy returned {actual_modes}"
        )

    crosscorr = calcCrossCorr(anm)

    if crosscorr.shape != (n_ca, n_ca):
        raise RuntimeError(
            f"Unexpected cross-correlation shape {crosscorr.shape}; "
            f"expected ({n_ca}, {n_ca})"
        )

    diag = np.diag(crosscorr)
    max_asymmetry = float(np.max(np.abs(crosscorr - crosscorr.T)))
    diag_deviation = float(np.max(np.abs(diag - 1.0)))
    min_cc = float(np.min(crosscorr))
    max_cc = float(np.max(crosscorr))

    print(f"Matrix shape       : {crosscorr.shape}")
    print(f"CC range           : {min_cc:.6f} to {max_cc:.6f}")
    print(f"Max asymmetry      : {max_asymmetry:.3e}")
    print(f"Max |diag - 1|     : {diag_deviation:.3e}")

    np.savetxt(outfile, crosscorr, delimiter=",", fmt="%.10f")
    print(f"Cross-corr written : {outfile}")

    residue_file = outfile.replace("_crosscorr.csv", "_residues.csv")
    resnums = ca.getResnums()
    resnames = ca.getResnames()
    chids = ca.getChids()

    with open(residue_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Matrix_index", "Chain", "Residue_number", "Residue_name"])
        for i, (chid, resnum, resname) in enumerate(
            zip(chids, resnums, resnames), start=1
        ):
            writer.writerow([i, chid, int(resnum), resname])

    print(f"Residue list written: {residue_file}")


def main():
    ap = argparse.ArgumentParser(
        description="Calculate ProDy ANM cross-correlation matrices for C subunit."
    )
    ap.add_argument(
        "--bos",
        default="1TCO_C_Subunit.pdb",
        help="1TCO C-subunit PDB (default: 1TCO_C_Subunit.pdb)",
    )
    ap.add_argument(
        "--candida",
        default="6TZ6_C_Subunit.pdb",
        help="6TZ6 C-subunit PDB (default: 6TZ6_C_Subunit.pdb)",
    )
    ap.add_argument(
        "--chain",
        default="C",
        help="Protein chain ID to analyze (default: C)",
    )
    ap.add_argument("--cutoff", type=float, default=15.0)
    ap.add_argument("--gamma", type=float, default=1.0)
    ap.add_argument("--modes", type=int, default=20)
    args = ap.parse_args()

    if args.modes <= 0:
        raise ValueError("--modes must be a positive integer")
    if args.cutoff <= 0:
        raise ValueError("--cutoff must be positive")
    if args.gamma <= 0:
        raise ValueError("--gamma must be positive")

    analyze_structure(
        args.bos, "1TCO_C", args.chain,
        args.cutoff, args.gamma, args.modes,
        "1TCO_C_crosscorr.csv",
    )
    analyze_structure(
        args.candida, "6TZ6_C", args.chain,
        args.cutoff, args.gamma, args.modes,
        "6TZ6_C_crosscorr.csv",
    )

    print("\n" + "=" * 72)
    print("FINISHED")
    print("=" * 72)
    print("Created:")
    print("  1TCO_C_crosscorr.csv")
    print("  6TZ6_C_crosscorr.csv")
    print("  1TCO_C_residues.csv")
    print("  6TZ6_C_residues.csv")
    print("\nNext step: mapping-based DeltaCC comparison and top/bottom 20.")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        sys.exit(1)
