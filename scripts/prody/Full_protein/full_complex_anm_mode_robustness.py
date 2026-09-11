# full_complex_anm_mode_robustness.py
#
# Robustness analysis for full-complex ANM comparison between
# Bos taurus 1TCO and Candida albicans 6TZ6.
#
# The script evaluates 10, 20, 30 and 50 lowest-frequency non-zero modes.
# It keeps CnA + CnB + FKBP12 in a single ANM network and excludes FK506,
# waters, metals and other heteroatoms by selecting only protein C-alpha atoms.
#
# Required files in the same folder:
#   1TCO_prepared.pdb
#   6TZ6_prepared.pdb
#
# Required packages:
#   python -m pip install prody biopython pandas numpy

import numpy as np
import pandas as pd

from prody import parsePDB, ANM, calcCrossCorr
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.SeqUtils import seq1


# ============================================================
# USER SETTINGS
# ============================================================

PDB_1TCO = "1TCO_prepared.pdb"
PDB_6TZ6 = "6TZ6_prepared.pdb"

MODE_COUNTS = [10, 20, 30, 50]
MAX_MODES = max(MODE_COUNTS)

CUTOFF = 15.0
GAMMA = 1.0

TOP_RESIDUES = 20
TOP_PAIRS = 100

CHAINS_1TCO = {
    "CnA": "A",
    "CnB": "B",
    "FKBP12": "C",
}

CHAINS_6TZ6 = {
    "CnA": "A",
    "CnB": "B",
    "FKBP12": "C",
}


# ============================================================
# RESIDUE NAME NORMALIZATION
# ============================================================

RESNAME_FIX = {
    "HID": "HIS",
    "HIE": "HIS",
    "HIP": "HIS",
    "ASH": "ASP",
    "GLH": "GLU",
    "LYN": "LYS",
    "CYX": "CYS",
    "CYM": "CYS",
}


def one_letter(resname):
    resname = RESNAME_FIX.get(resname.upper(), resname.upper())
    try:
        return seq1(resname)
    except Exception:
        return "X"


# ============================================================
# LOAD FULL PROTEIN COMPLEX
# ============================================================

def load_complex(filename, chains):
    ag = parsePDB(filename)
    chain_string = " ".join(chains.values())

    ca = ag.select(f"protein and name CA and chain {chain_string}")

    if ca is None:
        raise RuntimeError(
            f"No C-alpha atoms found in {filename}. Check the chain IDs."
        )

    print(f"\n{filename}")
    print(f"Total C-alpha atoms: {ca.numAtoms()}")

    for name, chain in chains.items():
        n = np.sum(ca.getChids() == chain)
        print(f"  {name} (chain {chain}): {n} residues")

    return ca


# ============================================================
# BUILD ANM ONCE UP TO 50 MODES
# ============================================================

def build_anm(ca, name):
    print(f"\nBuilding ANM for {name}...")
    anm = ANM(name)
    anm.buildHessian(ca, cutoff=CUTOFF, gamma=GAMMA)
    anm.calcModes(n_modes=MAX_MODES, zeros=False)

    if anm.numModes() < MAX_MODES:
        raise RuntimeError(
            f"{name} returned only {anm.numModes()} non-zero modes; "
            f"{MAX_MODES} are required."
        )

    print(f"{name}: calculated {anm.numModes()} non-zero modes")
    return anm


# ============================================================
# SEQUENCE ALIGNMENT
# ============================================================

def get_chain_positions(ca, chain):
    return np.where(ca.getChids() == chain)[0]


def get_sequence(ca, positions):
    resnames = ca.getResnames()[positions]
    return "".join(one_letter(res) for res in resnames)


def align_chains(ca1, chain1, ca2, chain2):
    pos1 = get_chain_positions(ca1, chain1)
    pos2 = get_chain_positions(ca2, chain2)

    seq1_string = get_sequence(ca1, pos1)
    seq2_string = get_sequence(ca2, pos2)

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    alignment = aligner.align(seq1_string, seq2_string)[0]
    blocks1, blocks2 = alignment.aligned

    mapped1 = []
    mapped2 = []

    for block1, block2 in zip(blocks1, blocks2):
        start1, end1 = block1
        start2, end2 = block2

        if (end1 - start1) != (end2 - start2):
            raise RuntimeError("Unexpected unequal aligned block lengths.")

        for i, j in zip(range(start1, end1), range(start2, end2)):
            mapped1.append(pos1[i])
            mapped2.append(pos2[j])

    return np.array(mapped1, dtype=int), np.array(mapped2, dtype=int)


# ============================================================
# FULL-COMPLEX HOMOLOGOUS MAPPING
# ============================================================

def build_mapping(ca1, ca2):
    idx1_all = []
    idx2_all = []
    annotations = []

    for subunit in ["CnA", "CnB", "FKBP12"]:
        ch1 = CHAINS_1TCO[subunit]
        ch2 = CHAINS_6TZ6[subunit]

        idx1, idx2 = align_chains(ca1, ch1, ca2, ch2)

        print(f"{subunit}: {len(idx1)} mapped residues")

        for i1, i2 in zip(idx1, idx2):
            annotations.append(
                {
                    "Subunit": subunit,
                    "1TCO_chain": ca1.getChids()[i1],
                    "1TCO_resnum": int(ca1.getResnums()[i1]),
                    "1TCO_resname": ca1.getResnames()[i1],
                    "6TZ6_chain": ca2.getChids()[i2],
                    "6TZ6_resnum": int(ca2.getResnums()[i2]),
                    "6TZ6_resname": ca2.getResnames()[i2],
                }
            )

        idx1_all.extend(idx1)
        idx2_all.extend(idx2)

    ann = pd.DataFrame(annotations)

    # Stable residue-pair identifier
    ann["Residue_pair"] = (
        ann["Subunit"]
        + ":"
        + ann["1TCO_resname"].astype(str)
        + ann["1TCO_resnum"].astype(str)
        + "/"
        + ann["6TZ6_resname"].astype(str)
        + ann["6TZ6_resnum"].astype(str)
    )

    return (
        np.array(idx1_all, dtype=int),
        np.array(idx2_all, dtype=int),
        ann,
    )


# ============================================================
# INTERFACE HELPERS
# ============================================================

def block_indices(annotation_df, subunit):
    return np.where(annotation_df["Subunit"].values == subunit)[0]


def interface_mean_abs_delta(delta_matrix, annotation_df, unit1, unit2):
    i = block_indices(annotation_df, unit1)
    j = block_indices(annotation_df, unit2)
    block = delta_matrix[np.ix_(i, j)]
    return float(np.mean(np.abs(block)))


# ============================================================
# TOP INTER-SUBUNIT PAIRS
# ============================================================

def get_top_inter_subunit_pairs(cc1, cc6, annotation_df, mode_count, top_n=100):
    records = []
    n = len(annotation_df)

    for i in range(n):
        for j in range(i + 1, n):
            sub_i = annotation_df.iloc[i]["Subunit"]
            sub_j = annotation_df.iloc[j]["Subunit"]

            if sub_i == sub_j:
                continue

            diff = cc6[i, j] - cc1[i, j]

            pair_label_1tco = (
                f"{sub_i}:{annotation_df.iloc[i]['1TCO_resname']}"
                f"{annotation_df.iloc[i]['1TCO_resnum']}"
                f"--{sub_j}:{annotation_df.iloc[j]['1TCO_resname']}"
                f"{annotation_df.iloc[j]['1TCO_resnum']}"
            )

            pair_label_6tz6 = (
                f"{sub_i}:{annotation_df.iloc[i]['6TZ6_resname']}"
                f"{annotation_df.iloc[i]['6TZ6_resnum']}"
                f"--{sub_j}:{annotation_df.iloc[j]['6TZ6_resname']}"
                f"{annotation_df.iloc[j]['6TZ6_resnum']}"
            )

            # Order-independent stable identifier based on mapped-row indices
            pair_id = f"{min(i,j)}__{max(i,j)}"

            records.append(
                {
                    "Mode_count": mode_count,
                    "Pair_ID": pair_id,
                    "Subunit_1": sub_i,
                    "Subunit_2": sub_j,
                    "1TCO_pair": pair_label_1tco,
                    "6TZ6_pair": pair_label_6tz6,
                    "CC_1TCO": cc1[i, j],
                    "CC_6TZ6": cc6[i, j],
                    "Difference_6TZ6_minus_1TCO": diff,
                    "Absolute_difference": abs(diff),
                }
            )

    df = pd.DataFrame(records)
    df = df.sort_values("Absolute_difference", ascending=False).head(top_n).copy()
    df["Rank_within_mode"] = np.arange(1, len(df) + 1)
    return df


# ============================================================
# OPTIONAL KEY RESIDUE PANEL
# ============================================================
# These positions are biologically relevant to the present manuscript.
# Rows are included only if the exact 1TCO residue is present in the mapping.

KEY_1TCO_RESIDUES = [
    ("CnA", 315),
    ("CnB", 115),
    ("CnB", 118),
    ("CnB", 119),
    ("CnB", 169),
    ("FKBP12", 55),
    ("FKBP12", 56),
    ("FKBP12", 59),
    ("FKBP12", 61),
    ("FKBP12", 62),
    ("FKBP12", 63),
    ("FKBP12", 64),
    ("FKBP12", 65),
]


# ============================================================
# MAIN
# ============================================================

print("\n===================================================")
print("FULL-COMPLEX ANM MODE ROBUSTNESS ANALYSIS")
print("Modes:", MODE_COUNTS)
print("===================================================")

# Load structures
ca1 = load_complex(PDB_1TCO, CHAINS_1TCO)
ca6 = load_complex(PDB_6TZ6, CHAINS_6TZ6)

# Build ANMs once
anm1 = build_anm(ca1, "1TCO_full_complex")
anm6 = build_anm(ca6, "6TZ6_full_complex")

# Build homologous mapping once
idx1, idx6, annotation_ordered = build_mapping(ca1, ca6)

interfaces = [
    ("CnA", "CnB"),
    ("CnA", "FKBP12"),
    ("CnB", "FKBP12"),
]

all_interface_rows = []
all_residue_rows = []
all_top_pair_rows = []
all_key_rows = []

for n_modes in MODE_COUNTS:
    print(f"\n--- Processing {n_modes} modes ---")

    # Cross-correlation from the first n lowest-frequency non-zero modes
    cc1_full = calcCrossCorr(anm1[:n_modes])
    cc6_full = calcCrossCorr(anm6[:n_modes])

    # Reduce to homologous mapped residues in identical order
    cc1 = cc1_full[np.ix_(idx1, idx1)]
    cc6 = cc6_full[np.ix_(idx6, idx6)]

    delta_signed = cc6 - cc1
    delta_abs = np.abs(delta_signed)

    # --------------------------------------------------------
    # Interface summary
    # --------------------------------------------------------
    interface_values = []

    for u1, u2 in interfaces:
        value = interface_mean_abs_delta(
            delta_signed,
            annotation_ordered,
            u1,
            u2
        )

        interface_values.append((f"{u1}-{u2}", value))

    # Rank interfaces within each mode count
    interface_values_sorted = sorted(
        interface_values,
        key=lambda x: x[1],
        reverse=True
    )

    interface_rank = {
        name: rank
        for rank, (name, _) in enumerate(interface_values_sorted, start=1)
    }

    for name, value in interface_values:
        all_interface_rows.append(
            {
                "Mode_count": n_modes,
                "Interface": name,
                "Mean_absolute_CC_difference": value,
                "Interface_rank": interface_rank[name],
            }
        )

    # --------------------------------------------------------
    # Residue-wise ΔCC
    # --------------------------------------------------------
    residue_delta = np.mean(delta_abs, axis=1)

    res_df = annotation_ordered.copy()
    res_df["Mode_count"] = n_modes
    res_df["DeltaCC_full_complex"] = residue_delta
    res_df["Rank_within_mode"] = (
        pd.Series(residue_delta)
        .rank(method="min", ascending=False)
        .astype(int)
        .values
    )
    res_df["Top20"] = res_df["Rank_within_mode"] <= TOP_RESIDUES

    all_residue_rows.append(res_df)

    # --------------------------------------------------------
    # Key residues
    # --------------------------------------------------------
    for subunit, resnum in KEY_1TCO_RESIDUES:
        hit = res_df[
            (res_df["Subunit"] == subunit)
            & (res_df["1TCO_resnum"] == resnum)
        ]

        if len(hit) == 1:
            row = hit.iloc[0]
            all_key_rows.append(
                {
                    "Mode_count": n_modes,
                    "Subunit": subunit,
                    "1TCO_residue": f"{row['1TCO_resname']}{row['1TCO_resnum']}",
                    "6TZ6_residue": f"{row['6TZ6_resname']}{row['6TZ6_resnum']}",
                    "DeltaCC_full_complex": row["DeltaCC_full_complex"],
                    "Rank_within_mode": int(row["Rank_within_mode"]),
                    "Top20": bool(row["Top20"]),
                }
            )

    # --------------------------------------------------------
    # Top inter-subunit residue-pair differences
    # --------------------------------------------------------
    top_pairs_df = get_top_inter_subunit_pairs(
        cc1,
        cc6,
        annotation_ordered,
        mode_count=n_modes,
        top_n=TOP_PAIRS
    )
    all_top_pair_rows.append(top_pairs_df)


# ============================================================
# COMBINE OUTPUTS
# ============================================================

interface_df = pd.DataFrame(all_interface_rows)

residue_df = pd.concat(all_residue_rows, ignore_index=True)

key_df = pd.DataFrame(all_key_rows)

top_pairs_df = pd.concat(all_top_pair_rows, ignore_index=True)


# ============================================================
# RESIDUE STABILITY ACROSS MODE COUNTS
# ============================================================

residue_stability = (
    residue_df.groupby(
        [
            "Residue_pair",
            "Subunit",
            "1TCO_resname",
            "1TCO_resnum",
            "6TZ6_resname",
            "6TZ6_resnum",
        ],
        as_index=False
    )
    .agg(
        Mean_DeltaCC=("DeltaCC_full_complex", "mean"),
        SD_DeltaCC=("DeltaCC_full_complex", "std"),
        Min_DeltaCC=("DeltaCC_full_complex", "min"),
        Max_DeltaCC=("DeltaCC_full_complex", "max"),
        Mean_rank=("Rank_within_mode", "mean"),
        Best_rank=("Rank_within_mode", "min"),
        Worst_rank=("Rank_within_mode", "max"),
        Top20_count=("Top20", "sum"),
    )
)

residue_stability["Top20_fraction"] = (
    residue_stability["Top20_count"] / len(MODE_COUNTS)
)

residue_stability = residue_stability.sort_values(
    ["Top20_count", "Mean_rank", "Mean_DeltaCC"],
    ascending=[False, True, False]
)


# ============================================================
# INTER-SUBUNIT PAIR STABILITY ACROSS MODE COUNTS
# ============================================================

pair_stability = (
    top_pairs_df.groupby(
        ["Pair_ID", "1TCO_pair", "6TZ6_pair"],
        as_index=False
    )
    .agg(
        Modes_in_top100=("Mode_count", "nunique"),
        Mean_absolute_difference=("Absolute_difference", "mean"),
        Min_absolute_difference=("Absolute_difference", "min"),
        Max_absolute_difference=("Absolute_difference", "max"),
        Mean_rank=("Rank_within_mode", "mean"),
        Best_rank=("Rank_within_mode", "min"),
        Worst_rank=("Rank_within_mode", "max"),
    )
)

pair_stability = pair_stability.sort_values(
    ["Modes_in_top100", "Mean_rank", "Mean_absolute_difference"],
    ascending=[False, True, False]
)


# ============================================================
# INTERFACE ROBUSTNESS PIVOT
# ============================================================

interface_pivot = interface_df.pivot(
    index="Interface",
    columns="Mode_count",
    values="Mean_absolute_CC_difference"
)

interface_pivot.columns = [f"{c}_modes" for c in interface_pivot.columns]
interface_pivot = interface_pivot.reset_index()

rank_pivot = interface_df.pivot(
    index="Interface",
    columns="Mode_count",
    values="Interface_rank"
)

rank_pivot.columns = [f"Rank_{c}_modes" for c in rank_pivot.columns]
rank_pivot = rank_pivot.reset_index()

interface_robustness = interface_pivot.merge(rank_pivot, on="Interface")


# ============================================================
# SAVE OUTPUT FILES
# ============================================================

interface_df.to_csv(
    "mode_robustness_interface_summary_long.csv",
    index=False
)

interface_robustness.to_csv(
    "mode_robustness_interface_summary.csv",
    index=False
)

residue_df.to_csv(
    "mode_robustness_all_residue_rankings.csv",
    index=False
)

residue_stability.to_csv(
    "mode_robustness_residue_stability.csv",
    index=False
)

key_df.to_csv(
    "mode_robustness_key_residues.csv",
    index=False
)

top_pairs_df.to_csv(
    "mode_robustness_top100_inter_subunit_pairs_all_modes.csv",
    index=False
)

pair_stability.to_csv(
    "mode_robustness_inter_subunit_pair_stability.csv",
    index=False
)


# ============================================================
# PRINT COMPACT SUMMARY
# ============================================================

print("\n===================================================")
print("INTERFACE SUMMARY")
print("===================================================")

for n_modes in MODE_COUNTS:
    x = interface_df[interface_df["Mode_count"] == n_modes].sort_values(
        "Mean_absolute_CC_difference",
        ascending=False
    )
    print(f"\n{n_modes} modes:")
    for _, row in x.iterrows():
        print(
            f"  {row['Interface']:12s} "
            f"{row['Mean_absolute_CC_difference']:.6f} "
            f"(rank {int(row['Interface_rank'])})"
        )

print("\n===================================================")
print("TOP RESIDUES STABLE ACROSS MODES")
print("===================================================")

stable_top = residue_stability[
    residue_stability["Top20_count"] == len(MODE_COUNTS)
].head(20)

if len(stable_top) == 0:
    print("No residue remained in the top 20 for all mode counts.")
else:
    for _, row in stable_top.iterrows():
        print(
            f"  {row['Residue_pair']:28s} "
            f"Top20 {int(row['Top20_count'])}/{len(MODE_COUNTS)} modes, "
            f"mean rank {row['Mean_rank']:.1f}"
        )

print("\n===================================================")
print("ANALYSIS COMPLETE")
print("===================================================")

print("\nMain files to send for interpretation:")
print("1. mode_robustness_interface_summary.csv")
print("2. mode_robustness_residue_stability.csv")
print("3. mode_robustness_key_residues.csv")
print("4. mode_robustness_inter_subunit_pair_stability.csv")
print("\nAdditional detailed files are also saved.")
