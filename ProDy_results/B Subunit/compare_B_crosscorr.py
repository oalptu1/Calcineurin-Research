import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser


# ============================================================
# DOSYA İSİMLERİ
# ============================================================

BOS_CORR = "1TCO_B_crosscorr.csv"
CANDIDA_CORR = "6TZ6_B_crosscorr.csv"

BOS_PDB = "1TCO_B.pdb"
CANDIDA_PDB = "6TZ6_B.pdb"

MAPPING_FILE = "mapping_B2.csv"


# ============================================================
# 1. CROSS-CORRELATION MATRİSLERİNİ OKU
# ============================================================

print("Cross-correlation matrices loading...")

bos_corr = np.loadtxt(BOS_CORR, delimiter=",")
candida_corr = np.loadtxt(CANDIDA_CORR, delimiter=",")

print("1TCO matrix :", bos_corr.shape)
print("6TZ6 matrix :", candida_corr.shape)


# ============================================================
# 2. PDB DOSYALARINDAN RESIDUE NUMARALARINI AL
# ============================================================

parser = PDBParser(QUIET=True)


def get_residue_numbers(pdb_file):

    structure = parser.get_structure("structure", pdb_file)

    residue_numbers = []

    for model in structure:

        for chain in model:

            for residue in chain:

                # yalnızca protein residue'leri
                if "CA" in residue:

                    residue_numbers.append(
                        residue.id[1]
                    )

        # sadece ilk model
        break

    return residue_numbers


print("\nReading PDB residue numbering...")

bos_residues = get_residue_numbers(BOS_PDB)
candida_residues = get_residue_numbers(CANDIDA_PDB)

print(
    "1TCO residues:",
    len(bos_residues)
)

print(
    "6TZ6 residues:",
    len(candida_residues)
)


# ============================================================
# 3. RESIDUE NUMARASI → MATRİS İNDEKSİ
# ============================================================

bos_index = {
    residue: i
    for i, residue in enumerate(bos_residues)
}

candida_index = {
    residue: i
    for i, residue in enumerate(candida_residues)
}


# ============================================================
# 4. MAPPING DOSYASINI OKU
# ============================================================

print("\nReading residue mapping...")

mapping = pd.read_csv(MAPPING_FILE)

print(
    "Mapping columns:",
    mapping.columns.tolist()
)


# ============================================================
# 5. ORTAK RESIDUE'LERİ BUL
# ============================================================

bos_indices = []
candida_indices = []

valid_rows = []

for _, row in mapping.iterrows():

    bos_res = int(row["Bos residue"])
    candida_res = int(row["Candida residue"])

    if (
        bos_res in bos_index
        and candida_res in candida_index
    ):

        bos_indices.append(
            bos_index[bos_res]
        )

        candida_indices.append(
            candida_index[candida_res]
        )

        valid_rows.append(row)


mapping_used = pd.DataFrame(valid_rows)


print(
    "\nMapped residues used:",
    len(mapping_used)
)


# ============================================================
# 6. MATRİSLERİ MAPPING'E GÖRE KÜÇÜLT
# ============================================================

print("\nAligning matrices using residue mapping...")

bos_aligned = bos_corr[
    np.ix_(
        bos_indices,
        bos_indices
    )
]

candida_aligned = candida_corr[
    np.ix_(
        candida_indices,
        candida_indices
    )
]


print(
    "Aligned 1TCO matrix:",
    bos_aligned.shape
)

print(
    "Aligned 6TZ6 matrix:",
    candida_aligned.shape
)


# ============================================================
# 7. CROSS-CORRELATION FARKI
# ============================================================

difference = np.abs(
    candida_aligned - bos_aligned
)


# ============================================================
# 8. HER RESIDUE İÇİN ORTALAMA FARK
# ============================================================

residue_difference = difference.mean(axis=1)


# ============================================================
# 9. SONUÇ TABLOSU
# ============================================================

results = mapping_used.copy()

results["Dynamic_Difference"] = residue_difference


# ============================================================
# 10. EN FAZLA FARKLILAŞAN RESIDUE'LER
# ============================================================

most_different = results.sort_values(
    "Dynamic_Difference",
    ascending=False
)


print("\n")
print("==============================================")
print("TOP 20 MOST DIFFERENT RESIDUES")
print("==============================================")

print(
    most_different.head(20).to_string(
        index=False
    )
)


# ============================================================
# 11. EN AZ FARKLILAŞAN RESIDUE'LER
# ============================================================

least_different = results.sort_values(
    "Dynamic_Difference",
    ascending=True
)


print("\n")
print("==============================================")
print("TOP 20 LEAST DIFFERENT RESIDUES")
print("==============================================")

print(
    least_different.head(20).to_string(
        index=False
    )
)


# ============================================================
# 12. DOSYALARI KAYDET
# ============================================================

results.to_csv(
    "B_crosscorr_residue_difference.csv",
    index=False
)

results.to_excel(
    "B_crosscorr_residue_difference.xlsx",
    index=False
)

most_different.head(20).to_csv(
    "B_top20_most_different.csv",
    index=False
)

least_different.head(20).to_csv(
    "B_top20_least_different.csv",
    index=False
)


# ============================================================
# 13. FARK MATRİSİNİ KAYDET
# ============================================================

np.savetxt(
    "B_crosscorr_difference_matrix.csv",
    difference,
    delimiter=","
)


# ============================================================
# 14. HEATMAP
# ============================================================

plt.figure(
    figsize=(10, 8)
)

plt.imshow(
    difference,
    cmap="hot",
    origin="lower",
    interpolation="nearest"
)

plt.colorbar(
    label="Absolute cross-correlation difference"
)

plt.xlabel(
    "Mapped residue index"
)

plt.ylabel(
    "Mapped residue index"
)

plt.title(
    "1TCO-B vs 6TZ6-B Cross-Correlation Difference"
)

plt.tight_layout()

plt.savefig(
    "B_crosscorr_difference_heatmap.png",
    dpi=600,
    bbox_inches="tight"
)

plt.close()


# ============================================================
# 15. BİTTİ
# ============================================================

print("\n")
print("==============================================")
print("ANALYSIS COMPLETED SUCCESSFULLY")
print("==============================================")

print(
    "\nOutput files:"
)

print(
    "B_crosscorr_residue_difference.csv"
)

print(
    "B_crosscorr_residue_difference.xlsx"
)

print(
    "B_top20_most_different.csv"
)

print(
    "B_top20_least_different.csv"
)

print(
    "B_crosscorr_difference_matrix.csv"
)

print(
    "B_crosscorr_difference_heatmap.png"
)