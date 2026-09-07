#!/usr/bin/env bash
set -Eeuo pipefail

# 6TZ6 FK506/derivatives: one-PDB -> topology -> single-snapshot MM/GBSA pipeline
# Usage:
#   bash run_analysis_6TZ6_from_PDB.sh [input.pdb] [prefix]
# Example:
#   bash run_analysis_6TZ6_from_PDB.sh 6TZ6_Derivative09.pdb 6TZ6_Derivative09

INPUT_PDB="${1:-6TZ6_Derivative04.pdb}"
PREFIX="${2:-6TZ6_Derivative04}"
LIGAND="FK5"
LIGAND_CHARGE=0
VDW_LIMIT=100000
RUN_DIR="${PREFIX}_MMGBSA_$(date +%Y%m%d_%H%M%S)"

die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "'$1' bulunamadı. Önce: conda activate ambertools"; }

[[ -s "$INPUT_PDB" ]] || die "Girdi PDB bulunamadı: $INPUT_PDB"
for exe in python antechamber parmchk2 tleap cpptraj ante-MMPBSA.py MMPBSA.py; do need "$exe"; done

INPUT_PDB="$(python -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$INPUT_PDB")"
mkdir -p "$RUN_DIR"
cd "$RUN_DIR"
exec > >(tee pipeline.log) 2>&1

printf '\n6TZ6 PDB -> MM/GBSA\nInput : %s\nRun   : %s\n\n' "$INPUT_PDB" "$PWD"

# Clean the PDB, retain/recreate chain breaks, split protein and FK5, and run geometry QC.
python - "$INPUT_PDB" "$PREFIX" "$LIGAND" <<'PY'
import math, sys
from collections import defaultdict

src, prefix, ligand = sys.argv[1:]
expected_gaps = {('B', 32, 40), ('B', 84, 88), ('C', 75, 79)}
water = {'HOH','WAT','SOL','TIP3','TIP3P'}

def atom_info(line):
    return (line[12:16].strip(), line[16:17], line[17:20].strip(),
            line[21:22] or ' ', int(line[22:26]), line[76:78].strip().upper())

atoms=[]
lig_all=[]
with open(src, errors='replace') as fh:
    for line in fh:
        if not line.startswith(('ATOM  ','HETATM')): continue
        name, alt, res, chain, num, elem = atom_info(line)
        if alt not in (' ','A'): continue
        record=(line, name, res, chain, num, elem)
        # Preserve Maestro-assigned FK5 hydrogens for Antechamber/AM1-BCC.
        # Protein hydrogens are deliberately removed and rebuilt by TLEaP.
        if res == ligand:
            lig_all.append(record)
        if elem in ('H','D') or name.upper().startswith(('H','D')): continue
        if res in water: continue
        if line.startswith('HETATM') and res != ligand: continue
        atoms.append(record)

lig=[a for a in atoms if a[2] == ligand]
prot=[a for a in atoms if a[2] != ligand]
if not lig: raise SystemExit(f'ERROR: {ligand} bulunamadı.')
lig_h_count=len(lig_all)-len(lig)
if lig_h_count == 0:
    raise SystemExit('ERROR: Ligand hidrojeni bulunamadı. Maestro-prep edilmiş, '
                     'tüm hidrojenleri içeren PDB kullanılmalı.')
lig_res={(a[3],a[4]) for a in lig}
if len(lig_res) != 1: raise SystemExit(f'ERROR: Tek bir {ligand} bekleniyordu; bulundu: {sorted(lig_res)}')

res_order=[]
for a in prot:
    key=(a[3],a[4],a[2])
    if not res_order or res_order[-1] != key: res_order.append(key)
seen_gaps=set()
for (c1,n1,_),(c2,n2,_) in zip(res_order,res_order[1:]):
    if c1 == c2 and n2 != n1+1: seen_gaps.add((c1,n1,n2))
if seen_gaps != expected_gaps:
    raise SystemExit('ERROR: 6TZ6 zincir kesintileri beklenenden farklı.\n'
                     f'Beklenen: {sorted(expected_gaps)}\nBulunan: {sorted(seen_gaps)}')
if len(res_order) != 618:
    raise SystemExit(f'ERROR: 618 protein kalıntısı bekleniyordu; {len(res_order)} bulundu.')
if len(prot) != 4954:
    raise SystemExit(f'ERROR: Protein ağır atom sayısı beklenenden farklı: {len(prot)} (beklenen 4954).')
# Ligand ağır atom ve hidrojen sayıları türeve göre değişir. Makul olmayan bir
# atom sayısını yine de durdurarak yanlış/eksik ligand dosyalarını yakala.
if not 40 <= len(lig) <= 80:
    raise SystemExit(f'ERROR: Ligand ağır atom sayısı şüpheli: {len(lig)} (beklenen aralık 40-80).')

def normalized(a, ligand_mode=False):
    line,name,res,chain,num,elem=a
    line=line[:16]+' '+line[17:]
    if ligand_mode:
        line=line[:21]+'L'+f'{1:4d}'+line[26:]
    return line.rstrip()+'\n'

def write_protein(path):
    with open(path,'w') as out:
        prev=None
        for a in prot:
            key=(a[3],a[4],a[2])
            if prev and key != prev:
                if key[0] != prev[0] or (prev[0],prev[1],key[1]) in expected_gaps:
                    out.write('TER\n')
            out.write(normalized(a))
            prev=key
        out.write('TER\nEND\n')

write_protein(f'{prefix}_protein.pdb')
with open(f'{prefix}_{ligand}.pdb','w') as out:
    for a in lig_all: out.write(normalized(a, True))
    out.write('TER\nEND\n')
with open(f'{prefix}_clean.pdb','w') as out:
    with open(f'{prefix}_protein.pdb') as p: out.write(p.read().replace('END\n',''))
    with open(f'{prefix}_{ligand}.pdb') as l: out.write(l.read())

# Report sub-1 A inter-residue heavy-atom contacts in the supplied coordinates.
coords=[]
for a in atoms:
    line=a[0]; coords.append((a[3],a[4],a[2],a[1],float(line[30:38]),float(line[38:46]),float(line[46:54])))
bad=[]
for i,x in enumerate(coords):
    for y in coords[i+1:]:
        if x[:3] == y[:3]: continue
        d2=(x[4]-y[4])**2+(x[5]-y[5])**2+(x[6]-y[6])**2
        if d2 < 1.0: bad.append((math.sqrt(d2),x[:4],y[:4]))
with open('input_QC.txt','w') as q:
    q.write(f'protein_residues={len(res_order)}\nprotein_heavy_atoms={len(prot)}\nligand_heavy_atoms={len(lig)}\nligand_hydrogen_atoms={lig_h_count}\n')
    q.write(f'chain_gaps={sorted(seen_gaps)}\ninter_residue_heavy_pairs_below_1A={len(bad)}\n')
    for row in sorted(bad)[:20]: q.write(repr(row)+'\n')
if bad: raise SystemExit(f'ERROR: PDB içinde kalıntılar arası <1 A ağır atom çakışması var ({len(bad)} çift). input_QC.txt dosyasına bakın.')
print(f'PASS: 618 protein kalıntısı, 3 doğru TER kesintisi, protein ağır atom sayısı, '
      f'ligand={len(lig)} ağır atom/{lig_h_count} H ve <1 Å clash kontrolü.')
PY

antechamber -fi pdb -fo mol2 -i "${PREFIX}_${LIGAND}.pdb" \
  -o "${PREFIX}_${LIGAND}.mol2" -c bcc -nc "$LIGAND_CHARGE" -at gaff2 -rn "$LIGAND" -s 2
parmchk2 -i "${PREFIX}_${LIGAND}.mol2" -f mol2 -o "${PREFIX}_${LIGAND}.frcmod" -s 2

cat > tleap.in <<EOF
source leaprc.protein.ff14SB
source leaprc.gaff2
loadamberparams ${PREFIX}_${LIGAND}.frcmod
protein = loadpdb ${PREFIX}_protein.pdb
ligand = loadmol2 ${PREFIX}_${LIGAND}.mol2
complex = combine { protein ligand }
check complex
saveamberparm complex ${PREFIX}.prmtop ${PREFIX}.inpcrd
savepdb complex ${PREFIX}_tleap.pdb
quit
EOF
tleap -f tleap.in | tee tleap.log
[[ -s "${PREFIX}.prmtop" && -s "${PREFIX}.inpcrd" ]] || die "TLEaP topolojiyi oluşturamadı. tleap.log dosyasına bakın."

# Detect the unique FK5 topology residue automatically.
LIGAND_RESIDUE="$(python - "${PREFIX}.prmtop" "$LIGAND" <<'PY'
import sys
pdb, target=sys.argv[1:]
lines=open(pdb).read().splitlines(); labels=[]
for i,line in enumerate(lines):
    if line.startswith('%FLAG RESIDUE_LABEL'):
        j=i+2
        while j<len(lines) and not lines[j].startswith('%FLAG'):
            labels += [lines[j][k:k+4].strip() for k in range(0,len(lines[j]),4)]
            j+=1
        break
hits=[i+1 for i,x in enumerate(labels) if x == target]
if len(hits)!=1: raise SystemExit(f'ERROR: Topolojide tek {target} bekleniyordu; indeksler={hits}')
print(hits[0])
PY
)"
[[ "$LIGAND_RESIDUE" == "619" ]] || die "FK5 topoloji kalıntısı 619 olmalıydı; bulunan: $LIGAND_RESIDUE"
printf 'PASS: FK5 topology residue = %s\n' "$LIGAND_RESIDUE" | tee -a input_QC.txt

cat > single_frame.in <<EOF
parm ${PREFIX}.prmtop
trajin ${PREFIX}.inpcrd
trajout ${PREFIX}.nc netcdf
run
EOF
cpptraj -i single_frame.in > cpptraj.log

ante-MMPBSA.py -p "${PREFIX}.prmtop" -c "${PREFIX}.prmtop" \
  -r "${PREFIX}_receptor.prmtop" -l "${PREFIX}_ligand.prmtop" -n ":${LIGAND_RESIDUE}"

cat > mmpbsa.in <<'EOF'
&general
  startframe=1, endframe=1, interval=1,
  verbose=1, keep_files=2,
/
&gb
  igb=5, saltcon=0.150,
/
&decomp
  idecomp=1, dec_verbose=1,
/
EOF

MMPBSA.py -O -i mmpbsa.in \
  -cp "${PREFIX}.prmtop" -rp "${PREFIX}_receptor.prmtop" -lp "${PREFIX}_ligand.prmtop" \
  -y "${PREFIX}.nc" \
  -o "FINAL_RESULTS_${PREFIX}_MMPBSA.dat" \
  -do "FINAL_RESULTS_${PREFIX}_DECOMP.dat"

# Abort on OVERFLOW/NaN or implausibly large VDW values in complex/receptor/ligand mdout files.
python - "$VDW_LIMIT" <<'PY'
import glob, math, re, sys
limit=float(sys.argv[1]); files=glob.glob('_MMPBSA_*gb.mdout*')
if not files: raise SystemExit('ERROR: MMPBSA mdout dosyaları bulunamadı; VDW kontrolü yapılamadı.')
bad=[]; rows=[]
for fn in sorted(files):
    text=open(fn,errors='replace').read()
    if re.search(r'OVERFLOW|\bNaN\b|\bInf\b',text,re.I): bad.append(f'{fn}: OVERFLOW/NaN/Inf')
    vals=[]
    for m in re.finditer(r'VDWAALS\s*=\s*([-+0-9.Ee]+)',text):
        try: vals.append(float(m.group(1)))
        except ValueError: pass
    if not vals: bad.append(f'{fn}: VDWAALS bulunamadı')
    else:
        rows.append((fn,min(vals),max(vals),max(map(abs,vals))))
        if max(map(abs,vals)) > limit: bad.append(f'{fn}: |VDWAALS| > {limit:g}')
with open('VDW_QC.txt','w') as out:
    out.write('file\tmin_VDWAALS\tmax_VDWAALS\tmax_abs_VDWAALS\n')
    for r in rows: out.write('%s\t%.6f\t%.6f\t%.6f\n'%r)
    out.write('STATUS\t'+('FAIL' if bad else 'PASS')+'\n')
    for x in bad: out.write(x+'\n')
if bad: raise SystemExit('ERROR: VDW kalite kontrolü başarısız. VDW_QC.txt dosyasına bakın.')
print('PASS: Kompleks/reseptör/ligand VDW taşma kontrolü.')
PY

python - "$PREFIX" <<'PY'
import glob, os, sys, zipfile
p=sys.argv[1]; out=f'{p}_MMGBSA_results.zip'
patterns=['FINAL_RESULTS_*','VDW_QC.txt','input_QC.txt','pipeline.log','tleap.log','mmpbsa.in',
          f'{p}.prmtop',f'{p}.inpcrd',f'{p}_receptor.prmtop',f'{p}_ligand.prmtop',f'{p}_tleap.pdb']
files=[]
for pat in patterns: files.extend(glob.glob(pat))
with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED) as z:
    for fn in sorted(set(files)):
        if os.path.isfile(fn): z.write(fn)
print(out)
PY

printf '\nTAMAMLANDI.\nSonuç: %s/FINAL_RESULTS_%s_MMPBSA.dat\nQC:    %s/VDW_QC.txt\nZIP:   %s/%s_MMGBSA_results.zip\n' \
  "$PWD" "$PREFIX" "$PWD" "$PWD" "$PREFIX"
