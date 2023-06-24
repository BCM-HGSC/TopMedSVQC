import sys
import pandas as pd
import truvari
import joblib

data = pd.read_csv("new.svafotate.bed.gz", sep='\t')
data['#CHROM'] = data['#CHROM'].astype(str)
keep_chrs = [str(_) for _ in range(1, 23)]
data = data[data["#CHROM"].isin(keep_chrs)].copy()
DOSINGLE, outname = sys.argv[1:]
DOSINGLE = int(DOSINGLE)
# 0 - 'all' 
# 1 - 'singletons'
# 2 - 'non-singletons'

sources = list(data["SOURCE"].unique())
svtypes = ["DEL", "DUP", "INV"]

if DOSINGLE == 1:
    a_data = data[(data['Het'] + data['HomAlt']) == 1].copy()
elif DOSINGLE == 2:
    a_data = data[(data['Het'] + data['HomAlt']) != 1].copy()
else:
    a_data = data

rows = []
for asrc in sources:
    for bsrc in sources:
        if asrc == bsrc:
            continue
        for svt in svtypes:
            print(asrc, bsrc, svt)
            adat = a_data[(a_data["SOURCE"] == asrc) & (a_data["SVTYPE"] == svt)]
            bdat = data[(data["SOURCE"] == bsrc) & (data["SVTYPE"] == svt)]
            adat[["#CHROM", "START", "END"]].to_csv("a.bed", sep='\t', index=False, header=False)
            bdat[["#CHROM", "START", "END"]].to_csv("b.bed", sep='\t', index=False, header=False)
            ret = truvari.cmd_exe(f"bedtools intersect -r -f 0.8 -u -a a.bed -b b.bed | wc -l", pipefail=True)
            rows.append([asrc, bsrc, svt, len(adat), len(bdat), int(ret.stdout.strip())])

for asrc in sources:
    for svt in svtypes:
        print(asrc, 'allother', svt)
        adat = a_data[(a_data["SOURCE"] == asrc) & (a_data["SVTYPE"] == svt)]
        bdat = data[(data["SOURCE"] != asrc) & (data["SVTYPE"] == svt)]
        adat[["#CHROM", "START", "END"]].to_csv("a.bed", sep='\t', index=False, header=False)
        bdat[["#CHROM", "START", "END"]].to_csv("b.bed", sep='\t', index=False, header=False)
        ret = truvari.cmd_exe(f"bedtools intersect -r -f 0.8 -u -a a.bed -b b.bed | wc -l", pipefail=True)
        rows.append([asrc, 'allother', svt, len(adat), len(bdat), int(ret.stdout.strip())])

rows = pd.DataFrame(rows, columns=["a_src", "b_src", "svtype", "a_cnt", "b_cnt", "ab_ovl"])
joblib.dump(rows, outname)
