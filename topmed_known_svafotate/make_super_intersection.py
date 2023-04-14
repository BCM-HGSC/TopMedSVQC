import pandas as pd
import truvari
import joblib

data = pd.read_csv("new.svafotate.bed.gz", sep='\t')
DOSINGLE=False

sources = list(data["SOURCE"].unique())
svtypes = ["DEL", "DUP", "INV"]

if DOSINGLE:
    data = data[(data['Het'] + data['HomAlt']) == 1]
rows = []
for asrc in sources:
    for bsrc in sources:
        if asrc == bsrc:
            continue
        for svt in svtypes:
            print(asrc, bsrc, svt)
            adat = data[(data["SOURCE"] == asrc) & (data["SVTYPE"] == svt)]
            bdat = data[(data["SOURCE"] == bsrc) & (data["SVTYPE"] == svt)]
            adat[["#CHROM", "START", "END"]].to_csv("a.bed", sep='\t', index=False, header=False)
            bdat[["#CHROM", "START", "END"]].to_csv("b.bed", sep='\t', index=False, header=False)
            ret = truvari.cmd_exe(f"bedtools intersect -r -f 0.8 -u -a a.bed -b b.bed | wc -l", pipefail=True)
            rows.append([asrc, bsrc, svt, len(adat), len(bdat), int(ret.stdout.strip())])

for asrc in sources:
    for svt in svtypes:
        print(asrc, 'allother', svt)
        adat = data[(data["SOURCE"] == asrc) & (data["SVTYPE"] == svt)]
        bdat = data[(data["SOURCE"] != asrc) & (data["SVTYPE"] == svt)]
        adat[["#CHROM", "START", "END"]].to_csv("a.bed", sep='\t', index=False, header=False)
        bdat[["#CHROM", "START", "END"]].to_csv("b.bed", sep='\t', index=False, header=False)
        ret = truvari.cmd_exe(f"bedtools intersect -r -f 0.8 -u -a a.bed -b b.bed | wc -l", pipefail=True)
        rows.append([asrc, 'allother', svt, len(adat), len(bdat), int(ret.stdout.strip())])

rows = pd.DataFrame(rows, columns=["a_src", "b_src", "svtype", "a_cnt", "b_cnt", "ab_ovl"])
if DOSINGLE:
    joblib.dump(rows, "biginter_single.jl")
else:
    joblib.dump(rows, "biginter.jl")


