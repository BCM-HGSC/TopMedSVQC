# Simple Karyoplot
import gzip
import joblib
import numpy as np
import seaborn as sb
import pandas as pd
import matplotlib.pyplot as plt
# I need a genome file...
genome = {}
with open("/users/u233287/scratch/insertion_ref/reference/grch38/GRCh38_1kg_mainchrs.fa.fai") as fh:
    for line in fh:
        data = line.strip().split('\t')
        genome[data[0]] = int(data[1])

lookup = {}
with gzip.open("topmed.sv.bed.gz") as fh:
    for line in fh:
        line = line.decode()
        data = line.strip().split('\t')
        if data[0] not in lookup:
            lookup[data[0]] = {}
        if data[3] not in lookup[data[0]]:
            lookup[data[0]][data[3]] = np.zeros(genome[data[0]])
        lookup[data[0]][data[3]][int(data[1]): int(data[2])] += 1

new_den = []
BIN = int(5e5)
for chrom in lookup:
    for svtype in lookup[chrom]:
        for i in range(0, len(lookup[chrom][svtype]), BIN):
            new_den.append([chrom, svtype, i, lookup[chrom][svtype][i:i+BIN].sum() / BIN])

new_den = pd.DataFrame(new_den, columns=["chrom", "svtype", "Pos", "Density"])
joblib.dump(new_den, "karyo.jl")
p = sb.FacetGrid(data=new_den, row="chrom", col="svtype")
p.map(sb.lineplot, "Pos", "Density")
plt.savefig("karyoplot.png")

