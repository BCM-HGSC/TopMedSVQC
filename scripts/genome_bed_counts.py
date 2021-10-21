import sys
import pysam
import numpy as np


INFILE="/stornext/snfs4/next-gen/scratch/english/round2/topmed_analysis/annosv/results/topmed.results.vcf.gz"
v = pysam.VariantFile(INFILE)
BINSIZE = 10000
overlap = {}
starting = {}
length_lookup = {}
for chrom in v.header.contigs.values():
    length_lookup[chrom.name] = chrom.length
    length = chrom.length // BINSIZE
    overlap[chrom.name] = np.zeros(length, dtype=np.int)
    starting[chrom.name] = np.zeros(length, dtype=np.int)


for entry in v:
    overlap[entry.chrom][entry.pos // BINSIZE : entry.stop // BINSIZE] += 1
    starting[entry.chrom][entry.pos // BINSIZE] += 1


for chrom in overlap:
    num_bins = len(overlap[chrom])
    for pos in range(num_bins):
        start = pos * BINSIZE
        end = start + BINSIZE if pos + 1 != num_bins else length_lookup[chrom]
        print(f"{chrom}\t{start}\t{end}\t{starting[chrom][pos]}\t{overlap[chrom][pos]}")
