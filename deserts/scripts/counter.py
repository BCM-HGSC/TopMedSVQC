import sys
from intervaltree import IntervalTree
from collections import Counter
import pandas as pd
import pysam
import truvari

genome = "/users/u233287/scratch/topmed_analysis/deserts/grch38.windows.bed.gz"
tree, cnt = truvari.build_anno_tree(genome)
out_fn = sys.argv[1]
counts = Counter()
for i in sys.argv[2:]:
    v = pysam.VariantFile(i)
    cnt = 0
    for entry in v:
        #if 'PASS' not in entry.filter:
        #    continue
        if 'INS' not in entry.info["SVTYPE"]:
            continue
        for intv in tree[entry.chrom].overlap(entry.start, entry.stop):
            cnt += 1
            counts[intv.data] += 1

data = []
keep_chrs = set()
for chrom in tree:
    for intv in tree[chrom]:
        cnt = counts[intv.data]
        data.append([chrom, intv.begin, intv.end, cnt])
        if cnt:
            keep_chrs.add(chrom)

data = pd.DataFrame(data, columns=['chrom', 'start', 'end', 'count'])
#data = data[data["chrom"].isin(list(keep_chrs))]
# Do the hotspot work
desc = data["count"].describe()
hs_threshold = desc["mean"] + (3 * desc["std"])

data["anno"] = None
data.loc[data["count"] == 0, "anno"] = "sparse"
data.loc[data["count"] > hs_threshold, "anno"] = "dense"
data.to_csv(out_fn, sep='\t', index=False, header=False)
