import sys
import gzip
import pysam
import truvari
from collections import defaultdict
import pandas as pd

id_fn = sys.argv[1]
in_fn = sys.argv[2]
#out is stdout

"""
ids = defaultdict(bool)
with gzip.open(id_fn) as fh:
    for line in fh:
        line = line.decode()
        if line.startswith("#"):
            continue
        dat = line.strip().split('\t')
        ids[dat[2]] = True
"""
ids = pd.read_csv(sys.argv[1])
ids.set_index("ID", inplace=True)
"""
ids = defaultdict(int)
with gzip.open(id_fn) as fh:
    fh.readline()
    for line in fh:
        line = line.decode()
        dat = line.strip().split('\t')
        score = min(100, int(float(dat[1]) * 100))
        ids[dat[0]] = score
"""

def other_filters(entry):
    """
    """
    ret = []
    unk, ref, het, hom = entry.info["GTCNT"]
    missing = unk / (unk + ref + het + hom)
    if missing < 0.15:
        ret.append("Missing")
    if entry.qual < 0.98:
        ret.append("LowQual")
    return ret
    #consol_dup["MissingPct"] = consol_dup["GTCNT"].apply(calc_missingpct)
    #sb.displot(data=consol_dup, x="MissingPct", hue="OnlyTP", multiple="dodge", binwidth=0.01)

fh = pysam.VariantFile(in_fn)
n_header = fh.header
n_header.add_line(('##FILTER=<ID=Collapse,'
                   'Description="Call collapsed due to high overlap and genotype concordance with other call">'))
n_header.add_line(('##FILTER=<ID=LowQual,'
                   'Description="Assessment of normalized depth indicates low variant quality">'))

out = pysam.VariantFile('/dev/stdout', 'w', header=n_header)
for entry in fh:
    if entry.alts[0] != '<DUP>': continue
    #entry = truvari.copy_entry(entry, n_header)
    if entry.id not in ids.index:
        entry.filter.add("Collapse")
    else:
        n_qual = ids.loc[entry.id]
        entry.qual = n_qual["qual"]
        if not n_qual["pass"]:
            entry.filter.add("LowQual")
    out.write(entry)
