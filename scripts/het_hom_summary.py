import sys
import h5py
import joblib
import allel
import pandas as pd
out_fn = sys.argv[1]
in_h5 = sys.argv[2:]
rows = []
for i in in_h5:
    h5 = h5py.File(i)
    meta = i.split('.')
    chrom = meta[-5]
    svtype = meta[-4]
    a = allel.GenotypeArray(h5["calldata"]["GT"][:])
    het = pd.Series(a.count_het(axis=0))
    het["chrom"] = chrom
    het["svtype"] = svtype
    het["gt"] = 'het'
    rows.append(het)

    hom = pd.Series(a.count_hom_alt(axis=0))
    hom["chrom"] = chrom
    hom["svtype"] = svtype
    hom["gt"] = 'hom'
    rows.append(hom)

x = pd.DataFrame(rows)
#x[x["gt"] == 'het'].sum().iloc[:-3] / x[x["gt"] == 'hom'].sum().iloc[:-3]
joblib.dump(x, out_fn)
    
    
