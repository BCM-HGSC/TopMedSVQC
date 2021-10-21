import sys
import pysam
import allel
import joblib
import numpy as np
np.random.seed(40)
import pandas as pd
from collections import defaultdict
import glob
import truvari
import itertools
import h5py
import random
random.seed(40)

ancestry_fn = sys.argv[1] # indexed_ancestries.jl
output_fn = sys.argv[2] # fst_data.jl
in_h5s = sys.argv[3:] # *.h5 - can reuse from the pca's vcf_to_df

ances = joblib.load(ancestry_fn)

lookup = {}
prog = 0

all_cd = []
for fn in in_h5s:
    data = h5py.File(fn)
    # Filter out singletons... wait.. I don't 
    flt = np.array(data["pca"]["flt"][:])
    # Get the variant alts
    alts = np.array(data["variants"]["ALT"])
    # Subset to Dels/DUPs not INV.
    keep = alts == b'<DEL>'
    # I don't need to do 'one'?
    print(data['calldata/GT'].shape)
    print(keep.shape)
    print(flt.shape)
    view = data["calldata/GT"][:][keep & flt]
    print(view.shape)

    all_cd.append(view)

# Concat and turn into genotypeArray
garr = allel.GenotypeArray(np.concatenate(all_cd))
print(garr.shape)
#For each of the ancestries, I need to run count_alleles...
# if subpop >= something, I'll want to random.sample it
labs = []
ac = []
for anc, dat in ances.groupby(["Ancestry"]):
    subpop = dat["index"].values
    if len(subpop) < 1000:
        continue
    if len(subpop) > 7000:
        subpop = np.random.choice(dat["index"].values, 7000)
        
    labs.append(anc)
    ac.append(garr.count_alleles(subpop=subpop))

#Then do the fst
to_save = {}
for i, j in itertools.combinations(list(range(len(labs))), 2):
    num, den = allel.hudson_fst(ac[i], ac[j])
    num = num[~np.isnan(num)]
    den = den[~np.isnan(den)]
    fst = np.sum(num) / np.sum(den)
    key = "%s-%s" % (labs[i], labs[j])
    to_save[key] = {}
    to_save[key]["num"] = num
    to_save[key]["den"] = den
    to_save[key]["fst"] = fst
    print(labs[i], labs[j], fst)

joblib.dump(to_save, output_fn)

"""
for part, chrom in enumerate(parts):
    chrom_name= chrom.split('/')[-1].split('.')[2]
    data = joblib.load(chrom)
    view = data["variants/ALT"][:, 0] != '__<DEL>'
    max_AF = data["variants/AF"][:,0] !=  -1210
    garr = allel.GenotypeArray(data["calldata/GT"][view & max_AF])
    prog += 1
    #bar.update(prog)
    ac = []
    labs = []
    for anc, dat in ances.groupby(["Ancestry"]):
        labs.append(anc)
        subpop = dat["index"].values
        ac.append(garr.count_alleles(subpop=subpop))
        continue
        ac = garr.count_alleles(subpop=subpop)
        if anc not in lookup:
            lookup[anc] = ac.values
        else:
            lookup[anc] = np.concatenate([lookup[anc], ac.values])
    # This is doing the PCA.
    for i, j in itertools.combinations(list(range(len(labs))), 2):
        print(labs[i], labs[j])
        num, den = allel.hudson_fst(ac[i], ac[j])
        #num, den = allel.patterson_fst(ac[i], ac[j])
        num = num[~np.isnan(num)]
        den = den[~np.isnan(den)]
        fst = np.sum(num) / np.sum(den)
        print(fst)
    
    prog += 1
    #bar.update(prog)

#bar.finish()
joblib.dump(lookup, "allele_counts.jl")
"""
