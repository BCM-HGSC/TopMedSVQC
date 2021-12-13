import sys
import h5py
import joblib

lookup = {}
for i in sys.argv[1:]:
    if "DEL" in i:
        key = "DEL"
    elif "INV" in i:
        key = "INV"
    elif "DUP" in i:
        key = "DUP"
    d = h5py.File(i)
    gts = d["calldata/GT"][:]
    cnts = (gts == 1).any(axis=2).sum(axis=0)
    if key not in lookup:
        lookup[key] = cnts
    else:
        lookup[key] += cnts
joblib.dump(lookup, "counts.jl")
