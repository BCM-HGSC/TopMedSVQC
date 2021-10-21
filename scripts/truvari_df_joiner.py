import os
import sys
import glob
import joblib
import pandas as pd

indir = sys.argv[1]
outname = sys.argv[2]

counts = []
for i in glob.glob(os.path.join(indir, "*")):
    def_fn = i + "/data.jl"
    data = joblib.load(def_fn)
    #data["is_multi"] = data.index.duplicated()
    data = data[~data.index.duplicated(keep="first")] # Only keep one
    view = data.groupby(["state", "szbin", "svtype"]).size()
    view.name = "count"
    view = view.reset_index()
    samp = os.path.basename(i)
    if samp.endswith("_up"):
        view['sample'] = samp[:-len("_up")]
        view['params'] = "maximum"
    else:
        view['sample'] = samp
        view['params'] = "default"
    counts.append(view)

counts = pd.concat(counts)
joblib.dump(counts, outname)



