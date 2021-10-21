import sys
import joblib
import random
import pandas as pd
random.seed(42)
import time
import numpy as np
np.random.seed(42)
import h5py
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')
import pandas
import allel
print('scikit-allel', allel.__version__)

#out_fig = sys.argv[1]
out_data = sys.argv[1]
in_h5s = sys.argv[2:]

full_gns = []
for fn in in_h5s:
    data = h5py.File(fn)
    one = data["pca"]["gn"][:]
    # For the ones filtered by singleton etc
    flt = np.array(data["pca"]["flt"][:])
    # Get the variant alts
    # alts = np.array(data["variants"]["ALT"])[flt]
    afs = np.array(data["variants"]["AF"])[flt]
    passing = np.array(data["variants"]["FILTER_PASS"])[flt]
    # Subset to Dels/DUPs not INV.
    # keep = alts == b'<DUP>'
    keep = (afs >= 0.05) & (passing)
    # PCA on those items not in LD and of svtype
    #new_one = one.compress(data["pca"]["loc_unlinked"], axis=0)
    new_one = one.compress(data["pca"]["loc_unlinked"][keep], axis=0)
    print('loaded', fn, new_one.shape)
    full_gns.append(new_one)

gnu = np.concatenate(full_gns)
print('final size', gnu.shape)
# Randomized PCA should be faster...
# I need to figure out how to subset by SVTYPE
coords5, model5 = allel.randomized_pca(gnu, n_components=10, scaler='patterson')
joblib.dump((coords5, model5), out_data)




