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
anc_fn = sys.argv[1]
dat_fn = sys.argv[2]
out_fig = sys.argv[3]
ancestries = pd.read_csv(anc_fn)
#ancestries = ancestries[~ancestries["NWDID"].isin(["NWD714003",
                                                   #"NWD170197",
                                                   #"NWD910621",
                                                   #"NWD320728",
                                                   #"NWD592708",
                                                   #"NWD859633",
                                                   #"NWD786098",
                                                   #"NWD480965",
                                                   #"NWD369359",
                                                   #"NWD719720",
                                                   #"NWD884712"
                                                  #])
                       #]
sample_population = ancestries["Ancestry"]
populations = sample_population.unique()

coords5, model5 = joblib.load(dat_fn)

d = pd.DataFrame(coords5, index=sample_population, columns=["PC%d" % x for x in range(1, 11)])
d = d.round(1)
d.to_csv("temp.csv")
