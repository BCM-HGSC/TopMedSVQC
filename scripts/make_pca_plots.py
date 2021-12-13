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
# Dropped samples NWD714003
ancestries = ancestries[~ancestries["NWDID"].isin(["NWD714003",
                                                   "NWD170197",
                                                   "NWD910621",
                                                   "NWD320728",
                                                   "NWD592708",
                                                   "NWD859633",
                                                   "NWD786098",
                                                   "NWD480965",
                                                   "NWD369359",
                                                   "NWD719720",
                                                   "NWD884712"
                                                  ])
                       ]
sample_population = ancestries["Ancestry"]
populations = sample_population.unique()

pop_colours = {
    'AFR': '#FF0000',
    'AMR': '#008000',
    'EAS': '#00FFFF',
    'EUR': '#90EE90',
    'MES': '#FFA500',
    'OCN': '#8B0000',
    'SAS': '#1E90FF',
    'KES': '#808080',
    'CMS': '#0000FF',
}


coords5, model5 = joblib.load(dat_fn)

def plot_pca_coords(coords, model, pc1, pc2, ax, sample_population):
    sns.despine(ax=ax, offset=5)
    x = coords[:, pc1]
    y = coords[:, pc2]
    for pop in populations:
        flt = (sample_population == pop)
        ax.plot(x[flt], y[flt], marker='o', linestyle=' ', color=pop_colours[pop], 
                label=pop, markersize=6, mec='k', mew=.5)
    ax.set_xlabel('PC%s (%.1f%%)' % (pc1+1, model.explained_variance_ratio_[pc1]*100))
    ax.set_ylabel('PC%s (%.1f%%)' % (pc2+1, model.explained_variance_ratio_[pc2]*100))

def fig_pca(coords, model, title, sample_population=None, out_name="test.pca.png"):
    if sample_population is None:
        sample_population = df_samples.population.values
    # plot coords for PCs 1 vs 2, 3 vs 4
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(1, 2, 1)
    plot_pca_coords(coords, model, 0, 1, ax, sample_population)
    #plt.xlim((-1, 10))
    #plt.ylim((-1, 100))
    ax = fig.add_subplot(1, 2, 2)
    plot_pca_coords(coords, model, 2, 3, ax, sample_population)
    #plt.xlim((-1, 100))
    ax.legend(bbox_to_anchor=(1, 1), loc='upper left')
    fig.suptitle(title, y=1.02)
    fig.tight_layout()
    plt.savefig(out_name, bbox_inches='tight')

fig_pca(coords5, model5, 'Randomized PCA.', sample_population, out_fig)
