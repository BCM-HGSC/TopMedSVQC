import sys
import random
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
import pysam
print('scikit-allel', allel.__version__)


in_name = sys.argv[1]
out_name = sys.argv[2]
fields = ["variants/CHROM", "variants/FILTER", "variants/POS", "variants/AF", "variants/ALT", "calldata/GT"]

# Filtering the weird sample for DUPS
v = pysam.VariantFile(in_name)
#samps = [_ for _ in v.header.samples if _ not in ["NWD976804"]]
allel.vcf_to_hdf5(in_name, out_name, fields=fields, alt_number=1)#, samples=samps)


callset = h5py.File(out_name, 'a')
g = allel.GenotypeChunkedArray(callset["calldata/GT"])
callset.create_group("pca")
ac = g.count_alleles()[:]
# DUP
#ac = ac[callset["variants/FILTER_PASS"]]

print("Num multiallelics", np.count_nonzero(ac.max_allele() > 1))
print("Bi-allelic singletons", np.count_nonzero((ac.max_allele() == 1) & ac.is_singleton(1)))

# Filtering
flt = (ac.max_allele() == 1) & (ac[:, :2].min(axis=1) > 1)
callset["pca"].create_dataset('flt', data=flt)
gf = g.compress(flt, axis=0)
# Input for the PCA
gn = gf.to_n_alt()

#chr_name = sys.argv[1].split('/')[-1].split('.')[2]

callset["pca"].create_dataset('gn', data=gn)

# LD prune
#def plot_ld(gn, title, out_name="test.png"):
#    m = allel.rogers_huff_r(gn) ** 2
#    ax = allel.plot_pairwise_ld(m)
#    ax.set_title(title)
#    plt.savefig(out_name)


def ld_prune(gn, size, step, threshold=.1, n_iter=1):
    for i in range(n_iter):
        loc_unlinked = allel.locate_unlinked(gn, size=size, step=step, threshold=threshold)
        n = np.count_nonzero(loc_unlinked)
        n_remove = gn.shape[0] - n
        print('iteration', i+1, 'retaining', n, 'removing', n_remove, 'variants')
    return loc_unlinked


#plot_ld(gn, 'Figure 1. Pairwise LD', "pca." + chr_name + '.preld.png')

loc_unlinked = ld_prune(gn, size=100, step=20, threshold=.1, n_iter=1)
# Make the ld plots (before/after)

#plot_ld(gnu, 'Figure 2. Pairwise LD after pruning', "pca." + chr_name + ".postld.png")
#callset["pca"].create_dataset('gn_prune', data=gnu.values)
callset["pca"].create_dataset('loc_unlinked', data=loc_unlinked)

