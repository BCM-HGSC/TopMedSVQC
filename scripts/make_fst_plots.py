
import sys
import joblib
import numpy as np
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt

data = joblib.load(sys.argv[1])
out_name = sys.argv[2]

ancestries = set()
for key in data:
    key = key.split('-')
    ancestries.add(key[0])
    ancestries.add(key[1])

ancestries = list(ancestries)
print(ancestries)
full_arr = np.zeros((len(ancestries), len(ancestries)))

for join_key in data:
    key = join_key.split('-')
    idx1 = ancestries.index(key[0])
    idx2 = ancestries.index(key[1])
    full_arr[idx1, idx2] = data[join_key]["fst"]
    full_arr[idx2, idx1] = data[join_key]["fst"]
data = pd.DataFrame(full_arr, index=ancestries, columns=ancestries)
p = sb.clustermap(data.fillna(0))
#p.legend(title="Fst")
#p.ax_heatmap.legend(title="Fst")
#p.fig.suptitle("Pairwise Fst")
plt.savefig(out_name, bbox_inches='tight')
