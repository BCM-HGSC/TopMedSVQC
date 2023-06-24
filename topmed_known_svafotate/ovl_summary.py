import sys
import pandas as pd
import joblib
import numpy as np

out_frame = {}
d = joblib.load(sys.argv[1])
oname = sys.argv[2]
aidx = list(d['a_src'].unique())
bidx = list(d['b_src'].unique())
d['a_idx'] = d['a_src'].apply(lambda x: aidx.index(x))
d['b_idx'] = d['b_src'].apply(lambda x: bidx.index(x))
d['pct'] = d['ab_ovl'] / d['a_cnt']

for svtype in ["DEL", "DUP", "INV"]:
    cnts = np.zeros((d['a_src'].nunique(), d['b_src'].nunique()))
    pcts = np.zeros((d['a_src'].nunique(), d['b_src'].nunique()))

    for _, row in d[d['svtype'] == svtype].iterrows():
        cnts[row['a_idx'], row['b_idx']] = row['ab_ovl']
        pcts[row['a_idx'], row['b_idx']] = row['pct']
    cnts = pd.DataFrame(cnts, index=aidx, columns=bidx)
    pcts = pd.DataFrame(pcts * 100, index=aidx, columns=bidx)
    
    out_frame[svtype] = {'cnts':cnts, 'pcts':pcts}
    print(('='* 4) + svtype + ('=' * 4))
    print("row intersect col counts")
    print(cnts.astype(int))
    print("row intersect col pct of row")
    print(pcts.round(1))

print("Total Counts")
print(d[['a_src', 'svtype', 'a_cnt']].sort_values(['a_src', 'svtype', 'a_cnt']).drop_duplicates())
joblib.dump(out_frame, oname)
