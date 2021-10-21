import sys
import joblib

in_fn = sys.argv[1]
out_fn = sys.argv[2]
d = joblib.load(in_fn)


keep = ['id', 'svtype', 'svlen', 'szbin', 'qual', 'filter', 'is_pass', 'AC', 'NS', 'AF', 'CALLRATE', 'END', 'SVLEN', 'SVTYPE',
        'Biallelic', "GD_AF", "GD_ID", "1000g_event", "GTCNT", "Gene_name"]

d = d[keep]
d["AF"] = d["AF"].apply(lambda x: x[0])
joblib.dump(d, out_fn)
