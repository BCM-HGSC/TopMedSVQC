import sys
import pandas as pd

in_file, out_file = sys.argv[1:]
data = pd.read_csv(in_file, sep='\t', header=None)
data.columns = ["chrom", "start", "end", "id", "svcount"]


desc = data["svcount"].describe()


print(desc)
hs_threshold = desc["mean"] + (3* desc["std"])
print(f"Setting HotSpot threshold at {hs_threshold:.2f}")

data["anno"] = None
data.loc[data["svcount"] == 0, "anno"] = "Des"
data.loc[data["svcount"] > hs_threshold, "anno"] = "Hot"

print("Counts")
print(data["anno"].value_counts())

data.to_csv(out_file, sep='\t', index=False)
