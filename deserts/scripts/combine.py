import joblib
import pandas as pd


snps = pd.read_csv("results/snps.bed", header=None, sep='\t')
snps.columns = ["chrom", "start", "end", "snp_count", "snp_density"]
snps['snp_density'] = snps['snp_density'].where(snps['snp_density'] != "None", None)
snps.set_index(['chrom', 'start', 'end'], inplace=True)

hg_del = pd.read_csv("results/hgsvc.DEL.bed", header=None, sep='\t')
hg_del.columns = ["chrom", "start", "end", "hg_del_count", "hg_del_density"]
hg_del.set_index(['chrom', 'start', 'end'], inplace=True)

hg_ins = pd.read_csv("results/hgsvc.INS.bed", header=None, sep='\t')
hg_ins.columns = ["chrom", "start", "end", "hg_ins_count", "hg_ins_density"]
hg_ins.set_index(['chrom', 'start', 'end'], inplace=True)

tp_del = pd.read_csv("results/topmed.DEL.bed", header=None, sep='\t')
tp_del.columns = ["chrom", "start", "end", "tp_del_count", "tp_del_density"]
tp_del.set_index(['chrom', 'start', 'end'], inplace=True)

tp_dup = pd.read_csv("results/topmed.DUP.bed", header=None, sep='\t')
tp_dup.columns = ["chrom", "start", "end", "tp_dup_count", "tp_dup_density"]
tp_dup.set_index(['chrom', 'start', 'end'], inplace=True)

tp_inv = pd.read_csv("results/topmed.INV.bed", header=None, sep='\t')
tp_inv.columns = ["chrom", "start", "end", "tp_inv_count", "tp_inv_density"]
tp_inv.set_index(['chrom', 'start', 'end'], inplace=True)

data = snps.join(hg_del, how='outer').join(hg_ins, how='outer').join(tp_del, how='outer').join(tp_dup, how='outer').join(tp_inv, how='outer')

# First time around
#data.reset_index()[["chrom", "start", "end"]].to_csv("full_intersection.bed", index=False, sep='\t', header=False)
#print(data)
# Second time around
annos = pd.read_csv("annosv/full_intersection.annotated.tsv", sep='\t')
annos.set_index(["SV chrom", "SV start", "SV end"], inplace=True)
annos = annos[annos["AnnotSV type"] == 'full']["Gene name"].to_frame()
annos.reset_index(inplace=True)
annos.columns = ["chrom", "start", "end", "gene"]
annos['chrom'] = annos['chrom'].apply(lambda x: 'chr' + str(x))
annos.set_index(["chrom", "start", "end"], inplace=True)
full = data.join(annos, how='left')
print(full)
joblib.dump(full, "full_intersection.jl")
