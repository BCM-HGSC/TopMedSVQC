import os
import sys
import glob
import pysam
import truvari
import pandas as pd
from collections import Counter, defaultdict

in_dir = sys.argv[1]
manifest_fn = sys.argv[2]
tt_mars_fn = sys.argv[3]
output_fn = sys.argv[4]

def make_consensus_gt(entry, samples):
    """
    Collecitons counter
    Also samples needs to be a list, so... hmm
    """
    cnt = Counter()
    for samp in samples:
        cnt[truvari.get_gt(entry.samples[samp]["GT"]).name] += 1
    most_common = cnt.most_common()[0][0]
    pre_cnt = cnt["HET"] + cnt["HOM"]
    con_cnt = cnt[most_common]
    return most_common, pre_cnt, con_cnt
        

# load manifest
manifest = pd.read_csv(manifest_fn, sep='\t')

subject_samples = {}
sample_to_subject = {}
for i, j in manifest.groupby(['study_subject_id']):
    subject_samples[i] = j["sample.id"].values
    for k in j["sample.id"].values:
        sample_to_subject[k] = i

info_order = []
info_order.append(("NA12878_BenchState", 
                   ('##INFO=<ID=NA12878_BenchState,Number=1,Type=String,'
                    'Description="NA12878 TP, FP, or .">')))
info_order.append(("NA12878_ConsensusGT", 
                   ('##INFO=<ID=NA12878_ConsensusGT,Number=1,Type=String,'
                    'Description="NA12878 most common genotype">')))
info_order.append(("NA19238_BenchState", 
                   ('##INFO=<ID=NA19238_BenchState,Number=1,Type=String,'
                    'Description="NA19238 TP, FP, or .">')))
info_order.append(("NA12878_ConsensusGT", 
                   ('##INFO=<ID=NA12878_ConsensusGT,Number=1,Type=String,'
                    'Description="NA12878 most common genotype">')))
info_order.append(("NA19238_ConsensusGT", 
                   ('##INFO=<ID=NA19238_ConsensusGT,Number=1,Type=String,'
                    'Description="NA19238 most common genotype">')))
info_order.append(("NA12878_PresentCount", 
                   ('##INFO=<ID=NA12878_PresentCount,Number=1,Type=Integer,'
                    'Description="NA12878 present replicates count">')))
info_order.append(("NA19238_PresentCount",
                   ('##INFO=<ID=NA19238_PresentCount,Number=1,Type=Integer,'
                    'Description="NA19238 present replicates count">')))
info_order.append(("NA12878_ConcordantCount",
                   ('##INFO=<ID=NA12878_ConcordantCount,Number=1,Type=Integer,'
                    'Description="NA12878 gt concordant replicates count">')))
info_order.append(("NA19238_ConcordantCount",
                   ('##INFO=<ID=NA19238_ConcordantCount,Number=1,Type=Integer,'
                    'Description="NA19238 gt concordant replicates count">')))
info_order.append(("TTMars_NA12878", 
                  ('##INFO=<ID=TTMars_NA12878,Number=1,Type=String,'
                   'Description="TT-Mars NA12878 results: TP, FP, NA or .">')))
info_order.append(("TTMars_NA19238",
                  ('##INFO=<ID=TTMars_NA19238,Number=1,Type=String,'
                   'Description="TT-Mars NA19238 results: TP, FP, NA or .">')))

tt_lookup = {}
v = pysam.VariantFile(tt_mars_fn)
for entry in v:
    tt_lookup[entry.id] = {"TTMars_NA12878": entry.info["TTMars_NA12878"],
                           "TTMars_NA19238": entry.info["TTMars_NA19238"]}

var_lookup = {}
big_header = None
for cur_dir in glob.glob(os.path.join(in_dir, "NWD*")):
    cur_sample = sample_to_subject[os.path.basename(cur_dir)]
    
    tps = pysam.VariantFile(os.path.join(cur_dir, 'tp-call.vcf.gz'))                                 
    if big_header is None:
        big_header = tps.header.copy()

    for entry in tps:
        var_id = entry.id
        # Adding a new TP
        if var_id not in var_lookup:
            var_lookup[var_id] = {}
            for i in info_order:
                var_lookup[var_id][i[0]] = None
            var_lookup[var_id]['variant'] = entry
            var_lookup[var_id]['score'] = entry.info["TruScore"]

        #If it is in lookup, but it's a FP, I need to overwrite the entry so I get the 'match score' stuff
        if var_lookup[var_id][cur_sample + '_BenchState'] != 'TP':
            var_lookup[var_id]['variant'] = entry
            var_lookup[var_id]['score'] = entry.info["TruScore"]
            

        # That match isn't the highest scoring occurance of the variant
        if var_lookup[var_id]['score'] < entry.info["TruScore"]:
            var_lookup[var_id]['variant'] = entry
            var_lookup[var_id]['score'] = entry.info["TruScore"]

        # Overwrite '.' or 'FP' to TP
        var_lookup[var_id][cur_sample + '_BenchState'] = 'TP'
        # probably redundant, but whatever
        con_gt, pre_cnt, con_cnt = make_consensus_gt(entry, subject_samples[cur_sample])
        var_lookup[var_id][cur_sample + '_ConsensusGT'] = con_gt
        var_lookup[var_id][cur_sample + '_PresentCount'] = pre_cnt
        var_lookup[var_id][cur_sample + '_ConcordantCount'] = con_cnt
        for key in tt_lookup[var_id]:
            var_lookup[var_id][key] = tt_lookup[var_id][key]

    fps = pysam.VariantFile(os.path.join(cur_dir, 'fp.vcf.gz'))
    for entry in fps:
        var_id = entry.id
        if var_id not in var_lookup:
            var_lookup[var_id] = {}
            for i in info_order:
                var_lookup[var_id][i[0]] = None
            var_lookup[var_id]["variant"] = entry

        if var_lookup[var_id][cur_sample + '_BenchState'] is None:
            # Can only overwrite '.' to FP
            var_lookup[var_id][cur_sample + '_BenchState'] = "FP"  

        # Just populate every time
        con_gt, pre_cnt, con_cnt = make_consensus_gt(entry, subject_samples[cur_sample])
        var_lookup[var_id][cur_sample + '_ConsensusGT'] = con_gt
        var_lookup[var_id][cur_sample + '_PresentCount'] = pre_cnt
        var_lookup[var_id][cur_sample + '_ConcordantCount'] = con_cnt
        for key in tt_lookup[var_id]:
            var_lookup[var_id][key] = tt_lookup[var_id][key]
## END
for i in info_order:
    big_header.add_line(i[1])

out = pysam.VariantFile(output_fn, 'w', header=big_header)
for key in var_lookup:
    m_dict = dict(var_lookup[key])
    entry = truvari.copy_entry(m_dict["variant"], big_header)
    for info in info_order:
        if m_dict[info[0]]:
            entry.info[info[0]] = m_dict[info[0]]
    out.write(entry)


