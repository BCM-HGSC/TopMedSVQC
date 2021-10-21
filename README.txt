SEPARATE BY SVTYPE
==================
There should be separate DEL, DUP, and INV files. DUP files are a special case and will be added later.
To start, the pipeline will focus on DEL/INV

bash scripts/mk_separate_by_type_cmds.sh make the jobs and then submit them

EXTRACT NEW DUPS
================

bash scripts/mk_subsetted_dups_cmds.sh 

ESTIMATING PPV OF DATA
======================
I need the 'bench' samples vcfs

1) I'll need to extract just the replicate samples so I have an easy VCF to work with
	scripts/mk_extract_bench_samples_cmds.sh make the jobs and then submit them
	I'm running on a per-chrom bases so I need to consolidate manually with vcf-concat and create
	the validation/topmed.replicate.DEL.vcf.gz
3) The script to does the validation..
	bash scripts/mk_truvari_cmds.sh make the jobs and then submit them
	NOTE! 
	the COMPVCF and OUTDIR variables need to be manually set inside of the script
4) I then need to make the dataframe and the FP/TP consolidation of the results
	Make the summary dataframe for analysis
	python scripts/truvari_df_joiner.py validation/results/deletions/ validation/results/deletions/summary.jl

5) Make the consolidation, annotated TP/FPs from the 'defaults'
	python scripts/validation_consolidator.py validation/results/deletions/ \
			metadata/validation_replicate.manifest.txt \
			validation/results/annotated.deletions.vcf

6) And then I need to clean up the notebook so I can run one simple thing on the summary.jl and consolidated data frame

DOING THE PCA
=====================
Note - this doesn't work for the DUPs with filter values. So I manually concat and filter down to only DUPs.
Then skip to step #2

1) bash scripts/mk_vcf2pca_cmds.sh make the jobs and then submit them

This will setup pca_data/ directory with the per svtype files

2) Run pca pipeline on the h5 files

	python scripts/make_pca_data.py pca_data/deletions.jl pca_data/*DEL*h5
	python scripts/make_pca_plots.py metadata/ancestries.txt pca_data/inversions.jl pca_data/inversions.png

PCA FILTERING
=============
I think there might be samples I need to remove because they're looking weird in the PCAs?
PC2 <= -20
    2505    NWD895484
    2529    NWD593511
    2632    NWD134290
    3006    NWD594428
    3008    NWD989091
PC4 >= 500
    162     NWD207905
    179     NWD112322
    2505    NWD895484
    2529    NWD593511
    2632    NWD134290
    3006    NWD594428
    3008    NWD989091
Should look at the QC metrics on these...
For DUPs
>>> ancestries[pc1 > 150]
       NWDID          Ancestry
48220  NWD976804      AFR

DOING THE FST
=============
python scripts/make_fst_data.py metadata/indexed_ancestries.jl fst_data/inversions.jl pca_data/*INV*h5
python scripts/make_fst_plots.py fst_data/inversions.jl fst_data/inversions.png
	Note. Inversions fail...
	Traceback (most recent call last):
	  File "scripts/make_fst_data.py", line 58, in <module>
	    ac.append(garr.count_alleles(subpop=subpop))
	  File "/stornext/snfs5/next-gen/scratch/english/round2/pyenv/lib/python3.8/site-packages/allel/model/ndarray.py", line
	1835, in count_alleles
	    max_allele = self.max()
	  File "/stornext/snfs5/next-gen/scratch/english/round2/pyenv/lib/python3.8/site-packages/numpy/core/_methods.py", line
	40, in _amax
	    return umr_maximum(a, axis, None, out, keepdims, initial, where)
	ValueError: zero-size array to reduction operation maximum which has no identity
DELs:
    AFR AMR 0.08346761091158061
    AFR EAS 0.10828191976566769
    AFR EUR 0.0731332436710205
    AFR SAS 0.07136415180653874
    AMR EAS 0.0570303261965993
    AMR EUR 0.03520831082656243
    AMR SAS 0.03736846401732368
    EAS EUR 0.08069643282330462
    EAS SAS 0.05827352777745373
    EUR SAS 0.019150763825462696

ANNOTATING VARIANTS
===================
scripts/mk_annosv_cmds.sh  awef.vcf.gz aoba.vcf.gz to make the jobs and then submit them

Populates annosv/ directories

Build the VCFs
python scripts/put_annosv_in_vcf.py annosv/topmed.INV.vcf.gz/tmp_topmed.INV.vcf.gz.annotated.tsv call_only_vcfs/topmed.INV.vcf.gz \
	| bgzip > annosv/topmed.INV.anno.vcf.gz

Which I then can vcf2df
truvari vcf2df -i annosv/topmed.INV.anno.vcf.gz annosv/topmed.INV.anno.jl

MAKING THE OVERALL VARIANT QC
=============================
After I have the final, filtered set of variants, I then need to do the concat and restriction to the first sample.
Then I can run through annosv and make a variant summary
would be fine to concatenate by type at that point, also

I got a notebook that will work just fine for this

