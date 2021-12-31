#!/bin/bash
#PBS -N anno
#PBS -o anno.out
#PBS -e anno.err
#PBS -l nodes=1:ppn=2,mem=50gb
#PBS -l walltime=44:00:00
#PBS -q analysis
#PBS -A proj-ng0013

##vcftools --gzvcf $reads --weir-fst-pop population_white.txt --weir-fst-pop population_black.txt --weir-fst-pop population_Hist.txt --maf 0.05 --max-maf 0.9  --out $reads'.5k.fst'  --max-indv 5000

vcftools --gzvcf $reads --weir-fst-pop population_white.txt --weir-fst-pop population_black.txt --weir-fst-pop population_Hist.txt --maf 0.010 --max-maf 0.99  --out $reads'.15k_AF01_.fst'  --max-indv 15000

vcftools --gzvcf $reads --weir-fst-pop population_white.txt --weir-fst-pop population_black.txt  --maf 0.010 --max-maf 0.99  --out $reads'.15k_AF01_WB.fst'  --max-indv 15000

vcftools --gzvcf $reads --weir-fst-pop population_white.txt  --weir-fst-pop population_Hist.txt --maf 0.010 --max-maf 0.99  --out $reads'.15k_AF01_WH.fst'  --max-indv 15000

vcftools --gzvcf $reads --weir-fst-pop population_black.txt --weir-fst-pop population_Hist.txt --maf 0.010 --max-maf 0.99  --out $reads'.15k_AF01_BH.fst'  --max-indv 15000


##zcat $reads | grep -v '<DUP>' | grep -v '<INV>' > $reads'_del.vcf'
##gzip $reads'_del.vcf'
