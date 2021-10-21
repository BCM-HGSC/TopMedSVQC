BASEDIR=/users/u233287/scratch/topmed_analysis
INDIR=$BASEDIR/freeze1
OUTDIR=$BASEDIR/freeze1.1
SRC=$BASEDIR/scripts/id_subsetter.py
SVMFN=$BASEDIR/metadata/duplications.svm.filter.10.14.txt.gz

for chr in $(seq 1 22)
do
    echo "python $SRC $SVMFN $INDIR/sv.freeze1.chr${chr}.gt.only.bcf \
        | bgzip > $OUTDIR/sv.freeze1.1.chr${chr}.DUP.gt.only.vcf.gz" > jobs/dup_idpull.chr${chr}.sh
done
