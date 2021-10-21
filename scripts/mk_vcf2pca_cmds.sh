BASEDIR=/users/u233287/scratch/topmed_analysis
INDIR=$BASEDIR/freeze1.1
OUTDIR=$BASEDIR/pca_data
mkdir -p $OUTDIR
for i in $INDIR/*DUP*.vcf.gz
do
    fname=$(basename $i)
    oname=$OUTDIR/${fname%.vcf.gz}.h5
    echo python /users/u233287/scratch/topmed_analysis/scripts/vcf_to_pca.py $i $oname > $BASEDIR/jobs/pcah5_${fname}.sh
done
