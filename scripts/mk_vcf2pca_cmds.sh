BASEDIR=/users/u233287/scratch/topmed_analysis
INDIR=$BASEDIR/freeze1.1
OUTDIR=$BASEDIR/pca_data
mkdir -p $OUTDIR
for i in $INDIR/*INV*.vcf.gz
do
    fname=$(basename $i)
    oname=$OUTDIR/${fname%.vcf.gz}.h5
    # Subsetting for removing samples
    #echo "bcftools view -S $BASEDIR/metadata/dup_passing_samples.txt -i \"FILTER == 'PASS'\" $i | bgzip > tmp/$fname" > $BASEDIR/jobs/pcah5_${fname}.sh
    echo "bcftools view -i \"FILTER == 'PASS'\" $i | bgzip > tmp/$fname" > $BASEDIR/jobs/pcah5_${fname}.sh
    echo python /users/u233287/scratch/topmed_analysis/scripts/vcf_to_pca.py tmp/$fname $oname >> $BASEDIR/jobs/pcah5_${fname}.sh
done
