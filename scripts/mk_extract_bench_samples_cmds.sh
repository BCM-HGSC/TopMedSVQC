mkdir -p jobs
BASEDIR=/users/u233287/scratch/topmed_analysis
INDIR=$BASEDIR/freeze1.1
OUTDIR=$BASEDIR/validation/raw_calls
mkdir -p $OUTDIR

sample_ids=$BASEDIR/metadata/validation_replicate.sampleids.txt
for i in $INDIR/*INV*.vcf.gz
do
	echo "bcftools view --force-samples -S ${sample_ids} $i | vcf-subset -c ${sample_ids} -e \
          | bgzip > $OUTDIR/$(basename $i)" > $BASEDIR/jobs/validation.subset.$(basename $i).sh
done
