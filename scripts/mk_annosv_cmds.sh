mkdir -p jobs/
BASEDIR=/users/u233287/scratch/topmed_analysis
OUTDIR=$BASEDIR/annosv
for i in "$@"
do
    jname=jobs/annosv_$(basename $i).sh
    inter=tmp_$(basename $i).vcf
    mkdir -p annosv/$(basename $i)
    echo "bcftools view $i | cut -f1-10 > \$TMPDIR/$inter" > $jname
    echo ~/scratch/misc_software/AnnotSV/bin/AnnotSV -genomeBuild GRCh38 -outputDir $OUTDIR/$(basename $i) -SVinputFile \$TMPDIR/$inter >> $jname
done
