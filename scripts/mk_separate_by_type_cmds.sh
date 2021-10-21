
BASEDIR=/users/u233287/scratch/topmed_analysis
#$TMPDIR/
INDIR=$BASEDIR/freeze1
OUTDIR=$BASEDIR/freeze1.1

for svtype in DEL INV
do
    for chr in $(seq 1 22)
    do
        echo bcftools view -i \"SVTYPE == \'${svtype}\'\" -O z \
            -o $OUTDIR/sv.freeze1.1.chr${chr}.${svtype}.gt.only.vcf.gz \
            $INDIR/sv.freeze1.chr${chr}.gt.only.bcf > jobs/by_type.chr${chr}.${svtype}.sh
    done
done
