BASEDIR=/stornext/snfs4/next-gen/scratch/english/round2/topmed_analysis
BASEVCFDIR=$BASEDIR/validation
MANIFEST=$BASEDIR/metadata/validation_replicate.manifest.txt
# INCLUDEBED=$BASEVCFDIR/HG002_SVs_Tier1_v0.6_liftover38.bed
INCLUDEBED=$BASEVCFDIR/hc_bed/

cat $MANIFEST | grep -v sample | while read sample subject center fn
do
    COMPVCF=$BASEDIR/validation/topmed.replicate.DUP.vcf.gz
    OUTDIR=$BASEDIR/validation/results/duplications
    PREFIX=DUP
    mkdir -p $OUTDIR

    job1=$BASEDIR/jobs/${sample}.${PREFIX}.truvari.sh
    echo '#!/bin/bash' > ${job1}
    echo truvari bench -b $BASEVCFDIR/${subject}.INS.base.vcf.gz -c $COMPVCF --cSample $sample \
            --no-ref c -P 0.50 -r 1000 -C 2000 -o ${OUTDIR}/${sample} --pctsim 0 \
            --passonly -S 10 --sizemax 100000000 --includebed $INCLUDEBED/${subject}.t1asmcov.bed >> ${job1}
    echo truvari vcf2df -i -b ${OUTDIR}/${sample} ${OUTDIR}/${sample}/data.jl >> ${job1}
    echo "for i in ${OUTDIR}/${sample}/*.vcf; do vcf_compress \$i; done" >> ${job1}
    echo "rm $OUTDIR/${sample}/*.vcf" >> ${job1}


    COMPVCF=$BASEDIR/validation/topmed.replicate.DEL.vcf.gz
    OUTDIR=$BASEDIR/validation/results/deletions
    PREFIX=DEL
    mkdir -p $OUTDIR

    job1=$BASEDIR/jobs/${sample}.${PREFIX}.truvari.sh
    echo '#!/bin/bash' > ${job1}
    echo truvari bench -b $BASEVCFDIR/${subject}.DEL.base.vcf.gz -c $COMPVCF --cSample $sample \
            --no-ref c -P 0.50 -r 1000 -C 2000 -o ${OUTDIR}/${sample} --pctsim 0 \
            --passonly -S 10 --sizemax 100000000 --includebed $INCLUDEBED/${subject}.t1asmcov.bed >> ${job1}
    echo truvari vcf2df -i -b ${OUTDIR}/${sample} ${OUTDIR}/${sample}/data.jl >> ${job1}
    echo "for i in ${OUTDIR}/${sample}/*.vcf; do vcf_compress \$i; done" >> ${job1}
    echo "rm $OUTDIR/${sample}/*.vcf" >> ${job1}
done
