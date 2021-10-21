# Runs Truvari With Defaults (no sequence similarity)


BASEDIR=/stornext/snfs4/next-gen/scratch/english/round2/topmed_analysis

COMPVCF=$BASEDIR/validation/topmed.replicate.DEL.vcf.gz
OUTDIR=$BASEDIR/validation/results/deletions_new

#COMPVCF=$BASEDIR/validation/topmed.replicate.DUP.vcf.gz
#OUTDIR=$BASEDIR/validation/results/duplications

BASEVCFDIR=$BASEDIR/validation
MANIFEST=$BASEDIR/metadata/validation_replicate.manifest.txt

mkdir -p $OUTDIR

cat $MANIFEST | grep -v sample | while read sample subject center fn
do
    job1=$BASEDIR/jobs/${sample}.truvari.sh
    echo truvari bench -b $BASEVCFDIR/${subject}.base.vcf.gz -c $COMPVCF --cSample $sample \
            --multimatch --no-ref c -o ${OUTDIR}/${sample} --pctsim 0 --passonly -s 10 -S 10 --sizemax 100000000 > ${job1}
    echo truvari vcf2df -i -b ${OUTDIR}/${sample} ${OUTDIR}/${sample}/data.jl >> ${job1}
    echo "for i in ${OUTDIR}/${sample}/*.vcf; do vcf_compress \$i; done" >> ${job1}
    echo "rm $OUTDIR/${sample}/*.vcf" >> ${job1}

    #job2=$BASEDIR/jobs/${sample}.truvari.up.sh
    #echo truvari bench -b $BASEVCFDIR/${subject}.base.vcf.gz -c $COMPVCF --cSample $sample -s 10 -S 10 --sizemax 100000000 \
         #--multimatch --no-ref c -o ${OUTDIR}/${sample}_up --pctsim 0 --pctsize 0 --refdist 1500 --passonly > ${job2}
    #echo truvari vcf2df -i -b ${OUTDIR}/${sample}_up ${OUTDIR}/${sample}_up/data.jl >> ${job2}
    #echo "for i in ${OUTDIR}/${sample}_up/*.vcf; do vcf_compress \$i; done" >> ${job2}
    #echo "rm $OUTDIR/${sample}_up/*.vcf" >> ${job2}
done
