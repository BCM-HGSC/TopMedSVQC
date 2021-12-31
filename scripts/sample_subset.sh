for i in freeze1.1/sv*DUP*.vcf.gz
do
    vcf-subset -c NWD976804 -u freeze1.1/sv.freeze1.1.chr1.DUP.gt.only.vcf.gz | bgzip > tmp/NWD976804.sv.freeze1.1.chr1.DUP.gt.only.vcf.gz
done 
