for i in pca_data/*.h5
do 
    bname=jobs/$(basename $i).het_hom.sh
    oname=pca_data/temp_het_hom/$(basename $i).summary.jl
    echo "python scripts/het_hom_summary.py $oname $i" > $bname
done
