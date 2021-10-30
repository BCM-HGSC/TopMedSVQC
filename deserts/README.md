# Deserts and Hotspots analysis

# Steps
##  Turn VCFs into BED

bcftools query -i "SVLEN >= 50 || SVLEN <= -50" -f "%CHROM\t%POS\t%END\n" ../call_only_vcfs/topmed.DEL.vcf.gz \
	| bedtools sort | bgzip > topmed.DEL.bed.gz

2) Make the control data. Until I get HGSV or whatever, I'll use my msru
bcftools query -i "SVLEN >= 50 || SVLEN <= -50" -i "SVTYPE == 'DEL'" -f "%CHROM\t%POS\t%END\n" \
~/scratch/insertion_ref/msru/data/inter_merge/grch38/strict/strict.vcf.gz | bedtools sort | bgzip > control.DEL.bed.gz

3) Here's the notes of what we're trying to do
the method we used in our paper for your reference:
  1. split human reference genome into 100 kb windows by "bedtools makewindows"
	$ bedtools makewindows -w 100000 -g grch38.genome.txt -i winnum > grch38.100kbwin.bed

###NOTE! I removed sex chromosomes manually from the genome bed because topmed does not have calls on it
# Also, the following steps are for documentation only. I don't run any of these commands

  2. calculate SV number in each window using "bedtools intersect -a genome_window.bed -b sv.bed -c"
	bedtools intersect -a grch38.100kbwin.bed -b control.DEL.bed.gz -c > control.del_svperwindow.bed

  3. a window containging > XX (mean + 3* standard deviation) SVs is defined as a hotspot

  4. we define a window that does not contain any SV as a desert


4) But first, I'm going to add a step to remove the telomeres/centromeres/gaps
I won't worry about expanding the regions too much, I reckon
Download gap and centromere mapping tracks for grch38 from UCSC Table Browser
Concatenate, and sort

cat <(cut -f1,2,3 grch38.centromeres.bed) grch38.gap.bed | bedtools sort | bgzip > grch38.exclude_regions.bed.gz

Then subtract from the windows

bedtools subtract -A -a grch38.100kbwin.bed -b grch38.exclude_regions.bed.gz | bedtools sort | bgzip > grch38.100kbwin.regions_excluded.bed.gz

4) Now I need to run the window intersection with the SVs
$ bedtools intersect -a grch38.100kbwin.regions_excluded.bed.gz  -b control.DEL.bed.gz -c > control.del_svperwindow.bed

$ bedtools intersect -a grch38.100kbwin.regions_excluded.bed.gz  -b topmed.DEL.bed.gz -c > topmed.del_svperwindow.bed

5) I have a custom sccript to annotate this bed file (steps 3.3 and 3.4 above)
    $ python hotspotDesert.py control.del_svperwindow.bed control.del_svperwindow_anno.bed
    $ python hotspotDesert.py topmed.del_svperwindow.bed topmed.del_svperwindow_anno.bed

### Control
count    26792.000000
mean         5.421954
std          9.081182
min          0.000000
25%          1.000000
50%          3.000000
75%          6.000000
max        220.000000
Name: svcount, dtype: float64
Setting HotSpot threshold at 32.67
Counts
Des    4690
Hot     428
Name: anno, dtype: int64

### Topmed
count    26792.000000
mean        10.640714
std         12.425114
min          0.000000
25%          6.000000
50%          9.000000
75%         13.000000
max        571.000000
Name: svcount, dtype: float64
Setting HotSpot threshold at 47.92
Counts
Des    502
Hot    188
Name: anno, dtype: int64


6) Quick Summary 
$ paste control.del_svperwindow_anno.bed  topmed.del_svperwindow_anno.bed | cut -f6,12 | sort | uniq -c
	ctrl 	tpmd
  21071
    440         Des
    163         Hot
   4613 Des
     54 Des     Des
     23 Des     Hot
    418 Hot
      8 Hot     Des
      2 Hot     Hot
      1 anno    anno


# Note here's the summary WITHOUT gap/cent filtering
  	ctrl	tpmd
  21166
   1838         Des
    174         Hot
   4626 Des
   2592 Des     Des
     25 Des     Hot
    430 Hot
     41 Hot     Des
      2 Hot     Hot
      1 anno    anno

And look at those interesting regions viao
paste control.del_svperwindow_anno.bed  topmed.del_svperwindow_anno.bed | awk '{if ($6 == "Des" && $12 == "Des") print $0}'
