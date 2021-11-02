# Deserts and Hotspots analysis

Goal is to find regions of the genome that have Fewer or More SVs than would be expected.

# Method description from Shaohua Fan
1. split human reference genome into 100 kb windows by "bedtools makewindows"
2. calculate SV number in each window using "bedtools intersect -a genome_window.bed -b sv.bed -c"
3. a window containging `> XX (mean + 3* standard deviation)` SVs is defined as a hotspot
4. we define a window that does not contain any SV as a desert

# Steps

##  Turn VCFs into BED
TopMed Data
```bash
bcftools query -i "SVLEN >= 50 || SVLEN <= -50" -f "%CHROM\t%POS\t%END\n" ../call_only_vcfs/topmed.DEL.vcf.gz \
	| bedtools sort | bgzip > topmed.DEL.bed.gz
```

Until I get HGSV or whatever, I'll use my msru as control data
```bash
bcftools query -i "SVLEN >= 50 || SVLEN <= -50" -i "SVTYPE == 'DEL'" -f "%CHROM\t%POS\t%END\n" \
    ~/scratch/insertion_ref/msru/data/inter_merge/grch38/strict/strict.vcf.gz \
    | bedtools sort | bgzip > control.DEL.bed.gz
```

And now that I have the HGSVC:
```bash
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/release/v2.0/integrated_callset/variants_freeze4_sv_insdel_alt.vcf.gz
bcftools query -i "SVLEN >= 50 || SVLEN <= -50" -i "SVTYPE == 'DEL'" -f "%CHROM\t%POS\t%END\n" \
    variants_freeze4_sv_insdel_alt.vcf.gz \
    | bedtools sort | bgzip > control.DEL.bed.gz
```

## Get genome file
This is just a `reference.fa.fai`. Ours is `grch38.genome.txt` in notes below

## Make Windows
```bash
bedtools makewindows -w 100000 -g grch38.genome.txt -i winnum > grch38.100kbwin.bed
```

**NOTE!** I removed sex chromosomes manually from the genome bed because topmed does not have calls there.

## Remove the telomeres/centromeres/gaps
Download gap and centromere mapping tracks for grch38 from UCSC Table Browser, concatenate, and sort

```bash
cat <(cut -f1,2,3 grch38.centromeres.bed) grch38.gap.bed \
    | bedtools sort \
    | bgzip > grch38.exclude_regions.bed.gz
```

Then subtract from the windows

```bash
bedtools subtract -A -a grch38.100kbwin.bed -b grch38.exclude_regions.bed.gz \
    | bedtools sort \
    | bgzip > grch38.100kbwin.regions_excluded.bed.gz
```

## Run intersection

```bash
bedtools intersect -a grch38.100kbwin.regions_excluded.bed.gz  -b control.DEL.bed.gz -c > control.del_svperwindow.bed
bedtools intersect -a grch38.100kbwin.regions_excluded.bed.gz  -b topmed.DEL.bed.gz -c > topmed.del_svperwindow.bed
```

## Annotate the regions 
```bash
python hotspotDesert.py control.del_svperwindow.bed control.del_svperwindow_anno.bed
python hotspotDesert.py topmed.del_svperwindow.bed topmed.del_svperwindow_anno.bed
```
Results:
```
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

### HGSVC
count    26792.000000
mean         1.456965
std          2.835152
min          0.000000
25%          0.000000
50%          1.000000
75%          2.000000
max         53.000000
Name: svcount, dtype: float64
Setting HotSpot threshold at 9.96
Counts
Des    11334
Hot      529
Name: anno, dtype: int64
```

## Quick Summary 
```bash
paste control.del_svperwindow_anno.bed  topmed.del_svperwindow_anno.bed | cut -f6,12 | sort | uniq -c
```

Results:
```
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
```

Note here's the summary WITHOUT gap/cent filtering
```
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
```

## hgsvc summary
```
	ctrl	hgsv
  14590
   6842         Des
    242         Hot
    210 Des
   4480 Des     Des
    129 Hot
     12 Hot     Des
    287 Hot     Hot
```

## Three Project Summary
```
        ctrl    tpmd    hgsv
  14223
   6620                 Des
    228                 Hot
    252         Des
    176         Des     Des
     12         Des     Hot
    115         Hot
     46         Hot     Des
      2         Hot     Hot
    200 Des
   4413 Des             Des
      6 Des     Des
     48 Des     Des     Des
      4 Des     Hot
     19 Des     Hot     Des
    125 Hot
     10 Hot             Des
    283 Hot             Hot
      4 Hot     Des
      2 Hot     Des     Des
      2 Hot     Des     Hot
      2 Hot     Hot     Hot
```

## Finding interesting regions
```bash
paste control.del_svperwindow_anno.bed topmed.del_svperwindow_anno.bed \
    | awk '{if ($6 == "Des" && $12 == "Des") print $0}'  > dry_candidates.bed
```

## Annotating Candidate Regions
```bash
~/scratch/misc_software/AnnotSV/bin/AnnotSV -genomeBuild GRCh38 -outputDir annosv -SVinputFile dry_candidates.bed
```


of the 54 dry candidates 39 hit gene(s)
