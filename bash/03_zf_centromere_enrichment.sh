################################################################################
### CALCULATING NON-B DNA MOTIF ENRICHMENT IN ZEBRA FINCH CENTROMERES
### CODE WRITTEN BY LINNÉA SMEDS

##### TABLE OF CONTENTS
# ENRICHMENT IN CENTROMERES
#   -GC CONTENT FOR CENTROMERES
#   -ENRICHMENT IN WINDOWS FOR STATISTICAL TEST



######################## ENRICHMENT IN CENTROMERES  ############################
# 1 July 2025

# Convert centromere gff to bed file
awk -v OFS="\t" '{start=$4-1; print $1,start,$5}' ref/bTaeGut7v0.4_MT_rDNA.centromere_detector.v0.1.gff >ref/bTaeGut7v0.4_MT_rDNA.centromere_detector.v0.1.bed

# Enrichment in centromeres compared to the genome-wide density and to chrom-
# specific density
prefix="bTaeGut7v0.4_MT_rDNA"
cen_bed_file="ref/${prefix}.centromere_detector.v0.1.bed"
module load bedtools/2.31.0
cat $cen_bed_file |while read -r chr start end;
do
  cenlen=`echo "$end-$start" |bc`
  echo cenlen is $cenlen for chr $chr
  echo '#!/bin/bash
  module load bedtools/2.31.0
  echo Working on '$chr'
  rm -f centromeres/'$chr'.tmp-centromere.txt
  cat coverage/'${prefix}'.nonB_genome_wide.txt |grep -v "all" | while read -r nb tot dens;
  do
    chrdens=`grep '$chr' coverage/'${prefix}'.nonB_per_chrom.txt |awk -v n=$nb '"'"'($2==n){print $4}'"'"'`
    echo "Coverage for '$chr' and $nb is $chrdens"
    stats=`grep '$chr' '$cen_bed_file' |sort -k1,1 -k2,2n |intersectBed -a - -b final_nonB/'${prefix}'.$nb.merged.bed | sort -k1,1 -k2,2n | mergeBed -i - |\
    awk -v l='$cenlen' -v dtot=$dens -v dchr=$chrdens '"'"'{sum+=$3-$2}END{d=sum/l; fracgw=d/dtot; fracchr=d/dchr; print d,fracgw,fracchr}'"'"'`
    echo '$chr'" "$nb" "$stats |sed "s/ /\t/g" >>centromeres/'$chr'.tmp-centromere.txt
  done
  ' | sbatch -J $chr --ntasks=1 --cpus-per-task=1 --mem-per-cpu=4G --out slurm/job.centromeric_enrichment.$chr.%j.out
done
#Merge
echo "Chr nonB Density Enrichment_GW Enrichment_Chr" |sed "s/ /\t/g" >centromeres/${prefix}.enrichment.txt
cat $cen_bed_file |cut -f1 |while read -r chr;
do
  cat centromeres/$chr.tmp-centromere.txt >> centromeres/${prefix}.enrichment.txt
done
# Remove temporary files 
rm centromeres/*.tmp-centromere*


# ~~~~~~~~~~~~~~~~~~~~~~~~~ GC CONTENT FOR CENTROMERES ~~~~~~~~~~~~~~~~~~~~~~~~~
# Check GC content in centromeres (this was not used in the paper)
prefix="bTaeGut7v0.4_MT_rDNA"
cen_bed_file="ref/${prefix}.centromere_detector.v0.1.bed"
module load bedtools/2.31.0
bedtools nuc -fi ref/$prefix.fa -bed $cen_bed_file >centromeres/cen_GC_fullInfo.nuc
echo "Chr GCcont" |sed 's/ /\t/g' >centromeres/cen_GC.tsv
awk '(NR>1){print $1"\t"$5}' centromeres/cen_GC_fullInfo.nuc >> centromeres/cen_GC.tsv


# ~~~~~~~~~~~~~~~~~ ENRICHMENT IN WINDOWS FOR STATISTICAL TEST ~~~~~~~~~~~~~~~~~
# To have something to compare the centromere to, we divide the rest of the
# chromosomes into segments with the same length as each centromere (excluding
# the actual centromere)

module load bedtools/2.31.0
mkdir centromere
mkdir -p centromeres/windows/

# Generating 100 windows of the same size as the centromere, not overlapping with
# the real centromere. Randomly selected if there are more than 100 windows,
# overlapping if there are less. For very large centromeres, removing overlap
# with the real centromere and non-complete windows at the end of the chromosome
# can result in less than 100 windows. To make sure there always are more than
# 100, I make the overlap so there are 120 windows to choose from, and allow end
# windows to be as small as 1M (if this is smaller than the centromere size).
prefix="bTaeGut7v0.4_MT_rDNA"
for chr in $(cut -f1 ref/$prefix.centromere_detector.v0.1.bed |uniq)
do
    echo '#!/bin/bash
    echo "Looking at '$chr'"
    module load bedtools/2.31.0
    c_len=`awk -v c='$chr' '"'"'($1==c){sum+=$3-$2}END{print sum}'"'"' ref/'$prefix'.centromere_detector.v0.1.bed`
    chrlen=`awk -v c='$chr' '"'"'($1==c){print $2}'"'"' ref/'$prefix'.fa.fai`
    poswin=`echo $chrlen"/"$c_len |bc`
    minsize=100000
    if [ "$poswin" -lt "120" ]
    then
      echo "For '$chr', there will be less than 120 nonovl windows"
      slide=`echo "("$chrlen"-"$c_len")/120"|bc`
      bedtools makewindows -g <(awk -v c='$chr' '"'"'(c==$1){print}'"'"' ref/'$prefix'.fa.fai) -s $slide -w $c_len |subtractBed -a - -b ref/'$prefix'.centromere_detector.v0.1.bed |awk -v min=$minsize -v clen=$c_len '"'"'($3-$2>min || $3-$2==clen){print}'"'"' |shuf -n 100 |sort -k2,2n >centromeres/windows/'$chr'.exclCen.bed
    else
      echo "For '$chr', there will be more than 120 nonovl windows"
      bedtools makewindows -g <(awk -v c='$chr' '"'"'(c==$1){print}'"'"' ref/'$prefix'.fa.fai) -w $c_len |subtractBed -a - -b ref/'$prefix'.centromere_detector.v0.1.bed |awk -v min=$minsize -v clen=$c_len '"'"'($3-$2>min || $3-$2==clen){print}'"'"' |shuf -n 100 |sort -k2,2n >centromeres/windows/'$chr'.exclCen.bed
    fi
    bedtools nuc -fi ref/'$prefix'.fa -bed ref/'$prefix'.centromere_detector.v0.1.bed >centromeres/windows/'$chr'.exclCen.with_GC.bed
    rm -f centromeres/windows/'$chr'.100rand.enrichment.tsv
    tmp="'$chr'"
    while read line
    do
      cat coverage/'${prefix}'.nonB_genome_wide.txt | while read -r non_b tot dens;
      do
          echo $line |sed "s/ /\t/g" >tmp.$SLURM_JOB_ID.bed
          d=`intersectBed -a tmp.$SLURM_JOB_ID.bed -b final_nonB/'$prefix'.${non_b}.merged.bed -wo |awk -v l=$c_len -v dtot=$dens '"'"'{sum+=$7}END{d=sum/l; frac=d/dtot; print frac}'"'"'`
          tmp=`echo $tmp" "$d`
          echo $tmp >tmp.$SLURM_JOB_ID
      done
      cat tmp.$SLURM_JOB_ID |sed "s/ /\t/g" >>centromeres/windows/'$chr'.100rand.enrichment.tsv
    done <centromeres/windows/'$chr'.exclCen.bed
    '| sbatch -J $chr --ntasks=1 --cpus-per-task=1 --time=5:00:00 --out slurm/job.centromeric_enrichment.$chr.%j.out
done
# Check if there are more or less than 120 windows
cat slurm/job.centromeric_enrichment.chr*_*at.*.out |grep "there will be"

# Merge background densities and GC
echo "Window Chr G4 APR DR DRfilt IR IRall MRall TRI STR Z Zgfa ZDNAm1 Zseeker G4quadron Any" |sed 's/ /\t/g'  >centromeres/${prefix}.background.enrichment.tsv
for file in $(ls centromeres/windows/chr*.100rand.enrichment.tsv)
do
  awk '{print NR,$0}' $file |sed 's/ /\t/g' >>centromeres/${prefix}.background.enrichment.tsv
done
echo "Window Chr GCcont" |sed 's/ /\t/g'  >centromeres/${prefix}.background.GC.tsv
for file in $(ls centromeres/windows/chr*.exclCen.with_GC.bed)
do
  awk '(NR>1){rownum=NR-1; print rownum,$1,$5}' $file |sed 's/ /\t/g' >>centromeres/${prefix}.background.GC.tsv
done

# The same but compared to the per chromosome densities instead of genome-wide
prefix="bTaeGut7v0.4_MT_rDNA"
for chr in $(cut -f1 ref/$prefix.centromere_detector.v0.1.bed |uniq |tail -n+2)
do
    echo '#!/bin/bash
    echo "Looking at '$chr'"
    module load bedtools/2.31.0
    c_len=`awk -v c='$chr' '"'"'($1==c){sum+=$3-$2}END{print sum}'"'"' ref/'$prefix'.centromere_detector.v0.1.bed`
    tmp="'$chr'"
    rm -f centromeres/windows/'$chr'.100rand.enrichmentCHR.tsv
    while read line
    do
      echo $line |sed "s/ /\t/g" >tmp.$SLURM_JOB_ID.bed
      cat coverage/'${prefix}'.nonB_per_chrom.txt |grep '$chr' |grep -v "ALL" | while read -r c non_b tot dens;
      do
          d=`intersectBed -a tmp.$SLURM_JOB_ID.bed -b final_nonB/'$prefix'.${non_b}.merged.bed -wo |awk -v l=$c_len -v dtot=$dens '"'"'{sum+=$7}END{d=sum/l; frac=d/dtot; print frac}'"'"'`
          tmp=`echo $tmp" "$d`
          echo $tmp >tmp.$SLURM_JOB_ID
      done
      cat tmp.$SLURM_JOB_ID |sed "s/ /\t/g" >>centromeres/windows/'$chr'.100rand.enrichmentCHR.tsv
    done <centromeres/windows/'$chr'.exclCen.bed
    '| sbatch -J $chr --ntasks=1 --cpus-per-task=1 --partition=open --time=5:00:00 --out slurm/job.centromeric_enrichmentCHR.$chr.%j.out
done
# Merge
echo "Window Chr G4 APR DR DRfilt IR IRall MRall TRI STR Z Zgfa ZDNAm1 Zseeker G4quadron Any" |sed 's/ /\t/g'  >centromeres/${prefix}.background.enrichment.CHR.tsv
for file in $(ls centromeres/windows/chr*.100rand.enrichmentCHR.tsv)
do
  awk '{print NR,$0}' $file |sed 's/ /\t/g' >>centromeres/${prefix}.background.enrichment.CHR.tsv
done


# Make an ordered chromosome file to be used for plotting
echo "Chr Length" |sed 's/ /\t/g'  >centromeres/chr_len.ordered.txt
for i in '1' '1A' '2' '3' '4' '4A' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' '23' '24' '25' '26' '27' '28' '29' '30' '31' '32' '33' '34' '35' '36' '37'
do
  for hap in "mat" "pat"
  do
    grep "chr${i}_$hap" ref/bTaeGut7v0.4_MT_rDNA.fa.fai |cut -f1,2 >>centromeres/chr_len.ordered.txt
  done
done
grep "chrZ" ref/bTaeGut7v0.4_MT_rDNA.fa.fai |cut -f1,2 >>centromeres/chr_len.ordered.txt
grep "chrW" ref/bTaeGut7v0.4_MT_rDNA.fa.fai |cut -f1,2 >>centromeres/chr_len.ordered.txt


