################################################################################
### CALCULATING NON-B DNA MOTIF ENRICHMENT IN ZEBRA FINCH REPEATS
### CODE WRITTEN BY LINNÉA SMEDS

##### TABLE OF CONTENTS
# PREPROCESSING THE REPEAT ANNOTATION
#   -FIND ALL OVERLAPS WITH NON-B DNA MOTIFS
#   -MAKE A REPEAT LIST WITH LENGTHS
# CALCULATE ENRICHMENT IN THE REPEAT TYPES
#   -MANUAL CHECK OF INTERESTING REPEATS
#   -CALCULATE ENRICHMENT PER CHROMOSOME TYPE


################### PREPROCESSING THE REPEAT ANNOTATION ########################
# 1 July 2025
mkdir annotation
# First, the three different repeat annotations are processed separately, then
# we merge them into one big annotation file.

# EDTA2 file (Extensive de-novo TE Annotator version 2)
zcat ref/bTaeGut7v0.4_MT_rDNA.EDTA2.v0.2.gtf.gz |grep -v "##" |\
awk -v OFS="\t" -F'\t' '{split($9, a, ";");
    name = ""; class = "";
    for (i in a) {
    #    if (a[i] ~ /^Name=/) {
    #        split(a[i], b, "=");
    #        name = b[2];
    #    }
        if (a[i] ~ /^classification=/) {
            split(a[i], c, "=");
            class = c[2];
        }
      }
      s=$4-1; print $1,s,$5,class,$3,$7}' >annotation/EDTA2.v0.2.bed

#There are both "LTR/unknown" and "LTR/Unknown" in this file! Merge them into
# one:
 sed 's/LTR\/unknown/LTR\/Unknown/g' annotation/EDTA2.v0.2.bed >edta.tmp
 mv edta.tmp annotation/EDTA2.v0.2.bed

# TRFs; add how many basepairs each repeat unit is, and group into categories
 cat ref/bTaeGut7v0.4_MT_rDNA.trf.sorted.v0.1.bed |\
  awk -v OFS="\t" '{l=length($7);
      if(l<5){g="1-4mer"}
      else if(l>=5 && l<11){g="5-10mer"}
      else if(l>=11 && l<51){g="11-50mer"}
      else if(l>=51 && l<101){g="51-100mer"}
      else if(l>=101){g="100+mer"};
      print $1,$2,$3,g,l"-mer"}' >annotation/TRF_withMers.bed

# The Satellites are already in bedformat with satellite name in the 4th column,
# but I'll make a new file just for consistency. 
cut -f1-4 ref/bTaeGut7v0.4_MT_rDNA.satellome.v0.1.bed >functional/annotation/Satellites.bed

# How many bp are each file?
module load bedtools/2.31.0
for type in "EDTA2.v0.2" "TRF_withMers" "Satellites"
do
  bp=`mergeBed -i annotation/$type.bed |awk '{sum+=$3-$2}END{print sum}'`
  echo $type" "$bp
done
#EDTA2.v0.2 322876394
#TRF_withMers 114340899
#Satellites 61074612

# And how much overlap are there:
intersectBed -a <(mergeBed -i annotation/EDTA2.v0.2.bed) -b <(mergeBed -i annotation/TRF_withMers.bed) |awk '{sum+=$3-$2}END{print sum}'
#62145892
intersectBed -a <(mergeBed -i annotation/EDTA2.v0.2.bed) -b <(mergeBed -i annotation/Satellites.bed) |awk '{sum+=$3-$2}END{print sum}'
#57934276
intersectBed -a <(mergeBed -i annotation/TRF_withMers.bed) -b <(mergeBed -i annotation/Satellites.bed) |awk '{sum+=$3-$2}END{print sum}'
#39932430
# Almost all satellites are in the EDTA file, and more than half are found in
# the TRF file.
# And more than half of the TRFs are in the EDTA2 file..

# Despite the overlap, I run the enrichment analysis for all these regions

# Also merge to get rid of overlaps 
echo '#!/bin/bash
module load bedtools/2.31.0
cat annotation/EDTA2.v0.2.bed annotation/TRF_withMers.bed annotation/Satellites.bed |cut -f1-3 |sort -k1,1 -k2,2n |mergeBed -i - >annotation/all_repeats.bed
' | sbatch -J repeats --ntasks=1 --cpus-per-task=2 --mem-per-cpu=4G --out slurm/job.merge-repeats.%j.out


# ~~~~~~~~~~~~~~~~~~ FIND ALL OVERLAPS WITH NON-B DNA MOTIFS ~~~~~~~~~~~~~~~~~~~
mkdir -p repeats/overlap
prefix="bTaeGut7v0.4_MT_rDNA"
module load bedtools/2.31.0
for non_b in  "G4" "APR" "DR" "DRfilt" "IR" "IRall" "MRall" "TRI" "STR" "Z" "Zgfa" "ZDNAm1" "Zseeker" "G4quadron" "Any"
do
  echo '#!/bin/bash
  intersectBed -a annotation/'${prefix}'.EDTA2.v0.2.bed -b <(sort -k1,1 -k2,2n final_nonB/'${prefix}'.'${non_b}'.bed |mergeBed -i - ) >repeats/overlap/'${non_b}'.TEs.bed
  intersectBed -a annotation/'${prefix}'.TRF_withMers.bed -b <(sort -k1,1 -k2,2n final_nonB/'${prefix}'.'${non_b}'.bed |mergeBed -i - ) >repeats/overlap/'${non_b}'.TRF.bed
  intersectBed -a annotation/'${prefix}'.Satellites.bed -b <(sort -k1,1 -k2,2n final_nonB/'${prefix}'.'${non_b}'.bed | mergeBed -i - ) >repeats/overlap/'${non_b}'.SAT.bed
  '| sbatch -J $non_b --ntasks=1 --cpus-per-task=1 --mem-per-cpu=4G --time=1:00:00 --partition=open --out slurm/job.overlap-repeats.$non_b.%j.out
done


# ~~~~~~~~~~~~~~~~~~~~~~ MAKE A REPEAT LIST WITH LENGTHS ~~~~~~~~~~~~~~~~~~~~~~~
# Using a script from the non-B in Apes paper, origin:
# https://raw.githubusercontent.com/makovalab-psu/T2T_primate_nonB/refs/heads/main/python/repeat_summary.py

module load python/3.11.2
cat functional/annotation/EDTA2.v0.2.bed functional/annotation/TRF_withMers.bed functional/annotation/Satellites.bed |python3 python/repeat_summary.py >repeats/TE_TRF_SAT_lengths.txt


# #################### CALCULATE ENRICHMENT IN THE REPEATS #####################
# Summarize the enrichment of each non-B motif in all repeat classes
# (Print NA if there are no such repeats in the genome)
# Run all repeat types simultaneously and merge afterwards
prefix="bTaeGut7v0.4_MT_rDNA"
cat repeats/TE_TRF_SAT_lengths.txt | while read -r rep replen;
do
  name=`echo $rep |sed 's/\//-/g'`
  echo $name
  echo '#!/bin/bash
  rm -f tmp.'$name'
  cat coverage/'${prefix}'.nonB_genome_wide.txt |grep -v "all" | while read -r non_b tot dens;
  do
    d=`cat repeats/overlap/${non_b}.*.bed | awk -v r='$rep' -v l='$replen' -v dtot=$dens '"'"'($4==r){sum+=$3-$2}END{if(l==0 || dtot==0){print "NA"} else{d=sum/l; frac=d/dtot; print d,frac}}'"'"'`
      echo '$rep'" "$non_b" "$d >>tmp.'$name'
  done
  '| sbatch -J $rep --ntasks=1 --cpus-per-task=1 --time=1:00:00 --partition=open --out slurm/repeat-enrich.$name.%j.out
done
# Merge
echo "Repeat nonB Coverage Enrichment_gw" |sed "s/ /\t/g" >repeats/${prefix}.enrichment.tsv
cat repeats/TE_TRF_SAT_lengths.txt | while read -r rep replen;
do
  name=`echo $rep |sed 's/\//-/g'`
  cat tmp.$name |sed "s/ /\t/g" >>repeats/${prefix}.enrichment.tsv
done

# The length file was ordered by hand with row numbers added to be used for
# plotting. It can be found in helpfiles/TE_TRF_SAT_lengths.ordered.txt


# Add a number column to the length file after reorder it by hand
awk '{print NR"\t"$0}' repeats/TE_TRF_SAT_lengths.txt >repeats/TE_TRF_SAT_lengths.ordered.txt



# ~~~~~~~~~~~~~~~~~~~~ MANUAL CHECK OF INTERESTING REPEATS ~~~~~~~~~~~~~~~~~~~~~
# CHECK NGARO ELEMENTS, STRONGLY ENRICHED FOR Z DNA
grep Ngaro repeats/overlap/Z.TEs.bed |awk '{if(NR==1){chr=$1; l=$3-$2; sum=l}else{if($1==chr){l=$3-$2; sum+=l}else{print chr,sum; chr=$1; l=$3-$2; sum=l}}}END{print chr,sum}'

# Check number of repeats per chromosome
awk '{if(NR==1){chr=$1; sum=$3-$2}else{if($1==chr){sum+=$3-$2}else{print chr,sum; chr=$1; sum=$3-$2}}}END{print chr,sum;}' annotation/all_repeats.bed |sort -k1,1 |join -1 1 -2 1 - <(cat ref/bTaeGut7v0.4_MT_rDNA.fa.fai |sort -k1,1| cut -f1,2) |awk '{d=$2/$3; print $0,d}' |sed 's/ /\t/g' >repeats/density.per.chr.txt

# Get some stats on different repeat types
# Satellites
grep "Satellite/Satellite" annotation/EDTA2.v0.2.bed |sort -k1,1 -k2,2n |mergeBed -i - |awk '{if(NR==1){chr=$1; sum=$3-$2}else{if($1==chr){sum+=$3-$2}else{print chr,sum; chr=$1; sum=$3-$2}}}END{print chr,sum;}' |sort -k1,1 |join -1 1 -2 1 - <(cat ref/bTaeGut7v0.4_MT_rDNA.fa.fai |sort -k1,1| cut -f1,2) |awk '{d=$2/$3; print $0,d}' |sed 's/ /\t/g' >repeats/satellite.per.chr.txt
# Simple repeats
grep "Simple_repeat" annotation/EDTA2.v0.2.bed |sort -k1,1 -k2,2n |mergeBed -i - |awk '{if(NR==1){chr=$1; sum=$3-$2}else{if($1==chr){sum+=$3-$2}else{print chr,sum; chr=$1; sum=$3-$2}}}END{print chr,sum;}' |sort -k1,1 |join -1 1 -2 1 - <(cat ref/bTaeGut7v0.4_MT_rDNA.fa.fai |sort -k1,1| cut -f1,2) |awk '{d=$2/$3; print $0,d}' |sed 's/ /\t/g' >repeats/simple.per.chr.txt
# TEs
grep -v "Satellite/Satellite" annotation/EDTA2.v0.2.bed |grep -v "Simple_repeat" |grep -v "rRNA" |grep -v "unknown" |sort -k1,1 -k2,2n |mergeBed -i - |awk '{if(NR==1){chr=$1; sum=$3-$2}else{if($1==chr){sum+=$3-$2}else{print chr,sum; chr=$1; sum=$3-$2}}}END{print chr,sum;}' |sort -k1,1 |join -1 1 -2 1 - <(cat ref/bTaeGut7v0.4_MT_rDNA.fa.fai |sort -k1,1| cut -f1,2) |awk '{d=$2/$3; print $0,d}' |sed 's/ /\t/g'  >repeats/TEs.per.chr.txt
# UNKNOWN
grep "unknown" annotation/EDTA2.v0.2.bed |sort -k1,1 -k2,2n |mergeBed -i - |awk '{if(NR==1){chr=$1; sum=$3-$2}else{if($1==chr){sum+=$3-$2}else{print chr,sum; chr=$1; sum=$3-$2}}}END{print chr,sum;}' |sort -k1,1 |join -1 1 -2 1 - <(cat ref/bTaeGut7v0.4_MT_rDNA.fa.fai |sort -k1,1| cut -f1,2) |awk '{d=$2/$3; print $0,d}' |sed 's/ /\t/g' >repeats/unknown.per.chr.txt


# ~~~~~~~~~~~~~~~~~~ CALCULATE ENRICHMENT PER CHROMOSOME TYPE ~~~~~~~~~~~~~~~~~~
# Added September 5

# First, we need length per type
module load python/3.11.2
for group in "macro" "micro" "microdot"
do
  cat annotation/${prefix}.EDTA2.v0.2.bed annotation/${prefix}.TRF_withMers.bed annotation/${prefix}.Satellites.bed |grep -f helpfiles/$group.txt - |python3 python/repeat_summary.py >repeats/$group.TE_TRF_SAT_lengths.txt
done

# Run all repeat types simultaneously and merge afterwards
prefix="bTaeGut7v0.4_MT_rDNA"
for group in "macro" # "microdot" "micro" "macro" 
do
  cat repeats/$group.TE_TRF_SAT_lengths.txt | while read -r rep replen;
  do
    name=`echo $rep |sed 's/\//-/g'`
    echo $name
    echo '#!/bin/bash
    rm -f tmp.'$group'.'$name'
    cat coverage/'${prefix}'.nonB_genome_wide.txt |grep -v "all" | while read -r non_b tot dens;
    do
      d=`grep -f '$group.txt' repeats/overlap/${non_b}.*.bed | awk -v r='$rep' -v l='$replen' -v dtot=$dens '"'"'($4==r){sum+=$3-$2}END{if(l==0 || dtot==0){print "NA"} else{d=sum/l; frac=d/dtot; print d,frac}}'"'"'`
        echo '$group'" "'$rep'" "$non_b" "$d >>tmp.'$group'.'$name'
    done
    '| sbatch -J $rep --ntasks=1 --cpus-per-task=1 --time=1:00:00 --partition=open --out slurm/repeat-enrich.$name.%j.out
  done
done

# Merge (here we want to include all repeats! If file is missing, add NA!)
echo "Group Repeat nonB Coverage Enrichment_gw" |sed "s/ /\t/g" >repeats/${prefix}.enrichment.group.tsv
for group in "macro" "micro" "microdot"
do
  cat repeats/TE_TRF_SAT_lengths.txt | while read -r rep replen;
  do
    name=`echo $rep |sed 's/\//-/g'`
    if test -e tmp.$group.$name
    then
      cat tmp.$group.$name |sed "s/ /\t/g" >>repeats/${prefix}.enrichment.group.tsv
    else
      echo "$group $rep APR NA NA
$group $rep DR NA NA
$group $rep STR NA NA
$group $rep IR NA NA
$group $rep MR NA NA
$group $rep TRI NA NA
$group $rep G4 NA NA
$group $rep Z NA NA
$group $rep DRfilt NA NA
$group $rep G4quadron NA NA
$group $rep IRall NA NA
$group $rep MRall NA NA
$group $rep Zseeker NA NA
$group $rep ZDNAm1 NA NA  
$group $rep Zgfa NA NA  
$group $rep Any NA NA" | sed "s/ /\t/g" >>repeats/${prefix}.enrichment.group.tsv
    fi
  done
done
