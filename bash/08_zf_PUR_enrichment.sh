################################################################################
### COMPARING PREVIOUSLY UNASSEMBLED REGIONS (PURs) TO PREVIOUSLY ASSEMBLED,
### IN TERMS OF NON-B DNA DENSITY
### CODE WRITTEN BY LINNÉA SMEDS


# THIS CODE WAS NOT USED IN THE NON-B DNA PAPER, BUT NUMBERS WERE INCLUDED 
# IN FORMENTI ET AL. 

# Take complement for non-PURs
module load bedtools/2.31.0
complementBed -i ref/bTaeGut7v0.4_MT_rDNA_MATonly_vs_GCF_000151805.1.PUR.fastga.v0.1.bed -g <(grep -E "_mat|chrZ_pat" ref/bTaeGut7v0.4_MT_rDNA.genome) >PUR/bTaeGut7v0.4_MT_rDNA_MATonly_vs_GCF_000151805.nonPUR.bed

# Check lengths
awk '{sum+=$3-$2}END{print sum}' ref/bTaeGut7v0.4_MT_rDNA_MATonly_vs_GCF_000151805.1.PUR.fastga.v0.1.bed
#152538735
awk '{sum+=$3-$2}END{print sum}' PUR/bTaeGut7v0.4_MT_rDNA_MATonly_vs_GCF_000151805.nonPUR.bed
#988609905

prefix="bTaeGut7v0.4_MT_rDNA"
pur_bed="ref/${prefix}_MATonly_vs_GCF_000151805.1.PUR.fastga.v0.1.bed"
nonpur_bed="PUR/${prefix}_MATonly_vs_GCF_000151805.nonPUR.bed"

module load bedtools/2.31.0

# Calculate non-B density inside and outside PUR
# FULL GENOME
echo "NonB PURTotBp PURNonB PURDens OldTotBp OldNonB OldDens" >PUR/genome_summary.MATonly.txt
for nb in "APR" "DR" "STR" "IR" "MR" "TRI" "G4" "Z" "ALL"
do
  echo '#!/bin/bash
  module load bedtools/2.31.0
  new=`intersectBed -a '$pur_bed' -b final_nonB/'${prefix}'.'${nb}'.merged.bed -wao |cut -f1,2,3,7| awk -v OFS="\t" '"'"'{if(NR==0){chr=$1; s=$2; e=$3; sum=$4}else{if($1==chr && $2==s){sum+=$4}else{print chr,s,e,sum; chr=$1; s=$2; e=$3; sum=$4}}}END{print chr,s,e,sum}'"'"' | sed "/^\s*$/d" |awk '"'"'{sum_l+=$3-$2; sum_nb+=$4}END{d=sum_nb/sum_l; print sum_l,sum_nb,d}'"'"'`
  old=`intersectBed -a '$nonpur_bed' -b final_nonB/'${prefix}'.'${nb}'.merged.bed -wao |cut -f1,2,3,7| awk -v OFS="\t" '"'"'{if(NR==0){chr=$1; s=$2; e=$3; sum=$4}else{if($1==chr && $2==s){sum+=$4}else{print chr,s,e,sum; chr=$1; s=$2; e=$3; sum=$4}}}END{print chr,s,e,sum}'"'"' | sed "/^\s*$/d" |awk '"'"'{sum_l+=$3-$2; sum_nb+=$4}END{d=sum_nb/sum_l; print sum_l,sum_nb,d}'"'"'`
  echo '$nb'" "$new" "$old >>PUR/genome_summary.MATonly.txt
 ' |sbatch -J $nb --ntasks=1 --cpus-per-task=1 --time=2:00:00 --partition=open --mem-per-cpu=4G --out slurm/job.PUR.$nb.%j.out
done

# Per chromosome
cat ref/${prefix}*.fa.fai |cut -f1 |grep -v "chrM" |grep -v "rDNA" |while read -r chr;
do
  echo '#!/bin/bash
  module load bedtools/2.31.0
  rm -f PUR/'$chr'_summary.MATonly.txt
  for nb in "APR" "DR" "STR" "IR" "MR" "TRI" "G4" "Z" "ALL"
  do
    new=`grep '$chr' '$pur_bed' | intersectBed -a - -b -b final_nonB/'${prefix}'.${nb}.merged.bed -wao |cut -f1,2,3,7| awk -v OFS="\t" '"'"'{if(NR==0){chr=$1; s=$2; e=$3; sum=$4}else{if($1==chr && $2==s){sum+=$4}else{print chr,s,e,sum; chr=$1; s=$2; e=$3; sum=$4}}}END{print chr,s,e,sum}'"'"' | sed "/^\s*$/d" |awk '"'"'{sum_l+=$3-$2; sum_nb+=$4}END{d=sum_nb/sum_l; print sum_l,sum_nb,d}'"'"'`
    old=`grep '$chr' '$nonpur_bed' | intersectBed -a - -b -b final_nonB/'${prefix}'.${nb}.merged.bed -wao |cut -f1,2,3,7| awk -v OFS="\t" '"'"'{if(NR==0){chr=$1; s=$2; e=$3; sum=$4}else{if($1==chr && $2==s){sum+=$4}else{print chr,s,e,sum; chr=$1; s=$2; e=$3; sum=$4}}}END{print chr,s,e,sum}'"'"' | sed "/^\s*$/d" |awk '"'"'{sum_l+=$3-$2; sum_nb+=$4}END{d=sum_nb/sum_l; print sum_l,sum_nb,d}'"'"'`
    echo '$chr'" "$nb" "$new" "$old >>PUR/'$chr'_summary.MATonly.txt
  done
 ' |sbatch -J $chr --ntasks=1 --cpus-per-task=1 --time=1:00:00 --partition=open --mem-per-cpu=4G --out slurm/job.PUR.$chr.%j.out
done
# Merge into one file
echo "#Chr NonB PURTotBp PURNonB PURDens OldTotBp OldNonB OldDens" >PUR/chromosome_summary.MATonly.txt
cat ref/${prefix}*.fa.fai |cut -f1 |grep -v "chrM" |grep -v "rDNA" |while read -r chr;
do
  cat PUR/${chr}_summary.txt >>PUR/chromosome_summary.MATonly.txt
done
