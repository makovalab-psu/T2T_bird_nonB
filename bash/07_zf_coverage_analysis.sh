################################################################################
### ANALYSING THE HIFI AND ONT COVERAGE IN RELATION TO NON-B DNA MOTIFS
### CODE WRITTEN BY LINNÉA SMEDS

########################### COVERAGE ANALYSIS ##################################
# Download HiFi and ONT coverage files
cd ref
wget https://genomeark.s3.amazonaws.com/species/Taeniopygia_guttata/bTaeGut7/assembly_verkko_0.1/manual_curation/bTaeGut7v0.4/mapping/v0.4_dip_hifi/v0.4_dip_hifi.pri.cov.wig
wget https://genomeark.s3.amazonaws.com/species/Taeniopygia_guttata/bTaeGut7/assembly_verkko_0.1/manual_curation/bTaeGut7v0.4/mapping/v0.4_dip_ont/v0.4_dip_ont.pri.cov.wig
cd ..
convert2bed -i wig <ref/v0.4_dip_hifi.pri.cov.wig >annotation/v0.4_dip_hifi.pri.cov.bed
convert2bed -i wig <ref/v0.4_dip_ont.pri.cov.wig >annotation/v0.4_dip_ont.pri.cov.bed

# Overlap coverage and non-B DNA motifs
prefix="bTaeGut7v0.4_MT_rDNA"
for seq in "hifi" "ont"
do
  for n in "G4" "APR" "DR" "IR" "MR" "TRI" "STR" "Z" "ALL"
  do
   ls ref/v0.4_dip_$seq.pri.cov.bed
   echo '#!/bin/bash
   module load bedtools/2.31.0
   intersectBed -wao -a <(cut -f1,2,3 ref/v0.4_dip_'$seq'.pri.cov.bed) -b <(cut -f1,2,3 final_nonB/'$prefix'.'${n}'.merged.bed) | cut -f1,2,3,7 |awk -v OFS="\t" '"'"'{if(NR==0){chr=$1; s=$2; e=$3; sum=$4}else{if($1==chr && $2==s){sum+=$4}else{print chr,s,e,sum; chr=$1; s=$2; e=$3; sum=$4}}}END{print chr,s,e,sum}'"'"' | sed "/^\s*$/d" >densities/'${prefix}'.'${n}'.'$seq'.bed
   '| sbatch -J $n.$seq --ntasks=1 --cpus-per-task=1 --time=5:00:00 --mem-per-cpu=8G --out slurm/density.$n.$seq.%j.out
  done
done

# Merge files to have seq coverage and nonB coverage in the same files
prefix="bTaeGut7v0.4_MT_rDNA"
for seq in "hifi" "ont"
do
  echo "Chr SeqCov NonB NonBCov" | sed 's/ /\t/g' >densities/$prefix.$seq.nonB_and_seqCov.bed
  for n in "G4" "APR" "DR" "IR" "MR" "TRI" "STR" "Z" "ALL"
  do
    intersectBed -a ref/v0.4_dip_$seq.pri.cov.bed -b densities/$prefix.$n.$seq.bed -wo |awk -v OFS="\t" -v n=$n '{dens=$9/($3-$2); print $1,$5,n,dens}' >>densities/$prefix.$seq.nonB_and_seqCov.bed
  done
done

# Combine with dot chromosome compartment information!
for seq in "hifi" "ont"
do
  for n in "G4" "APR" "DR" "IR" "MR" "TRI" "STR" "Z" "ALL"
  do
    rm -f tmp.compartment.cov.$seq.$n.txt
    echo '#!/bin/bash
    module load bedtools/2.31.0
    intersectBed -a ref/v0.4_dip_'$seq'.pri.cov.bed -b densities/'$prefix'.'$n'.'$seq'.bed -wo |grep -f microdot.txt |cut -f1-5,9 |intersectBed -a - -b compart/'$prefix'.dot.AB.10kb.bed -wao |cut -f1,5,6,10 |awk -v n='$n' '"'"'{print $1,n,$2,$3,$4}'"'"' >>tmp.compartment.cov.'$seq'.'$n'.txt
    '| sbatch -J $n.$seq --ntasks=1 --cpus-per-task=1 --time=5:00:00 --mem-per-cpu=8G --out slurm/density.compartment.$n.$seq.%j.out
  done
done
# Merge
for seq in "hifi" "ont"
do
  echo "Chr NonB SeqCov NonBCov Comp" | sed 's/ /\t/g' >densities/$prefix.$seq.nonB_and_seqCov.dotComp.tsv
  cat tmp.compartment.cov.$seq.*.txt | sed 's/ /\t/g' >>densities/$prefix.$seq.nonB_and_seqCov.dotComp.tsv
done

