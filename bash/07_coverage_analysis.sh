################################################################################
### ANALYSING THE HIFI AND ONT COVERAGE IN RELATION TO NON-B DNA MOTIFS
### CODE WRITTEN BY LINNÉA SMEDS

########################### COVERAGE ANALYSIS ##################################
# Download HiFi and ONT coverage files
cd ref
wget https://genomeark.s3.amazonaws.com/species/Taeniopygia_guttata/bTaeGut7/assembly_verkko_0.1/manual_curation/bTaeGut7v0.4/mapping/v0.4_dip_hifi/v0.4_dip_hifi.pri.cov.wig
wget https://genomeark.s3.amazonaws.com/species/Taeniopygia_guttata/bTaeGut7/assembly_verkko_0.1/manual_curation/bTaeGut7v0.4/mapping/v0.4_dip_ont/v0.4_dip_ont.pri.cov.wig
cd ..
convert2bed -i wig <ref/v0.4_dip_hifi.pri.cov.wig >ref/v0.4_dip_hifi.pri.cov.bed
convert2bed -i wig <ref/v0.4_dip_ont.pri.cov.wig >ref/v0.4_dip_ont.pri.cov.bed

# Overlap coverage and non-B DNA motifs
prefix="bTaeGut7v0.4_MT_rDNA"
for seq in "hifi" "ont"
do
  for n in "G4" "APR" "DR" "IR" "TRI" "STR" "Z" 
  do
   ls ref/v0.4_dip_$seq.pri.cov.bed
   echo '#!/bin/bash
   module load bedtools/2.31.0
   intersectBed -wao -a <(cut -f1,2,3 ref/v0.4_dip_'$seq'.pri.cov.bed) -b <(cut -f1,2,3 final_nonB/'$prefix'.'${n}'.merged.bed) | cut -f1,2,3,7 |\
   awk -v OFS="\t" '"'"'{if(NR==0){chr=$1; s=$2; e=$3; sum=$4}else{if($1==chr && $2==s){sum+=$4}else{print chr,s,e,sum; chr=$1; s=$2; e=$3; sum=$4}}}END{print chr,s,e,sum}'"'"' |\
    sed "/^\s*$/d" >coverage/'${prefix}'.'${n}'.'$seq'.bed
   '| sbatch -J $n.$seq --ntasks=1 --cpus-per-task=1 --time=5:00:00 --mem-per-cpu=6G --out slurm/density.$n.$seq.%j.out
  done
done

# Overlap coverage and repeats
prefix="bTaeGut7v0.4_MT_rDNA"
for seq in "hifi" "ont"
do
   ls ref/v0.4_dip_$seq.pri.cov.bed
   echo '#!/bin/bash
   module load bedtools/2.31.0
   intersectBed -wao -a <(cut -f1,2,3 ref/v0.4_dip_'$seq'.pri.cov.bed) -b <(cut -f1,2,3 annotation/all_repeats.bed)  -nonamechheck| cut -f1,2,3,7 |\
   awk -v OFS="\t" '"'"'{if(NR==0){chr=$1; s=$2; e=$3; sum=$4}else{if($1==chr && $2==s){sum+=$4}else{print chr,s,e,sum; chr=$1; s=$2; e=$3; sum=$4}}}END{print chr,s,e,sum}'"'"' |\
    sed "/^\s*$/d" >coverage/'${prefix}'.repeats.'$seq'.bed
   '| sbatch -J rep.$seq --ntasks=1 --cpus-per-task=1 --time=5:00:00 --mem-per-cpu=6G --out slurm/job.density.rep.$seq.%j.out
done

# GC content for each window with coverage (use intersectBed after nuc since the order might not be the same)
prefix="bTaeGut7v0.4_MT_rDNA"
for seq in "hifi" "ont"
do
  echo '#!/bin/bash
  module load bedtools/2.31.0
  bedtools nuc -fi ref/'$prefix'.fa -bed ref/v0.4_dip_'$seq'.pri.cov.bed |cut -f1-3,7 |tail -n+2 >coverage/'${prefix}'.GC.'$seq'.bed
 '| sbatch -J GC --ntasks=1 --cpus-per-task=1 --mem-per-cpu=6G --out slurm/run.coverage.GC.%j.out
done 

# For NonB and repeats, make one version with fraction instead of absolute number 
# Also noticed that the coverage files sometimes have windows that end 
# outside of the chromosome end, make sure to remove such windows because the 
# bedtools nuc will not report them
prefix="bTaeGut7v0.4_MT_rDNA"
cut -f1,2 ref/$prefix.fa.fai |awk 'BEGIN{OFS="\t"}{print $1, 0, $2}' >ref/$prefix.bed
for seq in "hifi" "ont"
do
  for type in "All" "repeats" 
  do
    awk -v OFS="\t" '{frac=$4/($3-$2); print $1,$2,$3,frac}' coverage/${prefix}.$type.$seq.bed | \
     intersectBed -a - -b ref/$prefix.bed -f 1.0 -nonamecheck >coverage/${prefix}.$type.frac.$seq.bed
  done 
done
# Also make sure to only keep full windows for the original Seq files 
for seq in "hifi" "ont"
do
  intersectBed -a ref/v0.4_dip_$seq.pri.cov.bed -b ref/$prefix.bed -f 1.0 -nonamecheck >coverage/${prefix}.SeqCov.$seq.bed
done 

# Get compartment information
prefix="bTaeGut7v0.4_MT_rDNA"
for seq in "hifi" #"ont"
do
   ls coverage/${prefix}.SeqCov.$seq.bed)
   echo '#!/bin/bash
   module load bedtools/2.31.0
   intersectBed -wao -a <(cut -f1,2,3 coverage/'${prefix}'.SeqCov.'$seq'.bed) -b <(cut -f1,2,3,4 ref/'${prefix}'.Cooltools.E1.200kbp.flipped.dip.collated.v0.1.bed) |\
   awk '"'"'BEGIN{OFS="\t"}{
    key = $1"\t"$2"\t"$3
    if (key != prev && prev != "") {
        print prevchr, prevstart, prevend, (n==1 ? comp : "NA")
        n=0
    }
    prev=key; prevchr=$1; prevstart=$2; prevend=$3
    comp=$7
    if(comp==".") comp="NA"
    n++
  }
  END{
      if (prev != "") print prevchr, prevstart, prevend, (n==1 ? comp : "NA")
  }'"'"' >coverage/'${prefix}'.compartment.'$seq'.bed
    '| sbatch -J comp.$seq --ntasks=1 --cpus-per-task=1 --time=5:00:00 --mem-per-cpu=6G --out slurm/job.compartment.$seq.%j.out 
done 
   
# For the linear model I want to have a big table with all info. Since all files are identical, then can be pasted together. 
for seq in "hifi" "ont"
do
  echo "Make sure lengths are the same for $seq"
  wc -l coverage/${prefix}.SeqCov.$seq.bed coverage/${prefix}.All.frac.$seq.bed \
    coverage/${prefix}.repeats.frac.$seq.bed coverage/${prefix}.GC.$seq.bed \
    coverage/${prefix}.compartment.$seq.bed
  echo -e "Chr\tStart\tEnd\tSeqCov\tNonBFrac\tRepFrac\tGCFrac\tComp" >coverage/${prefix}.mergedInfo.$seq.bed
  paste coverage/${prefix}.SeqCov.$seq.bed coverage/${prefix}.All.frac.$seq.bed \
    coverage/${prefix}.repeats.frac.$seq.bed coverage/${prefix}.GC.$seq.bed \
    coverage/${prefix}.compartment.$seq.bed |cut -f1,2,3,5,9,13,17,21 |grep -v rDNA |grep -v "MT" >>coverage/${prefix}.mergedInfo.$seq.bed
done 


# Merge files to have seq coverage and nonB in the same files
module load bedtools/2.31.0
prefix="bTaeGut7v0.4_MT_rDNA"
for seq in "hifi" "ont"
do
  echo "Chr SeqCov Element Coverage" | sed 's/ /\t/g' >coverage/$prefix.$seq.nonB_and_seqCov.bed
  for n in "G4" "APR" "DR" "IR" "TRI" "STR" "Z" "All"
  do
    intersectBed -a ref/v0.4_dip_$seq.pri.cov.bed -b coverage/$prefix.$n.$seq.bed -wo -nonamecheck | \
    awk -v OFS="\t" -v n=$n '{dens=$9/($3-$2); print $1,$5,n,dens}' >>coverage/$prefix.$seq.nonB_and_seqCov.bed
  done
  #    intersectBed -a ref/v0.4_dip_$seq.pri.cov.bed -b coverage/$prefix.repeats.$seq.bed -wo -nonamecheck |\
  #    awk -v OFS="\t" '{dens=$9/($3-$2); print $1,$5,"REP",dens}' >>coverage/$prefix.$seq.nonB_and_seqCov.bed
done

# Combine with dot chromosome compartment information!
for seq in "hifi" "ont"
do
  for n in "G4" "APR" "DR" "IR" "TRI" "STR" "Z" "Any"
  do
    rm -f tmp.compartment.cov.$seq.$n.txt
    echo '#!/bin/bash
    module load bedtools/2.31.0
    intersectBed -a ref/v0.4_dip_'$seq'.pri.cov.bed -b coverage/'$prefix'.'$n'.'$seq'.bed -wo -nonamecheck|\
    grep -f <(grep dot helpfiles/'$prefix'.groups.txt |cut -f1) - |cut -f1-5,9 |\
    intersectBed -a - -b compart/'$prefix'.dot.AB.200kb.bed -wao |cut -f1,5,6,10 |\
    awk -v n='$n' '"'"'{print $1,n,$2,$3,$4}'"'"' >>tmp.compartment.cov.'$seq'.'$n'.txt
    '| sbatch -J $n.$seq --ntasks=1 --cpus-per-task=1 --time=5:00:00 --mem-per-cpu=8G --out slurm/coverage.compartment.$n.$seq.%j.out
  done
done
# Merge
for seq in "hifi" "ont"
do
  echo "Chr NonB SeqCov NonBCov Comp" | sed 's/ /\t/g' >coverage/$prefix.$seq.nonB_and_seqCov.dotComp.tsv
  cat tmp.compartment.cov.$seq.*.txt | sed 's/ /\t/g' >>coverage/$prefix.$seq.nonB_and_seqCov.dotComp.tsv
done

# And GC content to each window with coverage:
prefix="bTaeGut7v0.4_MT_rDNA"
for seq in "hifi" "ont"
do
  module load bedtools/2.31.0
  bedtools nuc -fi ref/$prefix.fa -bed ref/v0.4_dip_$seq.pri.cov.bed >coverage/tmp.$seq.GC_fullInfo.nuc
  # Extract only dot chromosomes
  grep dot helpfiles/$prefix.groups.txt |grep -f - coverage/tmp.$seq.GC_fullInfo.nuc |awk '(NR>1){print $1"\t"$2"\t"$3"\t"$5,"\t"$7*100}' >coverage/$prefix.$seq.GC.dot.bed
done 
  
# Merge this with nonB content 
module load bedtools/2.31.0
for seq in "hifi" "ont"
do
  echo -e "Chr\tStart\tStop\tNonB\tSeqCov\tGC\tNonBCov\tComp" >coverage/$prefix.$seq.GC.nonB.dot.tsv
  for n in "G4" "APR" "DR" "IR" "TRI" "STR" "Z" "Any"
  do
    intersectBed -a coverage/$prefix.$seq.GC.dot.bed -b coverage/$prefix.$n.$seq.bed -wo |\
    cut -f1,2,3,4,5,9 |intersectBed -a - -b compart/$prefix.dot.AB.200kb.bed -wao |\
    awk -v OFS="\t" -v n=$n '{print $1,$2,$3,n,$4,$5,$6,$10}' >>coverage/$prefix.$seq.GC.nonB.dot.tsv
  done 
done

########################### LINEAR MIXED MODEL ##################################
# The linear mixed model is run in R, code found in:
R/LMM_for_seq_coverage.R

