################################################################################
### CALCULATING NON-B DNA MOTIF ENRICHMENT IN DOT CHROMOSOME A AND B COMPARTMENTS
### CODE WRITTEN BY LINNÉA SMEDS

##### TABLE OF CONTENTS
# CALCULATE ENRICHMENT
#    -FUNCTIONAL ENRICHMENT IN A & B COMPARTMENTS
#    -CALCULATE STATS
# INTRON TRFs IN DOT CHROMOSOME COMPARTMENTS

############################ CALCULATE ENRICHMENT ##############################
# Compartments are only available for the zebra finch T2T assembly
mkdir compart

# A and B compartment annotation for 10 kb windows processed by Simona Secomandi
# is found in: annotation/bTaeGut7v0.4_MT_rDNA.Cooltools.v0.2.E1.10Kb.flipped.dip.collated.AB.bed

# Extract only A and B compartments from dot chromosomes
grep dot helpfiles/bTaeGut7v0.4_MT_rDNA.groups.txt |cut -f1 |\
grep -f - annotation/bTaeGut7v0.4_MT_rDNA.Cooltools.v0.2.E1.10Kb.flipped.dip.collated.AB.bed >compart/bTaeGut7v0.4_MT_rDNA.dot.AB.10kb.bed
# Extract all windows from the new file from Simonas
grep dot helpfiles/bTaeGut7v0.4_MT_rDNA.groups.txt |cut -f1 |\
grep -f - ref/bTaeGut7v0.4_MT_rDNA.Cooltools.v0.2.E1.10Kb.flipped.dip.collated.AB.noNA.REP.GC.GA.MET.GENES.SAT.HiFicov.ONTcov.bed \
|cut -f1,2,3,5 >compart/bTaeGut7v0.4_MT_rDNA.dot.AB.10kb.windows.bed


prefix="bTaeGut7v0.4_MT_rDNA"
AB_bed="compart/${prefix}.dot.AB.10kb.windows.bed"
module load bedtools/2.31.0

# Per chromosome
for chr in `grep dot helpfiles/bTaeGut7v0.4_MT_rDNA.groups.txt |cut -f1`;
do
  echo '#!/bin/bash
  module load bedtools/2.31.0
  rm -f compart/'$chr'_summary.10kb.txt
  for nb in "APR" "DR" "STR" "IR" "TRI" "G4" "Z" "Any"
  do
    totdens=`grep $nb coverage/'${prefix}'.per_genome.tsv |cut -f3 -d" "`
    for region in "A" "B"
    do
      summary=`grep '$chr' '$AB_bed' |awk -v OFS="\t" -v r=$region '"'"'($4==r){print $1,$2,$3}'"'"'|\
       intersectBed -a - -b final_nonB/'${prefix}'.${nb}.merged.bed -wao |cut -f1,2,3,7|\
        awk -v OFS="\t" '"'"'{if(NR==0){chr=$1; s=$2; e=$3; sum=$4}else{if($1==chr && $2==s){sum+=$4}else{print chr,s,e,sum; chr=$1; s=$2; e=$3; sum=$4}}}END{print chr,s,e,sum}'"'"' |\
         sed "/^\s*$/d" |awk -v gwd=$totdens '"'"'{sum_l+=$3-$2; sum_nb+=$4}END{d=sum_nb/sum_l; enr=d/gwd; print sum_l,sum_nb,d,enr}'"'"'`
      echo '$chr'" "$nb" "$region" "$summary >>compart/'$chr'_summary.10kb.txt
    done
  done
 ' |sbatch -J dot.$chr --ntasks=1 --cpus-per-task=1 --time=1:00:00 --mem-per-cpu=4G --account=kdm16_sc_default --partition=sla-prio --out slurm/job.AB.$chr.%j.out
done
# Merge into one file
echo "Chr NonB Compartment CompLen NonBLen Density Enrichment" |sed 's/ /\t/g' >compart/dot_summary.10kb.txt
for chr in `grep dot helpfiles/bTaeGut7v0.4_MT_rDNA.groups.txt |cut -f1`;
do

  cat compart/${chr}_summary.10kb.txt |sed 's/ /\t/g' >>compart/dot_summary.10kb.txt
done


# ~~~~~~~~~~~~~~~~~ FUNCTIONAL ENRICHMENT IN A & B COMPARTMENTS ~~~~~~~~~~~~~~~~

# Make bedfiles with A&B microdot regions
awk -v OFS="\t" '($4=="A"){print}' compart/bTaeGut7v0.4_MT_rDNA.dot.AB.10kb.windows.bed |sort -k1,1 -k2,2n |mergeBed -i - >compart/dot.A.kb.bed
awk -v OFS="\t" '($4=="B"){print}' compart/bTaeGut7v0.4_MT_rDNA.dot.AB.10kb.windows.bed |sort -k1,1 -k2,2n |mergeBed -i - >compart/dot.B.kb.bed

prefix="bTaeGut7v0.4_MT_rDNA"
for comp in "A" "B"
do
  for class in "promoter" "intergenic" "introns" "CDS" "UTR5" "UTR3" "lncrna"
  do
      echo '#!/bin/bash
      module load bedtools/2.31.0
      len=`intersectBed -a annotation/'${prefix}'.'$class'.bed -b compart/dot.'$comp'.kb.bed |sort -k1,1 -k2,2n |mergeBed -i - |awk '"'"'{sum+=$3-$2}END{print sum}'"'"'`
      echo "Length of '$class' in '$comp' compartments is $len"
      rm -f tmp.dotanalysis.'$comp'.'$class'
      cat coverage/'${prefix}'.per_genome.tsv | while read -r non_b tot dens;
      do
          echo "looking at $non_b"
          d=`intersectBed -a <(intersectBed -a annotation/'${prefix}'.'$class'.bed -b compart/dot.'$comp'.kb.bed |cut -f1-3 |sort -k1,1 -k2,2n |mergeBed -i -) -b final_nonB/'${prefix}'.${non_b}.merged.bed -wo -nonamecheck |awk -v l=$len -v dtot=$dens '"'"'{sum+=$7}END{d=sum/l; frac=d/dtot; print d,frac}'"'"'`
          echo '$comp'" "'$class'" "$non_b" "$d >>tmp.dotanalysis.'$comp'.'$class'
      done
      ' | sbatch -J dot.$class --ntasks=1 --cpus-per-task=1 --mem-per-cpu=6G --time=5:00:00 --out slurm/functional.$comp.$class.%j.out
  done
done
# Merge the tmp files
echo "Compartment Class nonB Density Enrichment_gw" |sed "s/ /\t/g" >compart/$prefix.enrichment.dot.AB.tsv
for comp in "A" "B"
do
  for class in "promoter" "intergenic" "introns" "CDS" "UTR5" "UTR3" "lncrna" 
  do
    cat tmp.dotanalysis.$comp.$class |sed "s/ /\t/g" >>compart/$prefix.enrichment.dot.AB.tsv
  done
done
#rm tmp.dotanalysis*

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CALCULATE STATS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Some calculations for paper.

# Number of genes on dot chromosomes:
# Non protein coding:
grep dot helpfiles/$prefix.groups.txt |cut -f1 |grep -f - annotation/$prefix.lncrna.bed |cut -f4 |sort |uniq |wc -l
3927
# Protein coding:
grep dot helpfiles/$prefix.groups.txt |cut -f1 |grep -f - annotation/$prefix.CDS.bed |cut -f4 |sort |uniq |wc -l
4177
# All coding genes overlapping with any of A or B compartments (dot chrom only)
grep dot helpfiles/$prefix.groups.txt |cut -f1 |grep -f - annotation/$prefix.CDS.bed |intersectBed -a - -b compart/bTaeGut7v0.4_MT_rDNA.dot.AB.10kb.bed |cut -f4 |sort |uniq |wc -l
3876
# Only A
grep dot helpfiles/$prefix.groups.txt |cut -f1 |grep -f - annotation/$prefix.CDS.bed |intersectBed -a - -b compart/dot.A.kb.bed |cut -f4 |sort |uniq |wc -l
2629
# Only B
grep dot helpfiles/$prefix.groups.txt |cut -f1 |grep -f - annotation/$prefix.CDS.bed |intersectBed -a - -b compart/dot.B.kb.bed |cut -f4 |sort |uniq |wc -l
1582
# Save gene names
grep dot helpfiles/$prefix.groups.txt |cut -f1 |grep -f - annotation/$prefix.CDS.bed |intersectBed -a - -b compart/dot.B.kb.bed |cut -f4 |sort |uniq >tmp.B.genes
grep dot helpfiles/$prefix.groups.txt |cut -f1 |grep -f - annotation/$prefix.CDS.bed |intersectBed -a - -b compart/dot.A.kb.bed |cut -f4 |sort |uniq >tmp.A.genes
# How many genes are overlapping:
join tmp.A.genes tmp.B.genes |wc -l
335



################## INTRON TRFs IN DOT CHROMOSOME COMPARTMENTS ##################


prefix="bTaeGut7v0.4_MT_rDNA"
for comp in "A" "B"
do
    echo '#!/bin/bash
    module load bedtools/2.31.0
    rm -f tmp.'$comp'.intronTRF
    for i in "TRF" "noTRF"
    do
      len=`intersectBed -a repeats/intronTRF/'$prefix'.introns_${i}.bed -b compart/dot.'$comp'.kb.bed |\
      sort -k1,1 -k2,2n |mergeBed -i - |awk '"'"'{sum+=$3-$2}END{print sum}'"'"'`
      echo "Length of $i in '$comp' compartments is $len"
      cat coverage/'${prefix}'.per_genome.tsv | while read -r non_b tot dens;
      do
        echo "looking at $non_b"
        d=`intersectBed -a <(intersectBed -a repeats/intronTRF/'$prefix'.introns_${i}.bed -b compart/dot.'$comp'.kb.bed |\
        sort -k1,1 -k2,2n |mergeBed -i -) -b final_nonB/'${prefix}'.${non_b}.merged.bed -wo |\
        awk -v l=$len -v dtot=$dens '"'"'{sum+=$7}END{d=sum/l; frac=d/dtot; print d,frac}'"'"'`
        echo '$comp' $i $non_b $d >>tmp.'$comp'.intronTRF
      done
    done
    ' | sbatch -J $comp --ntasks=1 --cpus-per-task=1 --mem-per-cpu=6G --time=5:00:00 --out slurm/job.IntronTR.$comp.%j.out
done
# Merge the tmp files
echo "Compartment Subset NonB Density Enrichment_gw" |sed "s/ /\t/g" >repeats/$prefix.intron_enrichment.compartment.tsv
for comp in "A" "B"
do
  cat tmp.$comp.intronTRF |sed "s/ /\t/g"
done >>repeats/$prefix.intron_enrichment.compartment.tsv
#rm tmp.*

# Checking number of TRF repeats in introns in A and B compartments
for comp in "A" "B"
do
 intersectBed -a repeats/intronTRF/$prefix.introns_TRF.bed -b compart/dot.$comp.kb.bed | cut -f4 |\
 sort |uniq -c |awk -v OFS="\t" -v g=$comp '{print g,$2,$1}' >>repeats/$prefix.introns_TR.repeatsummary.tsv
done

# And exact length of repeats
echo -e "Comp\tRepUnitLen" >repeats/$prefix.introns_TR.dot.compart.lengths.tsv
for comp in "A" "B"
do
  intersectBed -a repeats/intronTRF/$prefix.introns_TRF.bed -b compart/dot.$comp.kb.bed |\
  awk -v c=$comp '{gsub(/-mer/, ""); print c"\t"$5}' >>repeats/$prefix.introns_TR.dot.compart.lengths.tsv
done 

# Add compartment lengths to the file initiated in 04_functional_enrichment.sh
for comp in "A" "B"
do
  grep "Length of" slurm/job.IntronTR.$comp.487*.out |tail -n2 |awk -v OFS="\t" '{print $5,$3,$8}' |uniq >>repeats/$prefix.introns_TR.lengthsummary.tsv
done
