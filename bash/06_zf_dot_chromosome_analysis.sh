################################################################################
### CALCULATING NON-B DNA MOTIF ENRICHMENT IN DOT CHROMOSOME A AND B COMPARTMENTS
### CODE WRITTEN BY LINNÉA SMEDS

##### TABLE OF CONTENTS
# CALCULATE ENRICHMENT
#    -FUNCTIONAL ENRICHMENT IN A & B COMPARTMENTS
#    -CALCULATE STATS
# INTRON TRFs IN DOT CHROMOSOME COMPARTMENTS

############################ CALCULATE ENRICHMENT ##############################
# 16 July 2025

mkdir compart

# A and B compartment annotation for 10 kb windows processed by Simona Secomandi
# is found in: annotation/bTaeGut7v0.4_MT_rDNA.Cooltools.v0.2.E1.10Kb.flipped.dip.collated.AB.bed

# Extract only A and B compartments from dot chromosomes
grep -f helpfiles/microdot.txt annotation/bTaeGut7v0.4_MT_rDNA.Cooltools.v0.2.E1.10Kb.flipped.dip.collated.AB.bed >compart/bTaeGut7v0.4_MT_rDNA.dot.AB.10kb.bed
# Extract all windows from the new file from Simonas
grep -f microdot.txt ref/bTaeGut7v0.4_MT_rDNA.Cooltools.v0.2.E1.10Kb.flipped.dip.collated.AB.noNA.REP.GC.GA.MET.GENES.SAT.HiFicov.ONTcov.bed |cut -f1,2,3,5 >compart/bTaeGut7v0.4_MT_rDNA.dot.AB.10kb.windows.bed


prefix="bTaeGut7v0.4_MT_rDNA"
AB_bed="compart/${prefix}.dot.AB.10kb.windows.bed"
module load bedtools/2.31.0

# Per chromosome
cat helpfiles/microdot.txt |while read -r chr;
do
  echo '#!/bin/bash
  module load bedtools/2.31.0
  rm -f compart/'$chr'_summary.10kb.txt
  for nb in "APR" "DR" "STR" "IR" "MR" "TRI" "G4" "Z" "ALL"
  do
    totdens=`grep $nb densities/'${prefix}'.nonB_genome_wide.txt |cut -f3 -d" "`
    for region in "A" "B"
    do
      summary=`grep '$chr' '$AB_bed' |awk -v OFS="\t" -v r=$region '"'"'($4==r){print $1,$2,$3}'"'"'| intersectBed -a - -b final_nonB/'${prefix}'.${nb}.merged.bed -wao |cut -f1,2,3,7| awk -v OFS="\t" '"'"'{if(NR==0){chr=$1; s=$2; e=$3; sum=$4}else{if($1==chr && $2==s){sum+=$4}else{print chr,s,e,sum; chr=$1; s=$2; e=$3; sum=$4}}}END{print chr,s,e,sum}'"'"' | sed "/^\s*$/d" |awk -v gwd=$totdens '"'"'{sum_l+=$3-$2; sum_nb+=$4}END{d=sum_nb/sum_l; enr=d/gwd; print sum_l,sum_nb,d,enr}'"'"'`
      echo '$chr'" "$nb" "$region" "$summary >>compart/'$chr'_summary.10kb.txt
    done
  done
 ' |sbatch -J $chr --ntasks=1 --cpus-per-task=1 --time=1:00:00 --mem-per-cpu=4G --out slurm/job.AB.$chr.%j.out
done
# Merge into one file
echo "Chr NonB Compartment CompLen NonBLen Density Enrichment" |sed 's/ /\t/g' >compart/microdot_summary.10kb.txt
cat microdot.txt |while read -r chr;
do
  cat compart/${chr}_summary.10kb.txt |sed 's/ /\t/g' >>compart/microdot_summary.10kb.txt
done


# ~~~~~~~~~~~~~~~~~ FUNCTIONAL ENRICHMENT IN A & B COMPARTMENTS ~~~~~~~~~~~~~~~~
# July 25, 2025

# Make bedfiles with A&B microdot regions
awk -v OFS="\t" '($5=="A"){print $1,$3,$4,$5}' compart/microdot_regions.10kb.txt |sort -k1,1 -k2,2n |uniq >compart/microdot.A.kb.bed
awk -v OFS="\t" '($5=="B"){print $1,$3,$4,$5}' compart/microdot_regions.10kb.txt |sort -k1,1 -k2,2n |uniq >compart/microdot.B.kb.bed

prefix="bTaeGut7v0.4_MT_rDNA"
for comp in "A" "B"
do
  for class in "promoter"  "intergenic"  "intronic" "CDS" "UTR5" "UTR3" "nonprotcoding" 
  do
      echo '#!/bin/bash
      module load bedtools/2.31.0
      len=`intersectBed -a annotation/'$class'.v3.bed -b compart/microdot.'$comp'.kb.bed |sort -k1,1 -k2,2n |mergeBed -i - |awk '"'"'{sum+=$3-$2}END{print sum}'"'"'`
      echo "Length of '$class' in '$comp' compartments is $len"
      rm -f tmp.'$comp'.'$class'
      cat densities/'${prefix}'.nonB_genome_wide.txt | while read -r non_b tot dens;
      do
          echo "looking at $non_b"
          d=`intersectBed -a <(intersectBed -a annotation/'$class'.bed -b compart/microdot.'$comp'.kb.bed |cut -f1-3 |sort -k1,1 -k2,2n |mergeBed -i -) -b final_nonB/'${prefix}'.${non_b}.merged.bed -wo |awk -v l=$len -v dtot=$dens '"'"'{sum+=$7}END{d=sum/l; frac=d/dtot; print d,frac}'"'"'`
          echo '$comp'" "'$class'" "$non_b" "$d >>tmp.'$comp'.'$class'
      done
      ' | sbatch -J $class --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8G --time=5:00:00 --partition=open --out slurm/functional.$comp.$class.%j.out
  done
done
# Merge the tmp files
echo "Compartment Class nonB Density Enrichment_gw" |sed "s/ /\t/g" >functional/enrichment_microdot.AB.v4.tsv
for comp in "A" "B"
do
  for class in "intronic" "promoter"  "intergenic" "CDS" "UTR5" "UTR3"  "nonprotcoding" #"all_repeats" #"pseudogenes"
  do
    cat tmp.$comp.$class |sed "s/ /\t/g" >>functional/enrichment_microdot.AB.v4.tsv
  done
done
rm tmp.*

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CALCULATE STATS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Some calculations for paper.

# Number of genes on dot chromosomes:
# Non protein coding:
grep -f microdot.txt annotation/nonprotcoding.bed |cut -f4 |sort |uniq |wc -l
4041
# Protein coding:
grep -f microdot.txt annotation/CDS.bed |cut -f4 |sort |uniq |wc -l
4193
# All coding genes overlapping with any of A or B compartments (dot chrom only)
grep -f microdot.txt annotation/CDS.bed |intersectBed -a - -b compart/bTaeGut7v0.4_MT_rDNA.dot.AB.10kb.bed |cut -f4 |sort |uniq |wc -l
3885
# Only A
grep -f microdot.txt annotation/CDS.bed |intersectBed -a - -b compart/microdot.A.kb.bed |cut -f4 |sort |uniq |wc -l
2533
# Only B
grep -f microdot.txt annotation/CDS.bed |intersectBed -a - -b compart/microdot.B.kb.bed |cut -f4 |sort |uniq |wc -l
1513
# Save gene names
grep -f microdot.txt annotation/CDS.bed |intersectBed -a - -b compart/microdot.B.kb.bed |cut -f4 |sort |uniq >tmp.B.genes
grep -f microdot.txt annotation/CDS.bed |intersectBed -a - -b compart/microdot.A.kb.bed |cut -f4 |sort |uniq >tmp.A.genes
# How many genes are overlapping:
join tmp.A.genes tmp.B.genes |wc
    328     328   14761


################## INTRON TRFs IN DOT CHROMOSOME COMPARTMENTS ##################


prefix="bTaeGut7v0.4_MT_rDNA"
for comp in "A" "B"
do
    echo '#!/bin/bash
    module load bedtools/2.31.0
    rm -f tmp.'$comp'.intronTRF
    for i in "TRF" "noTRF"
    do
      len=`intersectBed -a repeats/intronTRF/introns_${i}.bed -b compart/microdot.'$comp'.kb.bed |sort -k1,1 -k2,2n |mergeBed -i - |awk '"'"'{sum+=$3-$2}END{print sum}'"'"'`
      echo "Length of $i in '$comp' compartments is $len"
      cat densities/'${prefix}'.nonB_genome_wide.txt | while read -r non_b tot dens;
      do
        echo "looking at $non_b"
        d=`intersectBed -a <(intersectBed -a repeats/intronTRF/introns_${i}.bed -b compart/microdot.'$comp'.kb.bed |sort -k1,1 -k2,2n |mergeBed -i -) -b final_nonB/'${prefix}'.${non_b}.merged.bed -wo |awk -v l=$len -v dtot=$dens '"'"'{sum+=$7}END{d=sum/l; frac=d/dtot; print d,frac}'"'"'`
        echo '$comp' $i $non_b $d >>tmp.'$comp'.intronTRF
      done
    done
    ' | sbatch -J $comp --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8G --time=5:00:00 --partition=open --out slurm/job.intron-enrich.$comp.%j.out
done
# Merge the tmp files
echo "Compartment Subset NonB Density Enrichment_gw" |sed "s/ /\t/g" >repeats/compartment_intron_enrichment.tsv
for comp in "A" "B"
do
  cat tmp.$comp.intronTRF |sed "s/ /\t/g"
done >> repeats/compartment_intron_enrichment.tsv
rm tmp.*

# Checking number of TRF repeats in introns in A and B compartments
for comp in "A" "B"
do
 intersectBed -a repeats/intronTRF/introns_TRF.bed -b compart/microdot.$comp.kb.bed | cut -f4 |sort |uniq -c |awk -v OFS="\t" -v g=$comp '{print g,$2,$1}' >>repeats/intron_TR_repeatsummary.tsv
done
