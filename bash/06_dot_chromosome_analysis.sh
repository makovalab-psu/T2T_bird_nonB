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

# A and B compartment annotation for 200 kb windows is downloaded from 
wget https://genomeark.s3.amazonaws.com/species/Taeniopygia_guttata/bTaeGut7/manuscript/annotations/3D/bTaeGut7v0.4_MT_rDNA.Cooltools.E1.200kbp.flipped.dip.collated.v0.1.bed
# (place in ref/) 

# Extract only A and B compartments from dot chromosomes
grep dot helpfiles/bTaeGut7v0.4_MT_rDNA.groups.txt |cut -f1 |\
grep -f - ref/bTaeGut7v0.4_MT_rDNA.Cooltools.E1.200kbp.flipped.dip.collated.v0.1.bed >compart/bTaeGut7v0.4_MT_rDNA.dot.AB.200kb.bed

# Extract all windows from the new file from Simonas
#grep dot helpfiles/bTaeGut7v0.4_MT_rDNA.groups.txt |cut -f1 |\
#grep -f - ref/bTaeGut7v0.4_MT_rDNA.Cooltools.v0.2.E1.10Kb.flipped.dip.collated.AB.noNA.REP.GC.GA.MET.GENES.SAT.HiFicov.ONTcov.bed \
#|cut -f1,2,3,5 >compart/bTaeGut7v0.4_MT_rDNA.dot.AB.10kb.windows.bed


prefix="bTaeGut7v0.4_MT_rDNA"
AB_bed="compart/${prefix}.dot.AB.200kb.bed"
module load bedtools/2.31.0

# Per chromosome
for chr in `grep dot helpfiles/bTaeGut7v0.4_MT_rDNA.groups.txt |cut -f1`;
do
  echo '#!/bin/bash
  module load bedtools/2.31.0
  rm -f compart/'$chr'_summary.200kb.txt
  for nb in "APR" "DR" "STR" "IR" "TRI" "G4" "Z" "Any"
  do
    totdens=`grep $nb coverage/'${prefix}'.per_genome.tsv |cut -f3 -d" "`
    for region in "A" "B"
    do
      summary=`grep '$chr' '$AB_bed' |awk -v OFS="\t" -v r=$region '"'"'($4==r){print $1,$2,$3}'"'"'|\
       intersectBed -a - -b final_nonB/'${prefix}'.${nb}.merged.bed -wao |cut -f1,2,3,7|\
        awk -v OFS="\t" '"'"'{if(NR==0){chr=$1; s=$2; e=$3; sum=$4}else{if($1==chr && $2==s){sum+=$4}else{print chr,s,e,sum; chr=$1; s=$2; e=$3; sum=$4}}}END{print chr,s,e,sum}'"'"' |\
         sed "/^\s*$/d" |awk -v gwd=$totdens '"'"'{sum_l+=$3-$2; sum_nb+=$4}END{d=sum_nb/sum_l; enr=d/gwd; print sum_l,sum_nb,d,enr}'"'"'`
      echo '$chr'" "$nb" "$region" "$summary >>compart/'$chr'_summary.200kb.txt
    done
  done
 ' |sbatch -J dot.$chr --ntasks=1 --cpus-per-task=1 --time=1:00:00 --mem-per-cpu=4G --out slurm/job.AB.$chr.%j.out
done
# Merge into one file
echo "Chr NonB Compartment CompLen NonBLen Density Enrichment" |sed 's/ /\t/g' >compart/dot_summary.200kb.txt
for chr in `grep dot helpfiles/bTaeGut7v0.4_MT_rDNA.groups.txt |cut -f1`;
do

  cat compart/${chr}_summary.200kb.txt |sed 's/ /\t/g' >>compart/dot_summary.200kb.txt
done


# ~~~~~~~~~~~~~~~~~ FUNCTIONAL ENRICHMENT IN A & B COMPARTMENTS ~~~~~~~~~~~~~~~~

# Make bedfiles with A&B microdot regions
awk -v OFS="\t" '($4=="A"){print}' compart/bTaeGut7v0.4_MT_rDNA.dot.AB.200kb.bed |sort -k1,1 -k2,2n |mergeBed -i - >compart/dot.A.200kb.bed
awk -v OFS="\t" '($4=="B"){print}' compart/bTaeGut7v0.4_MT_rDNA.dot.AB.200kb.bed |sort -k1,1 -k2,2n |mergeBed -i - >compart/dot.B.200kb.bed

prefix="bTaeGut7v0.4_MT_rDNA"
for comp in "A" "B"
do
  for class in "promoter" "intergenic" "introns" "CDS" "UTR5" "UTR3" "lncrna"
  do
      echo '#!/bin/bash
      module load bedtools/2.31.0
      len=`intersectBed -a annotation/'${prefix}'.'$class'.bed -b compart/dot.'$comp'.200kb.bed |sort -k1,1 -k2,2n |mergeBed -i - |awk '"'"'{sum+=$3-$2}END{print sum}'"'"'`
      echo "Length of '$class' in '$comp' compartments is $len"
      rm -f tmp.dotanalysis.'$comp'.'$class'
      cat coverage/'${prefix}'.per_genome.tsv | while read -r non_b tot dens;
      do
          echo "looking at $non_b"
          d=`intersectBed -a <(intersectBed -a annotation/'${prefix}'.'$class'.bed -b compart/dot.'$comp'.200kb.bed |cut -f1-3 |sort -k1,1 -k2,2n |mergeBed -i -) -b final_nonB/'${prefix}'.${non_b}.merged.bed -wo -nonamecheck |awk -v l=$len -v dtot=$dens '"'"'{sum+=$7}END{d=sum/l; frac=d/dtot; print d,frac}'"'"'`
          echo '$comp'" "'$class'" "$non_b" "$d >>tmp.dotanalysis.'$comp'.'$class'
      done
      ' | sbatch -J dot.$class --ntasks=1 --cpus-per-task=1 --mem-per-cpu=6G --time=5:00:00 --out slurm/functional.$comp.$class.%j.out
  done
done
# Merge the tmp files
echo "Compartment Class nonB Density Enrichment_gw" |sed "s/ /\t/g" >compart/$prefix.enrichment.dot.AB.200kb.tsv
for comp in "A" "B"
do
  for class in "promoter" "intergenic" "introns" "CDS" "UTR5" "UTR3" "lncrna" 
  do
    cat tmp.dotanalysis.$comp.$class |sed "s/ /\t/g" >>compart/$prefix.enrichment.dot.AB.200kb.tsv
  done
done
#rm tmp.dotanalysis*

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CALCULATE STATS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Some calculations for paper.

# Number of genes on dot chromosomes:
# Non protein coding:
grep dot helpfiles/$prefix.groups.txt |cut -f1 |grep -f - annotation/$prefix.lncrna.bed |cut -f4 |sort |uniq |wc -l
#3927
# Protein coding:
grep dot helpfiles/$prefix.groups.txt |cut -f1 |grep -f - annotation/$prefix.CDS.bed |cut -f4 |sort |uniq |wc -l
#4177
# All coding genes overlapping with any of A or B compartments (dot chrom only)
grep dot helpfiles/$prefix.groups.txt |cut -f1 |grep -f - annotation/$prefix.CDS.bed |intersectBed -a - -b compart/bTaeGut7v0.4_MT_rDNA.dot.AB.200kb.bed |cut -f4 |sort |uniq |wc -l
#4010
# Only A
grep dot helpfiles/$prefix.groups.txt |cut -f1 |grep -f - annotation/$prefix.CDS.bed |intersectBed -a - -b compart/dot.A.200kb.bed |cut -f4 |sort |uniq |wc -l
#2414
# Only B
grep dot helpfiles/$prefix.groups.txt |cut -f1 |grep -f - annotation/$prefix.CDS.bed |intersectBed -a - -b compart/dot.B.200kb.bed |cut -f4 |sort |uniq |wc -l
#1614
# Save gene names
grep dot helpfiles/$prefix.groups.txt |cut -f1 |grep -f - annotation/$prefix.CDS.bed |intersectBed -a - -b compart/dot.B.200kb.bed |cut -f4 |sort |uniq >tmp.B.genes
grep dot helpfiles/$prefix.groups.txt |cut -f1 |grep -f - annotation/$prefix.CDS.bed |intersectBed -a - -b compart/dot.A.200kb.bed |cut -f4 |sort |uniq >tmp.A.genes
# How many genes are overlapping:
join tmp.A.genes tmp.B.genes |wc -l
#18

# Size of A and B compartments:
awk '{sum+=$3-$2}END{print sum}' compart/dot.A.200kb.bed
#43200000
awk '{sum+=$3-$2}END{print sum}' compart/dot.B.200kb.bed
#67600000

# Number of genes per compartment bp (density)
#2414/43200000=5.59e-5
#1614/67600000=2.39e-5

# "exon coverage" (how much of the compartment is exonic)
awk '{sum+=$3-$2}END{print sum}' <(intersectBed -a annotation/$prefix.CDS.bed -b compart/dot.A.200kb.bed |sort -k1,1 -k2,2n |mergeBed -i -)
#3966444
awk '{sum+=$3-$2}END{print sum}' <(intersectBed -a annotation/$prefix.CDS.bed -b compart/dot.B.200kb.bed |sort -k1,1 -k2,2n |mergeBed -i -)
#2113269
#
#3966444/43200000=0.0918
#2113269/67600000=0.0313


# "gene coverage" (how much of the compartment is genic)
awk -F"\t" '($3=="transcript" && $9~/transcript_biotype "mRNA"/){s=$4-1; print $1"\t"s"\t"$5}' ref/bTaeGut7v0.4_MT_rDNA.longest_only.gff |intersectBed -a - -b compart/dot.A.200kb.bed |sort -k1,1 -k2,2n  |mergeBed -i - |awk '{sum+=$3-$2}END{print sum}'
#31434230
awk -F"\t" '($3=="transcript" && $9~/transcript_biotype "mRNA"/){s=$4-1; print $1"\t"s"\t"$5}' ref/bTaeGut7v0.4_MT_rDNA.longest_only.gff |intersectBed -a - -b compart/dot.B.200kb.bed |sort -k1,1 -k2,2n  |mergeBed -i - |awk '{sum+=$3-$2}END{print sum}'
#13045867

#31434230/43200000=0.727
#13045867/67600000=0.193



################## INTRON TRFs IN DOT CHROMOSOME COMPARTMENTS ##################


prefix="bTaeGut7v0.4_MT_rDNA"
for comp in "A" "B"
do
    echo '#!/bin/bash
    module load bedtools/2.31.0
    rm -f tmp.'$comp'.intronTRF
    for i in "TRF" "noTRF"
    do
      len=`intersectBed -a repeats/intronTRF/'$prefix'.introns_${i}.bed -b compart/dot.'$comp'.200kb.bed |\
      sort -k1,1 -k2,2n |mergeBed -i - |awk '"'"'{sum+=$3-$2}END{print sum}'"'"'`
      echo "Length of $i in '$comp' compartments is $len"
      cat coverage/'${prefix}'.per_genome.tsv | while read -r non_b tot dens;
      do
        echo "looking at $non_b"
        d=`intersectBed -a <(intersectBed -a repeats/intronTRF/'$prefix'.introns_${i}.bed -b compart/dot.'$comp'.200kb.bed |\
        sort -k1,1 -k2,2n |mergeBed -i -) -b final_nonB/'${prefix}'.${non_b}.merged.bed -wo |\
        awk -v l=$len -v dtot=$dens '"'"'{sum+=$7}END{d=sum/l; frac=d/dtot; print d,frac}'"'"'`
        echo '$comp' $i $non_b $d >>tmp.'$comp'.intronTRF
      done
    done
    ' | sbatch -J $comp --ntasks=1 --cpus-per-task=1 --mem-per-cpu=6G --time=5:00:00 --out slurm/job.IntronTR.$comp.%j.out
done
# Merge the tmp files
echo "Compartment Subset NonB Coverage Enrichment_gw" |sed "s/ /\t/g" >repeats/$prefix.intron_enrichment.compartment.tsv
for comp in "A" "B"
do
  cat tmp.$comp.intronTRF |sed "s/ /\t/g"
done >>repeats/$prefix.intron_enrichment.compartment.tsv
#rm tmp.*

# Checking number of TRF repeats in introns in A and B compartments
for comp in "A" "B"
do
 intersectBed -a repeats/intronTRF/$prefix.introns_TRF.bed -b compart/dot.$comp.200kb.bed | cut -f4 |\
 sort |uniq -c |awk -v OFS="\t" -v g=$comp '{print g,$2,$1}' >>repeats/$prefix.introns_TR.repeatsummary.tsv
done

# And exact length of repeats
echo -e "Comp\tRepUnitLen" >repeats/$prefix.introns_TR.dot.compart.lengths.tsv
for comp in "A" "B"
do
  intersectBed -a repeats/intronTRF/$prefix.introns_TRF.bed -b compart/dot.$comp.200kb.bed |\
  awk -v c=$comp '{gsub(/-mer/, ""); print c"\t"$5}' >>repeats/$prefix.introns_TR.dot.compart.lengths.tsv
done 

# Add compartment lengths to the file initiated in 04_functional_enrichment.sh
for comp in "A" "B"
do
  grep "Length of" slurm/job.IntronTR.$comp.498*.out |tail -n2 |awk -v OFS="\t" '{print $5,$3,$8}' |uniq >>repeats/$prefix.introns_TR.lengthsummary.tsv
done


# Dot chromosome fraction of the genome 
awk '{sum+=$2}END{print sum}' helpfiles/bTaeGut7v0.4_MT_rDNA.groups.txt 
2162893666
awk '($3=="dot"){sum+=$2}END{print sum}' helpfiles/bTaeGut7v0.4_MT_rDNA.groups.txt 
131816572