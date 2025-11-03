################################################################################
### ANNOTATION OF NON-B DNA IN THE ZEBRA FINCH T2T GENOME, AND OVERLAPS 
### BETWEEN THE DIFFERENT MOTIFS. 
### CODE WRITTEN BY LINNÉA SMEDS

# This file contains the first code, including downloading the zebra finch 
# genome and running the annotation software. 

# All commands are assumed to be run from within the downloaded github repo.
# Time and/or memory consuming jobs were run on a slurm cluster. The code 
# within the slurm scripts can be copied and run on any machine with enough
# RAM/cores.

##### TABLE OF CONTENTS
# REQUIREMENTS 
# ANNOTATION WITH GFA
# ANNOTATION WITH QUADRON
# CREATE A JOINT DIRECTORY WITH ANNOTATIONS
# OVERLAP BETWEEN NON-B TYPES
# CHECKING NCBI FOR BIRD GENOMES 

################################# REQUIREMENTS #################################
# A LIST OF ALL SOFTWARE USED:
# bedtools/2.31.0
# gfa (downloaded Feb 4, 2025 from https://github.com/abcsFrederick/non-B_gfa)
# Quadron (dockerized version 1.0.0 from https://hub.docker.com/r/kxk302/quadron)
# circos/0.69-9
# python/3.11.2
# convert2bed/2.4.41

# Download Zebra Finch genome, centromere annotation, repeat annotation, genes,
#  from GenomeArk
mkdir ref
cd ref
wget https://genomeark.s3.amazonaws.com/species/Taeniopygia_guttata/bTaeGut7/assembly_verkko_0.1/manual_curation/bTaeGut7v0.4/bTaeGut7v0.4_MT_rDNA.fa.gz
wget https://genomeark.s3.amazonaws.com/species/Taeniopygia_guttata/bTaeGut7/manuscript/annotations/centromeres/bTaeGut7v0.4_MT_rDNA.centromere_detector.v0.1.gff
wget https://genomeark.s3.amazonaws.com/species/Taeniopygia_guttata/bTaeGut7/manuscript/annotations/repeats/bTaeGut7v0.4_MT_rDNA.EDTA2.v0.2.gtf.gz
wget https://genomeark.s3.amazonaws.com/species/Taeniopygia_guttata/bTaeGut7/manuscript/annotations/repeats/bTaeGut7v0.4_MT_rDNA.trf.sorted.v0.1.bed
wget https://genomeark.s3.amazonaws.com/species/Taeniopygia_guttata/bTaeGut7/manuscript/annotations/repeats/bTaeGut7v0.4_MT_rDNA.satellome.v0.1.bed
wget https://genomeark.s3.amazonaws.com/species/Taeniopygia_guttata/bTaeGut7/manuscript/annotations/genes/bTaeGut7v0.4_MT_rDNA.EGAPx.v0.1.gtf.gz
wget https://genomeark.s3.amazonaws.com/species/Taeniopygia_guttata/bTaeGut7/manuscript/annotations/genes/bTaeGut7v0.4.v0.1.annotation.modified.longest_isoform.only_exon.gtf
wget https://genomeark.s3.amazonaws.com/species/Taeniopygia_guttata/bTaeGut7/manuscript/annotations/genes/bTaeGut7v0.4.v0.1.annotation.modified.longest_isoform.only_gene.gtf
wget https://genomeark.s3.amazonaws.com/species/Taeniopygia_guttata/bTaeGut7/manuscript/annotations/genes/bTaeGut7v0.4.v0.1.annotation.modified.longest_isoform.only_intron.gtf
wget https://genomeark.s3.amazonaws.com/species/Taeniopygia_guttata/bTaeGut7/manuscript/annotations/methylation/bTaeGut7v0.4_MT_rDNA.PBmethylation.v0.1.bw
wget https://genomeark.s3.amazonaws.com/species/Taeniopygia_guttata/bTaeGut7/manuscript/annotations/3D/bTaeGut7v0.4_MT_rDNA.Cooltools.E1.200kbp.flipped.dip.collated.v0.1.bed
wget https://genomeark.s3.amazonaws.com/species/Taeniopygia_guttata/bTaeGut7/manuscript/annotations/3D/bTaeGut7v0.4_MT_rDNA.Cooltools.E1.10kbp.flipped.dip.collated.v0.1.bw
wget https://genomeark.s3.amazonaws.com/species/Taeniopygia_guttata/bTaeGut7/manuscript/annotations/PUR/bTaeGut7v0.4_MT_rDNA_MATonly_vs_GCF_000151805.1.PUR.fastga.v0.1.bed
wget https://genomeark.s3.amazonaws.com/species/Taeniopygia_guttata/bTaeGut7/assembly_verkko_0.1/manual_curation/bTaeGut7v0.4/mapping/v0.4_dip_hifi/v0.4_dip_hifi.pri.cov.wig
wget https://genomeark.s3.amazonaws.com/species/Taeniopygia_guttata/bTaeGut7/assembly_verkko_0.1/manual_curation/bTaeGut7v0.4/mapping/v0.4_dip_ont/v0.4_dip_ont.pri.cov.wig

gunzip bTaeGut7v0.4_MT_rDNA.fa.gz
cd ..


############################## ANNOTATION WITH GFA #############################
# Annnotating A-Phased repeats (APR), Direct Repeats (DR), Inverted Repeats (IR),
# Mirror Repeats (MR), Short Tandem Repeats (STR) and Z-DNA (Z) motifs with GFA
# June 12, 2025
mkdir slurm
mkdir gfa_annotation

# GFA is assumed to be installed in: 
# ~/software/non-B_gfa/gfa

# Set parameters
fasta="ref/bTaeGut7v0.4_MT_rDNA.fa"
prefix="bTaeGut7v0.4_MT_rDNA"
# Run gfa on local scratch, convert to bed and copy files back to storage
echo '#!/bin/bash
echo "##### Start gfa"
~/software/non-B_gfa/gfa -seq '$fasta' -out gfa_annotation/'$prefix'  -skipGQ -skipWGET
echo "##### translate to bed"
for type in "APR" "DR" "IR" "MR" "STR" "Z"
do
  echo "#####....translating $type...."
  awk -v OFS="\t" '"'"'(NR>1){s=$4-1; print $1,s,$5,$9":"$10":"$11,".",$8}'"'"' gfa_annotation/'${prefix}'_${type}.tsv >gfa_annotation'${prefix}'.$type.bed
done
awk -v OFS="\t" '"'"'($12==1){s=$4-1; print $1,s,$5,$9":"$10":"$11,".",$8}'"'"' gfa_annotation/'${prefix}'_MR.tsv  >/gfa_annotation/'${prefix}'.TRI.bed
' |sbatch -J gfa.$prefix --ntasks=1 --cpus-per-task=1 --mem-per-cpu=2G --time=2-00:00:00 -o slurm/job.gfa.zf.%j.%N.out


########################### ANNOTATION WITH QUADRON ############################
# Annotating G-quadruplex (G4) motifs with Quadron
# June 12, 2025

# The dockerized version of Quadron was converted into a singularity container 
# using 
singularity build quadron.sif docker://kxk302/quadron:1.0.0
mv quadron.sif software/quadron.sif

mkdir Quadron_annotation
# Divide fasta into chromosomes
mkdir Quadron_annotation/split_fasta
echo '#!/bin/bash
faidx --split-files ref/bTaeGut7v0.4_MT_rDNA.fa
mv *.fa Quadron_annotation/split_fasta/
' |sbatch -J $chr --ntasks=1 --cpus-per-task=4 --time=5:00:00 -o slurm/job.splitfasta.%j.out

# Run Quadron on each chromosome separately
# Quadron singularity image is assumed to be in ~/software/quadron.sif

#Note that Quadron is requiring a full path to the fasta file 
for fa in $(ls $PWD/Quadron_annotation/split_fasta/*.fa)
do
  chr=`basename $fa |sed 's/.fa//'`
  echo $chr
  out="$PWD/Quadron_annotation/$chr.out"
  echo '#!/bin/bash
  singularity exec software/quadron.sif sh -c "cd /scripts/; Rscript run_quadron.R '$fa' '$out' 2 1000000"
  ' |sbatch -J $chr --ntasks=1 --cpus-per-task=2 --time=1:00:00 -o slurm/job.quadron.$chr.%j.out
done

# If Quadron get's too little memory, it can sometimes stop calculating the 
# score, but still output the motif. Check how many "NA" scores there are 
# per file (a couple close to the ends are ok - Quadron can't calculate a 
# score if the flanking sequence is shorter than 50bp)
for file in $(ls Quadron_annotation/chr*.out)
do
  na=`grep "^DATA:" $file |awk '($5=="NA"){n++}END{print n}'`
  echo $file $na
done
# Looks ok!

# Convert to bed
prefix="bTaeGut7v0.4_MT_rDNA"
rm -f Quadron_annotation/$prefix.G4.bed
echo $prefix
for chr in $(cut -f1 ref/$prefix.fa.fai)
do
  awk -v chr=$chr -v OFS="\t" '(/^DATA/){if($5=="NA"){}else{score=$5; if(score<0){score=0}; s=$2-1; e=s+$4; print chr,s,e,$5,sprintf("%.f", score),$3}}' Quadron_annotation/$chr.out >>Quadron_annotation/$prefix.G4.bed
done


################## CREATE A JOINT DIRECTORY WITH ANNOTATIONS ###################
mkdir final_nonB
prefix="bTaeGut7v0.4_MT_rDNA"
cp gfa_annotation/bTaeGut7v0.4_MT_rDNA.*.bed final_nonB/
cp Quadron_annotation/bTaeGut7v0.4_MT_rDNA.G4.bed final_nonB/

# AND CREATE NEW NONB FILES WITH MERGED REGIONS (THAT CONTAIN NO OVERLAPS)
for n in "G4" "APR" "DR" "IR" "MR" "TRI" "STR" "Z" "ALL"
do
 echo '#!/bin/bash
 module load bedtools/2.31.0
 sort -k1,1 -k2,2n final_nonB/'$prefix'.'${n}'.bed |mergeBed -i - >final_nonB/'$prefix'.'${n}'.merged.bed
 '| sbatch -J $n --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8G --out slurm/mergeBed.$n.%j.out
 done


######################### OVERLAP BETWEEN NON-B TYPES ##########################
# Check overlap between all non-B types to be visulized with an upset plot

mkdir overlap
prefix="bTaeGut7v0.4_MT_rDNA"
# Check overlap for each chromosome separately
cat ref/${prefix}*.fa.fai |cut -f1 |tail -n+2 |while read -r chr;
do
  echo '#!/bin/bash
  echo "Running upset for '$chr'"
  module load python/3.11.2
  python3 python/upset_summary.py -b final_nonB/'$prefix'. -c '$chr' -o overlap/summary.'$chr'.txt
  ' | sbatch -J $chr --ntasks=1 --mem-per-cpu=6G --cpus-per-task=1 --time=1:00:00
done
# The longest chr2 took 3Gb of RAM and less than a minute to run.

# Merge autosomes
cat overlap/summary.chr[0-9]*.txt |sort | awk -v OFS="\t" '{if(NR==1){type=$1; sum=$2}else{if($1==type){sum+=$2}else{print type,sum; type=$1; sum=$2}}}END{print type,sum}' > overlap/summary.autosomes.txt
# And rDNA
cat overlap/summary.rDNA*.txt |sort | awk -v OFS="\t" '{if(NR==1){type=$1; sum=$2}else{if($1==type){sum+=$2}else{print type,sum; type=$1; sum=$2}}}END{print type,sum}' > overlap/summary.rDNA.txt
# And chromosome types 
for type in "macro" "micro" "microdot"
do
  for chr in $(cat helpfiles/$type.txt)
  do
    cat overlap/summary.$chr.txt
  done | sort | awk -v OFS="\t" '{if(NR==1){type=$1; sum=$2}else{if($1==type){sum+=$2}else{print type,sum; type=$1; sum=$2}}}END{print type,sum}' > overlap/summary.$type.txt
done

# PAIRWISE OVERLAP WITH FRACTION
# To get the fraction I first need the totals for all types
mkdir stats/nonB_totals
for file in $(ls overlap/summary.*.txt)
do
  chr=`echo $file | cut -f2 -d"."`
  echo $chr
  echo "#NonB bp" |sed "s/ /\t/g" >stats/nonB_totals/annotated.${chr}.txt
  for nb in "APR" "DR" "G4" "IR" "MR" "TRI" "STR" "Z"
  do
    grep $nb $file |awk -v nb=$nb '{sum+=$2}END{print nb"\t"sum}' >>stats/nonB_totals/annotated.${chr}.txt
  done
done

# Make bed files for each chr type (this is just to simplify the next step)
mkdir final_nonB/merged_per_chrom/
cat ref/${prefix}*.fa.fai |cut -f1 |while read -r chr;
do
  echo $chr
  for nb in "APR" "DR" "G4" "IR" "MR" "TRI" "STR" "Z"
  do
    grep $chr final_nonB/$prefix.$nb.bed |sort -k1,1 -k2,2n |mergeBed -i - >final_nonB/merged_per_chrom/$chr.$nb.bed
  done
done
# Autosomes & rDNA
for nb in "APR" "DR" "G4" "IR" "MR" "TRI" "STR" "Z"
do
  cat final_nonB/merged_per_chrom/chr[0-9]*.$nb.bed > final_nonB/merged_per_chrom/autosomes.$nb.bed
  cat final_nonB/merged_per_chrom/rDNA_*.$nb.bed > final_nonB/merged_per_chrom/rDNA.$nb.bed
done
# Macro, micro, dot
for type in "macro" "micro" "microdot"
do
  for nb in "APR" "DR" "G4" "IR" "MR" "TRI" "STR" "Z"
  do
    for chr in $(cat $type.txt)
    do
      cat final_nonB/merged_per_chrom/$chr.$nb.bed
    done >final_nonB/merged_per_chrom/$type.$nb.bed
  done
done

# Then we can combine all the types
for file in $(ls stats/nonB_totals/annotated.*.txt )
do
  chr=`echo $file | cut -f2 -d"."`
  echo "Chr NB1 NB2 Ovl Frac" | sed "s/ /\t/g" >overlap/pairwise.$chr.txt
  echo '#!/bin/bash
  module load bedtools/2.31.0
# Go through all combinations
  for nb in "APR" "DR" "G4" "IR" "MR" "TRI" "STR" "Z"
  do
    nbbp=`grep $nb '$file' |cut -f2`
    echo "$nb has $nbbp number of bp!"
    for nb2 in "APR" "DR" "G4" "IR" "MR" "TRI" "STR" "Z"
    do
        if [[ $nb != $nb2 ]]
        then
          intersectBed -a final_nonB/merged_per_chrom/'$chr'.$nb.bed -b final_nonB/merged_per_chrom/'$chr'.$nb2.bed -wo | awk -v tot=$nbbp -v chr='$chr' -v sum=0 -v nb1=$nb -v nb2=$nb2 -v OFS="\t" '"'"'{sum+=$7}END{f=sum/tot; print chr,nb1,nb2,sum,f}'"'"' >>overlap/pairwise.'$chr'.txt
        fi
    done
  done
  ' | sbatch -J $chr --ntasks=1 --cpus-per-task=1 --mem-per-cpu=4G --time=1:00:00 --out slurm/job.pairwise-overlap.$chr.%j.%N.out
done

# Merging the groups
echo "Region NonB Overlap" |sed 's/ /\t/g' >overlap/merged.summary.txt
echo "Region NB1 NB2 Ovl Frac" |sed 's/ /\t/g' >overlap/merged.pairwise.txt
for type in "macro" "micro" "microdot" "autosomes" "chrW_mat" "chrZ_pat"
do
  awk -v t=$type '{print t"\t"$0}' overlap/summary.$type.txt >>overlap/merged.summary.txt
  awk -v t=$type '(NR>1){print $0}' overlap/pairwise.$type.txt >>overlap/merged.pairwise.txt
done


# Some stats
# Number of non-B bases
grep cro overlap/merged.summary.txt |awk '{sum+=$3}END{print sum}'
#245759702
# Number of bases in overlap
grep cro overlap/merged.summary.txt |grep "-" |awk '{sum+=$3}END{print sum}'
#48347053
# Number of bp in at least three types
grep cro overlap/merged.summary.txt |grep "-" |awk -F'\t' 'gsub(/-/, "&", $2) >= 2' |awk '{sum+=$3}END{print sum}'
#17180744

# Number of bases with only one vs more than one non-B annotation
for set in "autosomes" "chrZ_pat" "chrW_mat" "macro" "micro" "microdot"
do
  u=`grep -v "-" overlap/summary.$set.txt |awk '{sum+=$2}END{print sum}'`
  o=`grep "-" overlap/summary.$set.txt |awk '{sum+=$2}END{print sum}'`
  echo $set $u $o
done 

# Overlap between G4 and DR, genome wide
less overlap/merged.pairwise.txt |grep "DR" |grep "G4" |grep "cro" |sort -k4n |cut -f4 |uniq |awk '{sum+=$1}END{print sum}'
#16567688
# Total number of non-B bases
grep ALL densities/bTaeGut7v0.4_MT_rDNA.nonB_genome_wide.txt
#ALL 245781668 0.113636


######################## CHECKING NCBI FOR BIRD GENOMES ########################
# July 28, 2025

# Search for "Aves" On NCBI, download table with the following columns:
# Assembly Accession      Assembly Name   Organism Name   Organism Taxonomic ID   Organism Infraspecific Names Breed      Organism Infraspecific Names Strain     Organism Infraspecific Names Cultivar      Organism Infraspecific Names Ecotype    Organism Infraspecific Names Isolate    Organism Infraspecific Names Sex        Annotation Name Assembly Level  Assembly Release Da

# Save as ~/Downloads/ncbi_dataset.tsv

# GET UNIQUE TAXON ID:
cut -f4 ~/Downloads/ncbi_dataset.tsv |sort |uniq |wc
    1570    1572   10664

# Download one version with less info and extract completeness info:
#Assembly Accession      Assembly Name   Organism Taxonomic ID   Assembly Stats Total Number of Chromosomes      Assembly Level  WGS project accession
cut -f3,5 ~/Downloads/ncbi_dataset-4.tsv |uniq |sort |uniq |grep Complete >NCBI_complete.txt
cut -f3,5 ~/Downloads/ncbi_dataset-4.tsv |uniq |sort |uniq |grep Chromosome >NCBI_chromosome.txt
cut -f3,5 ~/Downloads/ncbi_dataset-4.tsv |uniq |sort |uniq |grep Contig >NCBI_contig.txt
cut -f3,5 ~/Downloads/ncbi_dataset-4.tsv |uniq |sort |uniq |grep Scaffold >NCBI_scaffold.txt

# Number complete:
wc -l NCBI_complete.txt
#       1 NCBI_complete.txt
# Number chromosome (remove complete)
join -v1 -1 1 -2 1 <(cut -f1 NCBI_chromosome.txt) <(cut -f1 NCBI_complete.txt) |wc
#     231     231    1485
# Number of scaffold level (remove chromosome and complete)
join -v1 -1 1 -2 1 <(cut -f1 NCBI_scaffold.txt) <(cut -f1 NCBI_chromosome.txt) |wc
#    1302    1302    8916
