################################################################################
### ANNOTATION OF NON-B DNA IN THE ZEBRA FINCH T2T GENOME
### CODE WRITTEN BY LINNÉA SMEDS

# This file contains the first code, including downloading the genome and 
# running the annotation software. 

# All commands are assumed to be run from within the downloaded github repo.
# Time and/or memory consuming jobs were run on a slurm cluster. The code 
# within the slurm scripts can be copied and run on any machine with enough
# RAM/cores.

##### LIST OF CONTENTS 
# REQUIREMENTS 
# ANNOTATION WITH GFA
# ANNOTATION WITH QUADRON
# CREATE A JOINT DIRECTORY WITH ANNOTATIONS



################################# REQUIREMENTS #################################
# A LIST OF ALL SOFTWARE USED:
# bedtools/2.31.0
# gfa (downloaded Feb 4, 2025 from https://github.com/abcsFrederick/non-B_gfa)
# Quadron (dockerized version 1.0.0 from https://hub.docker.com/r/kxk302/quadron)


# Download Zebra Finch genome from GenomeArk
mkdir ref
cd ref
wget https://genomeark.s3.amazonaws.com/species/Taeniopygia_guttata/bTaeGut7/assembly_verkko_0.1/manual_curation/bTaeGut7v0.4/bTaeGut7v0.4_MT_rDNA.fa.gz
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

