################################################################################
### CODE FOR REPEATING PREVIOUS ANALYSES ON OTHER GENOMES:
### THE CHICKEN GENOME (T2T EXCEPT FOR THE W-CHROMOSOME)
### ANNAS HUMMINGBIRD (NOT T2T)
### WRITTEN BY LINNÉA SMEDS

# TABLE OF CONTENTS 
# DOWNLOAD AND PREPARE GENOMES 
# ANNOTATION WITH GFA
# ANNOTATION WITH QUADRON
# CALCULATE DENSITY
#     -PLOT DENSITIES
# SORT CHROMOSOMES INTO GROUPS
# PREPARE ENRICHMENT CALCULATIONS
#     -PRINT DENSITIES FOR TABLE
#     -GC CONTENT PER CHROM
# ENRICHMENT IN CENTROMERES
#     -GC CONTENT FOR CENTROMERES
#    -STATISTICAL TEST 


######################### DOWNLOAD AND PREPARE GENOME ##########################
# Download chicken genome and annotation from dropbox
# https://www.dropbox.com/scl/fo/plq2tm2w9lzlk0ua1rzph/h?rlkey=l6z3rgmjs7ec9azun8nundnzl&e=1&dl=0
# Zebra finch for homologies 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/771/995/GCF_048771995.1_bTaeGut7.mat/GCF_048771995.1_bTaeGut7.mat_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/771/995/GCF_048771995.1_bTaeGut7.mat/GCF_048771995.1_bTaeGut7.mat_assembly_report.txt
# Download Anna's hummingbird from NCBI: 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/957/555/GCF_003957555.1_bCalAnn1_v1.p/GCF_003957555.1_bCalAnn1_v1.p_assembly_report.txt
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/957/555/GCF_003957555.1_bCalAnn1_v1.p/GCF_003957555.1_bCalAnn1_v1.p_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/957/555/GCF_003957555.1_bCalAnn1_v1.p/GCF_003957555.1_bCalAnn1_v1.p_genomic.gff.gz
# Mallard 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/746/955/GCA_008746955.3_CAU-Wild1.1/GCA_008746955.3_CAU-Wild1.1_assembly_report.txt
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/746/955/GCA_008746955.3_CAU-Wild1.1/GCA_008746955.3_CAU-Wild1.1_genomic.fna.gz
#annotation from https://zenodo.org/records/12721248
# Great Bustard 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/026/413/225/GCA_026413225.1_OTswu/GCA_026413225.1_OTswu_assembly_report.txt 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/026/413/225/GCA_026413225.1_OTswu/GCA_026413225.1_OTswu_genomic.fna.gz 
# annotation from https://figshare.com/articles/dataset/Figure_source_data/23650419, 
# had to add "chr" and remove blank lines using "|grep -v ^$ |awk -v OFS="\t" '{$1="chr"$1; print}'"
# Ural owl 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/047/716/275/GCF_047716275.1_bStrUra1/GCF_047716275.1_bStrUra1_assembly_report.txt
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/047/716/275/GCF_047716275.1_bStrUra1/GCF_047716275.1_bStrUra1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/047/716/275/GCF_047716275.1_bStrUra1/GCF_047716275.1_bStrUra1_genomic.gff.gz
# place files in ref/ and unzip 

# Make a list of species 
echo "zebra_finch bTaeGut7v0.4_MT_rDNA bTaeGut7v0.4_MT_rDNA
chicken chicken.v23 chicken.v23
annas_hummingbird GCF_003957555.1_bCalAnn1_v1.p bCalAnn1_v1.p
mallard GCA_008746955.3_CAU-Wild1.1 CAU-Wild1.1
great_bustard GCA_026413225.1_OTswu OTswu
ural_owl GCF_047716275.1_bStrUra1 bStrUra1" > helpfiles/species_list.txt

# Change chromosome names to something shorter
awk '{if(/^>/){if(/chromosome/){gsub(/,/, "", $7); print ">chr"$7}else{print $1}}else{print $0}}' ref/GCF_003957555.1_bCalAnn1_v1.p_genomic.fna >ref/bCalAnn1_v1.p.fa

# Need to change the annotations files from NCBI to common chromosome names! 
module load python/3.11.2 
module load samtools/1.19.2
cat helpfiles/species_list.txt |grep bCalAnn1 |while read -r sp longname prefix
do
    echo "Looking at "$sp
    grep -v "#" ref/${longname}_assembly_report.txt |cut -f1,3,4,5,7 |awk '{if($3=="Chromosome"){new="chr"$2}else{new=$1}; if($5=="na"){old=$4}else{old=$5}; print old"\t"new}' >ref/$prefix.conversion-table.txt
    # If a gff file exists, change chromosome names there too
    if [ -f ref/${longname}_genomic.gff ]; then
        echo "Translate gff file for $sp"
        python3 T2T_bird_nonB/python/replace_chromosome_names.py -g ref/${longname}_genomic.gff -t ref/$prefix.conversion-table.txt -o ref/$prefix.raw.gff
        # remove the "region" lines in the gff since they represent full chromosomes
        awk '($3 != "region")' ref/$prefix.raw.gff > ref/${prefix}.gff  
    fi
    # Change name in fasta file
 #   awk 'NR==FNR{a[$1]=$2; next} /^>/ {key = substr($1, 2); $0 = ">"a[key]}1' ref/$prefix.conversion-table.txt ref/${longname}_genomic.fna > ref/$prefix.fa
    # Index genome 
 #   samtools faidx ref/$prefix.fa
done 


# For finding homologies, use RefSeq annotation for zebra finch 
longname="GCF_048771995.1_bTaeGut7.mat"
prefix="bTaeGut7.mat"
grep -v "#" ref/${longname}_assembly_report.txt |cut -f1,3,4,5,7 |awk '{new=$1; if($5=="na"){old=$4}else{old=$5}; print old"\t"new}' >ref/$prefix.conversion-table.txt
python3 T2T_bird_nonB/python/replace_chromosome_names.py -g ref/${longname}_genomic.gff -t ref/$prefix.conversion-table.txt -o ref/$prefix.raw.gff
awk '($3 != "region")' ref/$prefix.raw.gff > ref/${prefix}.gff  
cp ref/${prefix}.gff  ref/bTaeGut7v0.4_MT_rDNA.gff

################################# FIND HOMOLOGY ################################

# Get protein and cds sequences 
cat helpfiles/species_list.txt |grep mallard |while read -r sp longname prefix
do
  ~/software/gffread-0.12.7.Linux_x86_64/gffread ref/$prefix.gff \
  -g ref/$prefix.fa \
  -y ref/$prefix.proteins.faa \
  -x ref/$prefix.cds.fna \
  -F
done 
# Something is wrong with mallard, it seems like the annotation file comes from 
# another version of the mallard genome (some genes are outside chromosomes)

# Move to orthofinder and check there are no dots 
mkdir orthofinder/
cat helpfiles/species_list.txt |while read -r sp longname prefix
do
  awk '{if(/>/){print $0}else{gsub(/\./,"",$0); print $0}}' ref/$prefix.proteins.faa >orthofinder/$prefix.noDots.proteins.faa
done

# I accidentally removed dots also in the headers so I need a map file to go back to original names
cat helpfiles/species_list.txt |while read -r sp longname prefix
do 
  awk '{if(/>/){ gsub(/>/,"",$0); orig=$1; gsub(/\./,"",$1); print orig"\t"$1}}' ref/$prefix.proteins.faa >orthofinder/$prefix.tmp.id_map.tsv
done

#Orthodinder is installed in a virtual environment created in ~/software/, with
module load anaconda/2023.09
cd ~/software/
conda create -n of3_env python=3.10
conda activate of3_env
conda install -c bioconda orthofinder

# Run orthofinder
cd /storage/group/kdm16/default/lbs5874/ZebraFinch
echo '#!/bin/bash
orthofinder -f orthofinder/ -M msa -T fasttree
' |sbatch -J orthofinder --ntasks=1 --cpus-per-task=8 --time=1-00:00:00 -o slurm/job.orthofinder.%j.out


# Get map files between gene ID - chromosome 
cat helpfiles/species_list.txt |while read -r sp longname prefix
do
  awk '$3!="gene"{
    match($9,/Parent=([^;]+)/,a);
    print a[1] "\t" $1
  }' ref/${prefix}.gff  |sort |uniq > orthofinder/${prefix}_gene2chr.tsv
done

# Tmp map files without dots 
cat helpfiles/species_list.txt |while read -r sp longname prefix
do
  awk '$3!="gene"{
    match($9,/Parent=([^;]+)/,a);
    gsub(/\./,"",a[1]);
    print a[1] "\t" $1
  }' ref/${prefix}.gff |sort |uniq > orthofinder/${prefix}.tmp.gene2chr.tsv
done

# Use chatGPT generated python script to link chromosome between species. 
module load python/3.11.2 
pip install pandas
python3 T2T_bird_nonB/python/orthofinder_chr_homology.py
# The input files are hardcoded in this script. The output is saved to 
# orthofinder/sp1_vs_sp2_chr_homology.tsv  



######################## ALIGN MALLARD TO OTHER SPECIES ########################
# Because the mallard annotation was problematic, I decided to align it to the 
# chicken genome with lastz, using fast settings, and then just count the 
# number or length of alignments between the chromosome pairs. 

# Runing lastz on each chicken chromosome 
mkdir lastz
for chr in $(cut -f1 ref/chicken.v23.fa.fai |grep -v "mito")
do
  echo '#!/bin/bash
  gunzip -c ref/split_fasta/chicken.v23/'$chr'.fa.gz >/scratch/lbs5874/'$chr'.fa
  /storage/home/lbs5874/lastz-distrib/bin/lastz /scratch/lbs5874/'$chr'.fa ref/CAU-Wild1.1.fa \
      --nochain --gapped --gfextend \
      --format=general:name1,start1,end1,length1,name2,start2,end2,strand2,score >lastz/chicken.v23.'$chr'.vs.mallard.lastz.out
  ' |sbatch -J lastz.$chr --ntasks=1 --cpus-per-task=1 --time=1-00:00:00 -o slurm/job.lastz.%j.out
done 
# Reciprocally, runing lastz on each mallard chromosome
for chr in $(cut -f1 ref/CAU-Wild1.1.fa.fai |grep "chr")
do
  echo '#!/bin/bash
  gunzip -c ref/split_fasta/CAU-Wild1.1/'$chr'.fa.gz >/scratch/lbs5874/mallard.'$chr'.fa
  /storage/home/lbs5874/lastz-distrib/bin/lastz /scratch/lbs5874/mallard.'$chr'.fa ref/chicken.v23.fa \
      --notransition --step=20 --nogapped \
      --format=general:name1,start1,end1,length1,name2,start2,end2,strand2,score >lastz/mallard.'$chr'.vs.chicken.lastz.out
  ' |sbatch -J lastz.$chr --ntasks=1 --cpus-per-task=1 --time=1-00:00:00 -o slurm/job.lastz.%j.out
done 

echo "Mallard chromosome aligned to chicken:"
for chr in $(cut -f1 ref/CAU-Wild1.1.fa.fai |grep "chr")
do
  echo "################################################"
  less lastz/mallard.$chr.vs.chicken.lastz.out |awk '($9>10000){print}' |cut -f1,5 |sort |uniq -c |sort -k1,1nr |head -n5
done

# Mallard on bustard
for chr in $(cut -f1 ref/CAU-Wild1.1.fa.fai |grep "chr")
do
  echo '#!/bin/bash
  #gunzip -c ref/split_fasta/CAU-Wild1.1/'$chr'.fa.gz >/scratch/lbs5874/mallard.'$chr'.fa
  /storage/home/lbs5874/lastz-distrib/bin/lastz /scratch/lbs5874/mallard.'$chr'.fa ref/OTswu.fa \
      --notransition --step=20 --nogapped \
      --format=general:name1,start1,end1,length1,name2,start2,end2,strand2,score >lastz/mallard.'$chr'.vs.bustard.lastz.out
  ' |sbatch -J lastz.$chr --ntasks=1 --cpus-per-task=1 --time=1:00:00 -o slurm/job.lastz.%j.out
done 

echo "Mallard chromosome aligned to bustard:"
for chr in $(cut -f1 ref/CAU-Wild1.1.fa.fai |grep "chr11")
do
  echo "################################################"
  less lastz/mallard.$chr.vs.bustard.lastz.out |awk '($9>10000){print}' |cut -f1,5 |sort |uniq -c |sort -k1,1nr |head -n5
done

# Mallard on zebra finch
for chr in $(cut -f1 ref/CAU-Wild1.1.fa.fai |grep "chr")
do
  echo '#!/bin/bash
  #gunzip -c ref/split_fasta/CAU-Wild1.1/'$chr'.fa.gz >/scratch/lbs5874/mallard.'$chr'.fa
  /storage/home/lbs5874/lastz-distrib/bin/lastz /scratch/lbs5874/mallard.'$chr'.fa ref/bTaeGut7v0.4_MT_rDNA.fa \
      --notransition --step=20 --nogapped \
      --format=general:name1,start1,end1,length1,name2,start2,end2,strand2,score >lastz/mallard.'$chr'.vs.zebrafinch.lastz.out
  ' |sbatch -J lastz.$chr --ntasks=1 --cpus-per-task=1 --time=1:00:00 -o slurm/job.lastz.%j.out
done 

echo "Mallard chromosome aligned to zebra finch:"
for chr in $(cut -f1 ref/CAU-Wild1.1.fa.fai |grep "chr")
do
  echo "################################################"
  less lastz/mallard.$chr.vs.zebrafinch.lastz.out |awk '($9>10000){print}' |cut -f1,5 |sort |uniq -c |sort -k1,1nr |head -n5
done



############################## ANNOTATION WITH GFA #############################
# GFA is assumed to be installed in: 
# ~/software/non-B_gfa/gfa
# Set parameters
fasta="chicken.v23.fa"
#fasta="bCalAnn1_v1.p.fa"
prefix=`echo $fasta |sed 's/.fa//'`
# Run gfa and convert to bed
cat helpfiles/species_list.txt |while read -r sp longname prefix
do
#  echo "start gfa for $sp"
  echo '#!/bin/bash
  echo "##### Start gfa"
  ~/software/non-B_gfa/gfa -seq ref/'$prefix'.fa -out gfa_annotation/'$prefix'  -skipGQ -skipWGET
  echo "##### translate to bed"
  for type in "APR" "DR" "IR" "MR" "STR" "Z"
  do
    echo "#####....translating $type...."
    awk -v OFS="\t" '"'"'(NR>1){s=$4-1; print $1,s,$5,$9":"$10":"$11,".",$8}'"'"' gfa_annotation/'${prefix}'_${type}.tsv >final_nonB/'${prefix}'.$type.bed
  done
  ' |sbatch -J gfa.$prefix  --ntasks=1 --cpus-per-task=1 --mem-per-cpu=2G --time=2-00:00:00 -o slurm/job.gfa.$prefix.%j.%N.out
done 
# Change name of original IR, Z and MR files as they will not be used for 
# the main analysis 
cat helpfiles/species_list.txt |grep anna |while read -r sp longname prefix
do
  mv final_nonB/${prefix}.IR.bed final_nonB/${prefix}.IRgfa.bed
  mv final_nonB/${prefix}.Z.bed final_nonB/${prefix}.Zgfa.bed
  mv final_nonB/${prefix}.MR.bed final_nonB/${prefix}.MRgfa.bed
  #Subset DR (not sure we will use this), IR and MR
  awk -v OFS="\t" '($9<=50 && $10<=5){s=$4-1; print $1,s,$5,$9":"$10":"$11,".",$8}' gfa_annotation/${prefix}_DR.tsv  >final_nonB/${prefix}.DRfilt.bed
  awk -v OFS="\t" '($9<=30 && $10<=10){s=$4-1; print $1,s,$5,$9":"$10":"$11,".",$8}' gfa_annotation/${prefix}_IR.tsv  >final_nonB/${prefix}.IR.bed
  # max 10bp loop, and 100% purine or pyrimidine
  awk -v OFS="\t" '($10<=10){split($13, parts, "/"); 
      for (i=1; i<=4; i++) {
        cnt = parts[i]
        gsub(/[A-Z]/, "", cnt)  # remove A,C,G,T to get the numeric part
        vals[i] = cnt + 0
      }
      A = vals[1]; C = vals[2]; G = vals[3]; T = vals[4];
      total = A + C + G + T;
      pur = A + G;
      frac=pur/total;
    if(frac==0 || frac==1){s=$4-1; print $1,s,$5,$9":"$10":"$11":"frac,".",$8}}' gfa_annotation/${prefix}_MR.tsv  >final_nonB/${prefix}.TRI.bed
done 

# This was just run for chicken that I already had results for 
cp old.final_nonB/${prefix}.IR.bed final_nonB/${prefix}.IRall.bed
cp old.final_nonB/${prefix}.Z.bed final_nonB/${prefix}.Zgfa.bed
cp old.final_nonB/${prefix}.MR.bed final_nonB/${prefix}.MRall.bed
cp old.final_nonB/${prefix}.APR.bed final_nonB/${prefix}.APR.bed
cp old.final_nonB/${prefix}.STR.bed final_nonB/${prefix}.STR.bed
cp old.final_nonB/${prefix}.DR.bed final_nonB/${prefix}.DR.bed

############################# G4 WITH G4DISCOVERY ##############################

# Split input files 
cat helpfiles/species_list.txt |tail -n3 |while read -r sp longname prefix
do
  fasta=../..ref/$prefix.fa
    echo "Split fasta for $sp"
    mkdir -p ref/split_fasta/$prefix
    cd ref/split_fasta/$prefix
    echo '#!/bin/bash
    faidx --split-files ../../'$prefix'.fa
    ' |sbatch -J split --ntasks=1 --cpus-per-task=1 --time=5:00:00 -o job.splitfasta.$prefix.%j.out
    cd ../../..
done 

# Zip input files 
cat helpfiles/species_list.txt |tail -n3 |while read -r sp longname prefix
do
  echo "Compress fasta for $sp"
  echo '#!/bin/bash
  for file in ref/split_fasta/'$prefix'/*fa
  do 
      ls $file 
      gzip $file 
  done 
    ' |sbatch -J zip --ntasks=1 --cpus-per-task=1 --time=5:00:00 -o job.zipfasta.$prefix.%j.out
done 


# Start g4Discovery for all chromosome files 
cat helpfiles/species_list.txt |grep "chicken" |while read -r sp longname prefix
do
  echo "Start g4Discovery for $sp"
  mkdir -p g4discovery/$prefix
  for file in $(ls ref/split_fasta/$prefix/chr16.fa.gz)
  do
    # If file is bigger than 30Mb, use 4 cores
    if [ $(awk -v fa=$file 'BEGIN{cmd="stat -c%s "fa; cmd | getline size; close(cmd); print size}') -gt 30000000 ]; then
        mem="4G"
        cpus=4
    else
        mem="8G" #"4G"
        cpus=4  #"1"
    fi
    chr=$(basename ${file%.fa.gz})
    out=g4discovery/$prefix/$chr.bed.gz
    if [ -f $out ]; then
        #echo $out" exists, skipping"
        continue
    fi  
    echo "start g4Discovery for "$chr
    echo '#!/bin/bash
    module load python/3.11.2
    module load r/4.5.0
    python3 ~/software/g4Discovery.PanSN/src/g4Discovery.py -fa '$file' -o g4discovery/'$prefix'/'$chr'.bed 
      ' |sbatch -J g4D.$chr  --ntasks=1 --cpus-per-task=$cpus --mem-per-cpu=$mem --time=1-00:00:00 -o slurm/job.g4D.$chr.$prefix.%j.%N.out
  done
done 
# used 1 cores for all files smaller than 30Mb for which I used 4 cores. 
# For chicken, 4 chromosomes kept failing (chr 16, chr29, chr36, chr38)


# Combine results into one bed file NEED TO RERUN CHICKEN
cat helpfiles/species_list.txt |grep chicken |while read -r sp longname prefix
do
  echo "Combine g4Discovery results for $sp"
  rm -f final_nonB/$prefix.G4.bed
  for chr in $(cat ref/$prefix.fa.fai |cut -f1)
  do 
      zcat g4discovery/$prefix/$chr.bed.gz >>final_nonB/$prefix.G4.bed
  done 
done


##################### ANNOTATION OF Z DNA WITH ZDNA HUNTER  ####################

# ZDNA Hunter runs on a webserver located here: https://bioinformatics.ibp.cz/#/analyse/zdna
# We uploaded the full genomes and downloaded the results as bed graph files. 
# ZDNA Hunter outputs positions continusly for the whole genome, so we need to 
# convert the positions to proper chromosome coordinates.

# Place output in the directory ZDNAHunter
#fasta="chicken.v23.fa"
#fasta="bCalAnn1_v1.p.fa"
#prefix=`echo $fasta |sed 's/.fa//'`
module load python/3.11.2
cat helpfiles/species_list.txt |tail -n4 |while read -r sp longname prefix
do
  #python3 T2T_bird_nonB/python/remap_zdna.py -b ZDNAHunter/$longname.model1.bedgraph -f ref/$prefix.fa.fai -o final_nonB/$prefix.ZDNAm1.bed
  python3 T2T_bird_nonB/python/remap_zdna.py -b ZDNAHunter/$longname.model2.bedgraph -f ref/$prefix.fa.fai -o final_nonB/$prefix.Z.bed
done 
# For now, we will use model 2 so I name that just "Z"

######################## ANNOTATION OF G4s WITH QUADRON ########################
# Assumes the quadron singularity is here: ~/software/quadron.sif 
# Divide fasta into chromosomes
#fasta="chicken.v23.fa"
fasta="bCalAnn1_v1.p.fa"
prefix=`echo $fasta |sed 's/.fa//'`
mkdir ref/split_fasta/$prefix
echo '#!/bin/bash
faidx --split-files ref/'$fasta'
mv *.fa ref/split_fasta/'$prefix'/
' |sbatch -J split --ntasks=1 --cpus-per-task=4 --time=5:00:00 -o slurm/job.splitfasta.$prefix.%j.out

# Run Quadron on full chromosomes (separately)
# Note that Quadron needs a full path to the input files
mkdir Quadron_annotation/$prefix
#for fa in $(ls [insert_full_path_here]/ref/split_fasta/$prefix/chr*.fa)
for fa in $(ls /storage/group/kdm16/default/lbs5874/ZebraFinch/ref/split_fasta/$prefix/chr*.fa |grep -v "chr1.fa" |grep -v "chr2.")
do
  ls $fa
  chr=`basename $fa |sed 's/.fa//'`
  echo $chr
  out="/storage/group/kdm16/default/lbs5874/ZebraFinch/Quadron_annotation/'$prefix'/$chr.out"
  echo '#!/bin/bash
  singularity exec ~/software/quadron.sif sh -c "cd /scripts/; Rscript run_quadron.R '$fa' '$out' 4 1000000"
  ' |sbatch -J $chr --ntasks=1 --cpus-per-task=2 --mem-per-cpu=4G --time=1:00:00 -o slurm/job.quadron.$chr.%j.out
done
# Chr1-2 needed 4 cores, the others worked with 2 cores 

# Check that there aren't too many NA motifs
for file in $(ls Quadron_annotation/$prefix/chr*.out)
do
  na=`grep "^DATA:" $file |awk '($5=="NA"){n++}END{print n}'`
  echo $file $na
done
# Looks ok, most have 2 NA G4s, some have 1 or 3

# Convert to bed
prefix="chicken.v23"
rm -f Quadron_annotation/$prefix.G4.bed
echo $prefix
for chr in $(cut -f1 ref/$prefix.fa.fai)
do
  awk -v chr=$chr -v OFS="\t" '(/^DATA/){if($5=="NA"){}else{score=$5; if(score<0){score=0}; s=$2-1; e=s+$4; print chr,s,e,$5,sprintf("%.f", score),$3}}' Quadron_annotation/$prefix/$chr.out >>final_nonB/$prefix.G4quadron.bed
done


####################### CALCULATE COVERAGE IN WINDOWS ##########################
prefix="chicken.v23"
# Make sure all non-B bed files are in the final folder, and make merged versions 
cat helpfiles/species_list.txt |tail -n4 |while read -r sp longname prefix
do
  for n in "G4" "APR" "DR" "DRfilt" "IR" "IRall" "MRall" "TRI" "STR" "Z" "Zgfa" #"G4quadron" "ZDNAm1" "Z"
  do
  echo '#!/bin/bash
  module load bedtools/2.31.0
  sort -k1,1 -k2,2n final_nonB/'$prefix'.'${n}'.bed |mergeBed -i - >final_nonB/'$prefix'.'${n}'.merged.bed
  '| sbatch -J $n --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8G --out slurm/mergeBed.$n.%j.out
  done
done

# And a joint bed file called "Any" (Just using the chosen settings for each type)
cat helpfiles/species_list.txt |tail -n5 |while read -r sp longname prefix
do
  echo "Make Any bed for $sp"
  sort -m final_nonB/$prefix.APR.merged.bed final_nonB/$prefix.DR.merged.bed final_nonB/$prefix.IR.merged.bed final_nonB/$prefix.TRI.merged.bed final_nonB/$prefix.STR.merged.bed final_nonB/$prefix.G4.merged.bed final_nonB/$prefix.Z.merged.bed >final_nonB/$prefix.Any.bed
  sort -k1,1 -k2,2n final_nonB/$prefix.Any.bed |mergeBed -i - >final_nonB/$prefix.Any.merged.bed
done 

# Create genomic windows and make length file
module load bedtools/2.31.0
cat helpfiles/species_list.txt |tail -n4 |while read -r sp longname prefix
do
  bedtools makewindows -g ref/$prefix.fa.fai -w 100000 >windows/$prefix.100kb_windows.bed;
  bedtools makewindows -g ref/$prefix.fa.fai -w 10000 >windows/$prefix.10kb_windows.bed;
done 

cat helpfiles/species_list.txt |grep zebra |while read -r sp longname prefix
do
  for wind in "100kb" #"10kb" #
  do
    for n in  "G4" #"APR" "DR" "DRfilt" "IR" "IRall" "MRall" "TRI" "STR" "Z" "ZDNAgfa" #"G4quadron"
    do
      echo '#!/bin/bash
      module load bedtools/2.31.0
      intersectBed -wao -a windows/'$prefix'.'${wind}'_windows.bed -b final_nonB/'$prefix'.'${n}'.merged.bed | cut -f1,2,3,7 |awk -v OFS="\t" '"'"'{if(NR==0){chr=$1; s=$2; e=$3; sum=$4}else{if($1==chr && $2==s){sum+=$4}else{print chr,s,e,sum; chr=$1; s=$2; e=$3; sum=$4}}}END{print chr,s,e,sum}'"'"' | sed "/^\s*$/d" >coverage/'${prefix}'.'${n}'.'$wind'.bed
      '| sbatch -J $n --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8G --out slurm/coverage.$n.%j.out
      done
  done 
done 

# For plotting, merge the data:
cat helpfiles/species_list.txt |grep zebra |while read -r sp longname prefix
do
  echo "NonB Chr Start Stop Dens" |sed 's/ /\t/g' >coverage/$prefix.merged.100kb.txt
  for n in "APR" "DR" "G4" "IR" "TRI" "STR" "Z"
  do
    awk -v OFS="\t" -v nb=$n '{print nb,$0}' coverage/${prefix}.${n}.100kb.bed >>coverage/$prefix.merged.100kb.txt
  done
done 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PLOT COVERAGE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Assume circos is installed in ~/software/circos-0.69-9/
# The circos configuration files are found in circos/
mkdir -p circos/chicken/

# Centromeres for most chicken chromosome were kindly privided by Luohao Xu. I've
# placed them in helpfiles/chicken.CEN.bed. But the scripts below assumes it is
# placed in the ref/ directory.
mv helpfiles/chicken.CEN.bed ref/

# Create circos karyotype, nonB density and centromere band files
# BUT skip chrW
dir="circos/chicken/"
prefix="chicken.v23"
grep -v "mito" ref/$prefix.fa.fai |grep -v "chrW" | awk '{gsub(/chr/,""); print "chr - nn"$1,$1,"0",$2,"vvlgrey"}' |sort sort -k4,4n >$dir/karyotype.txt
# change order so chrZ is last
grep -v "nnZ" $dir/karyotype.txt >tmp
grep "nnZ" $dir/karyotype.txt >>tmp
mv tmp $dir/karyotype.txt
# Centromeres
cat ref/chicken.CEN.bed |grep -v "chrW" |awk '{gsub(/chr/,""); print "nn"$1,$2,$3}' >$dir/data/highlight_CEN.txt
for nb in "TRI" "APR" "DR" "G4" "IR" "STR" "Z"
do
  grep -v "mito" coverage/$prefix.${nb}.100kb.bed |grep -v "chrW" |sed '/^[[:space:]]*$/d' |awk '{gsub(/chr/, "nn"); $3=$3-1; print $1,$2,$3,$4}' >$dir/data/${nb}.100kb.txt
  max=`cut -f4 -d" " $dir/data/${nb}.100kb.txt |sort -nr |head -n1`
  echo "will replace max = $nb with max = $max"
  sed -i'' -e "s/max = $nb/max = $max/" $dir/circos_halfInnerCircle.conf
done

# Plot
cd $dir
~/software/circos-0.69-9/bin/circos -conf circos_halfInnerCircle.conf
mv circos.png chicken_circos.png
cd ../..


######################### SORT CHROMOSOMES INTO GROUPS #########################
# Make a file with lengths and the chromosome types (macro, micro and microdot)
# The latter is added by hand, taken from the PNAS paper by Huang et al 2023.
cat helpfiles/species_list.txt |while read -r sp longname prefix
do
  echo "Chr Length Group" |sed 's/ /\t/g' >by.hand.$prefix.groups.txt
  cat ref/$prefix.fa.fai | awk -v OFS="\t" '{print $1,$2}' >by.hand.$prefix.groups.txt
done
# Added group information by hand to the files, based on synteny with Zebra finch. 
# Then add "unplaced" to the unplaced scaffolds 
cat helpfiles/species_list.txt |while read -r sp longname prefix
do
  awk '{if(NF==3){print}else{print $0"\tunplaced"}}' by.hand.$prefix.groups.txt >helpfiles/$prefix.groups.txt
done



###################### PREPARE ENRICHMENT CALCULATIONS #########################
# Calculate genome-wide coverage to use for normalization
cat helpfiles/species_list.txt |while read -r sp longname prefix
do
  echo '#!/bin/bash
  module load bedtools/2.31.0
  rm -f coverage/'${prefix}'.per_genome.tsv
  totlen=`cat ref/'${prefix}'*.fa.fai | awk '"'"'{sum+=$2}END{print sum}'"'"'`
  echo "Totlen: "$totlen
  for non_b in "APR" "DR" "STR" "IR" "TRI" "G4" "Z" "Any"
  do
    awk -v n=$non_b -v tot=$totlen '"'"'{sum+=$3-$2}END{d=sum/tot; print n, sum, d}'"'"' final_nonB/'${prefix}'.${non_b}.merged.bed  >>coverage/'${prefix}'.per_genome.tsv
  done
  ' |sbatch -J enrichment --ntasks=1 --cpus-per-task=2 --mem-per-cpu=8G --out slurm/job.enrichment.$prefix.%j.out
done 


# Per group (not only enrichment, but per chromosome coverage)
cat helpfiles/species_list.txt |grep zebra |while read -r sp longname prefix
do
  for group in "macro" "micro" "dot" "unplaced"
  do
    echo '#!/bin/bash
    totlen=`grep '$group' helpfiles/'$prefix'.groups.txt | \
    awk '"'"'NR==FNR {keep[$1]; next} $1 in keep'"'"' /dev/stdin ref/'$prefix'.fa.fai | \
    awk '"'"'{sum+=$2}END{print sum}'"'"'`
    module load bedtools/2.31.0
    rm -f tmp.'$prefix'.'$group'.txt
    for non_b in "APR" "DR" "G4" "IR" "TRI" "STR" "Z" "Any"
    do
      grep '$group' helpfiles/'$prefix'.groups.txt | \
      awk '"'"'NR==FNR { keep[$1]; next } $1 in keep'"'"' /dev/stdin final_nonB/'${prefix}'.${non_b}.merged.bed | \
      awk -v OFS="\t" -v n=$non_b -v tot=$totlen -v g='$group' '"'"'{sum+=$3-$2}END{d=sum/tot; print g, n, sum, d}'"'"' >>tmp.'$prefix'.'$group'.txt
    done
  '| sbatch -J $group.$sp --ntasks=1 --cpus-per-task=1 --time=1:00:00 --partition=open
  done
done 

#Merge output 
echo -e "Group\tNonB\tBp\tCoverage" >coverage/$prefix.per_group.tsv
cat helpfiles/species_list.txt |grep zebra |while read -r sp longname prefix
do
  for group in "macro" "micro" "dot" "unplaced"
  do
    cat tmp.$prefix.$group.txt >>coverage/$prefix.per_group.tsv
  done 
done


# Per chromosome (not only enrichment, but per chromosome coverage)
cat helpfiles/species_list.txt |while read -r sp longname prefix
do
  echo '#!/bin/bash
  module load bedtools/2.31.0
  echo "Chr NonB Bp Coverage" |sed "s/ /\t/g" >coverage/'${prefix}'.per_chrom.tsv
  cat ref/'${prefix}'*.fa.fai |grep "chr" |cut -f1,2 | while read -r chr len;
  do
      for non_b in "APR" "DR" "G4" "IR" "TRI" "STR" "Z" "Any"
      do
          cat final_nonB/'${prefix}'.${non_b}.merged.bed |awk -v c=$chr '"'"'($1==c){print}'"'"' | awk -v OFS="\t" -v n=$non_b -v tot=$len '"'"'{sum+=$3-$2; c=$1}END{d=sum/tot; print c, n, sum, d}'"'"' >>coverage/'${prefix}'.per_chrom.tsv
      done
  done
  '| sbatch -J per_chrom --ntasks=1 --cpus-per-task=1 --mem-per-cpu=4G --time=1:00:00 --partition=open
done



# ~~~~~~~~~~~~~~~~~~~~~~~~~ PRINT DENSITIES FOR TABLE ~~~~~~~~~~~~~~~~~~~~~~~~~~
# Write it out with the nonB types as columns instead
# Whole genome
echo "Sp Chromosome APR DR G4 IR TRI STR Z ALL"
cat helpfiles/species_list.txt |while read -r sp longname prefix
do
  tmp="$sp genome"
  for non_b in "APR" "DR" "G4" "IR" "TRI" "STR" "Z" "Any"
  do
    n=`awk -v n=$non_b '($1==n){print $3}' coverage/${prefix}.nonB_genome_wide.txt`
    tmp=$tmp" "$n
  done
  echo $tmp
done 

# Per chromosome
cat ref/${prefix}*.fa.fai |grep "chr" |while read -r chr;
do
    tmp="$chr"
    for non_b in "APR" "DR" "G4" "IR" "TRI" "STR" "Z" "Any"
    do
      n=`awk -v c=$chr '($1==c){print}' coverage/${prefix}.nonB_per_chrom.txt |grep " $non_b " |cut -f4 -d" "`
      tmp=$tmp" "$n
    done
    echo $tmp
done


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ GC CONTENT PER CHROM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
prefix="chicken.v23"

# Find GC and AT of full genomes
echo '#!/bin/bash
cat ref/'$prefix'.fa | awk '"'"'(!/^>/){gc+=gsub(/[gGcC]/,""); at+=gsub(/[aAtT]/,"")}END{print gc"\t"at}'"'"' >stats/'$prefix'.nuc
'| sbatch -J GC --ntasks=1 --cpus-per-task=1 --time=1:00:00 -o slurm/calculate.gc.$prefix.%j.out

# And for each chromosome separately
# (used --cpus-per-task=2 --mem-per-cpu=4G for chr1 and chr2)
cat ref/${prefix}*.fa.fai |cut -f1  |tail -n+3 |while read -r chr;
do
  echo $chr
  echo '#!/bin/bash
  faidx ref/'$prefix'.fa '$chr' | awk '"'"'(!/^>/){gc+=gsub(/[gGcC]/,""); at+=gsub(/[aAtT]/,"")}END{print gc"\t"at}'"'"' >stats/'$chr'.nuc
  '| sbatch -J $chr --ntasks=1 --cpus-per-task=1 --mem-per-cpu=4G --time=15:00 -o slurm/calculate.gc.$chr.%j.out
done

# Merge all and calculate GC content
echo "Chr GC_num AT_num GC_cont AT_cont" |sed 's/ /\t/g' >stats/$prefix.GC_per_chrom.txt
cat ref/${prefix}*.fa.fai |cut -f1 |while read -r chr;
do
 cat stats/$chr.nuc |awk -v c=$chr -v OFS="\t" '{sum=$1+$2; gc=$1/sum; at=$2/sum; print c,$1,$2,gc,at}'
done >>stats/$prefix.GC_per_chrom.txt



########################## ENRICHMENT IN CENTROMERES ###########################
# Chicken positions from Luohao Xu.
sed 's/ /\t/g' centromeres_chicken.txt >ref/chicken.CEN.bed

# He actually sent another set of centromeres later that were wider - I don't
# which are most correct, but I will test this file as well.

# Enrichment in centromeres compared to genome-wide AND to chromosome level
prefix="chicken.v23"
cen_bed_file="ref/chicken.CEN.bed"
module load bedtools/2.31.0
cat $cen_bed_file |cut -f1-3 |while read -r chr start end;
do
  cenlen=`echo "$end-$start" |bc`
  echo cenlen is $cenlen for chr $chr
  echo '#!/bin/bash
  module load bedtools/2.31.0
  echo Working on '$chr'
  rm -f centromeres/'$prefix'.'$chr'.tmp-centromere.txt
  cat coverage/'${prefix}'.nonB_genome_wide.txt |grep -v "ALL" | while read -r nb tot dens;
  do
    chrdens=`grep "'$chr' " coverage/'${prefix}'.nonB_per_chrom.txt |grep " $nb " |cut -f4 -d" "`
    echo "Density for '$chr' and $nb is $chrdens"
    stats=`echo '$chr' '$start' '$end' |sed "s/ /\t/g" |intersectBed -a - -b final_nonB/'${prefix}'.$nb.merged.bed | sort -k1,1 -k2,2n | mergeBed -i - |\
    awk -v l='$cenlen' -v dtot=$dens -v dchr=$chrdens '"'"'{sum+=$3-$2}END{d=sum/l; fracgw=d/dtot; fracchr=d/dchr; print d,fracgw,fracchr}'"'"'`
    echo '$chr'" "$nb" "$stats |sed "s/ /\t/g" >>centromeres/'$prefix'.'$chr'.tmp-centromere.txt
  done
  ' | sbatch -J $chr --ntasks=1 --cpus-per-task=1 --mem-per-cpu=4G --out slurm/job.centromeric_enrichment.$chr.%j.out
done
#Merge
echo "Chr nonB Coverage Enrichment_GW Enrichment_Chr" |sed "s/ /\t/g" >centromeres/${prefix}.enrichment.tsv
cat $cen_bed_file |cut -f1 |while read -r chr;
do
  cat centromeres/$prefix.$chr.tmp-centromere.txt >> centromeres/${prefix}.enrichment.tsv
done


# ~~~~~~~~~~~~~~~~~~~~~~~~~ GC CONTENT FOR CENTROMERES ~~~~~~~~~~~~~~~~~~~~~~~~~
prefix="chicken.v23"
cen_bed_file="ref/chicken.CEN.bed"
module load bedtools/2.31.0
bedtools nuc -fi ref/$prefix.fa -bed $cen_bed_file >centromeres/$prefix.cen_GC_fullInfo.nuc
echo "Chr GCcont" |sed 's/ /\t/g' >centromeres/$prefix.cen_GC.tsv
awk '(NR>1){print $1"\t"$5}' centromeres/$prefix.cen_GC_fullInfo.nuc >> centromeres/$prefix.cen_GC.tsv


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ STATISTICAL TEST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
prefix="chicken.v23"
cen_bed_file="ref/chicken.CEN.bed"
for chr in $(cut -f1 $cen_file |uniq)
do
    echo '#!/bin/bash
    echo "Looking at '$chr'"
    module load bedtools/2.31.0
    c_len=`awk -v c='$chr' '"'"'($1==c){sum+=$3-$2}END{print sum}'"'"' '$cen_file'`
    chrlen=`awk -v c='$chr' '"'"'($1==c){print $2}'"'"' ref/'$prefix'.fa.fai`
 #   poswin=`echo $chrlen"/"$c_len |bc`
 #   minsize=100000
#    if [ "$poswin" -lt "120" ]
 #   then
 #     echo "For '$chr', there will be less than 120 nonovl windows"
 #     slide=`echo "("$chrlen"-"$c_len")/120"|bc`
 #     bedtools makewindows -g <(awk -v c='$chr' '"'"'(c==$1){print}'"'"' ref/'$prefix'.fa.fai) -s $slide -w $c_len |subtractBed -a - -b '$cen_file' |awk -v min=$minsize -v clen=$c_len '"'"'($3-$2>min || $3-$2==clen){print}'"'"' |shuf -n 100 |sort -k2,2n >centromeres/windows/'$prefix'.'$chr'.exclCen.bed
 #   else
 #     echo "For '$chr', there will be more than 120 nonovl windows"
 #     bedtools makewindows -g <(awk -v c='$chr' '"'"'(c==$1){print}'"'"' ref/'$prefix'.fa.fai) -w $c_len |subtractBed -a - -b '$cen_file' |awk -v min=$minsize -v clen=$c_len '"'"'($3-$2>min || $3-$2==clen){print}'"'"' |shuf -n 100 |sort -k2,2n >centromeres/windows/'$prefix'.'$chr'.exclCen.bed
 #   fi
 #   bedtools nuc -fi ref/'$prefix'.fa -bed '$cen_file' >centromeres/windows/'$prefix'.'$chr'.exclCen.with_GC.bed
    rm -f centromeres/windows/'$prefix'.'$chr'.100rand.enrichment.tsv
    tmp="'$chr'"
    while read line
    do
      cat coverage/'${prefix}'.nonB_genome_wide.txt  | while read -r non_b tot dens;
      do
          echo $line |sed "s/ /\t/g" >tmp.$SLURM_JOB_ID.bed
          d=`intersectBed -a tmp.$SLURM_JOB_ID.bed -b final_nonB/'$prefix'.${non_b}.merged.bed -wo |awk -v l=$c_len -v dtot=$dens '"'"'{sum+=$7}END{d=sum/l; frac=d/dtot; print frac}'"'"'`
          tmp=`echo $tmp" "$d`
          echo $tmp >tmp.$SLURM_JOB_ID
      done
      cat tmp.$SLURM_JOB_ID |sed "s/ /\t/g" >>centromeres/windows/'$prefix'.'$chr'.100rand.enrichment.tsv
    done <centromeres/windows/'$prefix'.'$chr'.exclCen.bed
    '| sbatch -J $schr --ntasks=1 --cpus-per-task=1 --time=5:00:00 --out slurm/job.centromeric_enrichment.$prefix.$chr.%j.out
done
# Check if there are more or less than 120 windows
cat slurm/job.centromeric_enrichment.chr*.*.out |grep "there will be"

# Merge background densities and GC
echo "Window Chr APR DR STR IR TRI G4 Z Any" |sed 's/ /\t/g'  >centromeres/$prefix.background_enrichment_merged.tsv
for file in $(ls centromeres/windows/$prefix.chr*.100rand.enrichment.tsv)
do
  awk '{print NR,$0}' $file |sed 's/ /\t/g' >>centromeres/$prefix.background_enrichment_merged.tsv
done
echo "Window Chr GCcont" |sed 's/ /\t/g'  >centromeres/$prefix.background_GC.tsv
for file in $(ls centromeres/windows/$prefix.chr*.exclCen.with_GC.bed)
do
  awk '(NR>1){rownum=NR-1; print rownum,$1,$5}' $file |sed 's/ /\t/g' >>centromeres/$prefix.background_GC.tsv
done

# The same but compared to the per chromosome densities instead of genome-wide
prefix="chicken.v23"
cen_file="ref/chicken.CEN.bed"
for chr in $(cut -f1 $cen_file |uniq)
do
    echo '#!/bin/bash
    echo "Looking at '$chr'"
    module load bedtools/2.31.0
    c_len=`awk -v c='$chr' '"'"'($1==c){sum+=$3-$2}END{print sum}'"'"' '$cen_file'`
    echo "length is $c_len"
    tmp="'$chr'"
    rm -f centromeres/windows/'$prefix'.'$chr'.100rand.enrichmentCHR.tsv
    while read line
    do
      echo "looking at line $line"
      echo $line |sed "s/ /\t/g" >tmp.$SLURM_JOB_ID.bed
      cat densities/'${prefix}'.nonB_per_chrom.txt |awk -v c='$chr' '"'"'(c==$1){print}'"'"' |grep -v "ALL" | while read -r c non_b tot dens;
      do
          d=`intersectBed -a tmp.$SLURM_JOB_ID.bed -b final_nonB/'$prefix'.${non_b}.merged.bed -wo |awk -v l=$c_len -v dtot=$dens '"'"'{sum+=$7}END{d=sum/l; frac=d/dtot; print frac}'"'"'`
          tmp=`echo $tmp" "$d`
          echo $tmp >tmp.$SLURM_JOB_ID
      done
      cat tmp.$SLURM_JOB_ID |sed "s/ /\t/g" >>centromeres/windows/'$prefix'.'$chr'.100rand.enrichmentCHR.tsv
    done <centromeres/windows/'$prefix'.'$chr'.exclCen.bed
    '| sbatch -J $chr --ntasks=1 --cpus-per-task=1 --time=5:00:00 --out slurm/job.centromeric_enrichmentCHR.$chr.%j.out
done
# Merge
echo "Window Chr APR DR G4 IR MR TRI STR Z" |sed 's/ /\t/g'  >centromeres/$prefix.background_enrichmentCHR_merged.tsv
for file in $(ls centromeres/windows/$prefix.chr*.100rand.enrichmentCHR.tsv)
do
  awk '{print NR,$0}' $file |sed 's/ /\t/g' >>centromeres/$prefix.background_enrichmentCHR_merged.tsv
done


########################### ENRICHMENT IN GENES ###############################
# September 4, 2025

# Gene annotation in input/
# Couldn't repeat the longest transcript analysis I did for zebra finch because
# the nomenclature kept changing over the chromosomes for chicken, and it was
# too tricky to include all different namings.

# Instead, I just extract all different categories, and remove overlaps.
module load bedtools/2.31.0

# CDS FIRST:
prefix="chicken.v23"
cat input/$prefix.gff |\
  awk -F'\t' '($3=="CDS"){match($9, /Parent=[^;]+/, type);
  gsub(/Parent=/, "", type[0]); s=$4-1; print $1,s,$5,type[0],".",$7}' |sed 's/ /\t/g' >functional/annotation/$prefix.CDS.raw.bed
# UTR5
cat input/$prefix.gff |\
awk '($3=="five_prime_UTR"){match($9, /Parent=[^;]+/, type);
gsub(/Parent=/, "", type[0]); s=$4-1; print $1,s,$5,type[0],".",$7}' |sed 's/ /\t/g' |subtractBed -a - -b functional/annotation/$prefix.CDS.raw.bed >functional/annotation/$prefix.UTR5.raw.bed
# UTR3
cat input/$prefix.gff |\
awk '($3=="three_prime_UTR"){match($9, /Parent=[^;]+/, type);
gsub(/Parent=/, "", type[0]); s=$4-1; print $1,s,$5,type[0],".",$7}' |sed 's/ /\t/g' |subtractBed -a - -b functional/annotation/$prefix.CDS.raw.bed >functional/annotation/$prefix.UTR3.raw.bed
# PROMOTERS
# Take 1k upstream from 5'UTRs)
cat functional/annotation/$prefix.UTR5.raw.bed |awk -v OFS="\t" '{if(NR==1){chr=$1; tr=$4; str=$6; if($6=="+"){start=$2}else if($6=="-"){start=$3}else{start="error"}}else{if($4==tr){if($6=="-"){start=$3}}else{if(str=="+"){e=start; s=e-1000}else{s=start; e=s+1000}; print chr,s,e,tr,".",str; chr=$1; tr=$4; str=$6; if($6=="+"){start=$2}else if($6=="-"){start=$3}else{start="error"}}}}END{if(str=="+"){e=start; s=e-1000}else{s=start; e=s+1000}; print chr,s,e,tr,".",str;}' |subtractBed -a - -b functional/annotation/$prefix.CDS.raw.bed >functional/annotation/$prefix.promoter.raw.bed
# INTRONS
cat input/$prefix.gff |\
awk '($3=="exon"){match($9, /Parent=[^;]+/, type);
gsub(/Parent=/, "", type[0]); s=$4-1; print $1,s,$5,type[0],".",$7}' |awk -v OFS="\t" '{if(NR==1){chr=$1; end=$3; name=$4; strand=$6}else{if($1==chr && $4==name){print chr,end,$2,name,".",strand; end=$3}else{chr=$1; end=$3; name=$4;strand=$6}}}' |subtractBed -a - -b functional/annotation/$prefix.CDS.raw.bed >functional/annotation/$prefix.introns.raw.bed
# INTERGENIC - take everything in between genes and remove overlapping exons and
# introns and promoters
awk -v OFS='\t' '{if($3=="gene"){s=$4-1; print $1,s,$5}}' input/$prefix.gff |sort -k1,1 -k2,2n |mergeBed -i -|awk -v OFS='\t' '{if(NR==1){chr=$1; end=$3}else{if($1==chr){print chr,end,$2; end=$3}else{chr=$1; end=$3}}}' |subtractBed -a - -b functional/annotation/$prefix.promoter.raw.bed |subtractBed -a - -b functional/annotation/$prefix.CDS.raw.bed |subtractBed -a - -b functional/annotation/$prefix.UTR5.raw.bed |subtractBed -a - -b functional/annotation/$prefix.UTR3.raw.bed |subtractBed -a - -b functional/annotation/$prefix.introns.raw.bed >functional/annotation/$prefix.intergenic.raw.bed



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~ CALCULATE ENRICHMENT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
prefix="chicken.v23"
for class in "promoter" "CDS" "UTR5" "UTR3" "intergenic" "introns"
do
    echo '#!/bin/bash
    module load bedtools/2.31.0
    len=`sort -k1,1 -k2,2n functional/annotation/'$prefix'.'$class'.raw.bed |mergeBed -i - |awk '"'"'{sum+=$3-$2}END{print sum}'"'"'`
    echo "Length of '$class' is $len"
    rm -f tmp.'$prefix'.'$class'
    cat coverage/'${prefix}'.nonB_genome_wide.txt | while read -r non_b tot dens;
    do
        echo "looking at $non_b"
        d=`intersectBed -a <(cut -f1-3 functional/annotation/'$prefix'.'$class'.raw.bed |sort -k1,1 -k2,2n |mergeBed -i -) -b final_nonB/'${prefix}'.${non_b}.merged.bed -wo |awk -v l=$len -v dtot=$dens '"'"'{sum+=$7}END{d=sum/l; frac=d/dtot; print d,frac}'"'"'`
        echo '$class'" "$non_b" "$d >>tmp.'$prefix'.'$class'
    done
    ' | sbatch -J $class --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8G --time=5:00:00 --partition=open --out slurm/functional.$class.%j.out
done
# Merge after
echo "Class nonB Coverage Enrichment_gw" |sed "s/ /\t/g" >functional/$prefix.enrichment.tsv
for class in "CDS" "UTR5" "UTR3" "promoter" "introns" "intergenic"
do
    cat tmp.$prefix.$class |sed "s/ /\t/g" >>functional/$prefix.enrichment.tsv
done

# ~~~~~~~~~~~~~~~~~ CALCULATE ENRICHMENT FOR CHROMOSOME TYPES ~~~~~~~~~~~~~~~~~~
# September 5, 2025

# Chicken chromosomes doesn't have unique names, so we can not use grep!
# Make annotation files specific for each group
prefix="chicken.v23"
for group in "macro" "micro" "microdot"
do
  for class in "promoter" "intergenic" "introns" "CDS" "UTR5" "UTR3"
  do
    for chr in $(awk -v g=$group '($3==g){print $1}' $prefix.groups.txt)
    do
      awk -v c=$chr '($1==c){print}' functional/annotation/$prefix.$class.raw.bed
    done >functional/annotation/$prefix.$class.$group.raw.bed
  done
done

prefix="chicken.v23"
for group in "macro" "micro" "microdot"
do
  for class in "promoter" "intergenic" "introns" "CDS" "UTR5" "UTR3"
  do
      echo '#!/bin/bash
      module load bedtools/2.31.0
      len=`sort -k1,1 -k2,2n functional/annotation/'$prefix'.'$class'.'$group'.raw.bed |mergeBed -i - |awk '"'"'{sum+=$3-$2}END{print sum}'"'"'`
      echo "Length of '$class' is $len"
      rm -f tmp.'$prefix'.'$group'.'$class'
      cat coverage/'${prefix}'.nonB_genome_wide.txt | while read -r non_b tot dens;
      do
          echo "looking at $non_b, with $tot bp, density is $dens"
          d=`intersectBed -a <(sort -k1,1 -k2,2n functional/annotation/'$prefix'.'$class'.'$group'.raw.bed |cut -f1-3 |mergeBed -i -) -b final_nonB/'${prefix}'.${non_b}.merged.bed -wo |awk -v l=$len -v dtot=$dens '"'"'{sum+=$7}END{d=sum/l; frac=d/dtot; print d,frac}'"'"'`
          echo '$group'" "'$class'" "$non_b" "$d >>tmp.'$prefix'.'$group'.'$class'
      done
      ' | sbatch -J $class --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8G --time=5:00:00 --partition=open --out slurm/functional.$prefix.$group.$class.%j.out
  done
done

# Merge the tmp files
echo "Group Class nonB Coverage Enrichment_gw" |sed "s/ /\t/g" >functional/$prefix.enrichment.groups.tsv
for group in "macro" "micro" "microdot"
do
  for class in  "promoter" "intergenic" "introns" "CDS" "UTR5" "UTR3"
  do
    cat tmp.$prefix.$group.$class |sed "s/ /\t/g" >>functional/$prefix.enrichment.groups.tsv
  done
done



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CHECK SIGNIFICANCE OF THE ENRICHMENT BARS
# Spet 12, 2025
# Subsample the functional regions to get a distribution of the enrichment
# values. If distribution overlaps with 1 (genome average), we deem this
# comparison NOT to be significally different from the genome average.

# Subsample 50 percent of the regions 100 times
prefix="chicken.v23"
mkdir -p functional/resample/
perc=50
for group in "macro" "micro" "microdot"
do
  for class in "promoter" "intergenic" "introns" "CDS" "UTR5" "UTR3"
  do
    num=`cat functional/annotation/$prefix.$class.$group.raw.bed |awk -v p=$perc '{}END{num=int(NR*(p/100)); print num}' `
    echo "we will extract $num rows from $class on $group"
    # Subsample the sequences 100 times
    echo '#!/bin/bash
    module load bedtools/2.31.0
    for i in {1..100}
    do
        echo "working on $i..."
        cat functional/annotation/'$prefix'.'$class'.'$group'.raw.bed |shuf -n '$num' |sort -k1,1 -k2,2n |mergeBed -i -  >functional/resample/'$prefix'.Resamp.'$perc'perc.'$class'.'$group'.$i.bed
        l=`wc -l functional/resample/'$prefix'.Resamp.'$perc'perc.'$class'.'$group'.$i.bed`
        echo "...extracting $l lines"
    done
    ' |sbatch -J $class.$group -o slurm/job.resamp.$prefix.${perc}perc.$class.$group.%j.out --requeue --ntasks=1 --mem-per-cpu=1G --cpus-per-task=1 --time=1:00:00 --partition=open
  done
done

# Calculate the enrichment for each of the subsamples
prefix="chicken.v23"
for subset in "50perc"
do
  for group in "macro" "micro" "microdot"
  do
    for class in "promoter" "intergenic" "introns" "CDS" "UTR5" "UTR3"
    do
        echo '#!/bin/bash
        module load bedtools/2.31.0
        for i in {1..100}
        do
            rm -f functional/resample/tmp.'$prefix'.'$subset'.'$class'.'$group'.$i.txt
            len=`awk '"'"'{sum+=$3-$2}END{print sum}'"'"' functional/resample/'$prefix'.Resamp.'$subset'.'$class'.'$group'.$i.bed`
            cat coverage/'${prefix}'.nonB_genome_wide.txt |grep -v "ALL" | while read -r non_b tot dens;
            do
                d=`intersectBed -a functional/resample/'$prefix'.Resamp.'$subset'.'$class'.'$group'.$i.bed -b final_nonB/'$prefix'.${non_b}.merged.bed -wo |awk -v l=$len -v dtot=$dens '"'"'{sum+=$7}END{d=sum/l; frac=d/dtot; print frac}'"'"'`
                echo $i '$group' '$class' $non_b $d >>functional/resample/tmp.'$prefix'.'$subset'.'$class'.'$group'.$i.txt
            done
        done
        ' |sbatch -J $class.$group -o slurm/job.resamp-enrich.$prefix.$subset.$class.$group.%j.out --requeue --ntasks=1 --cpus-per-task=1 --mem-per-cpu=4G --time=10:00:00 --partition=open
    done
  done
done


# Merge the results
rep="100rep"
subset="50perc" #subset="10perc"
echo "Rep Group Class nonB Enrichment" |sed "s/ /\t/g" >functional/summary.$prefix.$subset.$rep.txt
cat functional/resample/tmp.$prefix.$subset.*.*.{1..100}.txt |sed "s/ /\t/g" >>functional/summary.$prefix.$subset.$rep.txt

# We want to find the min and max for each category, but remove the ~5% most
# extreme (meaning we remove the top and bottom two values, saving a 96% CI)
rep="100rep"
subset="50perc" #subset="10perc"
echo "Group Class nonB Min Max" |sed "s/ /\t/g" >functional/$prefix.$subset.$rep.minmax.96CI.txt
for group in "macro" "micro" "microdot"
do
  for class in "promoter" "intergenic" "introns" "CDS" "UTR5" "UTR3"
  do
      for nonb in APR DR G4 IR TRI STR Z
      do
          min=`grep $class functional/summary.$prefix.$subset.$rep.txt |awk -v nb=$nonb -v g=$group -v OFS="\t" '($4==nb && $2==g){print $5}' |sort -n |head -n3 |tail -n1`
          max=`grep $class functional/summary.$prefix.$subset.$rep.txt |awk -v nb=$nonb -v g=$group -v OFS="\t" '($4==nb && $2==g){print $5}' |sort -n |tail -n3 |head -n1`
          echo $group $class $nonb $min $max |sed "s/ /\t/g" >>functional/$prefix.$subset.$rep.minmax.96CI.txt
      done
    done
done


############################## INVESTIGATING G4s ###############################

# Extract the sequences of the predicted G4 regions 
prefix="chicken.v23"
module load bedtools/2.31.0
bedtools getfasta -fi ref/$prefix.fa -bed final_nonB/$prefix.G4.bed -fo g4discovery/$prefix.G4.fasta

# Reverse complement sequences on the rev strand and print a list with
# position and sequence on one line 
prefix="chicken.v23"
file=g4discovery/$prefix.G4.fasta
  awk '/^>/ {h=substr($0,2); next}; {seq=toupper($0); if (substr(seq,1,1)=="C") {
      # reverse complement
      rc=""
      for (i=length(seq); i>0; i--) {
          base=substr(seq,i,1)
          if (base=="A") b="T"
          else if (base=="T") b="A"
          else if (base=="C") b="G"
          else if (base=="G") b="C"
          else b=base
          rc = rc b
      }
      seq=rc
  }
  print h "\t" seq
  }' $file  >g4discovery/$prefix.revcomp.txt

# Sort and count the sequence 
cut -f2 g4discovery/$prefix.revcomp.txt| sort | uniq -c | sort -k1,1nr |awk '{print $2"\t"$1}' >g4discovery/$prefix.g4.sequence.counts.txt

# Check how many of the 100 most common G4s are also found in zebra finch
head -n100 g4discovery/$prefix.g4.sequence.counts.txt |cut -f1 |sort |join -1 1 -2 1 - <(sort -k1,1 test.common.g4s.txt) |wc
#     47      94    1166



# HAVE NOT RUN THE BELOW CODE 

# For each common G4, find positions in the genome and chose one to extract the
# scores from the original g4Discovery file 
for chr in $(cat m*cro.txt |grep "chr2_mat")
do
    echo "Processing chromosome: "$chr
    cat test.common.g4s.txt |head -n30 |while read line
    do 
        seq=$(echo $line |cut -f1 -d" ")
        #echo "Processing sequence: "$seq
        no=`grep $seq RESULTS/g4discovery/$chr.g4discovery.revcomp.txt |wc -l`
        start=`grep $seq RESULTS/g4discovery/$chr.g4discovery.revcomp.txt |awk -F':|-|\t' -v seq=$seq '($4==seq){print $2}' - |shuf |head -n1`
        info=`zcat RESULTS/g4discovery/$chr.bed.gz | awk -v st=$start '($2==st){print $0}'`
        echo "  "$info"  "$no" "$seq
    done 
done 

# And make files with all positions for each common G4, to overlap with regions 
mkdir RESULTS/g4discovery/common_g4s/
cat test.common.g4s.txt |head -n30 |while read line
do 
    seq=$(echo $line |cut -f1 -d" ")
    rm -f RESULTS/g4discovery/common_g4s/$seq.bed
    for chr in $(cat ref/bTaeGut7v0.4_MT_rDNA.fa.fai |cut -f1)
    do 
        grep $seq RESULTS/g4discovery/$chr.g4discovery.revcomp.txt |awk -F':|-|\t' '{print $1"\t"$2"\t"$3"\t"$4}' >>RESULTS/g4discovery/common_g4s/$seq.bed
    done 
done
