################################################################################
### ANNOTATION OF NON-B DNA IN THE ZEBRA FINCH T2T GENOME AND OTHER BIRD GENOMES
### FOR COMPARISON. THIS FILE CONTAINS THE FIRST CODE, INCLUDING DOWNLOADING AND 
### SOFTWARE REQUIREMENTS.
### CODE WRITTEN BY LINNÉA SMEDS

# All commands are assumed to be run from within the downloaded github repo.
# Time and/or memory consuming jobs were run on a slurm cluster. The code 
# within the slurm scripts can be copied and run on any machine with enough
# RAM/cores.

##### TABLE OF CONTENTS
# REQUIREMENTS 
# DOWNLOAD FILES 
# PROCESS FILES 
# CHECKING NCBI FOR BIRD GENOMES 

################################# REQUIREMENTS #################################
# A LIST OF ALL SOFTWARE USED:
# bedtools/2.31.0
# gfa (downloaded Feb 4, 2025 from https://github.com/abcsFrederick/non-B_gfa)
# Z-DNA Hunter (webserver: https://bioinformatics.ibp.cz/#/analyse/zdna)
# g4discovery (downloaded Nov 20 2025 from https://github.com/saswat-km/g4Discovery.PanSN.git,
#   also needs R packages seqinr, BiocManager, rtracklayer and Biostrings)
# Quadron (dockerized version 1.0.0 from https://hub.docker.com/r/kxk302/quadron)
# circos/0.69-9
# python/3.11.2
# samtools/1.19.2
# convert2bed/2.4.41
# gffread-0.12.7
# Orthofinder (run through anaconda/2023.09)
# R/4.5.1 with tidyverse/2.0.0 and patchwork_1.3.2


################################ DOWNLOAD FILES #################################
# Download Zebra Finch genome and annotation files from GenomeArk
mkdir ref
cd ref
wget https://genomeark.s3.amazonaws.com/species/Taeniopygia_guttata/bTaeGut7/assembly_verkko_0.1/manual_curation/bTaeGut7v0.4/bTaeGut7v0.4_MT_rDNA.fa.gz
wget https://genomeark.s3.amazonaws.com/species/Taeniopygia_guttata/bTaeGut7/manuscript/annotations/centromeres/bTaeGut7v0.4_MT_rDNA.centromere_detector.v0.1.gff
wget https://genomeark.s3.amazonaws.com/species/Taeniopygia_guttata/bTaeGut7/manuscript/annotations/repeats/bTaeGut7v0.4_MT_rDNA.EDTA2.v0.2.gtf.gz
wget https://genomeark.s3.amazonaws.com/species/Taeniopygia_guttata/bTaeGut7/manuscript/annotations/repeats/bTaeGut7v0.4_MT_rDNA.trf.sorted.v0.1.bed
wget https://genomeark.s3.amazonaws.com/species/Taeniopygia_guttata/bTaeGut7/manuscript/annotations/repeats/bTaeGut7v0.4_MT_rDNA.satellome.v0.1.bed
wget https://genomeark.s3.amazonaws.com/species/Taeniopygia_guttata/bTaeGut7/manuscript/annotations/genes/bTaeGut7v0.4_MT_rDNA.EGAPx.v0.1.gtf.gz
wget https://genomeark.s3.amazonaws.com/species/Taeniopygia_guttata/bTaeGut7/manuscript/annotations/methylation/bTaeGut7v0.4_MT_rDNA.PBmethylation.v0.1.bw
wget https://genomeark.s3.amazonaws.com/species/Taeniopygia_guttata/bTaeGut7/manuscript/annotations/3D/bTaeGut7v0.4_MT_rDNA.Cooltools.E1.200kbp.flipped.dip.collated.v0.1.bed
wget https://genomeark.s3.amazonaws.com/species/Taeniopygia_guttata/bTaeGut7/manuscript/annotations/3D/bTaeGut7v0.4_MT_rDNA.Cooltools.E1.10kbp.flipped.dip.collated.v0.1.bw
wget https://genomeark.s3.amazonaws.com/species/Taeniopygia_guttata/bTaeGut7/manuscript/annotations/PUR/bTaeGut7v0.4_MT_rDNA_MATonly_vs_GCF_000151805.1.PUR.fastga.v0.1.bed
wget https://genomeark.s3.amazonaws.com/species/Taeniopygia_guttata/bTaeGut7/assembly_verkko_0.1/manual_curation/bTaeGut7v0.4/mapping/v0.4_dip_hifi/v0.4_dip_hifi.pri.cov.wig
wget https://genomeark.s3.amazonaws.com/species/Taeniopygia_guttata/bTaeGut7/assembly_verkko_0.1/manual_curation/bTaeGut7v0.4/mapping/v0.4_dip_ont/v0.4_dip_ont.pri.cov.wig
gunzip bTaeGut7v0.4_MT_rDNA.fa.gz


# Download other bird genomes for comparison:
# Download chicken genome and annotation from dropbox
# https://www.dropbox.com/scl/fo/plq2tm2w9lzlk0ua1rzph/h?rlkey=l6z3rgmjs7ec9azun8nundnzl&e=1&dl=0
# Zebra finch maternal+Z for gene homologies and functional annotation 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/771/995/GCF_048771995.1_bTaeGut7.mat/GCF_048771995.1_bTaeGut7.mat_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/771/995/GCF_048771995.1_bTaeGut7.mat/GCF_048771995.1_bTaeGut7.mat_assembly_report.txt
# Download Anna's hummingbird from NCBI: 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/957/555/GCF_003957555.1_bCalAnn1_v1.p/GCF_003957555.1_bCalAnn1_v1.p_assembly_report.txt
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/957/555/GCF_003957555.1_bCalAnn1_v1.p/GCF_003957555.1_bCalAnn1_v1.p_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/957/555/GCF_003957555.1_bCalAnn1_v1.p/GCF_003957555.1_bCalAnn1_v1.p_genomic.gff.gz
# Band-tailed pigeon
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/037/038/585/GCF_037038585.1_bPatFas1.hap1/GCF_037038585.1_bPatFas1.hap1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/037/038/585/GCF_037038585.1_bPatFas1.hap1/GCF_037038585.1_bPatFas1.hap1_assembly_report.txt
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/037/038/585/GCF_037038585.1_bPatFas1.hap1/GCF_037038585.1_bPatFas1.hap1_genomic.gff.gz
# Emu
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/370/855/GCF_036370855.1_bDroNov1.hap1/GCF_036370855.1_bDroNov1.hap1_assembly_stats.txt
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/370/855/GCF_036370855.1_bDroNov1.hap1/GCF_036370855.1_bDroNov1.hap1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/370/855/GCF_036370855.1_bDroNov1.hap1/GCF_036370855.1_bDroNov1.hap1_genomic.gff.gz
# Pekin Duck
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/047/663/525/GCF_047663525.1_IASCAAS_PekinDuck_T2T/GCF_047663525.1_IASCAAS_PekinDuck_T2T_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/047/663/525/GCF_047663525.1_IASCAAS_PekinDuck_T2T/GCF_047663525.1_IASCAAS_PekinDuck_T2T_assembly_report.txt
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/047/663/525/GCF_047663525.1_IASCAAS_PekinDuck_T2T/GCF_047663525.1_IASCAAS_PekinDuck_T2T_genomic.gff.gz
# Great Bustard 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/026/413/225/GCA_026413225.1_OTswu/GCA_026413225.1_OTswu_assembly_report.txt 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/026/413/225/GCA_026413225.1_OTswu/GCA_026413225.1_OTswu_genomic.fna.gz 
# annotation from https://figshare.com/articles/dataset/Figure_source_data/23650419, 
# had to add "chr" and remove blank lines using "|grep -v ^$ |awk -v OFS="\t" '{$1="chr"$1; print}'"
# Ural owl 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/047/716/275/GCF_047716275.1_bStrUra1/GCF_047716275.1_bStrUra1_assembly_report.txt
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/047/716/275/GCF_047716275.1_bStrUra1/GCF_047716275.1_bStrUra1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/047/716/275/GCF_047716275.1_bStrUra1/GCF_047716275.1_bStrUra1_genomic.gff.gz
cd ..
 
for file in ref/*.gz
do
  gunzip $file
done

################################ PROCESS FILES #################################

# Make a list of species and file name prefixes
echo "zebra_finch bTaeGut7v0.4_MT_rDNA bTaeGut7v0.4_MT_rDNA
chicken chicken.v23 chicken.v23
annas_hummingbird GCF_003957555.1_bCalAnn1_v1.p bCalAnn1_v1.p
great_bustard GCA_026413225.1_OTswu OTswu
ural_owl GCF_047716275.1_bStrUra1 bStrUra1
bandtailed_pigeon GCF_037038585.1_bPatFas1.hap1 bPatFas1.hap1
emu GCF_036370855.1_bDroNov1.hap1 bDroNov1.hap1
peking_duck GCF_047663525.1_IASCAAS_PekinDuck_T2T IASCAAS_PekinDuck_T2T" > helpfiles/species_list.txt

# Change the annotations files from NCBI to readable chromosome names
module load python/3.11.2 
module load samtools/1.19.2
cat helpfiles/species_list.txt |tail -n3 |while read -r sp longname prefix
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
    awk 'NR==FNR{a[$1]=$2; next} /^>/ {key = substr($1, 2); $0 = ">"a[key]}1' ref/$prefix.conversion-table.txt ref/${longname}_genomic.fna > ref/$prefix.fa
    # Index genome 
    samtools faidx ref/$prefix.fa
done 

# Use RefSeq annotation for zebra finch primary assembly (maternal+Z)
longname="GCF_048771995.1_bTaeGut7.mat"
prefix="bTaeGut7.mat"
grep -v "#" ref/${longname}_assembly_report.txt |cut -f1,3,4,5,7 |awk '{new=$1; if($5=="na"){old=$4}else{old=$5}; print old"\t"new}' >ref/$prefix.conversion-table.txt
python3 T2T_bird_nonB/python/replace_chromosome_names.py -g ref/${longname}_genomic.gff -t ref/$prefix.conversion-table.txt -o ref/$prefix.raw.gff
awk '($3 != "region")' ref/$prefix.raw.gff > ref/${prefix}.gff  
cp ref/${prefix}.gff  ref/bTaeGut7v0.4_MT_rDNA.gff



################################# FIND HOMOLOGY ################################

# Get protein and cds sequences 
cat helpfiles/species_list.txt |while read -r sp longname prefix
do
  ~/software/gffread-0.12.7.Linux_x86_64/gffread ref/$prefix.gff \
  -g ref/$prefix.fa \
  -y ref/$prefix.proteins.faa \
  -x ref/$prefix.cds.fna \
  -F
done 

# Move to orthofinder and check there are no dots (orthofinder cannot handle ".") 
mkdir orthofinder/
cat helpfiles/species_list.txt |while read -r sp longname prefix
do
  awk '{if(/>/){print $0}else{gsub(/\./,"",$0); print $0}}' ref/$prefix.proteins.faa >orthofinder/$prefix.noDots.proteins.faa
done

#Orthofinder is installed in a virtual environment created in ~/software/, with
module load anaconda/2023.09
cd ~/software/
conda create -n of3_env python=3.10
conda activate of3_env
conda install -c bioconda orthofinder

# Run orthofinder
cd /storage/group/kdm16/default/lbs5874/ZebraFinch
echo '#!/bin/bash
orthofinder -f orthofinder/ -M msa -T fasttree
' |sbatch -J orthofinder --ntasks=1 --cpus-per-task=6 --time=2-00:00:00 -o slurm/job.orthofinder.%j.out

# Get map files between gene ID - chromosome 
cat helpfiles/species_list.txt |tail -n3 |while read -r sp longname prefix
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

# Use python script to link chromosome between species. 
module load python/3.11.2 
pip install pandas
python3 T2T_bird_nonB/python/orthofinder_chr_homology.py
# The input files are hardcoded in this script. The output is saved to 
# orthofinder/sp1_vs_sp2_chr_homology.tsv  




######################## CHECKING NCBI FOR BIRD GENOMES ########################
# Jan 29, 2026

# Search for "Aves" On NCBI, download table with the following columns:
# Assembly Accession      Assembly Name   Organism Name   Organism Taxonomic ID   Organism Infraspecific Names Breed      Organism Infraspecific Names Strain     Organism Infraspecific Names Cultivar      Organism Infraspecific Names Ecotype    Organism Infraspecific Names Isolate    Organism Infraspecific Names Sex        Annotation Name Assembly Level  Assembly Release Da

# Saved as ~/Downloads/ncbi_dataset.tsv

# GET UNIQUE TAXON ID:
cut -f5 ~/Downloads/ncbi_dataset-7.tsv |sort |uniq |wc
    1766    3594   34751


# Download one version with less info and extract completeness info:
#Assembly Accession      Assembly Name   Organism Taxonomic ID   Assembly Stats Total Number of Chromosomes      Assembly Level  WGS project accession
cut -f5,15  ~/Downloads/ncbi_dataset-7.tsv |uniq |sort |uniq |grep Complete >NCBI_complete.txt
cut -f5,15  ~/Downloads/ncbi_dataset-7.tsv |uniq |sort |uniq |grep Chromosome >NCBI_chromosome.txt
cut -f5,15  ~/Downloads/ncbi_dataset-7.tsv |uniq |sort |uniq |grep Contig >NCBI_contig.txt
cut -f5,15  ~/Downloads/ncbi_dataset-7.tsv |uniq |sort |uniq |grep Scaffold >NCBI_scaffold.txt

# Number complete:
wc -l NCBI_complete.txt
#       3 NCBI_complete.txt
# Number chromosome (remove complete)
join -v1 -1 1 -2 1 <(cut -f1 NCBI_chromosome.txt) <(cut -f1 NCBI_complete.txt) |wc -l
# 269
# Number of scaffold level (remove chromosome and complete)
join -v1 -1 1 -2 1 <(cut -f1 NCBI_scaffold.txt) <(cut -f1 NCBI_chromosome.txt) |wc -l
#    1316    1302    8916
