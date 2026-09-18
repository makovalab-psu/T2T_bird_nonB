###############################################################################
### PDAL-SEQ PREP 
# This code just prep the data - the actual analysis was done by Jacob Sieg,
# link to code: https://github.com/makovalab-psu/ZebraFinch_PDAL-Seq


########################### PREP REFERENCE ASSEMBLY ############################
# For the PDAL-Seq analysis we need a haploid assembly. Choose mat+Z since it 
# seems to be the primary assembly on NCBI.
module load samtools/1.21
list=`tail -n+40 ref/bTaeGut7v0.4_MT_rDNA.fa.fai |cut -f1 |tr "\n" " "`
samtools faidx ref/bTaeGut7v0.4_MT_rDNA.fa $list >ref/bTaeGut7v0.4_MT_rDNA.matZ.fa
samtools faidx ref/bTaeGut7v0.4_MT_rDNA.matZ.fa

# Also remove paternal haplotype from annotation files:
awk '$1 !~ /_pat$/ || $1 == "chrZ_pat"' ref/bTaeGut7v0.4_MT_rDNA.CEN.bed > ref/bTaeGut7v0.4_MT_rDNA.matZ.CEN.bed
awk '$1 !~ /_pat$/ || $1 == "chrZ_pat"' annotation/bTaeGut7v0.4_MT_rDNA.PBmethylation.v0.1.bed >annotation/bTaeGut7v0.4_MT_rDNA.matZ.PBmethylation.v0.1.bed
awk '$1 !~ /_pat$/ || $1 == "chrZ_pat"' ref/bTaeGut7v0.4_MT_rDNA.gff >ref/bTaeGut7v0.4_MT_rDNA.matZ.gff
awk '$1 !~ /_pat$/ || $1 == "chrZ_pat"' annotation/bTaeGut7v0.4_MT_rDNA.EDTA2.v0.2.bed >annotation/bTaeGut7v0.4_MT_rDNA.matZ.EDTA2.v0.2.bed
awk '$1 !~ /_pat$/ || $1 == "chrZ_pat"' annotation/bTaeGut7v0.4_MT_rDNA.Satellites.bed >annotation/bTaeGut7v0.4_MT_rDNA.matZ.Satellites.bed
awk '$1 !~ /_pat$/ || $1 == "chrZ_pat"' annotation/bTaeGut7v0.4_MT_rDNA.TRF_withMers.bed >annotation/bTaeGut7v0.4_MT_rDNA.matZ.TRF_withMers.bed
awk '$1 !~ /_pat$/ || $1 == "chrZ_pat"' annotation/bTaeGut7v0.4_MT_rDNA.Cooltools.v0.2.E1.10Kb.flipped.dip.collated.AB.bed |cut -f1-4 >annotation/bTaeGut7v0.4_MT_rDNA.matZ.ABCompart.bed
# And all the separate annotation files:
for type in "promoter" "CDS" "intergenic" "introns" "UTR5" "UTR3" "lncrna"
do
 awk '$1 !~ /_pat$/ || $1 == "chrZ_pat"' annotation/bTaeGut7v0.4_MT_rDNA.$type.bed >annotation/bTaeGut7v0.4_MT_rDNA.matZ.$type.bed
done



############################# DOWNLOAD RNA-SEQ DATA ############################
# RNA-seq data from the original cell line publication 

cd DATA
echo '#!/bin/bash
~/software/sratoolkit.3.2.1-ubuntu64/bin/fasterq-dump --split-files SRR17849680
~/software/sratoolkit.3.2.1-ubuntu64/bin/fasterq-dump --split-files SRR17849681
' |sbatch -J sra  --ntasks=1 --cpus-per-task=1 --mem-per-cpu=2G --time=1-00:00:00 

# Zip files 
for fastq in DATA/SRR17849681*.fastq
do
    ls $fastq 
    echo '#!/bin/bash
    gzip '$fastq'
    ' |sbatch -J zip  --ntasks=1 --cpus-per-task=1 --mem-per-cpu=2G --time=1-00:00:00 
done 
