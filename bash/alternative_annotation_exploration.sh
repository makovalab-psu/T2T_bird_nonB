################################################################################
### TESTING OTHER ANNOTATION SOFTWARE FOR N0N-B DNA MOTIFS 
### CODE WRITTEN BY LINNÉA SMEDS

##### TABLE OF CONTENTS
# 

############################## Z-DNA WITH ZSEEKER ##############################
# November 18, 2025

# Using the ZSeeker python script due to problems with installation of the binary
mkdir -p RESULTS/zseeker
cat helpfiles/species_list.txt |grep -v "zebra_finch" |while read -r sp longname prefix
do
    fasta="ref/"$prefix".fa"
    echo '#!/bin/bash
    module load python/3.11.2
    python3 ~/software/ZSeeker-zseeker_cli_v1/zseeker/zseeker.py --fasta '$fasta' --output_dir zseeker
    ' | sbatch -J zseeker  --ntasks=1 --cpus-per-task=4 --mem-per-cpu=4G --time=1-00:00:00  -o slurm/job.zseeker.%j.%N.out
done 

#Positions are 0-based, need to add +1 to the ends to get bed format. 
awk -v FS=',' -v OFS="\t" '(NR>1){e=$3+1; print $1,int($2),e,$4}' RESULTS/zseeker/bTaeGut7v0.4_MT_rDNA_zdna_score.csv  >final_nonB/bTaeGut7v0.4_MT_rDNA.Zseeker.bed


############################## Z-DNA WITH Z-HUNTER #############################
# November 18, 2025

# Ran ZDNA-hunter on the webserver, and saved each chromosomes locally in the 
# subfolder ZDNAHunter/ as ${chr}.bedgraph

# Just realized the maternal and Z chromosomes are called "chr1mat" etc instead
# of chr1_mat.. need to rename each row in the bedgraph files
for file in ZDNAHunter/chr*_mat.bedgraph
do
    ls $file
    sed 's/mat/_mat/g' $file >temp.bedgraph
    mv temp.bedgraph $file
done    
sed 's/pat/_pat/g' ZDNAHunter/chrZ_pat.bedgraph >temp.bedgraph
mv temp.bedgraph ZDNAHunter/chrZ_pat.bedgraph

# Combine all chromosomes into one bed file (make sure to use TAB separaptors)
rm -f final_nonB/bTaeGut7v0.4_MT_rDNA.ZDNAHunter.bed
for chr in $(cat ref/bTaeGut7v0.4_MT_rDNA.fa.fai |cut -f1)
do 
    tail -n+3 ZDNAHunter/${chr}.bedgraph |awk -v OFS="\t" '{print $1,$2,$3,$4}' >>final_nonB/bTaeGut7v0.4_MT_rDNA.ZDNAHunter.bed
done


############################# G4 WITH G4DISCOVERY ##############################
# November 20, 2025

module load python/3.11.2
module load r/4.5.0
R 
install.packages("seqinr")
# whe nasked if you want to do create a custom location: press yes, choose mirror
install.packages("BiocManager")
BiocManager::install(c("pqsfinder", "rtracklayer", "Biostrings"))

cd ~/software/ 
git clone https://github.com/saswat-km/g4Discovery.PanSN.git

cd /storage/group/kdm16/default/lbs5874/ZebraFinch
prefix="bTaeGut7v0.4_MT_rDNA"
mkdir g4discovery/$prefix

# NOTES: The script assumes the g4discovery.PanSN folder is in your current 
# directory. I manually modified theg4DiscoveryFuncs.py script to add the full 
# path to the R script run_pqsfinder.R. The script also assumes zipped input.

# Prep input files 
mv Quadron_annotation/split_fasta/*.fa ref/split_fasta/bTaeGut7v0.4_MT_rDNA/
for file in $(ls ref/split_fasta/bTaeGut7v0.4_MT_rDNA/*fa); do 
    ls $file 
    gzip $file 
done 


# Start g4Discovery for all fasta files 
#for file in $(ls ref/split_fasta/bTaeGut7v0.4_MT_rDNA/chr*_*fa.gz |grep -v "chr1_" |grep -v "chr2_" |grep -v "chr3_") 
#or file in $(ls ref/split_fasta/bTaeGut7v0.4_MT_rDNA/chr{1..3}_*fa.gz) 
for file in ref/split_fasta/bTaeGut7v0.4_MT_rDNA/rDNA_*fa.gz
do 
    chr=$(basename ${file%.fa.gz})
    echo '#!/bin/bash
    module load python/3.11.2
    module load r/4.5.0
    python3 ~/software/g4Discovery.PanSN/src/g4Discovery.py -fa '$file' -o g4discovery/'$prefix'/'$chr'.bed
    ' |sbatch -J g4D.$chr  --ntasks=1 --cpus-per-task=1 --mem-per-cpu=2G --time=2-00:00:00 -o slurm/job.g4D.$chr.%j.%N.out
done
# The three biggest chromosomes OOM failed, I restarted them with 3 cores 
# (the others used 1 core)

# Combine results into one bed file
prefix="bTaeGut7v0.4_MT_rDNA"
rm -f final_nonB/$prefix.G4.bed
for chr in $(cat ref/$prefix.fa.fai |cut -f1)
do 
    zcat g4discovery/$prefix/$chr.bed.gz >>final_nonB/$prefix.G4.bed
done 

####################### A DEEPER LOOK INTO THE G4 SEQUENCES ####################

# Extract the sequences of the predicted G4 regions 
module load bedtools/2.31.0
for file in $(ls RESULTS/g4discovery/chrZ_*.bed.gz) 
do 
    prefix=$(basename ${file%.bed.gz})
    echo $prefix
    zcat $file | bedtools getfasta -fi ref/bTaeGut7v0.4_MT_rDNA.fa -bed - -fo RESULTS/g4discovery/$prefix.g4discovery.fasta
done

# Reverse complement sequences on the rev strand and print a list with
# position and sequence on one line 
for file in $(ls RESULTS/g4discovery/chr*.fasta) 
do 
    prefix=$(basename ${file%.fasta})
    echo $prefix
    awk '/^>/ {h=substr($0,2); next}; {seq=$0; if (substr(seq,1,1)=="C") {
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
    }' $file  >RESULTS/g4discovery/$prefix.revcomp.txt
done


# For each chromosome, sort and count the sequences 
for file in $(ls RESULTS/g4discovery/chrZ_*.g4discovery.revcomp.txt) 
do 
    prefix=$(basename ${file%.g4discovery.revcomp.txt})
    echo $prefix
    cut -f2 $file | sort | uniq -c | sort -k1,1nr |awk '($1>1){print}' >RESULTS/g4discovery/$prefix.g4.sequence.counts.txt
done 


# Merge all species and count the ocurrences of G4s in the diploid genome
cat RESULTS/g4discovery/*.g4.sequence.counts.txt |awk '{print $2"\t"$1}' |sort -k1,1 |awk '{if(NR==1){seq=$1; sum=$2}else{if($1==seq){sum+=$2}else{print seq"\t"sum; seq=$1; sum=$2}}}END{print seq"\t"sum}' |sort -k2,2nr >test.common.g4s.txt

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

# Overlap with repeat and gene annotations, but I want to separate the different
# repeat types first 
awk '($4=="rRNA" || $4=="Simple_repeat" || $4=="unknown"){print $0}' annotation/EDTA2.v0.2.bed >annotation/Other.bed
awk '($4!="rRNA" && $4!="Simple_repeat" && $4!="unknown" && $4!="Satellite/Satellite"){print $0}' annotation/EDTA2.v0.2.bed >annotation/TEs.bed

# Make table that can be copied into the supplementary material
module load bedtools/2.31.0
rm -f RESULTS/g4discovery/common_g4s/common_g4s.annotation.table.txt
cat test.common.g4s.txt |head -n30 |while read line
do 
    seq=$(echo $line |cut -f1 -d" ")
    file=RESULTS/g4discovery/common_g4s/$seq.bed
    len=`awk '{sum+=($3-$2)}END{print sum}' $file`
    line=""
    for type in "TEs" "Satellites" "TRF_withMers" "Other" "UTR5" "promoter" "CDS" "intronic" "UTR3" "nonprotcoding" "intergenic"
    do
        #summary=`intersectBed  -a $file -b annotation/$type.bed -wao | awk '{len+=$3-$2; overlap+=$NF}END{print len"\t"overlap"\t"(overlap/len)}'`
        summary=`intersectBed  -a $file -b annotation/$type.bed -wao | awk -v len=$len '{overlap+=$NF}END{print (overlap/len)}'`
        line=$line" "$summary
    done
    echo  -e $seq" "$line >>RESULTS/g4discovery/common_g4s/common_g4s.annotation.table.txt
done

# Should rerun the following, including rDNA! 

# Check more specific Satellites and TEs
module load bedtools/2.31.0
cat test.common.g4s.txt |head -n30 |while read -r seq no
do 
    file=RESULTS/g4discovery/common_g4s/$seq.bed
    len=`awk '{sum+=($3-$2)}END{print sum}' $file`
    line=""
    intersectBed -a $file -b annotation/TEs.bed -wao | awk '{print $8,$11}' |sort -k1,1 |awk -v l=$len -v seq=$seq '{if(NR==1){rep=$1;replen=$2}else{if(rep==$1){replen+=$2}else{frac=replen/l; if(frac>0.01){print seq"\t"rep"\t"frac}; rep=$1; replen=$2}}}END{frac=replen/l; if(frac>0.01){print seq"\t"rep"\t"frac}}'
done >RESULTS/g4discovery/common_g4s/common_g4s.detailed_TE.table.txt

# Satellites
cat test.common.g4s.txt |head -n30 |while read -r seq no
do 
    file=RESULTS/g4discovery/common_g4s/$seq.bed
    len=`awk '{sum+=($3-$2)}END{print sum}' $file`
    line=""
    intersectBed -a $file -b annotation/Satellites.bed -wao | awk '{print $8,$9}' |sort -k1,1 |awk -v l=$len -v seq=$seq '{if(NR==1){rep=$1;replen=$2}else{if(rep==$1){replen+=$2}else{frac=replen/l; if(frac>0.01){print seq"\t"rep"\t"frac}; rep=$1; replen=$2}}}END{frac=replen/l; if(frac>0.01){print seq"\t"rep"\t"frac}}'
done >RESULTS/g4discovery/common_g4s/common_g4s.detailed_Satellite.table.txt


# Check if the most common G4s have been validated experimentally
cat test.common.g4s.txt |head -n100 |while read -r seq no
do 
    echo "Searching for "$seq
    grep -i $seq helpfiles/experimentally_validated.txt
    echo "#######################################################"
done 


############## CHECKING G4 SCORES IN DIFFERENT FUNCTIONAL REGIONS ##############

mkdir score_analysis
module load bedtools/2.31.0
prefix="bTaeGut7v0.4_MT_rDNA"
echo -e "Region\tChromosome_type\tPQS\tG4Hunter" >score_analysis/$prefix.G4.tsv
for region in "UTR5" "promoter" "CDS" "intronic" "UTR3" "nonprotcoding" "intergenic"
do 
    for type in "macro" "micro" "dot"
    do 
        echo "Processing region: "$region" on "$type" chromosomes"
        intersectBed -a <(grep $type helpfiles/$prefix.chr_types.txt |cut -f1 |grep -f - final_nonB/$prefix.G4.bed) -b annotation/$prefix.$region.bed | awk -v t=$type -v r=$region '{print r"\t"t"\t"$4"\t"$7}' >>score_analysis/$prefix.G4.tsv
    done 
done 


# NEED TO PLOT THIS AND SEE IF THE DISTRIBUTIONS DIFFER 



