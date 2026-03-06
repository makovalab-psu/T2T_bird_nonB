###############################################################################
### ANNOTATION OF NON-B DNA IN THE ZEBRA FINCH AND OTHER BIRD GENOMES.
### CODE WRITTEN BY LINNÉA SMEDS

# This file contains the code for annotating non-B DNA in the zebra finch and
# other bird genomes, as well as finding overlap between the different motifs.

##### TABLE OF CONTENTS
# ANNOTATION WITH GFA
# ANNOTATION OF Z DNA WITH ZDNA HUNTER
# G4 WITH G4DISCOVERY
# ANNOTATION OF G4s WITH QUADRON
# CREATE A JOINT DIRECTORY WITH ANNOTATIONS
# OVERLAP BETWEEN NON-B TYPES
# EXTRACT SPACER AND REPEAT ARM INFORMATION



############################## ANNOTATION WITH GFA #############################
# GFA is assumed to be installed in:  ~/software/non-B_gfa/gfa
mkdir -p slurm
mkdir -p gfa_annotation
mkdir -p final_nonB

# Run gfa and convert to bed
cat helpfiles/species_list.txt |while read -r sp longname prefix
do
  echo "start gfa for $sp"
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
cat helpfiles/species_list.txt |tail -n3 |while read -r sp longname prefix
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TESTING SOME OTHER SUBSETS, NOT USED IN PAPER 
awk -v OFS="\t" '($10<=15){s=$4-1; print $1,s,$5,$9":"$10":"$11,".",$8}' gfa_annotation/${prefix}_IR.tsv  >final_nonB/${prefix}.IR_SM15.bed
# MR, max 15bp loop, and 100% purine or pyrimidine
awk -v OFS="\t" '($10<=15){split($13, parts, "/"); 
                for (i=1; i<=4; i++) {
                  cnt = parts[i]
                  gsub(/[A-Z]/, "", cnt)  # remove A,C,G,T to get the numeric part
                  vals[i] = cnt + 0
                }
                A = vals[1]; C = vals[2]; G = vals[3]; T = vals[4];
                total = A + C + G + T;
                pur = A + G;
                frac=pur/total;
              if(frac==0 || frac==1){s=$4-1; print $1,s,$5,$9":"$10":"$11":"frac,".",$8}}' gfa_annotation/${prefix}_MR.tsv  >final_nonB/${prefix}.MR0-15.bed
# Max 8bp loop, and min 90% purine or pyrimidine
awk -v OFS="\t" '($10<=8){split($13, parts, "/"); 
                for (i=1; i<=4; i++) {
                  cnt = parts[i]
                  gsub(/[A-Z]/, "", cnt)  # remove A,C,G,T to get the numeric part
                  vals[i] = cnt + 0
                }
                A = vals[1]; C = vals[2]; G = vals[3]; T = vals[4];
                total = A + C + G + T;
                pur = A + G;
                frac=pur/total;
              if(frac<0.1 || frac>0.9){s=$4-1; print $1,s,$5,$9":"$10":"$11":"frac,".",$8}}' gfa_annotation/${prefix}_MR.tsv  >final_nonB/${prefix}.MR01-8.bed
# This is what I think gfa meant when they defined the subset "TRI", but it 
# seems like they have set the purine/pyrimidine content till MORE than 10% 
# instead of LESS than 10%



##################### ANNOTATION OF Z DNA WITH ZDNA HUNTER  ####################

# ZDNA Hunter runs on a webserver located here: https://bioinformatics.ibp.cz/#/analyse/zdna
# We uploaded the full genomes and downloaded the results as bed graph files. 
# ZDNA Hunter outputs positions continusly for the whole genome, so we need to 
# convert the positions to proper chromosome coordinates. Since the zebra finch
# genome file is larger than 2Gb, I ran the primary and alternative haplotypes
# separately.

# Place output in the directory ZDNAHunter
module load python/3.11.2
cat helpfiles/species_list.txt |grep -v zebra_finch |while read -r sp longname prefix
do
  #python3 T2T_bird_nonB/python/remap_zdna.py -b ZDNAHunter/$longname.model1.bedgraph -f ref/$prefix.fa.fai -o final_nonB/$prefix.ZDNAm1.bed
  python3 T2T_bird_nonB/python/remap_zdna.py -b ZDNAHunter/$prefix.model2.bedgraph -f ref/$prefix.fa.fai -o final_nonB/$prefix.Z.bed
done 
# We will use model 2 so I name that just "Z"

# Zebra finch 
prefix="bTaeGut7v0.4_MT_rDNA"
for part in "ZF_mat_Z" "ZF_pat_autosomes"
do
  python3 T2T_bird_nonB/python/remap_zdna.py -b ZDNAHunter/full.$part.model1.bedgraph -f ref/$part.fa.fai -o ZDNAHunter/$part.ZDNA1.bed
  python3 T2T_bird_nonB/python/remap_zdna.py -b ZDNAHunter/full.$part.model2.bedgraph -f ref/$part.fa.fai -o ZDNAHunter/$part.ZDNA2.bed
done 
# Merge into one file and move to nonB directory:
cat ZDNAHunter/ZF_*.ZDNA1.bed |sort -k1,1 -k2,2n >final_nonB/$prefix.ZDNA1.bed
cat ZDNAHunter/ZF_*.ZDNA2.bed |sort -k1,1 -k2,2n >final_nonB/$prefix.ZDNA2.bed



############################ G4 WITH G4DISCOVERY ##############################
# Install R packages needed for g4Discovery
module load python/3.11.2
module load r/4.5.0
R 
install.packages("seqinr")
# when asked if you want to do create a custom location: press yes, choose mirror
install.packages("BiocManager")
BiocManager::install(c("pqsfinder", "rtracklayer", "Biostrings"))

git clone https://github.com/saswat-km/g4Discovery.PanSN.git
# NOTES: The script assumes the g4discovery.PanSN folder is in your current 
# directory. I manually modified theg4DiscoveryFuncs.py script to add the full 
# path to the R script run_pqsfinder.R. The script assumes zipped input fasta.

cat helpfiles/species_list.txt |while read -r sp longname prefix
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
cat helpfiles/species_list.txt |while read -r sp longname prefix
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
cat helpfiles/species_list.txt |while read -r sp longname prefix
do
  echo "Start g4Discovery for $sp"
  mkdir -p g4discovery/$prefix
  for file in $(ls ref/split_fasta/$prefix/chrW.fa.gz)
  do
    # If file is bigger than 30Mb, use 4 cores
    if [ $(awk -v fa=$file 'BEGIN{cmd="stat -c%s "fa; cmd | getline size; close(cmd); print size}') -gt 30000000 ]; then
        mem="6G"
        cpus=2
    else
        mem="4G" #"4G"
        cpus=1  #"1"
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
# last three species needed just above 8Gb RAM and 1 core (I had given it 2)
# 1 core/4Gb was not enough for pigeon chrW. rerun with 2 cores/6Gb

# Combine results into one bed file per species 
cat helpfiles/species_list.txt |while read -r sp longname prefix
do
  echo "Combine g4Discovery results for $sp"
  rm -f final_nonB/$prefix.G4.bed
  for chr in $(cat ref/$prefix.fa.fai |cut -f1)
  do 
      zcat g4discovery/$prefix/$chr.bed.gz >>final_nonB/$prefix.G4.bed
  done 
done


######################## ANNOTATION OF G4s WITH QUADRON ########################
# Annotating G-quadruplex (G4) motifs with Quadron
singularity build quadron.sif docker://kxk302/quadron:1.0.0
mv quadron.sif software/quadron.sif

# Use the split chromosomes from before, but assumes unzipped files 
#fasta="chicken.v23.fa"

# Run Quadron on full chromosomes (separately)
# Note that Quadron needs a full path to the input files
cat helpfiles/species_list.txt |while read -r sp longname prefix
do

  mkdir Quadron_annotation/$prefix
  #for fa in $(ls [insert_full_path_here]/ref/split_fasta/$prefix/chr*.fa)
  for fa in $(ls /storage/group/kdm16/default/lbs5874/ZebraFinch/ref/split_fasta/$prefix/chr*.fa)
  do
    ls $fa
    chr=`basename $fa |sed 's/.fa//'`
    echo $chr
    out="/storage/group/kdm16/default/lbs5874/ZebraFinch/Quadron_annotation/'$prefix'/$chr.out"
    echo '#!/bin/bash
    singularity exec ~/software/quadron.sif sh -c "cd /scripts/; Rscript run_quadron.R '$fa' '$out' 4 1000000"
    ' |sbatch -J $chr --ntasks=1 --cpus-per-task=2 --mem-per-cpu=4G --time=1:00:00 -o slurm/job.quadron.$chr.%j.out
  done
done 
# Chr1-2 needed 4 cores, the others worked with 2 cores 

# If Quadron get's too little memory, it can sometimes stop calculating the 
# score, but still output the motif. Check how many "NA" scores there are 
# per file (a couple close to the ends are ok - Quadron can't calculate a 
# score if the flanking sequence is shorter than 50bp)
cat helpfiles/species_list.txt |while read -r sp longname prefix
do
  for file in $(ls Quadron_annotation/$prefix/chr*.out)
  do
    na=`grep "^DATA:" $file |awk '($5=="NA"){n++}END{print n}'`
    echo $file $na
  done
done 
# Looks ok, most have 2 NA G4s, some have 1 or 3

# Convert to bed
cat helpfiles/species_list.txt |grep -v zebra_finch |while read -r sp longname prefix
do
  rm -f Quadron_annotation/$prefix.G4.bed
  echo $prefix
  for chr in $(cut -f1 ref/$prefix.fa.fai)
  do
    awk -v chr=$chr -v OFS="\t" '(/^DATA/){if($5=="NA"){}else{score=$5; if(score<0){score=0}; s=$2-1; e=s+$4; print chr,s,e,$5,sprintf("%.f", score),$3}}' Quadron_annotation/$prefix/$chr.out >>final_nonB/$prefix.G4quadron.bed
  done
done 


################## CREATE A JOINT DIRECTORY WITH ANNOTATIONS ###################
#Make sure all non-B bed files are in the final folder, and make merged versions 
cat helpfiles/species_list.txt |while read -r sp longname prefix
do
  for n in "G4" "APR" "DR" "IR" "TRI" "STR" "Z" #"IRgfa" "MRgfa" "Zgfa" "G4quadron" "ZDNAm1" 
  do
  echo '#!/bin/bash
  module load bedtools/2.31.0
  sort -k1,1 -k2,2n final_nonB/'$prefix'.'${n}'.bed |mergeBed -i - >final_nonB/'$prefix'.'${n}'.merged.bed
  '| sbatch -J $n --ntasks=1 --cpus-per-task=1 --mem-per-cpu=6G --out slurm/mergeBed.$n.%j.out
  done
done

# And a joint bed file called "Any" (Just using the chosen settings for each type)
module load bedtools/2.31.0
cat helpfiles/species_list.txt |grep Pat |while read -r sp longname prefix
do
  echo "Make Any bed for $sp"
  sort -m final_nonB/$prefix.APR.merged.bed final_nonB/$prefix.DR.merged.bed final_nonB/$prefix.IR.merged.bed final_nonB/$prefix.TRI.merged.bed final_nonB/$prefix.STR.merged.bed final_nonB/$prefix.G4.merged.bed final_nonB/$prefix.Z.merged.bed >final_nonB/$prefix.Any.bed
  sort -k1,1 -k2,2n final_nonB/$prefix.Any.bed |mergeBed -i - >final_nonB/$prefix.Any.merged.bed
done 



######################### OVERLAP BETWEEN NON-B TYPES ##########################
# Check overlap between all non-B types to be visulized with an upset plot

mkdir overlap

# Check overlap for each chromosome separately
cat helpfiles/species_list.txt |while read -r sp longname prefix
do
  echo '#!/bin/bash
  cat ref/'$prefix'.fa.fai |grep "chr" |cut -f1 |while read -r chr;
  do
    echo "Running upset for $chr"
    module load python/3.11.2
    python3 T2T_bird_nonB/python/upset_summary.py -b final_nonB/'$prefix'. -c $chr -o overlap/'$prefix'.summary.7types.$chr.txt
  done
   ' | sbatch -J $prefix --ntasks=1 --mem-per-cpu=8G --cpus-per-task=1 --time=1:00:00 --out slurm/job.overlap.$prefix.%j.out
done 
# The longest chr2 took 3Gb of RAM and less than a minute to run.
# When doing all 17 types, it used 5Gb and tool 1.5min

cat helpfiles/species_list.txt |grep zebra |while read -r sp longname prefix
do
  # Merge autosomes
  cat overlap/$prefix.summary.7types.chr[0-9]*.txt |sort | awk -v OFS="\t" '{if(NR==1){type=$1; sum=$2}else{if($1==type){sum+=$2}else{print type,sum; type=$1; sum=$2}}}END{print type,sum}' > overlap/$prefix.summary.7types.autosomes.txt
  # And chromosome types 
  for type in "macro" "micro" "dot"
  do
    for chr in $(grep $type helpfiles/$prefix.groups.txt |cut -f1)
    do
      cat overlap/$prefix.summary.7types.$chr.txt
    done | sort | awk -v OFS="\t" '{if(NR==1){type=$1; sum=$2}else{if($1==type){sum+=$2}else{print type,sum; type=$1; sum=$2}}}END{print type,sum}' > overlap/$prefix.summary.7types.$type.txt
  done
done

# PAIRWISE OVERLAP WITH FRACTION
# To get the fraction I first need the totals for all types
# This does not work when some of the names are overlapping 
mkdir stats/nonB_totals
cat helpfiles/species_list.txt |while read -r sp longname prefix
do
  for file in $(ls overlap/$prefix.summary.7types.*.txt)
  do
    chr="${file##*.7types.}"
    chr="${chr%.txt}"
    echo $chr
    echo "#NonB bp" |sed "s/ /\t/g" >stats/nonB_totals/$prefix.annotated.7types.${chr}.txt
    for nb in "APR" "DR" "G4" "IR" "TRI" "STR" "Z" 
    do
      grep $nb $file |awk -v nb=$nb '{sum+=$2}END{print nb"\t"sum}' >>stats/nonB_totals/$prefix.annotated.7types.${chr}.txt
    done
  done

  # Make bed files for each chr type (this is just to simplify the next step)
  module load bedtools/2.31.0
  mkdir -p final_nonB/merged_per_chrom/
  cat ref/${prefix}*.fa.fai |cut -f1 |grep chr |while read -r chr;
  do
    echo $chr
    for nb in  "APR" "DR" "G4" "IR" "TRI" "STR" "Z" 
    do
      awk -v c="$chr" '$1 == c' final_nonB/$prefix.$nb.merged.bed \
      > final_nonB/merged_per_chrom/$prefix.$chr.$nb.bed
    done
  done
  # Autosomes 
  for nb in "APR" "DR" "G4" "IR" "TRI" "STR" "Z" 
  do
    cat final_nonB/merged_per_chrom/$prefix.chr[0-9]*.$nb.bed > final_nonB/merged_per_chrom/$prefix.autosomes.$nb.bed
  done
  # Macro, micro, dot
  for type in "macro" "micro" "dot"
  do
    for nb in "APR" "DR" "G4" "IR" "TRI" "STR" "Z" 
    do
      for chr in $(grep $type helpfiles/$prefix.groups.txt |cut -f1)
      do
        cat final_nonB/merged_per_chrom/$prefix.$chr.$nb.bed
      done >final_nonB/merged_per_chrom/$prefix.$type.$nb.bed
    done
  done

  # Then we can combine all the types
  for file in $(ls stats/nonB_totals/$prefix.annotated.7types.*.txt )
  do
    chr="${file##*.7types.}"
    chr="${chr%.txt}"
    echo $chr
    echo "Chr NB1 NB2 Ovl Frac" | sed "s/ /\t/g" >overlap/$prefix.pairwise.7types.$chr.txt
    echo '#!/bin/bash
    module load bedtools/2.31.0
  # Go through all combinations
    for nb in "APR" "DR" "G4" "IR" "TRI" "STR" "Z"
    do
      nbbp=`grep $nb '$file' |cut -f2`
      echo "$nb has $nbbp number of bp!"
      for nb2 in "APR" "DR" "G4" "IR" "TRI" "STR" "Z"
      do
          if [[ $nb != $nb2 ]]
          then
            intersectBed -a final_nonB/merged_per_chrom/'$prefix'.'$chr'.$nb.bed -b final_nonB/merged_per_chrom/'$prefix'.'$chr'.$nb2.bed -wo | awk -v tot=$nbbp -v chr='$chr' -v sum=0 -v nb1=$nb -v nb2=$nb2 -v OFS="\t" '"'"'{sum+=$7}END{f=sum/tot; print chr,nb1,nb2,sum,f}'"'"' >>overlap/'$prefix'.pairwise.7types.'$chr'.txt
          fi
      done
    done
    ' | sbatch -J $chr --ntasks=1 --cpus-per-task=1 --mem-per-cpu=4G --time=1:00:00 --out slurm/job.pairwise-overlap.$prefix.7types.$chr.%j.%N.out
  done
done 

#When all jobs are done, merge  
cat helpfiles/species_list.txt |grep zebra |while read -r sp longname prefix
do
  echo "Region NonB Overlap" |sed 's/ /\t/g' >overlap/$prefix.merged.7types.summary.txt
  echo "Region NB1 NB2 Ovl Frac" |sed 's/ /\t/g' >overlap/$prefix.merged.7types.pairwise.txt
  for type in "macro" "micro" "dot" "autosomes" "chrW*" "chrZ*"
  do
    awk -v t=$type '{print t"\t"$0}' overlap/$prefix.summary.7types.$type.txt |sed 's/*//g' >>overlap/$prefix.merged.7types.summary.txt
    awk -v t=$type '(NR>1){print $0}' overlap/$prefix.pairwise.7types.$type.txt |sed 's/*//g' >>overlap/$prefix.merged.7types.pairwise.txt
  done
done 


# Some stats
cat helpfiles/species_list.txt |grep zebra |while read -r sp longname prefix
do
  echo "+++++++++++++++++++++++++++++++++++++++++++"  
  echo $sp
  # Number of non-B bases
  grep -E 'macro|micro|dot' overlap/$prefix.merged.7types.summary.txt |awk '{sum+=$3}END{print "Total non-B bases:", sum}'
  # Number of bases in overlap
  grep -E 'macro|micro|dot' overlap/$prefix.merged.7types.summary.txt |grep "-" |awk '{sum+=$3}END{print "Bases in overlap:", sum}'
  # Number of bp in at least three types
  grep -E 'macro|micro|dot' overlap/$prefix.merged.7types.summary.txt |grep "-" |awk -F'\t' 'gsub(/-/, "&", $2) >= 2' |awk '{sum+=$3}END{print "Bases in at least three types:", sum}'
done
+++++++++++++++++++++++++++++++++++++++++++
zebra_finch
Total non-B bases: 163492668
Bases in overlap: 24564720
Bases in at least three types: 8075289
+++++++++++++++++++++++++++++++++++++++++++
chicken
Total non-B bases: 83955375
Bases in overlap: 15495711
Bases in at least three types: 4627580
+++++++++++++++++++++++++++++++++++++++++++
annas_hummingbird
Total non-B bases: 62988302
Bases in overlap: 7290343
Bases in at least three types: 2542044
+++++++++++++++++++++++++++++++++++++++++++
great_bustard
Total non-B bases: 70403513
Bases in overlap: 6881955
Bases in at least three types: 1563257
+++++++++++++++++++++++++++++++++++++++++++
ural_owl
Total non-B bases: 80767149
Bases in overlap: 8275996
Bases in at least three types: 2505854
+++++++++++++++++++++++++++++++++++++++++++
bandtailed_pigeon
Total non-B bases: 106412589
Bases in overlap: 13666750
Bases in at least three types: 3708370
+++++++++++++++++++++++++++++++++++++++++++
emu
Total non-B bases: 84434601
Bases in overlap: 9370744
Bases in at least three types: 2766569
+++++++++++++++++++++++++++++++++++++++++++
peking_duck
Total non-B bases: 102946645
Bases in overlap: 17117681
Bases in at least three types: 6259198



# Number of bases with only one vs more than one non-B annotation
cat helpfiles/species_list.txt |while read -r sp longname prefix
do
  echo "+++++++++++++++++++++++++++++++++++++++++++"  
  echo $sp
  for set in "autosomes" "chrZ*" "chrW*" "macro" "micro" "dot"
  do
    u=`grep -v "-" overlap/$prefix.summary.7types.$set.txt |awk '{sum+=$2}END{print sum}'`
    o=`grep "-" overlap/$prefix.summary.7types.$set.txt |awk '{sum+=$2}END{print sum}'`
    echo $set $u $o
  done 
done 
+++++++++++++++++++++++++++++++++++++++++++
zebra_finch
autosomes 132298158 23506626
chrZ* 4531722 828026
chrW* 2098068 230068
macro 77394539 10053687
micro 35998971 6174009
dot 25534438 8337024
+++++++++++++++++++++++++++++++++++++++++++
chicken
autosomes 62109067 13213590
chrZ* 4463140 920508
chrW* 1887457 1361613
macro 40689642 6891024
micro 15716199 3435557
dot 12053823 5169130
+++++++++++++++++++++++++++++++++++++++++++
annas_hummingbird
autosomes 51057188 6607208
chrZ* 3597306 553129
chrW* 1043465 130006
macro 38107467 4959739
micro 17148600 2239798
dot 441892 90806
+++++++++++++++++++++++++++++++++++++++++++
great_bustard
autosomes 59328371 6335453
chrZ* 4193187 546502
grep: overlap/OTswu.summary.7types.chrW*.txt: No such file or directory
grep: overlap/OTswu.summary.7types.chrW*.txt: No such file or directory
chrW*
macro 43296427 3975233
micro 15870196 1453287
dot 4354935 1453435
+++++++++++++++++++++++++++++++++++++++++++
ural_owl
autosomes 67489565 7674845
chrZ* 5001588 601151
grep: overlap/bStrUra1.summary.7types.chrW*.txt: No such file or directory
grep: overlap/bStrUra1.summary.7types.chrW*.txt: No such file or directory
chrW*
macro 46027356 5188501
micro 22957664 2569957
dot 3506133 517538
+++++++++++++++++++++++++++++++++++++++++++
bandtailed_pigeon
autosomes 82755918 11068106
chrZ* 4286621 849059
chrW* 5703300 1749585
macro 54390356 8310794
micro 25528046 2246837
dot 12827437 3109119
+++++++++++++++++++++++++++++++++++++++++++
emu
autosomes 67103938 8477643
chrZ* 4176815 426693
chrW* 3783104 466408
macro 50285756 5626971
micro 20679312 2212469
dot 4098789 1531304
+++++++++++++++++++++++++++++++++++++++++++
peking_duck
autosomes 79953647 15230945
chrZ* 4924861 1676674
chrW* 950456 210062
macro 49998800 10812492
micro 24488558 4145728
dot 11341606 2159461






################## EXTRACT SPACER AND REPEAT ARM INFORMATION ###################
# For plotting to investigate the properties of DR, MR and IR motifs

for nb in "DR" "IR" "MR"
do
  cut -f9,10 gfa_annotation/bTaeGut7v0.4_MT_rDNA_$nb.tsv >stats/bTaeGut7v0.4_MT_rDNA_$nb.rep-space.txt
done



