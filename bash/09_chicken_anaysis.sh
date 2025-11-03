################################################################################
### CODE FOR REPEATING PREVIOUS ANALYSES ON THE CHICKEN GENOME (T2T EXCEPT FOR
### THE W-CHROMOSOME)
### WRITTEN BY LINNÉA SMEDS

# TABLE OF CONTENTS 
# DOWNLOAD AND PREPARE GENOME 
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
# Download genome and annotation from dropbox
# https://www.dropbox.com/scl/fo/plq2tm2w9lzlk0ua1rzph/h?rlkey=l6z3rgmjs7ec9azun8nundnzl&e=1&dl=0
# place files in ref/
# Index genome 
samtools faidx REF/chicken.v23.fa

############################## ANNOTATION WITH GFA #############################
# GFA is assumed to be installed in: 
# ~/software/non-B_gfa/gfa
# Set parameters
fasta="chicken.v23.fa"
prefix="chicken.v23"
# Run gfa and convert to bed
echo '#!/bin/bash
echo "##### Start gfa"
~/software/non-B_gfa/gfa -seq '$fasta' -out '$prefix'  -skipGQ -skipWGET
echo "##### translate to bed"
for type in "APR" "DR" "IR" "MR" "STR" "Z"
do
  echo "#####....translating $type...."
  awk -v OFS="\t" '"'"'(NR>1){s=$4-1; print $1,s,$5,$9":"$10":"$11,".",$8}'"'"' '${prefix}'_${type}.tsv >'${prefix}'.$type.bed
done
awk -v OFS="\t" '"'"'($12==1){s=$4-1; print $1,s,$5,$9":"$10":"$11,".",$8}'"'"' '${prefix}'_MR.tsv  >'${prefix}'.TRI.bed
' |sbatch -J gfa.$prefix  --ntasks=1 --cpus-per-task=1 --mem-per-cpu=2G --time=2-00:00:00 -o slurm/job.gfa.$prefix.%j.%N.out


########################### ANNOTATION WITH QUADRON ############################
# Assumes the quadron singularity is here: ~/software/quadron.sif 
# Divide fasta into chromosomes
fasta="chicken.v23.fa"
prefix="chicken.v23"
mkdir Quadron_annotation/split_fasta/$prefix
echo '#!/bin/bash
faidx --split-files ref/'$fasta'
mv *.fa Quadron_annotation/split_fasta/'$prefix'/
' |sbatch -J split --ntasks=1 --cpus-per-task=4 --time=5:00:00 -o slurm/job.splitfasta.$prefix.%j.out

# Run Quadron on full chromosomes (separately)
# Note that Quadron needs a full path to the input files
mkdir Quadron_annotation/$prefix
for fa in $(ls [insert_full_path_here]/Quadron_annotation/split_fasta/$prefix/chr*.fa)
do
  ls $fa
  chr=`basename $fa |sed 's/.fa//'`
  echo $chr
  out="Quadron_annotation/'$prefix'/$chr.out"
  echo '#!/bin/bash
  singularity exec ~/software/quadron.sif sh -c "cd /scripts/; Rscript run_quadron.R '$fa' '$out' 4 1000000"
  ' |sbatch -J $chr --ntasks=1 --cpus-per-task=4 --time=1:00:00 -o slurm/job.quadron.$chr.%j.out
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
  awk -v chr=$chr -v OFS="\t" '(/^DATA/){if($5=="NA"){}else{score=$5; if(score<0){score=0}; s=$2-1; e=s+$4; print chr,s,e,$5,sprintf("%.f", score),$3}}' Quadron_annotation/$prefix/$chr.out >>Quadron_annotation/$prefix.G4.bed
done


############################# CALCULATE DENSITY ################################
prefix="chicken.v23"
# Place all non-B bed files in the final folder
cp gfa_annotation/$prefix.*.bed final_nonB/
cp Quadron_annotation/$prefix.G4.bed final_nonB/

# Create genomic windows and make length file
module load bedtools/2.31.0
bedtools makewindows -g ref/$prefix.fa.fai -w 100000 >windows/$prefix.all_100k_windows.bed;

# Note - the following will exclude chrMT because it's shorter than 100 kb!
for n in "G4" "APR" "DR" "IR" "MR" "TRI" "STR" "Z"
do
 echo '#!/bin/bash
 module load bedtools/2.31.0
 intersectBed -wao -a windows/'$prefix'.all_100k_windows.bed -b <(cut -f1,2,3 final_nonB/'$prefix'.'${n}'.bed |sort -k1,1 -k2,2n |mergeBed -i -) | cut -f1,2,3,7 |awk -v OFS="\t" '"'"'{if(NR==0){chr=$1; s=$2; e=$3; sum=$4}else{if($1==chr && $2==s){sum+=$4}else{print chr,s,e,sum; chr=$1; s=$2; e=$3; sum=$4}}}END{print chr,s,e,sum}'"'"' | sed "/^\s*$/d" >densities/'${prefix}'.'${n}'.100kb.bed
 '| sbatch -J $n --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8G --out slurm/density.$n.%j.out
done

# For plotting, merge the data:
echo "NonB Chr Start Stop Dens" |sed 's/ /\t/g' >densities/$prefix.merged.100kb.txt
for n in "APR" "DR" "G4" "IR" "MR" "TRI" "STR" "Z"
do
  awk -v OFS="\t" -v nb=$n '{print nb,$0}' densities/${prefix}.${n}.100kb.bed >>densities/$prefix.merged.100kb.txt
done

# CREATE NEW NONB FILES WITH MERGED REGIONS (NO OVERLAPS!)
for n in "G4" "APR" "DR" "IR" "MR" "TRI" "STR" "Z" "ALL"
do
 echo '#!/bin/bash
 module load bedtools/2.31.0
 sort -k1,1 -k2,2n final_nonB/'$prefix'.'${n}'.bed |mergeBed -i - >final_nonB/'$prefix'.'${n}'.merged.bed
 '| sbatch -J $n --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8G --out slurm/mergeBed.$n.%j.out
 done

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PLOT DENSITIES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
for nb in "TRI" "APR" "DR" "G4" "IR" "MR" "STR" "Z"
do
  grep -v "mito" densities/$prefix.${nb}.100kb.bed |grep -v "chrW" |sed '/^[[:space:]]*$/d' |awk '{gsub(/chr/, "nn"); $3=$3-1; print $1,$2,$3,$4}' >$dir/data/${nb}.100kb.txt
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
echo "Chr Length Group" |sed 's/ /\t/g' >chicken.v23.groups.txt
cat ref/$prefix.fa.fai | awk -v OFS="\t" '{print $1,$2}' >>chicken.v23.groups.txt


###################### PREPARE ENRICHMENT CALCULATIONS #########################
# Calculate genome-wide density to use for normalization
prefix="chicken.v23"
echo '#!/bin/bash
module load bedtools/2.31.0
rm -f densities/'${prefix}'.nonB_genome_wide.txt
totlen=`cat ref/'${prefix}'*.fa.fai | awk '"'"'{sum+=$2}END{print sum}'"'"'`
echo "Totlen: "$totlen
for non_b in "APR" "DR" "STR" "IR" "MR" "TRI" "G4" "Z"
do
  sort -k1,1 -k2,2n final_nonB/'${prefix}'.${non_b}.bed |mergeBed -i - | awk -v n=$non_b -v tot=$totlen '"'"'{sum+=$3-$2}END{d=sum/tot; print n, sum, d}'"'"' >>densities/'${prefix}'.nonB_genome_wide.txt
done
' |sbatch -J enrichment --ntasks=1 --cpus-per-task=2 --mem-per-cpu=8G --out slurm/job.enrichment.$prefix.%j.out

# And all non-B combined
prefix="chicken.v23"
echo '#!/bin/bash
module load bedtools/2.31.0
totlen=`cat ref/'${prefix}'*.fa.fai | awk '"'"'{sum+=$2}END{print sum}'"'"'`
cat final_nonB/'${prefix}'.*.bed |sort -k1,1 -k2,2n |mergeBed -i - >final_nonB/'${prefix}'.ALL.bed
awk -v n="ALL" -v tot=$totlen '"'"'{sum+=$3-$2}END{d=sum/tot; print n, sum, d}'"'"' final_nonB/'${prefix}'.ALL.bed  >>densities/'${prefix}'.nonB_genome_wide.txt
' |sbatch -J enrichment --ntasks=1 --cpus-per-task=2 --mem-per-cpu=8G --out slurm/job.enrichment.$prefix.%j.out

# Per chromosome (not only enrichment, but per chromosome density)
echo '#!/bin/bash
module load bedtools/2.31.0
echo "Chr NonB Bp Density" |sed "s/ /\t/g" >densities/'${prefix}'.nonB_per_chrom.txt
cat ref/'${prefix}'*.fa.fai |cut -f1,2 | while read -r chr len;
do
    for non_b in "APR" "DR" "G4" "IR" "MR" "TRI" "STR" "Z" "ALL"
    do
        cat final_nonB/'${prefix}'.${non_b}.bed |awk -v c=$chr '"'"'($1==c){print}'"'"' |sort -k1,1 -k2,2n |mergeBed -i - | awk -v n=$non_b -v tot=$len '"'"'{sum+=$3-$2; c=$1}END{d=sum/tot; print c, n, sum, d}'"'"' >>densities/'${prefix}'.nonB_per_chrom.txt
    done
done
'| sbatch -J per_chrom --ntasks=1 --cpus-per-task=1 --mem-per-cpu=4G --time=1:00:00 --partition=open


# ~~~~~~~~~~~~~~~~~~~~~~~~~ PRINT DENSITIES FOR TABLE ~~~~~~~~~~~~~~~~~~~~~~~~~~
# Write it out with the nonB types as columns instead
# Whole genome
prefix="chicken.v23"
tmp="genome"
for non_b in "APR" "DR" "G4" "IR" "MR" "TRI" "STR" "Z" "ALL"
do
  n=`grep $non_b densities/${prefix}.nonB_genome_wide.txt |cut -f3 -d" "`
  tmp=$tmp" "$n
done
echo "Chromosome APR DR G4 IR MR TRI STR Z ALL"
echo $tmp

# Per chromosome
cat ref/${prefix}*.fa.fai |cut -f1 |while read -r chr;
do
    tmp="$chr"
    for non_b in "APR" "DR" "G4" "IR" "MR" "TRI" "STR" "Z" "ALL"
    do
      n=`awk -v c=$chr '($1==c){print}' densities/${prefix}.nonB_per_chrom.txt |grep " $non_b " |cut -f4 -d" "`
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
  rm -f enrichment/'$prefix'.'$chr'.tmp-centromere.txt
  cat densities/'${prefix}'.nonB_genome_wide.txt |grep -v "ALL" | while read -r nb tot dens;
  do
    chrdens=`grep "'$chr' " densities/'${prefix}'.nonB_per_chrom.txt |grep " $nb " |cut -f4 -d" "`
    echo "Density for '$chr' and $nb is $chrdens"
    stats=`echo '$chr' '$start' '$end' |sed "s/ /\t/g" |intersectBed -a - -b final_nonB/'${prefix}'.$nb.merged.bed | sort -k1,1 -k2,2n | mergeBed -i - |\
    awk -v l='$cenlen' -v dtot=$dens -v dchr=$chrdens '"'"'{sum+=$3-$2}END{d=sum/l; fracgw=d/dtot; fracchr=d/dchr; print d,fracgw,fracchr}'"'"'`
    echo '$chr'" "$nb" "$stats |sed "s/ /\t/g" >>enrichment/'$prefix'.'$chr'.tmp-centromere.txt
  done
  ' | sbatch -J $chr --ntasks=1 --cpus-per-task=1 --mem-per-cpu=4G --out slurm/job.centromeric_enrichment.$chr.%j.out
done
#Merge
echo "Chr nonB Density Enrichment_GW Enrichment_Chr" |sed "s/ /\t/g" >enrichment/${prefix}.centromere.txt
cat $cen_bed_file |cut -f1 |while read -r chr;
do
  cat enrichment/$prefix.$chr.tmp-centromere.txt >> enrichment/${prefix}.centromere.txt
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
    poswin=`echo $chrlen"/"$c_len |bc`
    minsize=100000
    if [ "$poswin" -lt "120" ]
    then
      echo "For '$chr', there will be less than 120 nonovl windows"
      slide=`echo "("$chrlen"-"$c_len")/120"|bc`
      bedtools makewindows -g <(awk -v c='$chr' '"'"'(c==$1){print}'"'"' ref/'$prefix'.fa.fai) -s $slide -w $c_len |subtractBed -a - -b '$cen_file' |awk -v min=$minsize -v clen=$c_len '"'"'($3-$2>min || $3-$2==clen){print}'"'"' |shuf -n 100 |sort -k2,2n >centromeres/windows/'$prefix'.'$chr'.exclCen.bed
    else
      echo "For '$chr', there will be more than 120 nonovl windows"
      bedtools makewindows -g <(awk -v c='$chr' '"'"'(c==$1){print}'"'"' ref/'$prefix'.fa.fai) -w $c_len |subtractBed -a - -b '$cen_file' |awk -v min=$minsize -v clen=$c_len '"'"'($3-$2>min || $3-$2==clen){print}'"'"' |shuf -n 100 |sort -k2,2n >centromeres/windows/'$prefix'.'$chr'.exclCen.bed
    fi
    bedtools nuc -fi ref/'$prefix'.fa -bed '$cen_file' >centromeres/windows/'$prefix'.'$chr'.exclCen.with_GC.bed
    rm -f centromeres/windows/'$prefix'.'$chr'.100rand.enrichment.tsv
    tmp="'$chr'"
    while read line
    do
      cat densities/'${prefix}'.nonB_genome_wide.txt |grep -v "ALL" | while read -r non_b tot dens;
      do
          echo $line |sed "s/ /\t/g" >tmp.$SLURM_JOB_ID.bed
          d=`intersectBed -a tmp.$SLURM_JOB_ID.bed -b final_nonB/'$prefix'.${non_b}.merged.bed -wo |awk -v l=$c_len -v dtot=$dens '"'"'{sum+=$7}END{d=sum/l; frac=d/dtot; print frac}'"'"'`
          tmp=`echo $tmp" "$d`
          echo $tmp >tmp.$SLURM_JOB_ID
      done
      cat tmp.$SLURM_JOB_ID |sed "s/ /\t/g" >>centromeres/windows/'$prefix'.'$chr'.100rand.enrichment.tsv
    done <centromeres/windows/'$prefix'.'$chr'.exclCen.bed
    '| sbatch -J $schr --ntasks=1 --cpus-per-task=1 --time=5:00:00 --out slurm/job.centromeric_enrichment.$chr.%j.out
done
# Check if there are more or less than 120 windows
cat slurm/job.centromeric_enrichment.chr*.*.out |grep "there will be"

# Merge background densities and GC
echo "Window Chr APR DR G4 IR MR TRI STR Z" |sed 's/ /\t/g'  >centromeres/$prefix.background_enrichment_merged.tsv
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

