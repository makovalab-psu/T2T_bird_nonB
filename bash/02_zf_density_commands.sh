################################################################################
### CALCULATING NON-B DNA MOTIF DENSITY IN THE ZEBRA FINCH T2T GENOME
### CODE WRITTEN BY LINNÉA SMEDS

##### TABLE OF CONTENTS
# CALCULATE DENTISTY
#   -PLOTTING DENSITIES USING CIRCOS
# PREPARE ENRICHMENT CALCULATIONS
#   -PRINT DENSITIES FOR TABLE
#   -GC CONTENT PER CHROMOSOME


############################# CALCULATE DENSITY ################################
# June 24, 2025
prefix="bTaeGut7v0.4_MT_rDNA"

# Create genomic windows and make length file, trying both 100kb and 10kb windows
mkdir -p windows/
bedtools makewindows -g ref/$prefix.fa.fai -w 100000 >windows/all_100kb_windows.bed;
bedtools makewindows -g ref/$prefix.fa.fai -w 10000 >windows/all_10kb_windows.bed;

mkdir densities
for wind in "10kb" "100kb"
do
  for n in "G4" "APR" "DR" "IR" "MR" "TRI" "STR" "Z"
  do
   echo '#!/bin/bash
   module load bedtools/2.31.0
   intersectBed -wao -a windows/all_'$wind'_windows.bed -b <(cut -f1,2,3 final_nonB/'$prefix'.'${n}'.merged.bed) | cut -f1,2,3,7 |awk -v OFS="\t" '"'"'{if(NR==0){chr=$1; s=$2; e=$3; sum=$4}else{if($1==chr && $2==s){sum+=$4}else{print chr,s,e,sum; chr=$1; s=$2; e=$3; sum=$4}}}END{print chr,s,e,sum}'"'"' | sed "/^\s*$/d" >densities/'${prefix}'.'${n}'.'$wind'.bed
   '| sbatch -J $n --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8G --out slurm/density.$n.%j.out
  done
done

# For plotting, merge the data:
for wind in "10kb" "100kb"
do
  echo "NonB Chr Start Stop Dens" |sed 's/ /\t/g' >densities/$prefix.merged.$wind.txt
  for n in "APR" "DR" "G4" "IR" "MR" "TRI" "STR" "Z"
  do
    awk -v OFS="\t" -v nb=$n '{print nb,$0}' densities/${prefix}.${n}.$wind.bed >>densities/$prefix.merged.$wind.txt
  done
done


# ~~~~~~~~~~~~~~~~~~~~~~~ PLOT DENSITIES USING CIRCOS ~~~~~~~~~~~~~~~~~~~~~~~~~~
# June 26 2025

# ASsume circos is installed in ~/software/circos-0.69-9/
# The circos configuration files are found in circos/

# Create density and centromere band files separately for MAT and PAT
prefix="bTaeGut7v0.4_MT_rDNA"
for hap in "mat" "pat"
do
    mkdir -p circos/$hap
    for nb in "TRI" "APR" "DR" "G4" "IR" "MR" "STR" "Z"
    do
        grep "_$hap" densities/$prefix.${nb}.noovl.100kb.bed |sed '/^[[:space:]]*$/d' |awk '{gsub(/chr/, "nn"); split($1,s,"_"); $3=$3-1; print s[1],$2,$3,$4}' >circos/$hap/data/${nb}.100kb.txt
        max=`cut -f4 -d" " circos/$hap/data/${nb}.100kb.txt |sort -nr |head -n1`
    done
    grep "_$hap" ref/bTaeGut7v0.4_MT_rDNA.centromere_detector.v0.1.gff |awk '{gsub(/chr/,""); split($1,s,"_"); start=$4-1; print "nn"s[1],start,$5}' >circos/$hap/data/highlight_CEN.txt
done

# Plot with
for hap in "mat" "pat"
do
  cd circos/${hap}/
  ~/software/circos-0.69-9/bin/circos -conf circos_halfInnerCircle.conf
  mv circos.png ${hap}_circos.png
  cd ../..
done
# This looks cool - a lot of enrichment on microchromosomes!



########################## ENRICHMENT CALCULATIONS #############################
# June 26 2025
# To calculate enrichment in a certain region, we need to have a base line to 
# compare the density to. We can use either genome-wide density, or the density
# in a certain chromosome group (macro, micro, microdot), or per chromosome. 

# Calculate genome-wide density to use for normalization
prefix="bTaeGut7v0.4_MT_rDNA"
echo '#!/bin/bash
module load bedtools/2.31.0
rm -f densities/'${prefix}'.nonB_genome_wide.txt
totlen=`cat ref/'${prefix}'*.fa.fai | awk '"'"'{sum+=$2}END{print sum}'"'"'`
echo "Totlen: "$totlen
for non_b in "APR" "DR" "STR" "IR" "MR" "TRI" "G4" "Z"
do
  cat final_nonB/'${prefix}'.${non_b}.merged.bed | awk -v n=$non_b -v tot=$totlen '"'"'{sum+=$3-$2}END{d=sum/tot; print n, sum, d}'"'"' >>densities/'${prefix}'.nonB_genome_wide.txt
done
' |sbatch -J enrichment --ntasks=1 --cpus-per-task=2 --mem-per-cpu=8G --out slurm/job.enrichment.$prefix.%j.out

# And all non-B combined
prefix="bTaeGut7v0.4_MT_rDNA"
echo '#!/bin/bash
module load bedtools/2.31.0
totlen=`cat ref/'${prefix}'*.fa.fai | awk '"'"'{sum+=$2}END{print sum}'"'"'`
cat final_nonB/'${prefix}'.*.bed |sort -k1,1 -k2,2n |mergeBed -i - >final_nonB/'${prefix}'.ALL.bed
awk -v n="ALL" -v tot=$totlen '"'"'{sum+=$3-$2}END{d=sum/tot; print n, sum, d}'"'"' final_nonB/'${prefix}'.ALL.bed  >>densities/'${prefix}'.nonB_genome_wide.txt
' |sbatch -J enrichment --ntasks=1 --cpus-per-task=2 --mem-per-cpu=8G --out slurm/job.enrichment.$prefix.%j.out

# Per chromosome density
echo '#!/bin/bash
module load bedtools/2.31.0
echo "Chr NonB Bp Density" |sed "s/ /\t/g" >densities/'${prefix}'.nonB_per_chrom.txt
cat ref/'${prefix}'*.fa.fai |cut -f1,2 |grep -v "chrM" |grep -v "rDNA" |while read -r chr len;
do
    for non_b in "APR" "DR" "G4" "IR" "MR" "TRI" "STR" "Z" "ALL"
    do
        cat final_nonB/'${prefix}'.${non_b}.bed |grep $chr |sort -k1,1 -k2,2n |mergeBed -i - | awk -v n=$non_b -v tot=$len '"'"'{sum+=$3-$2; c=$1}END{d=sum/tot; print c, n, sum, d}'"'"' >>densities/'${prefix}'.nonB_per_chrom.txt
    done
done
'| sbatch -J per_chrom --ntasks=1 --cpus-per-task=1 --mem-per-cpu=4G --time=1:00:00 --partition=open

# Per chromosome category density
for group in "macro" "micro" "microdot"
do
  totlen=`grep -f helpfiles/$group.txt ref/${prefix}*.fa.fai | awk '{sum+=$2}END{print sum}'`
  echo '#!/bin/bash
  module load bedtools/2.31.0
  rm -f tmp.density.'${group}'
  for non_b in "APR" "DR" "G4" "IR" "MR" "TRI" "STR" "Z" "ALL"
  do
      grep -f helpfiles/'$group'.txt final_nonB/'${prefix}'.${non_b}.merged.bed | awk -v g='$group '-v n=$non_b -v tot='$totlen' '"'"'{sum+=$3-$2}END{d=sum/tot; print g, n, sum, d}'"'"' >>tmp.density.'${group}'
  done
  '| sbatch -J per_chrom --ntasks=1 --cpus-per-task=1 --mem-per-cpu=4G --time=1:00:00 --partition=open
done

# Merge all categories together into one density file 
echo "Group NonB Bp Density" |sed "s/ /\t/g" >densities/${prefix}.nonB_per_group.txt
for group in "macro" "micro" "microdot"
do
  cat tmp.density.${group} | sed "s/ /\t/g" >>densities/${prefix}.nonB_per_group.txt
done

# ~~~~~~~~~~~~~~~~~~~~~~~~~ PRINT DENSITIES FOR TABLE ~~~~~~~~~~~~~~~~~~~~~~~~~~
# FOR SUPPLEMENTARY TABLE 1
# Whole genome
prefix="bTaeGut7v0.4_MT_rDNA"
tmp="genome"
for non_b in "APR" "DR" "G4" "IR" "MR" "TRI" "STR" "Z" "ALL"
do
  n=`grep $non_b densities/${prefix}.nonB_genome_wide.txt |cut -f3 -d" "`
  tmp=$tmp" "$n
done
echo "Chromosome APR DR G4 IR MR TRI STR Z ALL"
echo $tmp
# Per chromosome
cat ref/${prefix}*.fa.fai |cut -f1 |grep -v "chrM" |grep -v "rDNA" |while read -r chr;
do
    tmp="$chr"
    for non_b in "APR" "DR" "G4" "IR" "MR" "TRI" "STR" "Z" "ALL"
    do
      n=`grep $chr densities/${prefix}.nonB_per_chrom.txt |grep " $non_b " |cut -f4 -d" "`
      tmp=$tmp" "$n
    done
    echo $tmp
done
# Per group
for group in "macro" "micro" "microdot"
do
  tmp="$group"
  for non_b in "APR" "DR" "G4" "IR" "MR" "TRI" "STR" "Z" "ALL"
  do
    n=`awk -v g=$group '($1==g){print}' densities/${prefix}.nonB_per_group.txt |grep "$non_b" |cut -f4`
    tmp=$tmp" "$n
  done
    echo $tmp
done

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ GC CONTENT PER CHROM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mkdir stats/
prefix="bTaeGut7v0.4_MT_rDNA"

# Find GC and AT of full genomes
echo '#!/bin/bash
cat ref/'$prefix'.fa | awk '"'"'(!/^>/){gc+=gsub(/[gGcC]/,""); at+=gsub(/[aAtT]/,"")}END{print gc"\t"at}'"'"' >stats/'$prefix'.nuc
'| sbatch -J GC --ntasks=1 --cpus-per-task=1 --time=1:00:00 -o slurm/calculate.gc.$prefix.%j.out

# And for each chromosome separately
# (had to add --cpus-per-task=2 --mem-per-cpu=4G for chr2)
cat ref/${prefix}*.fa.fai |cut -f1 |grep "chr2_" |grep -v "chrM" |grep -v "rDNA" |while read -r chr;
do
  echo '#!/bin/bash
  faidx ref/'$prefix'.fa '$chr' | awk '"'"'(!/^>/){gc+=gsub(/[gGcC]/,""); at+=gsub(/[aAtT]/,"")}END{print gc"\t"at}'"'"' >stats/'$chr'.nuc
  '| sbatch -J GC --ntasks=1 --cpus-per-task=2 --mem-per-cpu=4G --time=15:00 -o slurm/calculate.gc.$chr.%j.out
done

# Merge all and calculate GC content
echo "Chr GC_num AT_num GC_cont AT_cont" |sed 's/ /\t/g' >stats/$prefix.GC_per_chrom.txt
cat ref/${prefix}*.fa.fai |cut -f1 |grep -v "chrM" |grep -v "rDNA" |while read -r chr;
do
 cat stats/$chr.nuc |awk -v c=$chr -v OFS="\t" '{sum=$1+$2; gc=$1/sum; at=$2/sum; print c,$1,$2,gc,at}'
done >>stats/$prefix.GC_per_chrom.txt