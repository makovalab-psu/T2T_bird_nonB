################################################################################
### CALCULATING NON-B DNA MOTIF DENSITY IN THE ZEBRA FINCH T2T GENOME
### CODE WRITTEN BY LINNÉA SMEDS

##### TABLE OF CONTENTS
# CALCULATE COVERAGE
#   -PLOTTING COVERAGE USING CIRCOS
# PREPARE ENRICHMENT CALCULATIONS
#   -PRINT DENSITIES FOR TABLE
#   -PRINT NUMBER OF MOTIFS AND BASEPAIRS
#   -SPACER LENGTH DISTRIBUTION FOR DR, IR AND MR
#   -GC CONTENT PER CHROMOSOME


######################## CALCULATE COVERAGE IN WINDOWS ##########################
# June 24, 2025
prefix="bTaeGut7v0.4_MT_rDNA"

# Create genomic windows and make length file, trying both 100kb and 10kb windows
mkdir -p windows/
bedtools makewindows -g ref/$prefix.fa.fai -w 100000 >windows/$prefix.100kb_windows.bed;
bedtools makewindows -g ref/$prefix.fa.fai -w 10000 >windows/$prefix.10kb_windows.bed;

mkdir coverage
for wind in "10kb" "100kb"
do
  for n in  "G4" "APR" "DR" "DRfilt" "IR" "IRall" "MRall" "TRI" "STR" "Z" "ZDNAm1" "Zseeker" "G4quadron"
  do
   echo '#!/bin/bash
   module load bedtools/2.31.0
   intersectBed -wao -a windows/$prefix.'${wind}'_windows.bed -b <(cut -f1,2,3 final_nonB/'$prefix'.'${n}'.merged.bed) | cut -f1,2,3,7 |awk -v OFS="\t" '"'"'{if(NR==0){chr=$1; s=$2; e=$3; sum=$4}else{if($1==chr && $2==s){sum+=$4}else{print chr,s,e,sum; chr=$1; s=$2; e=$3; sum=$4}}}END{print chr,s,e,sum}'"'"' | sed "/^\s*$/d" >coverage/'${prefix}'.'${n}'.'$wind'.bed
   '| sbatch -J $n --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8G --out slurm/density.$n.%j.out
  done
done

# For plotting, merge the data:
for wind in "10kb" "100kb"
do
  echo "NonB Chr Start Stop Coverage" |sed 's/ /\t/g' >coverage/$prefix.merged.$wind.txt
  for n in  "G4" "APR" "DR" "DRfilt" "IR" "IRall" "MRall" "TRI" "STR" "Z" "ZDNAm1" "Zseeker" "G4quadron"
  do
    awk -v OFS="\t" -v nb=$n '{print nb,$0}' coverage/${prefix}.${n}.$wind.bed >>coverage/$prefix.merged.$wind.txt
  done
done


# ~~~~~~~~~~~~~~~~~~~~~~~ PLOT COVERAGE USING CIRCOS ~~~~~~~~~~~~~~~~~~~~~~~~~~
# June 26 2025

# Assume circos is installed in ~/software/circos-0.69-9/
# The circos configuration files are found in circos/

# Create density and centromere band files separately for MAT and PAT
prefix="bTaeGut7v0.4_MT_rDNA"
for hap in "mat" "pat"
do
    mkdir -p circos/$hap
    for nb in  "APR" "DR" "G4" "IR" "TRI" "STR" "Z"
    do
        grep "_$hap" coverage/$prefix.${nb}.100kb.bed |sed '/^[[:space:]]*$/d' |awk '{gsub(/chr/, "nn"); split($1,s,"_"); $3=$3-1; print s[1],$2,$3,$4}' >circos/$hap/data/${nb}.100kb.txt
        max=`cut -f4 -d" " circos/$hap/data/${nb}.100kb.txt |sort -nr |head -n1`
        # Make sure the right max is set for the y scale 
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
# compare the coverage to. We can use either genome-wide coverage, or the coverage
# in a certain chromosome group (macro, micro, microdot), or per chromosome. 

# Calculate genome-wide density to use for normalization
prefix="bTaeGut7v0.4_MT_rDNA"
echo '#!/bin/bash
module load bedtools/2.31.0
rm -f coverage/'${prefix}'.nonB_genome_wide.txt
totlen=`cat ref/'${prefix}'*.fa.fai | awk '"'"'{sum+=$2}END{print sum}'"'"'`
echo "Totlen: "$totlen
for non_b in  "G4" "APR" "DR" "DRfilt" "IR" "IRall" "MRall" "TRI" "STR" "Z" "Zgfa" "ZDNAm1" "Zseeker" "G4quadron" "Any"
do
  cat final_nonB/'${prefix}'.${non_b}.merged.bed | awk -v n=$non_b -v tot=$totlen '"'"'{sum+=$3-$2}END{d=sum/tot; print n, sum, d}'"'"' >>coverage/'${prefix}'.nonB_genome_wide.txt
done
' |sbatch -J enrichment --ntasks=1 --cpus-per-task=2 --mem-per-cpu=8G --out slurm/job.enrichment.$prefix.%j.out


# Per chromosome coverage
cat ref/${prefix}.fa.fai |cut -f1,2 |grep -v "chrM" |grep -v "rDNA" |while read -r chr len;
do
  echo '#!/bin/bash
    module load bedtools/2.31.0
      rm -f tmp.coverage.'$chr'
    for non_b in  "G4" "APR" "DR" "DRfilt" "IR" "IRall" "MRall" "TRI" "STR" "Z" "Zgfa" "ZDNAm1" "Zseeker" "G4quadron" "Any"
    do
        cat final_nonB/'$prefix'.${non_b}.merged.bed | awk -v n=$non_b -v tot='$len' -v chr='$chr' '"'"'($1==chr){sum+=$3-$2; c=$1}END{d=sum/tot; print c, n, sum, d}'"'"' >>tmp.coverage.'$chr' 
    done
'| sbatch -J $chr --ntasks=1 --cpus-per-task=1 --mem-per-cpu=4G --time=1:00:00 --partition=open
done

# Merge all chromosomes together into one coverage file 
echo -e "Chr\tNonB\tBp\tCoverage" >coverage/${prefix}.nonB_per_chrom.txt
for chr in $(cat ref/${prefix}.fa.fai |cut -f1 |grep -v "chrM" |grep -v "rDNA") 
do
  cat tmp.coverage.${chr} | sed "s/ /\t/g" >>coverage/${prefix}.nonB_per_chrom.txt
done

# Per chromosome category coverage
for group in "macro" "micro" "microdot"
do
  totlen=`grep -f helpfiles/$group.txt ref/${prefix}*.fa.fai | awk '{sum+=$2}END{print sum}'`
  echo '#!/bin/bash
  module load bedtools/2.31.0
  rm -f tmp.coverage.'${group}'
  for non_b in  "G4" "APR" "DR" "DRfilt" "IR" "IRall" "MRall" "TRI" "STR" "Z" "Zgfa" "ZDNAm1" "Zseeker" "G4quadron" "Any"
  do
      grep -f helpfiles/'$group'.txt final_nonB/'${prefix}'.${non_b}.merged.bed | awk -v g='$group '-v n=$non_b -v tot='$totlen' '"'"'{sum+=$3-$2}END{d=sum/tot; print g, n, sum, d}'"'"' >>tmp.coverage.'${group}'
  done
  '| sbatch -J per_group --ntasks=1 --cpus-per-task=1 --mem-per-cpu=4G --time=1:00:00 --partition=open
done

# Merge all categories together into one coverage file 
echo "Group NonB Bp Coverage" |sed "s/ /\t/g" >coverage/${prefix}.nonB_per_group.txt
for group in "macro" "micro" "microdot"
do
  cat tmp.coverage.${group} | sed "s/ /\t/g" >>coverage/${prefix}.nonB_per_group.txt
done

# ~~~~~~~~~~~~~~~~~~~~~~~~~ PRINT COVERAGE FOR TABLE ~~~~~~~~~~~~~~~~~~~~~~~~~~
# FOR SUPPLEMENTARY TABLE 1
# Whole genome
prefix="bTaeGut7v0.4_MT_rDNA"
tmp="genome"
for non_b in "APR" "DR" "IR" "TRI" "STR"  "G4"  "Z" "Any" # "DRfilt" "IRall" "MRall" "G4quadron" "Zgfa" "ZDNAm1" "Zseeker" 
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


# ~~~~~~~~~~~~~~~~~~~~ PRINT NUMBER OF MOTIFS AND BASEPAIRS ~~~~~~~~~~~~~~~~~~~~
cd /storage/group/kdm16/default/lbs5874/ZebraFinch

# Write it out with the nonB types as columns
# Whole genome
prefix="bTaeGut7v0.4_MT_rDNA"
tmp="genome"
for non_b in "APR" "DR" "G4" "IR" "MR" "TRI" "STR" "Z" "ALL"
do
  n=`wc -l final_nonB/bTaeGut7v0.4_MT_rDNA.$non_b.bed |cut -f1 -d" "`
#  echo $n
  len=`grep $non_b densities/${prefix}.nonB_genome_wide.txt |cut -f2 -d" "`
#mb=`echo "scale=2; $len/1000000.00" |bc`
  tmp=$tmp" "$n" "$len
done
echo $tmp
# Per chromosome
cat ref/${prefix}*.fa.fai |cut -f1 |grep -v "chrM" |grep -v "rDNA" |while read -r chr;
do
    tmp="$chr"
    for non_b in "APR" "DR" "G4" "IR" "MR" "TRI" "STR" "Z" "ALL"
    do
      n=`grep $chr final_nonB/bTaeGut7v0.4_MT_rDNA.$non_b.bed |wc -l |cut -f1 -d" "`
      len=`grep $chr densities/${prefix}.nonB_per_chrom.txt |grep " $non_b " |cut -f3 -d" "`
      tmp=$tmp" "$n" "$len
    done
    echo $tmp
done
# Per group
for group in "macro" "micro" "microdot"
do
  tmp="$group"
  for non_b in "APR" "DR" "G4" "IR" "MR" "TRI" "STR" "Z" "ALL"
  do
      n=`grep -f $group.txt final_nonB/bTaeGut7v0.4_MT_rDNA.$non_b.bed |wc -l |cut -f1 -d" "`
    len=`awk -v g=$group '($1==g){print}' densities/${prefix}.nonB_per_group.txt |grep "$non_b" |cut -f3`
    tmp=$tmp" "$n" "$len
  done
  echo $tmp
done


# ~~~~~~~~~~~~~~~~~ CHECK THE SPACER LENGTH FOR DR, IR AND MR ~~~~~~~~~~~~~~~~~~
# Calculate the distribution of spacer lengths, put all in one file

prefix="bTaeGut7v0.4_MT_rDNA"
echo "Count Spacer NonB" |sed 's/ /\t/g' >stats/$prefix.spacers.txt
for nb in "DR" "IR" "MR"
do
 cut -f10 gfa_annotation/${prefix}_${nb}.tsv |tail -n+2 |sort |uniq -c |\
sort -k2,2n |awk -v OFS="\t" -v nb=$nb '{print $1,$2,nb}' >>stats/$prefix.spacers.txt
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