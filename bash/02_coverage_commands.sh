################################################################################
### CALCULATING NON-B DNA MOTIF COVERAGE IN THE BIRD GENOMES
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
prefix="bTaeGut7v0.4_MT_rDNA"

# Create genomic windows and make length file, trying both 100kb and 10kb windows
mkdir -p windows/
mkdir -p coverage
# Create genomic windows, trying both 100kb and 10kb windows
cat helpfiles/species_list.txt |tail -n3 |while read -r sp longname prefix
do
  bedtools makewindows -g ref/$prefix.fa.fai -w 100000 >windows/$prefix.100kb_windows.bed;
  bedtools makewindows -g ref/$prefix.fa.fai -w 10000 >windows/$prefix.10kb_windows.bed;
done 

cat helpfiles/species_list.txt |grep zebra |while read -r sp longname prefix
do
  for wind in "100kb" #"10kb" #
  do
    for n in  "Zgfa" #"G4" "APR" "DR" "IR" "TRI" "STR" "Z" #"Zgfa" #"G4" "APR" "DR" "DRfilt" "IR" "IRgfa" "MRgfa" "TRI" "STR" "Z" "Zgfa" #"G4quadron"
    do
      echo '#!/bin/bash
      module load bedtools/2.31.0
      intersectBed -wao -a windows/'$prefix'.'${wind}'_windows.bed -b final_nonB/'$prefix'.'${n}'.merged.bed | cut -f1,2,3,7 |awk -v OFS="\t" '"'"'{if(NR==0){chr=$1; s=$2; e=$3; sum=$4}else{if($1==chr && $2==s){sum+=$4}else{print chr,s,e,sum; chr=$1; s=$2; e=$3; sum=$4}}}END{print chr,s,e,sum}'"'"' | sed "/^\s*$/d" >coverage/'${prefix}'.'${n}'.'$wind'.bed
      '| sbatch -J $n --ntasks=1 --cpus-per-task=1 --mem-per-cpu=6G --out slurm/coverage.$n.%j.out
      done
  done 
done 

# For plotting, merge the data:
cat helpfiles/species_list.txt |while read -r sp longname prefix
do
  echo "NonB Chr Start Stop Dens" |sed 's/ /\t/g' >coverage/$prefix.merged.100kb.txt
  for n in "APR" "DR" "G4" "IR" "TRI" "STR" "Z"
  do
    awk -v OFS="\t" -v nb=$n '{print nb,$0}' coverage/${prefix}.${n}.100kb.bed >>coverage/$prefix.merged.100kb.txt
  done
done 

# For Supplementary Plot
prefix="bTaeGut7v0.4_MT_rDNA"
echo "NonB Chr Start Stop Dens" |sed 's/ /\t/g' >coverage/$prefix.different.merged.100kb.txt
for n in "G4quadron" "G4" "Z" "Zgfa" "ZDNAm1" "Zseeker"
do
  awk -v OFS="\t" -v nb=$n '{print nb,$0}' coverage/${prefix}.${n}.100kb.bed >>coverage/$prefix.different.merged.100kb.txt
done




# ~~~~~~~~~~~~~~~~~~~~~~~ PLOT COVERAGE USING CIRCOS ~~~~~~~~~~~~~~~~~~~~~~~~~~

# Assume circos is installed in ~/software/circos-0.69-9/
# The circos configuration files are found in circos/

# RUN ZEBRA FINCH SEPARATELY AS IT HAS TWO HAPLOTYPES
# Create density and centromere band files separately for MAT and PAT
prefix="bTaeGut7v0.4_MT_rDNA"
for hap in "mat" "pat"
do
    dir="circos/$hap/"
    mkdir -p $dir
    cp -r circos/template/* $dir/
    grep "_$hap" ref/$prefix.fa.fai | sed "s/_$hap//g" |awk '{gsub(/chr/,""); print "chr - nn"$1,$1,"0",$2,"vvlgrey"}' |sort -k4,4n >$dir/karyotype.txt 
     # change order so chrW and chrZ is last
    grep -v "nnZ" $dir/karyotype.txt |grep -v "nnW" >tmp
    grep "nnW" $dir/karyotype.txt >>tmp
    grep "nnZ" $dir/karyotype.txt >>tmp
    mv tmp $dir/karyotype.txt
    for nb in  "APR" "DR" "G4" "IR" "TRI" "STR" "Z"
    do
        grep "_$hap" coverage/$prefix.${nb}.100kb.bed |sed '/^[[:space:]]*$/d' |awk '{gsub(/chr/, "nn"); split($1,s,"_"); $3=$3-1; print s[1],$2,$3,$4}' >$dir/data/${nb}.100kb.txt
        max=`cut -f4 -d" " $dir/data/${nb}.100kb.txt |sort -nr |head -n1`
        echo "will replace max = $nb with max = $max"
        sed -i'' -e "s/max = $nb/max = $max/" $dir/circos_halfInnerCircle.conf
    done
    if [[ "$hap" == "mat" ]]; then 
      sed -i'' -e "s/nnZ/nnW/" $dir/ideogram.conf
    fi
    grep "_$hap" ref/bTaeGut7v0.4_MT_rDNA.centromere_detector.v0.1.gff |awk '{gsub(/chr/,""); split($1,s,"_"); start=$4-1; print "nn"s[1],start,$5}' >$dir/data/highlight_CEN.txt
done

# Plot with
for hap in "mat" "pat"
do
  cd circos/${hap}/
  ~/software/circos-0.69-9/bin/circos -conf circos_halfInnerCircle.conf
  mv circos.png ${hap}_circos.png
  cp ${hap}_circos.png ~/Downloads/
  cd ../..
done

# Centromeres for most chicken chromosome were kindly privided by Luohao Xu. 
mv helpfiles/chicken.CEN.bed ref/

# THE OTHER SPECIES (ONLY CHICKEN HAS CENTROMERE INFO)
cat helpfiles/species_list.txt |grep -v zebra |grep -v chicken |while read -r sp longname prefix
do
  dir="circos/$sp/"
  mkdir -p $dir
  cp -r circos/template/* $dir/
  # Create circos karyotype, nonB density and centromere band files
  # BUT skip chrW for chicken
  grep "chr" ref/$prefix.fa.fai | awk '{gsub(/chr/,""); print "chr - nn"$1,$1,"0",$2,"vvlgrey"}' |sort -k4,4n >$dir/karyotype.txt
  if [ "$sp" == "chicken" ]; then
    echo "Special for chicken"
    grep "chr" ref/$prefix.fa.fai | grep -v "chrW" | awk '{gsub(/chr/,""); print "chr - nn"$1,$1,"0",$2,"vvlgrey"}' |sort -k4,4n >$dir/karyotype.txt 
     # Centromeres
    cat ref/$prefix.CEN.bed |grep -v "chrW" |awk '{gsub(/chr/,""); print "nn"$1,$2,$3}' >$dir/data/highlight_CEN.txt
  else
    # Remove centromere highlights from all but chicken
     sed -i'' -e "s/show_highlights = yes/show_highlights = no/" $dir/circos_halfInnerCircle.conf
  fi
  # change order so chrW and chrZ is last
  grep -v "nnZ" $dir/karyotype.txt |grep -v "nnW" >tmp
  grep "nnW" $dir/karyotype.txt >>tmp
  grep "nnZ" $dir/karyotype.txt >>tmp
  mv tmp $dir/karyotype.txt
  for nb in "TRI" "APR" "DR" "G4" "IR" "STR" "Z"
  do
    grep "chr" coverage/$prefix.${nb}.100kb.bed |sed '/^[[:space:]]*$/d' |awk '{gsub(/chr/, "nn"); $3=$3-1; print $1,$2,$3,$4}' >$dir/data/${nb}.100kb.txt
    max=`cut -f4 -d" " $dir/data/${nb}.100kb.txt |sort -nr |head -n1`
    echo "will replace max = $nb with max = $max"
    sed -i'' -e "s/max = $nb/max = $max/" $dir/circos_halfInnerCircle.conf
  done

  # Plot
  cd $dir
  ~/software/circos-0.69-9/bin/circos -conf circos_halfInnerCircle.conf  
  mv circos.png ${sp}_circos.png
  cp ${sp}_circos.png ~/Downloads/
  cd ../..
done 



########################## ENRICHMENT CALCULATIONS #############################
# To calculate enrichment in a certain region, we need to have a base line to 
# compare the coverage to. We can use either genome-wide coverage, or the coverage
# in a certain chromosome group (macro, micro, microdot), or per chromosome. 

# Calculate genome-wide density to use for normalization
cat helpfiles/species_list.txt |grep pigeon |while read -r sp longname prefix
do
  echo '#!/bin/bash
  module load bedtools/2.31.0
  rm -f coverage/'${prefix}'.per_genome.tsv
  totlen=`cat ref/'${prefix}'*.fa.fai | awk '"'"'{sum+=$2}END{print sum}'"'"'`
  echo "Totlen: "$totlen
  for non_b in  "TRI" "APR" "DR" "G4" "IR" "STR" "Z" "Any"
  do
    cat final_nonB/'${prefix}'.${non_b}.merged.bed | awk -v n=$non_b -v tot=$totlen '"'"'{sum+=$3-$2}END{d=sum/tot; print n, sum, d}'"'"' >>coverage/'${prefix}'.per_genome.tsv
  done
  ' |sbatch -J enrichment --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8G --out slurm/job.enrichment.$prefix.%j.out
done

# different annotation algorithms
cat helpfiles/species_list.txt |grep zebra |while read -r sp longname prefix
do
  echo '#!/bin/bash
  module load bedtools/2.31.0
  rm -f coverage/'${prefix}'.different.per_genome.tsv
  totlen=`cat ref/'${prefix}'*.fa.fai | awk '"'"'{sum+=$2}END{print sum}'"'"'`
  echo "Totlen: "$totlen
  for non_b in  "Z" "Zgfa" "ZDNAm1" "Zseeker" "G4quadron" "IRall" "MRall"
  do
    cat final_nonB/'${prefix}'.${non_b}.merged.bed | \
    awk -v n=$non_b -v tot=$totlen '"'"'{sum+=$3-$2}END{d=sum/tot; print n, sum, d}'"'"' >>coverage/'${prefix}'.different.per_genome.tsv
  done
  ' |sbatch -J enrichment --ntasks=1 --cpus-per-task=1 --out slurm/job.enrichment.$prefix.%j.out
done


# Per chromosome coverage
cat helpfiles/species_list.txt |while read -r sp longname prefix
do
  cat ref/${prefix}.fa.fai |cut -f1,2 |grep -v "chrM" |grep -v "rDNA" |grep -v scaffold |grep -v SCAFFOLD |grep -v contig |grep -v Unplaced |while read -r chr len;
  do
    echo '#!/bin/bash
      module load bedtools/2.31.0
        rm -f tmp.coverage.'$prefix'.'$chr'
      for non_b in  "TRI" "APR" "DR" "G4" "IR" "STR" "Z" "Any"
      do
          cat final_nonB/'$prefix'.${non_b}.merged.bed | awk -v n=$non_b -v tot='$len' -v chr='$chr' '"'"'($1==chr){sum+=$3-$2; c=$1}END{d=sum/tot; print c, n, sum, d}'"'"' >>tmp.coverage.'$prefix'.'$chr' 
      done
  '| sbatch -J $chr --ntasks=1 --cpus-per-task=1 --mem-per-cpu=4G --time=1:00:00
  done
done 

# Merge all chromosomes together into one coverage file 
cat helpfiles/species_list.txt |while read -r sp longname prefix
do
  echo -e "Chr\tNonB\tBp\tCoverage" >coverage/${prefix}.per_chrom.tsv
  for chr in $(cat ref/$prefix.fa.fai |cut -f1 |grep "chr" |grep -v "chrM" ) 
  do
    cat tmp.coverage.$prefix.$chr | sed "s/ /\t/g" >>coverage/$prefix.per_chrom.tsv
  done
done 

# Z-DNA, different annotation algorithms, only ZebraFinch
prefix="bTaeGut7v0.4_MT_rDNA"
cat ref/${prefix}.fa.fai |cut -f1,2 |grep -v "chrM" |grep -v "rDNA" |while read -r chr len;
do
  echo '#!/bin/bash
    module load bedtools/2.31.0
      rm -f tmp.different.coverage.'$prefix'.'$chr'
    for non_b in  "Z" "Zgfa" "ZDNAm1" "Zseeker" "G4quadron" "IRall" "MRall"
    do
        cat final_nonB/'$prefix'.${non_b}.merged.bed | \
        awk -v n=$non_b -v tot='$len' -v chr='$chr' '"'"'($1==chr){sum+=$3-$2; c=$1}END{d=sum/tot; print c, n, sum, d}'"'"' >>tmp.different.coverage.'$prefix'.'$chr' 
    done
'| sbatch -J $chr --ntasks=1 --cpus-per-task=1 --mem-per-cpu=4G --time=1:00:00
done
echo -e "Chr\tNonB\tBp\tCoverage" >coverage/${prefix}.different.per_chrom.tsv
for chr in $(cat ref/${prefix}.fa.fai |cut -f1  |grep -v "chrM" |grep -v "rDNA") 
do
  cat tmp.different.coverage.$prefix.${chr} | sed "s/ /\t/g" >>coverage/${prefix}.different.per_chrom.tsv
done 


# Per chromosome category coverage
cat helpfiles/species_list.txt |grep zebra |while read -r sp longname prefix
do
  for group in "macro" "micro" "dot" "unplaced"
  do
    totlen=`grep $group helpfiles/${prefix}.groups.txt |awk 'NR==FNR{a[$1];next} ($1 in a){print}' - ref/${prefix}*.fa.fai | awk '{sum+=$2}END{print sum}'`
    echo '#!/bin/bash
    module load bedtools/2.31.0
    rm -f tmp.coverage.'$prefix'.'${group}'
    for non_b in "TRI" "APR" "DR" "G4" "IR" "STR" "Z" "Any"
    do
        grep '$group' helpfiles/'${prefix}'.groups.txt |awk '"'"'NR==FNR{a[$1];next} ($1 in a){print}'"'"' - final_nonB/'${prefix}'.${non_b}.merged.bed |\
         awk -v g='$group '-v n=$non_b -v tot='$totlen' '"'"'{sum+=$3-$2}END{d=sum/tot; print g, n, sum, d}'"'"' >>tmp.coverage.'$prefix'.'${group}'
    done
    '| sbatch -J per_group --ntasks=1 --cpus-per-task=1 --account=kdm16_cr_default --partition=standard  --mem-per-cpu=2G --time=1:00:00
  done
done 
# Merge all categories together into one coverage file 
cat helpfiles/species_list.txt |grep humming |while read -r sp longname prefix
do
  echo "Group NonB Bp Coverage" |sed "s/ /\t/g" >coverage/${prefix}.per_group.tsv
  for group in "macro" "micro" "dot"
  do
    cat tmp.coverage.$prefix.$group | sed "s/ /\t/g" >>coverage/${prefix}.per_group.tsv
  done
done 

# Different annotation algorithms 
prefix="bTaeGut7v0.4_MT_rDNA"
for group in "macro" "micro" "dot"
do
  totlen=`grep $group helpfiles/${prefix}.groups.txt |awk 'NR==FNR{a[$1];next} ($1 in a){print}' - ref/${prefix}*.fa.fai | awk '{sum+=$2}END{print sum}'`
  echo '#!/bin/bash
  module load bedtools/2.31.0
  rm -f tmp.coverage.different.'$prefix'.'${group}'
  for non_b in "IRall" "MRall" "G4quadron" "Zgfa" "ZDNAm1" "Zseeker" 
  do
      grep '$group' helpfiles/'${prefix}'.groups.txt |awk '"'"'NR==FNR{a[$1];next} ($1 in a){print}'"'"' - final_nonB/'${prefix}'.${non_b}.merged.bed |\
        awk -v g='$group '-v n=$non_b -v tot='$totlen' '"'"'{sum+=$3-$2}END{d=sum/tot; print g, n, sum, d}'"'"' >>tmp.coverage.different.'$prefix'.'${group}'
  done
  '| sbatch -J per_group --ntasks=1 --cpus-per-task=1 --account=kdm16_cr_default --partition=basic  --mem-per-cpu=2G --time=1:00:00
done
echo "Group NonB Bp Coverage" |sed "s/ /\t/g" >coverage/${prefix}.different.per_group.tsv
for group in "macro" "micro" "dot"
do
  cat tmp.coverage.different.$prefix.$group | sed "s/ /\t/g" >>coverage/${prefix}.different.per_group.tsv
done




Group   NonB    Bp      Coverage
macro   APR     6907942 0.00479662
macro   DR      19542669        0.0135697
macro   G4      8922817 0.00619567
macro   IR      45793751        0.0317975
macro   TRI     3448355 0.00239441
macro   STR     15924025        0.011057
macro   Z       1420599 0.000986411
macro   Any     87448226        0.0607208
micro   APR     2177569 0.0042814
micro   DR      14555122        0.0286173
micro   G4      9093253 0.0178786
micro   IR      16172971        0.0317983
micro   TRI     1148985 0.00225906
micro   STR     6957555 0.0136795
micro   Z       524749  0.00103173
micro   Any     42172980        0.0829178
dot     APR     3316279 0.0154941
dot     DR      18089157        0.0845149
dot     G4      11449333        0.0534927
dot     IR      7008734 0.0327457
dot     TRI     650909  0.00304113
dot     STR     3718729 0.0173744
dot     Z       186347  0.000870637
dot     Any     33871462        0.158252

# ~~~~~~~~~~~~~~~~~~~~~~~~~ PRINT COVERAGE FOR TABLE ~~~~~~~~~~~~~~~~~~~~~~~~~~
# FOR SUPPLEMENTARY TABLE 1
prefix="bTaeGut7v0.4_MT_rDNA"
# Per chromosome
cat ref/${prefix}.fa.fai |cut -f1 |grep "chr" |grep -v "chrM" |grep -v "rDNA" |while read -r chr;
do
    tmp="$chr"
    for non_b in "APR" "DR" "STR" "IR" "TRI" "G4" "Z" "Any" #IRall" "MRall" "G4quadron" "Zgfa" "ZDNAm1" "Zseeker" 
    do
      n=`awk -v chr=$chr -v n=$non_b '($1==chr && $2==n){print $4}' <(cat coverage/${prefix}.per_chrom.tsv coverage/${prefix}.different.per_chrom.tsv |sort |uniq)`
      tmp=$tmp" "$n
    done
    echo $tmp
done
# Whole genome
tmp="genome"
for non_b in "APR" "DR" "STR" "IR" "TRI" "G4" "Z" "Any" #IRall" "MRall" "G4quadron" "Zgfa" "ZDNAm1" "Zseeker" 
do
  n=`awk -v n=$non_b '($1==n){print $3}' <(cat coverage/${prefix}.per_genome.tsv coverage/${prefix}.different.per_genome.tsv |sort |uniq)`
  tmp=$tmp" "$n
done
echo $tmp
# Per group
for group in "macro" "micro" "dot"
do
  tmp="$group"
  for non_b in "APR" "DR" "STR" "IR" "TRI" "G4" "Z" "Any" "IRall" "MRall" "G4quadron" "Zgfa" "ZDNAm1" "Zseeker" 
  do
    n=`awk -v g=$group -v n=$non_b '($1==g && $2==n){print $4}' <(cat coverage/${prefix}.per_group.tsv coverage/${prefix}.different.per_group.tsv)`
    tmp=$tmp" "$n
  done
    echo $tmp
done


# SUPPLEMENTARY TABLE 3-9
cat helpfiles/species_list.txt |grep emu |while read -r sp longname prefix
do
  echo "+++++++++++++ PROCESSING $sp +++++++++++++++++++++"
  # Per chromosome
  cat ref/${prefix}.fa.fai |cut -f1 |grep "chr" |grep -v "chrM" |grep -v "rDNA" |while read -r chr;
  do
      tmp="$chr"
      for non_b in "APR" "DR" "STR" "IR" "TRI" "G4" "Z" "Any"
      do
        n=`awk -v chr=$chr -v n=$non_b '($1==chr && $2==n){print $4}' coverage/${prefix}.per_chrom.tsv`
        tmp=$tmp" "$n
      done
      echo $tmp
  done
  # Whole genome
  tmp="genome"
  for non_b in "APR" "DR" "STR" "IR" "TRI" "G4" "Z" "Any"
  do
    n=`awk -v n=$non_b '($1==n){print $3}' coverage/${prefix}.per_genome.tsv`
    tmp=$tmp" "$n
  done
  echo $tmp
  # Per group
  for group in "macro" "micro" "dot"
  do
    tmp="$group"
    for non_b in "APR" "DR" "STR" "IR" "TRI" "G4" "Z" "Any" 
    do
      n=`awk -v g=$group -v n=$non_b '($1==g && $2==n){print $4}' coverage/${prefix}.per_group.tsv`
      tmp=$tmp" "$n
    done
      echo $tmp
  done
done 



# ~~~~~~~~~~~~~~~~~~~~ PRINT NUMBER OF MOTIFS AND BASEPAIRS ~~~~~~~~~~~~~~~~~~~~
# SUPPLEMENTARY TABLE 2
# Write it out with the nonB types as columns

cat helpfiles/species_list.txt |grep zebra |while read -r sp longname prefix
do
  # Per chromosome
  cat ref/${prefix}.fa.fai |cut -f1 |grep -v "chrM" |grep -v "rDNA" |while read -r chr;
  do
      tmp="$chr"
      for non_b in "APR" "DR" "STR" "IR"  "TRI" "G4" "Z" "Any"
      do
        n=`grep $chr final_nonB/bTaeGut7v0.4_MT_rDNA.$non_b.bed |wc -l |cut -f1 -d" "`
        len=`awk -v chr=$chr -v n=$non_b '($1==chr && $2==n){print $3}' coverage/${prefix}.per_chrom.tsv`
        tmp=$tmp" "$n" "$len
      done
      echo $tmp
  done
  # Whole genome
  tmp="genome"
  for non_b in "APR" "DR" "STR" "IR" "TRI" "G4" "Z" "Any"
  do
    n=`wc -l final_nonB/bTaeGut7v0.4_MT_rDNA.$non_b.bed |cut -f1 -d" "`
  #  echo $n
    len=`awk -v n=$non_b '($1==n){print $2}' coverage/${prefix}.per_genome.tsv`
    tmp=$tmp" "$n" "$len
  done
  echo $tmp
  # Per group
  for group in "macro" "micro" "dot"
  do
    tmp="$group"
    for non_b in "APR" "DR" "STR" "IR" "TRI" "G4" "Z" "Any"
    do
        n=`grep $group helpfiles/$prefix.groups.txt |awk 'NR==FNR{a[$1];next} ($1 in a){print}' - final_nonB/$prefix.$non_b.bed |wc -l |cut -f1 -d" "`
      len=`awk -v g=$group -v n=$non_b '($1==g && $2==n){print $3}' coverage/$prefix.per_group.tsv`
      tmp=$tmp" "$n" "$len
    done
    echo $tmp
  done
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
cat helpfiles/species_list.txt |grep -v zebra |while read -r sp longname prefix
do
  echo '#!/bin/bash
  cat ref/'$prefix'.fa | awk '"'"'(!/^>/){gc+=gsub(/[gGcC]/,""); at+=gsub(/[aAtT]/,"")}END{print gc"\t"at}'"'"' >stats/'$prefix'.nuc
  '| sbatch -J GC --ntasks=1 --cpus-per-task=1 --time=1:00:00 -o slurm/calculate.gc.$prefix.%j.out
done 
# And for each chromosome separately
# (had to add --cpus-per-task=2 --mem-per-cpu=4G for chr2)
cat helpfiles/species_list.txt |tail -n4 |while read -r sp longname prefix
do
  cat ref/${prefix}*.fa.fai |cut -f1 |grep "chr" |grep -v "chrM" |grep -v "rDNA" |while read -r chr;
  do
    echo '#!/bin/bash
    faidx ref/'$prefix'.fa '$chr' | awk '"'"'(!/^>/){gc+=gsub(/[gGcC]/,""); at+=gsub(/[aAtT]/,"")}END{print gc"\t"at}'"'"' >stats/'$prefix'.'$chr'.nuc
    '| sbatch -J $sp.GC --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8G --time=15:00 -o slurm/calculate.gc.$prefix.$chr.%j.out
  done
done 

# Merge all and calculate GC content
cat helpfiles/species_list.txt  |grep -v zebra |while read -r sp longname prefix
do
  echo "Chr GC_num AT_num GC_cont AT_cont" |sed 's/ /\t/g' >stats/$prefix.GC_per_chrom.txt
  cat ref/${prefix}*.fa.fai |cut -f1 |grep "chr" |grep -v "chrM" |grep -v "rDNA" |while read -r chr;
  do
    cat stats/$prefix.$chr.nuc |awk -v c=$chr -v OFS="\t" '{sum=$1+$2; gc=$1/sum; at=$2/sum; print c,$1,$2,gc,at}'
  done >>stats/$prefix.GC_per_chrom.txt
done 