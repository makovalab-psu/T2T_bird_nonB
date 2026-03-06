################################################################################
### CALCULATING NON-B DNA MOTIF ENRICHMENT IN ZEBRA FINCH AND CHICKEN CENTROMERES
### CODE WRITTEN BY LINNÉA SMEDS

##### TABLE OF CONTENTS
# ENRICHMENT IN CENTROMERES
#   -GC CONTENT FOR CENTROMERES
#   -ENRICHMENT IN WINDOWS FOR STATISTICAL TEST



######################## ENRICHMENT IN CENTROMERES  ############################
mkdir -p centromeres
# Centromeres for chicken were kindly provided by Luhaou Xu (from Huang et al 2023)
# Convert Zebra finch centromere gff to bed file
awk -v OFS="\t" '{start=$4-1; print $1,start,$5}' ref/bTaeGut7v0.4_MT_rDNA.centromere_detector.v0.1.gff >ref/bTaeGut7v0.4_MT_rDNA.CEN.bed

# Enrichment in centromeres compared to the genome-wide density and to chrom-
# specific density
cat helpfiles/species_list.txt |grep chicken|while read -r sp longname prefix
do
  echo "Working on $sp centromeres"
  mkdir -p centromeres
  cen_bed_file="ref/${prefix}.CEN.bed"
  cat $cen_bed_file |cut -f1-3 |while read -r chr start end;
  do
    cenlen=`echo "$end-$start" |bc`
    echo cenlen is $cenlen for chr $chr
    echo '#!/bin/bash
    module load bedtools/2.31.0
    echo Working on '$chr'
    rm -f centromeres/'$prefix'.'$chr'.tmp-centromere.txt
    cat coverage/'${prefix}'.per_genome.tsv |grep -v "Any" | while read -r nb tot dens;
    do
      chrdens=`awk -v c='$chr' -v n=$nb '"'"'($2==n && $1==c){print $4}'"'"' coverage/'${prefix}'.per_chrom.tsv `
      echo "Coverage for '$chr' and $nb is $chrdens"
      stats=`echo '$chr' '$start' '$end' |sed "s/ /\t/g" |intersectBed -a - -b final_nonB/'${prefix}'.$nb.merged.bed | sort -k1,1 -k2,2n | mergeBed -i - |\
      awk -v l='$cenlen' -v dtot=$dens -v dchr=$chrdens '"'"'{sum+=$3-$2}END{d=sum/l; fracgw=d/dtot; fracchr=d/dchr; print d,fracgw,fracchr}'"'"'`
      echo '$chr'" "$nb" "$stats |sed "s/ /\t/g" >>centromeres/'$prefix'.'$chr'.tmp-centromere.txt
    done
    ' | sbatch -J $chr --ntasks=1 --cpus-per-task=1 --mem-per-cpu=4G --account=kdm16_sc_default --partition=sla-prio --out slurm/job.centromeric_enrichment.$prefix.$chr.%j.out
  done
done 


#Merge
cat helpfiles/species_list.txt |grep chicken |while read -r sp longname prefix
do
  cen_bed_file="ref/${prefix}.CEN.bed"
  echo "Chr NonB Coverage Enrichment_GW Enrichment_Chr" |sed "s/ /\t/g" >centromeres/${prefix}.enrichment.tsv
  cat $cen_bed_file |cut -f1 |while read -r chr;
  do
    cat centromeres/$prefix.$chr.tmp-centromere.txt >> centromeres/${prefix}.enrichment.tsv
  done
done 
# Remove temporary files 
#rm centromeres/*.tmp-centromere*


# ~~~~~~~~~~~~~~~~~ ENRICHMENT OF Z FOR DIFFERENT ANNOTATIONS ~~~~~~~~~~~~~~~~~~
# ONLY FOR ZEBRA FINCH 
prefix="bTaeGut7v0.4_MT_rDNA"
cen_bed_file="ref/${prefix}.CEN.bed"
cat $cen_bed_file |cut -f1-3 |while read -r chr start end;
do
  cenlen=`echo "$end-$start" |bc`
  echo cenlen is $cenlen for chr $chr
  echo '#!/bin/bash
  module load bedtools/2.31.0
  echo Working on '$chr'
  rm -f centromeres/'$prefix'.'$chr'.ZDNA.tmp-centromere.txt
  cat coverage/'${prefix}'.differentZ.per_genome.tsv | while read -r nb tot dens;
  do
    chrdens=`awk -v c='$chr' -v n=$nb '"'"'($2==n && $1==c){print $4}'"'"' coverage/'${prefix}'.differentZ.per_chrom.tsv`
    echo "Coverage for '$chr' and $nb is $chrdens"
    stats=`echo '$chr' '$start' '$end' |sed "s/ /\t/g" |intersectBed -a - -b final_nonB/'${prefix}'.$nb.merged.bed | sort -k1,1 -k2,2n | mergeBed -i - |\
    awk -v l='$cenlen' -v dtot=$dens -v dchr=$chrdens '"'"'{sum+=$3-$2}END{d=sum/l; fracgw=d/dtot; fracchr=d/dchr; print d,fracgw,fracchr}'"'"'`
    echo '$chr'" "$nb" "$stats |sed "s/ /\t/g" >>centromeres/'$prefix'.'$chr'.ZDNA.tmp-centromere.txt
  done
  ' | sbatch -J $chr --ntasks=1 --cpus-per-task=1 --mem-per-cpu=4G --account=kdm16_sc_default --partition=sla-prio --out slurm/job.diffZ.centromeric_enrichment.$prefix.$chr.%j.out
done 
# Merge
echo "Chr NonB Coverage Enrichment_GW Enrichment_Chr" |sed "s/ /\t/g" >centromeres/${prefix}.ZDNA.enrichment.tsv
cat $cen_bed_file |cut -f1 |while read -r chr;
do
  cat centromeres/$prefix.$chr.ZDNA.tmp-centromere.txt >> centromeres/${prefix}.ZDNA.enrichment.tsv
done 

# ~~~~~~~~~~~~~~~~~~~~~~~~~ GC CONTENT FOR CENTROMERES ~~~~~~~~~~~~~~~~~~~~~~~~~
# Check GC content in centromeres (this was not used in the paper)
prefix="bTaeGut7v0.4_MT_rDNA"
cen_bed_file="ref/${prefix}.centromere_detector.v0.1.bed"
module load bedtools/2.31.0
bedtools nuc -fi ref/$prefix.fa -bed $cen_bed_file >centromeres/cen_GC_fullInfo.nuc
echo "Chr GCcont" |sed 's/ /\t/g' >centromeres/cen_GC.tsv
awk '(NR>1){print $1"\t"$5}' centromeres/cen_GC_fullInfo.nuc >> centromeres/cen_GC.tsv


# ~~~~~~~~~~~~~~~~~ ENRICHMENT IN WINDOWS FOR STATISTICAL TEST ~~~~~~~~~~~~~~~~~
# To have something to compare the centromere to, we divide the rest of the
# chromosomes into segments with the same length as each centromere (excluding
# the actual centromere)

module load bedtools/2.31.0
mkdir -p centromeres/windows/

# Generating 100 windows of the same size as the centromere, not overlapping with
# the real centromere. Randomly selected if there are more than 100 windows,
# overlapping if there are less. For very large centromeres, removing overlap
# with the real centromere and non-complete windows at the end of the chromosome
# can result in less than 100 windows. To make sure there always are more than
# 100, I make the overlap so there are 120 windows to choose from, and allow end
# windows to be as small as 1M (if this is smaller than the centromere size).
cat helpfiles/species_list.txt |head -n2 |while read -r sp longname prefix
do
  echo "Working on $sp centromeric windows"
  cen_bed_file="ref/${prefix}.CEN.bed"
  for chr in $(cut -f1 $cen_bed_file |uniq)
  do
      echo '#!/bin/bash
      echo "Looking at '$chr'"
      module load bedtools/2.31.0
      c_len=`awk -v c='$chr' '"'"'($1==c){sum+=$3-$2}END{print sum}'"'"' '$cen_bed_file'`
 #     chrlen=`awk -v c='$chr' '"'"'($1==c){print $2}'"'"' ref/'$prefix'.fa.fai`
 #     poswin=`echo $chrlen"/"$c_len |bc`
 #     minsize=100000
 #     if [ "$poswin" -lt "120" ]
 #     then
 #       echo "For '$chr', there will be less than 120 nonovl windows"
 #       slide=`echo "("$chrlen"-"$c_len")/120"|bc`
 #       bedtools makewindows -g <(awk -v c='$chr' '"'"'(c==$1){print}'"'"' ref/'$prefix'.fa.fai) -s $slide -w $c_len |subtractBed -a - -b '$cen_bed_file' |awk -v min=$minsize -v clen=$c_len '"'"'($3-$2>min || $3-$2==clen){print}'"'"' |shuf -n 100 |sort -k2,2n >centromeres/windows/'$prefix'.'$chr'.exclCen.bed
 #     else
 #       echo "For '$chr', there will be more than 120 nonovl windows"
 #       bedtools makewindows -g <(awk -v c='$chr' '"'"'(c==$1){print}'"'"' ref/'$prefix'.fa.fai) -w $c_len |subtractBed -a - -b '$cen_bed_file' |awk -v min=$minsize -v clen=$c_len '"'"'($3-$2>min || $3-$2==clen){print}'"'"' |shuf -n 100 |sort -k2,2n >centromeres/windows/'$prefix'.'$chr'.exclCen.bed
 #     fi
 #     bedtools nuc -fi ref/'$prefix'.fa -bed '$cen_bed_file' >centromeres/windows/'$prefix'.'$chr'.exclCen.with_GC.bed
      rm -f centromeres/windows/'$prefix'.'$chr'.100rand.enrichment.tsv
      while read line
      do
        rm -f tmp.centro.$SLURM_JOB_ID.bed tmp.centro.$SLURM_JOB_ID
        cat coverage/'${prefix}'.per_genome.tsv | while read -r non_b tot dens;
        do
            chrcov=`awk -v c='$chr' -v n=$non_b '"'"'($2==n && $1==c){print $4}'"'"' coverage/'${prefix}'.per_chrom.tsv `
            echo $line |sed "s/ /\t/g" >tmp.centro.$SLURM_JOB_ID.bed
            d=`intersectBed -a tmp.centro.$SLURM_JOB_ID.bed -b final_nonB/'$prefix'.${non_b}.merged.bed -wo -nonamecheck |awk -v l=$c_len -v cov=$dens -v chrcov=$chrcov '"'"'{sum+=$7}END{d=sum/l; frac=d/cov; chfrac=d/chrcov; print frac, chfrac}'"'"'`
            echo '$chr'" "$non_b" "$d >>tmp.centro.$SLURM_JOB_ID
        done
        cat tmp.centro.$SLURM_JOB_ID |sed "s/ /\t/g" >>centromeres/windows/'$prefix'.'$chr'.100rand.enrichment.tsv
      done <centromeres/windows/'$prefix'.'$chr'.exclCen.bed
      '| sbatch -J $chr --ntasks=1 --cpus-per-task=1 --time=5:00:00 --account=kdm16_sc_default --partition=sla-prio  --out slurm/job.centromeric_windows.$chr.%j.out
  done
done 

# Check if there are more or less than 120 windows
cat slurm/job.centromeric_enrichment.chr*_*at.*.out |grep "there will be"

# Merge background densities and GC
cat helpfiles/species_list.txt |head -n2 |while read -r sp longname prefix
do
  echo -e "Chr\tNonB\tEnrichment_GW\tEnrichment_CHR" >centromeres/${prefix}.background.enrichment.tsv
  cat centromeres/windows/$prefix.chr*.100rand.enrichment.tsv  >>centromeres/${prefix}.background.enrichment.tsv
done 

# ZEBRAFINCH Z-DNA COMPARISON 
prefix="bTaeGut7v0.4_MT_rDNA"
cen_bed_file="ref/${prefix}.CEN.bed"
for chr in $(cut -f1 $cen_bed_file |uniq)
do
  echo '#!/bin/bash
    echo "Looking at '$chr'"
    c_len=`awk -v c='$chr' '"'"'($1==c){sum+=$3-$2}END{print sum}'"'"' '$cen_bed_file'`
    module load bedtools/2.31.0
    rm -f centromeres/windows/'$prefix'.'$chr'.100rand.ZDNA.enrichment.tsv
    while read line
    do
    rm -f tmp.ZDNA.centro.$SLURM_JOB_ID.bed tmp.ZDNA.centro.$SLURM_JOB_ID
    cat coverage/'${prefix}'.differentZ.per_genome.tsv | while read -r non_b tot dens;
    do
        chrcov=`awk -v c='$chr' -v n=$non_b '"'"'($2==n && $1==c){print $4}'"'"' coverage/'${prefix}'.differentZ.per_chrom.tsv`
        echo $line |sed "s/ /\t/g" >tmp.ZDNA.centro.$SLURM_JOB_ID.bed
        d=`intersectBed -a tmp.ZDNA.centro.$SLURM_JOB_ID.bed -b final_nonB/'$prefix'.${non_b}.merged.bed -wo -nonamecheck |\
        awk -v l=$c_len -v cov=$dens -v chrcov=$chrcov '"'"'{sum+=$NF}END{d=sum/l; frac=d/cov; chfrac=d/chrcov; print frac, chfrac}'"'"'`
        echo '$chr'" "$non_b" "$d >>tmp.ZDNA.centro.$SLURM_JOB_ID
    done
    cat tmp.ZDNA.centro.$SLURM_JOB_ID |sed "s/ /\t/g" >>centromeres/windows/'$prefix'.'$chr'.100rand.ZDNA.enrichment.tsv
  done <centromeres/windows/'$prefix'.'$chr'.exclCen.bed
  '| sbatch -J $chr --ntasks=1 --cpus-per-task=1 --time=5:00:00 --account=kdm16_sc_default --partition=sla-prio  --out slurm/job.ZDNA.centromeric_windows.$chr.%j.out
  done
# Merge
echo -e "Chr\tNonB\tEnrichment_GW\tEnrichment_CHR" >centromeres/${prefix}.background.ZDNA.enrichment.tsv
cat centromeres/windows/$prefix.chr*.100rand.ZDNA.enrichment.tsv  >>centromeres/${prefix}.background.ZDNA.enrichment.tsv



#echo "Window Chr GCcont" |sed 's/ /\t/g'  >centromeres/${prefix}.background.GC.tsv
#for file in $(ls centromeres/windows/chr*.exclCen.with_GC.bed)
#do
#  awk '(NR>1){rownum=NR-1; print rownum,$1,$5}' $file |sed 's/ /\t/g' >>centromeres/${prefix}.background.GC.tsv
#done

# Make an ordered chromosome file to be used for plotting
# ZEBRA FINCH 
echo "Chr Length" |sed 's/ /\t/g' >centromeres/bTaeGut7v0.4_MT_rDNA.chr_len.ordered.txt
for i in '1' '1A' '2' '3' '4' '4A' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' '23' '24' '25' '26' '27' '28' '29' '30' '31' '32' '33' '34' '35' '36' '37'
do
  for hap in "mat" "pat"
  do
    grep "chr${i}_$hap" ref/bTaeGut7v0.4_MT_rDNA.fa.fai |cut -f1,2 >>centromeres/bTaeGut7v0.4_MT_rDNA.chr_len.ordered.txt
  done
done
grep "chrZ" ref/bTaeGut7v0.4_MT_rDNA.fa.fai |cut -f1,2 >>centromeres/bTaeGut7v0.4_MT_rDNA.chr_len.ordered.txt
grep "chrW" ref/bTaeGut7v0.4_MT_rDNA.fa.fai |cut -f1,2 >>centromeres/bTaeGut7v0.4_MT_rDNA.chr_len.ordered.txt
# CHICKEN
echo "Chr Length" |sed 's/ /\t/g'  >centromeres/chicken.v23.chr_len.ordered.txt
for i in '1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' '23' '24' '25' '26' '27' '28' '29' '30' '31' '32' '34' '35' '36'
do
  awk -v c="chr${i}" '($1==c){print $1"\t"$2}' ref/chicken.v23.fa.fai >>centromeres/chicken.v23.chr_len.ordered.txt
done
grep "chrZ" ref/chicken.v23.fa.fai |cut -f1,2 >>centromeres/chicken.v23.chr_len.ordered.txt
grep "chrW" ref/chicken.v23.fa.fai |cut -f1,2 >>centromeres/chicken.v23.chr_len.ordered.txt

