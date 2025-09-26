################################################################################
### CALCULATING NON-B DNA MOTIF ENRICHMENT IN ZEBRA FINCH GENES, AND USING
### METHYLATION DATA TO ESTIMATE G4 FORMATION
### CODE WRITTEN BY LINNÉA SMEDS

##### TABLE OF CONTENTS
# PROCESS THE ANNOTATION FILE INTO GENE CATEGORIES
# CALCULATE ENRICHMENT
#   -CALCULATE ENRICHMENT FOR CHROMOSOME TYPES
#   -CALCULATE ENRICHMENT USING A CHROMOSOME CATEGORY-BASED AVG 
# SEPARATE G4S BASED ON STRANDNESS
# SUBSAMBLE REGIONS TO GET ERROR BARS
# COMBINATION OF INTRONS AND REPEATS
# METHYLATION

############ PROCESS THE ANNOTATION FILE INTO GENE CATEGORIES  #################
# 29 July 2025
mkdir -p functional
mkdir -p functional/annotation
prefix="bTaeGut7v0.4_MT_rDNA"
# Get the transcript-ids of all the longest transcripts (with exons):
awk -v FS=";" '{print $NF}' ref/bTaeGut7v0.4.v0.1.annotation.modified.longest_isoform.only_exon.gtf |sed 's/transcript_id=//g' |sort |uniq >functional/annotation/longest_transcript.txt

# Get the longest transcript CDSs:
awk '{print $1"\";"}' functional/annotation/longest_transcript.txt |grep -f - <(zcat ref/bTaeGut7v0.4_MT_rDNA.EGAPx.v0.1.gtf.gz) |awk '($3=="CDS"){s=$4-1; print $1,s,$5,$12,".",$7}' |sort -k1,1 -k2,2n |sed 's/ /\t/g' >annotation/$prefix.CDS.bed

# And the exons from those CDSs...
cut -f4 functional/annotation/CDS_longest_transcript.bed |grep -f - <(zcat ref/bTaeGut7v0.4_MT_rDNA.EGAPx.v0.1.gtf.gz) |awk '($3=="exon"){s=$4-1; print $1,s,$5,$12,".",$7}' |sort -k1,1 -k2,2n |sed 's/ /\t/g' >functional/annotation/exon_longest_transcript_with_CDS.bed

# ... to extract UTRs:
module load bedtools/2.31.0
subtractBed -s -a functional/annotation/exon_longest_transcript_with_CDS.bed -b functional/annotation/CDS_longest_transcript.bed >functional/annotation/UTR_longest_transcript.bed

# Get start codons:
cut -f4 functional/annotation/UTR_longest_transcript.bed |grep -f - <(zcat ref/bTaeGut7v0.4_MT_rDNA.EGAPx.v0.1.gtf.gz) | awk -v OFS='\t' '($3=="start_codon"){s=$4-1; print $1,s,$5,$12,".",$7}' |sort -k1,1 -k2,2n >functional/annotation/start_from_UTR_longest.bed
# and stop codons:
cut -f4 functional/annotation/UTR_longest_transcript.bed |grep -f - <(zcat ref/bTaeGut7v0.4_MT_rDNA.EGAPx.v0.1.gtf.gz) | awk -v OFS='\t' '($3=="stop_codon"){s=$4-1; print $1,s,$5,$12,".",$7}' |sort -k1,1 -k2,2n >functional/annotation/stop_from_UTR_longest.bed

# Based on this we can extract 5' UTRs
join -1 4 -2 4 <(sort -k4,4 functional/annotation/UTR_longest_transcript.bed) <(sort -k4,4 functional/annotation/start_from_UTR_longest.bed) |awk -v OFS="\t" '($2==$7){if($6=="+"){if($8>=$4){print $2,$3,$4,$1,$5,$6}}else if($6=="-"){if($9<=$3){print $2,$3,$4,$1,$5,$6}}}' |sort -k1,1 -k2,2n >annotation/$prefix.UTR5.bed
# And 3' UTR
join -1 4 -2 4 <(sort -k4,4 functional/annotation/UTR_longest_transcript.bed) <(sort -k4,4 functional/annotation/stop_from_UTR_longest.bed) |awk -v OFS="\t" '($2==$7){if($6=="+"){if($3>=$8){print $2,$3,$4,$1,$5,$6}}else if($6=="-"){if($4<=$9){print $2,$3,$4,$1,$5,$6}}}' |sort -k1,1 -k2,2n >annotation/$prefix.UTR3.bed


# Introns (remove any overlapping categories above, and only take introns from
# protein coding genes)
module load bedtools/2.31.0
awk '{print $1"\";"}' functional/annotation/longest_transcript.txt |grep -f - <(zcat ref/bTaeGut7v0.4_MT_rDNA.EGAPx.v0.1.gtf.gz) |awk '($3=="exon" && $0~/transcript_biotype "mRNA";/){s=$4-1; print $1,s,$5,$12,".",$7}' |sort -k4,4 -k2,2n | awk -v OFS="\t" '{if(NR==1){chr=$1; end=$3; name=$4; strand=$6}else{if($1==chr && $4==name){print chr,end,$2,name,".",strand; end=$3}else{chr=$1; end=$3; name=$4;strand=$6}}}' |sort -k1,1 -k2,2n | intersectBed -v -a - -b functional/annotation/exon_longest_transcript_with_CDS.bed >annotation/$prefix.intronic.bed

# Promoters (only protein coding genes)
awk '{print $1"\";"}' functional/annotation/longest_transcript.txt |grep -f - <(zcat ref/bTaeGut7v0.4_MT_rDNA.EGAPx.v0.1.gtf.gz) |awk '($3=="transcript" && $0~/transcript_biotype "mRNA";/){if($7=="+"){e=$4-1; s=e-1000}else if($7=="-"){s=$5; e=s+1000}else{}; print $1,s,e,$12,".",$7}' |sort -k1,1 -k2,2n |sed 's/ /\t/g' >annotation/$prefix.promoter.bed

# Intergenic - take everything in between genes and remove overlapping exons and
# introns and promoters
awk -v OFS='\t' '{s=$4-1; print $1,s,$5}' ref/bTaeGut7v0.4.v0.1.annotation.modified.longest_isoform.only_gene.gtf |sort -k1,1 -k2,2n |mergeBed -i - | awk -v OFS='\t' '{if(NR==1){chr=$1; end=$3}else{if($1==chr){print chr,end,$2; end=$3}else{chr=$1; end=$3}}}' |intersectBed -v -a - -b functional/annotation/mat_exon_longest_transcript_with_CDS.bed  | intersectBed -v -a - -b functional/annotation/exon_longest_transcript_with_CDS.bed | intersectBed -v -a - -b functional/annotation/intron_longest_transcript.bed |intersectBed -v -a - -b functional/annotation/promoter_longest_transcript.bed  >annotation/$prefix.intergenic.bed

# Non coding genes
grep -v "protein_coding" ref/bTaeGut7v0.4.v0.1.annotation.modified.longest_isoform.only_gene.gtf |grep -v "pseudogene" | sed 's/;/\t/g' |sed 's/gene_id=//g' |awk -v OFS='\t' '{s=$4-1; print $1,s,$5,$12,".",$7}' |sort -k1,1 -k2,2n > annotation/$prefix.nonprotcoding.bed


############################# CALCULATE ENRICHMENT #############################

prefix="bTaeGut7v0.4_MT_rDNA"
for class in  "promoter" "intergenic" "intronic" "CDS" "UTR5" "UTR3" "nonprotcoding"
do
    echo '#!/bin/bash
    module load bedtools/2.31.0
    len=`sort -k1,1 -k2,2n annotation/'$prefix'.'$class'.bed | mergeBed -i - |awk '"'"'{sum+=$3-$2}END{print sum}'"'"'`
    echo "Length of '$class' is $len"
    rm -f tmp.'$class'
    cat densities/'${prefix}'.nonB_genome_wide.txt | while read -r non_b tot dens;
    do
        echo "looking at $non_b"
        d=`intersectBed -a <(cut -f1-3 annotation/'$prefix'.'$class'.bed |sort -k1,1 -k2,2n |mergeBed -i -) -b <(sort -k1,1 -k2,2n final_nonB/'${prefix}'.${non_b}.bed |mergeBed -i -) -wo |awk -v l=$len -v dtot=$dens '"'"'{sum+=$7}END{d=sum/l; frac=d/dtot; print d,frac}'"'"'`
        echo '$class'" "$non_b" "$d >>tmp.'$class'
    done
    ' | sbatch -J $class --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8G --time=5:00:00 --out slurm/job.functional.$class.%j.out
done
# Merge after instead
echo "Class nonB Density Enrichment_gw" |sed "s/ /\t/g" >functional/$prefix.enrichment_fullgenome.tsv
for class in "promoter" "intergenic" "intronic" "CDS" "UTR5" "UTR3"  "nonprotcoding" 
do
    cat tmp.$class |sed "s/ /\t/g" >>functional/$prefix.enrichment_fullgenome.tsv
done
# Remove temp files 
rm tmp.*

# ~~~~~~~~~~~~~~~~~ CALCULATE ENRICHMENT FOR CHROMOSOME TYPES ~~~~~~~~~~~~~~~~~~
prefix="bTaeGut7v0.4_MT_rDNA"
for group in "macro" "micro" "microdot"
do
  for class in "promoter" "intergenic" "intronic" "CDS" "UTR5" "UTR3" "nonprotcoding" 
  do
      echo '#!/bin/bash
      module load bedtools/2.31.0
      len=`grep -f helpfiles/'$group'.txt annotation/'$prefix'.'$class'.bed |mergeBed -i - |awk '"'"'{sum+=$3-$2}END{print sum}'"'"'`
      echo "Length of '$class' is $len"
      rm -f tmp.'$group'.'$class'
      cat densities/'${prefix}'.nonB_genome_wide.txt | while read -r non_b tot dens;
      do
          echo "looking at $non_b"
          d=`intersectBed -a <( grep -f helpfiles/'$group'.txt annotation/'$prefix'.'$class'.bed |cut -f1-3 |mergeBed -i -) -b final_nonB/'${prefix}'.${non_b}.merged.bed -wo |awk -v l=$len -v dtot=$dens '"'"'{sum+=$7}END{d=sum/l; frac=d/dtot; print d,frac}'"'"'`
          echo '$group'" "'$class'" "$non_b" "$d >>tmp.'$group'.'$class'
      done
      ' | sbatch -J $class --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8G --time=5:00:00 --out slurm/job.functional.$group.$class.%j.out
  done
done

# Merge the tmp files
echo "Group Class nonB Density Enrichment_gw" |sed "s/ /\t/g" >functional/$prefix.enrichment_groups.tsv
for group in "macro" "micro" "microdot"
do
  for class in  "promoter"  "intergenic"  "intronic" "CDS" "UTR5" "UTR3" "nonprotcoding"
  do
    cat tmp.$group.$class |sed "s/ /\t/g" >>functional/$prefix.enrichment_groups.tsv
  done
done
rm tmp.*

# ~~~~~~~~ CALCULATE ENRICHMENT USING A CHROMOSOME CATEGORY-BASED AVG ~~~~~~~~~~
# Added September, this was not included in the paper as it doesn't make much
# sense for the dotchromosomes, where almost everything is coding. Another idea 
# would be to divide by the intergenic density. 
prefix="bTaeGut7v0.4_MT_rDNA"
for group in "macro" "micro" "microdot"
do
  for class in "promoter" "intergenic" "intronic" "CDS" "UTR5" "UTR3" "nonprotcoding"
  do
      echo '#!/bin/bash
      module load bedtools/2.31.0
      len=`grep -f helpfiles/'$group'.txt annotation/'$prefix'.'$class'.bed |mergeBed -i - |awk '"'"'{sum+=$3-$2}END{print sum}'"'"'`
      echo "Length of '$class' is $len"
      rm -f tmp.groupavg.'$group'.'$class'
      cat densities/'${prefix}'.nonB_per_group.txt |awk -v g='$group' '"'"'(g==$1){print}'"'"' | while read -r g non_b tot dens;
      do
          echo "looking at $non_b"
          d=`intersectBed -a <( grep -f helpfiles/'$group'.txt annotation/'$prefix'.'$class'.bed |cut -f1-3 |mergeBed -i -) -b final_nonB/'${prefix}'.${non_b}.merged.bed -wo |awk -v l=$len -v dtot=$dens '"'"'{sum+=$7}END{d=sum/l; frac=d/dtot; print d,frac}'"'"'`
          echo '$group'" "'$class'" "$non_b" "$d >>tmp.groupavg.'$group'.'$class'
      done
      ' | sbatch -J $class --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8G --time=5:00:00 --out slurm/job.functional.groupavg.$group.$class.%j.out
  done
done

# Merge the tmp files
echo "Group Class nonB Density Enrichment_gw" |sed "s/ /\t/g" >functional/$prefix.enrichment_groups.groupavg.tsv
for group in "macro" "micro" "microdot"
do
  for class in  "promoter"  "intergenic"  "intronic" "CDS" "UTR5" "UTR3" "nonprotcoding"
  do
    cat tmp.groupavg.$group.$class |sed "s/ /\t/g" >>functional/$prefix.enrichment_groups.groupavg.tsv
  done
done
rm tmp.groupavg.*


####################### SEPARATE G4s BASED ON STRANDNESS #######################
# We want to check if G4s occur more often on the coding or template strand

prefix="bTaeGut7v0.4_MT_rDNA"
for group in "macro" "micro" "microdot"
do
  for class in "promoter" "intronic" "CDS" "UTR5" "UTR3" "nonprotcoding"
  do
      echo '#!/bin/bash
      module load bedtools/2.31.0
      len=`grep -f helpfiles/'$group'.txt annotation/'$prefix'.'$class'.bed |mergeBed -i - |awk '"'"'{sum+=$3-$2}END{print sum}'"'"'`
      echo "Length of '$class' is $len"
      rm -f tmp.strand.'$group'.'$class'
      dsame=`intersectBed -a <(grep -f helpfiles/'$group'.txt annotation/'$prefix'.'$class'.bed) -b final_nonB/'${prefix}'.G4.bed -s |sort -k1,1 -k2,2n |mergeBed -i - |awk -v l=$len '"'"'{sum+=$3-$2}END{d=sum/l; print d}'"'"'`
      doppo=`intersectBed -a <(grep -f helpfiles/'$group'.txt annotation/'$prefix'.'$class'.bed) -b final_nonB/'${prefix}'.G4.bed -S |sort -k1,1 -k2,2n |mergeBed -i - |awk -v l=$len '"'"'{sum+=$3-$2}END{d=sum/l; print d}'"'"'`
      echo '$group'" "'$class'" Coding "$dsame >>tmp.strand.'$group'.'$class'
      echo '$group'" "'$class'" Template "$doppo >>tmp.strand.'$group'.'$class'
      ' | sbatch -J $class --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8G --time=5:00:00 --partition=open --out slurm/functional.$group.$class.%j.out
  done
done
#Merge
echo "Group Class Strand Density" |sed "s/ /\t/g" >functional/$prefix.G4_strand_density_groups.tsv
for group in "macro" "micro" "microdot"
do
  for class in  "promoter"  "intronic" "CDS" "UTR5" "UTR3" "nonprotcoding"
  do
    cat tmp.strand.$group.$class |sed "s/ /\t/g" >>functional/$prefix.G4_strand_density_groups.tsv
  done
done
rm tmp.strand.*

# To get distributions for Coding and Template mean bars, I would like to get an
# enrichment value for each gene to get a distribution
prefix="bTaeGut7v0.4_MT_rDNA"
for group in "macro" "micro" "microdot"
do
  for class in "promoter" "intronic" "CDS" "UTR5" "UTR3" "nonprotcoding"
  do
      echo '#!/bin/bash
      module load bedtools/2.31.0
      rm -f tmp.strand.perGene.'$group'.'$class'
      # SAME STRAND 
      intersectBed -a <(grep -f helpfiles/'$group'.txt annotation/'$prefix'.'$class'.bed) -b final_nonB/'${prefix}'.G4.bed -s -wao |\
      awk -v g='$group' -v c='$class' -v OFS="\t" '"'"'{if(NR==1){tr=$4; s=$2; tot=$3-$2; g4=$NF}
    else{if(tr==$4){g4+=$NF; if(s!=$2){tot+=$3-$2}; s=$2}
    else{d=g4/tot; print g,c,"Coding",tr,tot,g4,d; tr=$4; s=$2; tot=$3-$2; g4=$NF}
    }}END{d=g4/tot; print g,c,"Coding",tr,tot,g4,d}'"'"' >>tmp.strand.perGene.'$group'.'$class'
    # OPPOSITE STRAND 
      intersectBed -a <(grep -f helpfiles/'$group'.txt annotation/'$prefix'.'$class'.bed) -b final_nonB/'${prefix}'.G4.bed -S -wao |\
      awk -v g='$group' -v c='$class' -v OFS="\t" '"'"'{if(NR==1){tr=$4; s=$2; tot=$3-$2; g4=$NF}
    else{if(tr==$4){g4+=$NF; if(s!=$2){tot+=$3-$2}; s=$2}
    else{d=g4/tot; print g,c,"Template",tr,tot,g4,d; tr=$4; s=$2; tot=$3-$2; g4=$NF}
    }}END{d=g4/tot; print g,c,"Template",tr,tot,g4,d}'"'"' >>tmp.strand.perGene.'$group'.'$class'
      ' | sbatch -J $class --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8G --time=5:00:00 --out slurm/job.functional.perGene.$group.$class.%j.out
  done
done
# Merge
echo "Group Class Strand Gene GenLen G4Len Density" |sed "s/ /\t/g" >functional/$prefix.G4_strand_density_perGene.tsv
for group in "macro" "micro" "microdot"
do
  for class in  "promoter"  "intronic" "CDS" "UTR5" "UTR3" "nonprotcoding"
  do
    cat tmp.strand.perGene.$group.$class |sed 's/"//g' |sed 's/;//g' >>functional/$prefix.G4_strand_density_perGene.tsv
  done
done
rm tmp.strand.perGene.*


# #################### SUBSAMBLE REGIONS TO GET ERROR BARS #####################
# CHECK SIGNIFICANCE OF THE BARS IN FIGURE 2A

# Subsample the functional regions to get a distribution of the enrichment
# values. If distribution overlaps with 1 (genome average), we deem this
# comparison NOT to be significally different from the genome average.

# Subsample 50 percent of the regions 100 times
prefix="bTaeGut7v0.4_MT_rDNA"
mkdir -p functional/resample/
perc=50
for group in "macro" "micro" "microdot"
do
  for class in "promoter" "intergenic" "intronic" "CDS" "UTR5" "UTR3" "nonprotcoding"
  do
    num=`grep -f helpfiles/$group.txt annotation/'$prefix'.$class.bed |awk -v p=$perc '{}END{num=int(NR*(p/100)); print num}' `
    echo "we will extract $num rows from $class on $group"
    # Subsample the sequences 100 times
    echo '#!/bin/bash
    module load bedtools/2.31.0
    for i in {1..100}
    do
        echo "working on $i..."
        grep -f helpfiles/'$group'.txt annotation/'$prefix'.'$class'.bed |shuf -n '$num' |sort -k1,1 -k2,2n |mergeBed -i -  >functional/resample/Resamp.'$perc'perc.'$class'.'$group'.$i.bed
        l=`wc -l functional/resample/Resamp.'$perc'perc.'$class'.'$group'.$i.bed`
        echo "...extracting $l lines"
    done
    ' |sbatch -J $class.$group  --ntasks=1 --mem-per-cpu=1G --cpus-per-task=1 --time=1:00:00  -o slurm/job.resamp.${perc}perc.$class.$group.%j.out
  done
done

# Calculate the enrichment for each of the subsamples
prefix="bTaeGut7v0.4_MT_rDNA"
for subset in "50perc" 
do
  for group in "macro" "micro" "microdot"
  do
    for class in "promoter" "intergenic" "intronic" "CDS" "UTR5" "UTR3" "nonprotcoding"
    do
        echo '#!/bin/bash
        module load bedtools/2.31.0
        for i in {1..100}
        do
            rm -f functional/resample/tmp.'$subset'.'$class'.'$group'.$i.txt
            len=`awk '"'"'{sum+=$3-$2}END{print sum}'"'"' functional/resample/Resamp.'$subset'.'$class'.'$group'.$i.bed`
            cat densities/'${prefix}'.nonB_genome_wide.txt |grep -v "ALL" | while read -r non_b tot dens;
            do
                d=`intersectBed -a functional/resample/Resamp.'$subset'.'$class'.'$group'.$i.bed -b final_nonB/'$prefix'.${non_b}.merged.bed -wo |awk -v l=$len -v dtot=$dens '"'"'{sum+=$7}END{d=sum/l; frac=d/dtot; print frac}'"'"'`
                echo $i '$group' '$class' $non_b $d >>functional/resample/tmp.'$subset'.'$class'.'$group'.$i.txt
            done
        done
        ' |sbatch -J $class.$group -o slurm/job.resamp-enrich.$subset.$class.$group.%j.out --ntasks=1 --cpus-per-task=1 --mem-per-cpu=4G --time=10:00:00
    done
  done
done

# Merge the results
rep="100rep"
subset="50perc"
echo "Rep Group Class nonB Enrichment" |sed "s/ /\t/g" >functional/$prefix.summary.$subset.$rep.txt
cat functional/resample/tmp.$subset.*.*.{1..100}.txt |sed "s/ /\t/g" >>functional/$prefix.summary.$subset.$rep.txt

# We want to find the min and max for each category, but remove the ~5% most
# extreme (meaning we remove the top and bottom two values, saving a 96% 'CI')
rep="100rep"
subset="50perc" 
echo "Group Class nonB Min Max" |sed "s/ /\t/g" >functional/$prefix.minmax.$subset.$rep.96CI.txt
for group in "macro" "micro" "microdot"
do
  for class in "promoter" "intergenic" "intronic" "CDS" "UTR5" "UTR3" "nonprotcoding"
  do
      for nonb in APR DR G4 IR MR TRI STR Z
      do
          min=`grep $class functional/summary.$subset.$rep.txt |awk -v nb=$nonb -v g=$group -v OFS="\t" '($4==nb && $2==g){print $5}' |sort -n |head -n3 |tail -n1`
          max=`grep $class functional/summary.$subset.$rep.txt |awk -v nb=$nonb -v g=$group -v OFS="\t" '($4==nb && $2==g){print $5}' |sort -n |tail -n3 |head -n1`
          echo $group $class $nonb $min $max |sed "s/ /\t/g" >>functional/$prefix.minmax.$subset.$rep.96CI.txt
      done
    done
done


# ##################### COMBINATION OF INTRONS AND REPEATS #####################
# ADDED SEPT 10

# We want to see if it is the mini satellites in introns that are driving the
# high non-B density here.

# Overlap introns with the TR repeat file
module load bedtools/2.31.0
mkdir repeats/intronTRF
annotation/TRF_withMers.bed
prefix="bTaeGut7v0.4_MT_rDNA"
intersectBed -a annotation/TRF_withMers.bed -b annotation/$prefix.intronic.bed >repeats/intronTRF/introns_TRF.bed
subtractBed -a annotation/$prefix.intronic.bed -b repeats/intronTRF/introns_TRF.bed >repeats/intronTRF/introns_noTRF.bed

# For each chromosome category, check enrichment in these two
for group in "macro" "micro" "microdot"
do
    echo '#!/bin/bash
    module load bedtools/2.31.0
    rm -f repeats/intronTRF/enrichment.'$group'.txt
    for i in "TRF" "noTRF"
    do
        len=`grep -f helpfiles/'$group'.txt repeats/intronTRF/introns_${i}.bed |sort -k1,1 -k2,2n |mergeBed -i - |awk '"'"'{sum+=$3-$2}END{print sum}'"'"'`
        echo "length of $i: $len"
        cat densities/'${prefix}'.nonB_genome_wide.txt | while read -r non_b tot dens;
        do
            d=`intersectBed -a <(grep -f helpfiles/'$group'.txt repeats/intronTRF/introns_${i}.bed |sort -k1,1 -k2,2n |mergeBed -i -) -b final_nonB/'$prefix'.${non_b}.merged.bed -wo |awk -v l=$len -v dtot=$dens '"'"'{sum+=$7}END{d=sum/l; frac=d/dtot; print d, frac}'"'"'`
            echo '$group' $i $non_b $d >>repeats/intronTRF/enrichment.'$group'.txt
        done
    done
    ' |sbatch -J $group -o slurm/job.intron-enrich.$group.%j.out --requeue --ntasks=1 --cpus-per-task=1 --mem-per-cpu=4G --time=10:00:00 --partition=open
done
# Merge
echo "Group Subset NonB Density Enrichment_gw" |sed "s/ /\t/g" >repeats/$prefix.group_intron_enrichment.tsv
for group in "macro" "micro" "microdot"
do
  cat repeats/intronTRF/enrichment.$group.txt |sed "s/ /\t/g"
done >>repeats/$prefix.group_intron_enrichment.tsv


# ALSO MAKE FILES WITH THE NUMBER OF TRFs / noTRFs PER CATEGORY AND COMPARTMENT
echo "Group Subset Length" |sed "s/ /\t/g" >repeats/intron_TR_lengthsummary.tsv
for group in "macro" "micro" "microdot"
do
  grep "length of" slurm/job.intron-enrich.$group.*.out |tail -n2 |sed 's/://g' |awk -v OFS="\t" -v g=$group '{print g,$3,$4}' >>repeats/intron_TR_lengthsummary.tsv
done
for comp in "A" "B"
do
  grep "Length of" slurm/job.intron-enrich.$comp.*.out |tail -n2 |awk -v OFS="\t" '{print $5,$3,$8}' >>repeats/intron_TR_lengthsummary.tsv
done

# AND A FILE WITH TRF CLASSES
echo "Group Repeat Number" |sed "s/ /\t/g" >repeats/intron_TR_repeatsummary.tsv
for group in "macro" "micro" "microdot"
do
  grep -f $group.txt repeats/intronTRF/introns_TRF.bed |cut -f4 |sort |uniq -c |awk -v OFS="\t" -v g=$group '{print g,$2,$1}' >>repeats/intron_TR_repeatsummary.tsv
done

# And exact length of repeats
cut -f1,5 repeats/intronTRF/introns_TRF.bed |sed 's/-mer//g' >repeats/introns_TRF.lengths.txt


################################# METHYLATION ##################################
# July 30, 2025

mkdir -p methylation/functional

# Convert methylation bigwig file to bedgraph (assume UCSC bigWigToWig tool is
# installed in ~/software/)
~/software/bigWigToWig ref/bTaeGut7v0.4_MT_rDNA.PBmethylation.v0.1.bw ref/bTaeGut7v0.4_MT_rDNA.PBmethylation.v0.1.wig
# The Wig file is actually in BedGraph format, only need to remove the # lines
cat ref/bTaeGut7v0.4_MT_rDNA.PBmethylation.v0.1.wig |awk '($1!~/^#/){print}' >annotation/bTaeGut7v0.4_MT_rDNA.PBmethylation.v0.1.bed
# The file contains 22M methylated sites.

# Because the size of the file, I divide it into 3 parts before overlapping:
for group in "macro" "micro" "microdot"
do
  grep -f helpfiles/$group.txt annotation/$prefix.PBmethylation.v0.1.bed >methylation/$group.methylation.bed
done

# First I want to extract all G4s within functional regions, and the functional
# regions with G4s removed!
mkdir -p functional/annotation/G4s
module load bedtools/2.31.0
for class in "intronic" "nonprotcoding" "CDS" "UTR5" "UTR3" "promoter" "intergenic"
do
  subtractBed -a annotation/$prefix.$class.bed -b final_nonB/$prefix.G4.bed >functional/annotation/G4s/$class.exclG4s.bed
  intersectBed -a annotation/$prefix.$class.bed -b final_nonB/$prefix.G4.bed >functional/annotation/G4s/$class.onlyG4s.bed
  intersectBed -a annotation/$prefix.$class.bed -b final_nonB/$prefix.G4.bed -s >functional/annotation/G4s/$class.onlyG4s.samestrand.bed
  intersectBed -a annotation/$prefix.$class.bed -b final_nonB/$prefix.G4.bed -S >functional/annotation/G4s/$class.onlyG4s.oppositestrand.bed
done

# Make a list with all methylated sites within those regions!
mkdir -p methylation/functional/G4s/
prefix="bTaeGut7v0.4_MT_rDNA"
for group in "microdot" "macro" "micro"
do
  for class in "intronic" "nonprotcoding" "CDS" "UTR5" "UTR3" "promoter"
  do
      echo '#!/bin/bash
      module load bedtools/2.31.0
      intersectBed -a <(grep -f helpfiles/'$group'.txt functional/annotation/G4s/'$class'.exclG4s.bed) -b methylation/'$group'.methylation.bed -wao |awk -v OFS="\t" -v cl='$class' '"'"'($NF==1){print cl,"background",$10,$4}'"'"' >methylation/functional/G4s/'$group'.'$class'.txt
      intersectBed -a <(grep -f helpfiles/'$group'.txt functional/annotation/G4s/'$class'.onlyG4s.samestrand.bed) -b methylation/'$group'.methylation.bed -wao |awk -v OFS="\t" -v cl='$class' '"'"'($NF==1){print cl,"coding",$10,$4}'"'"' >>methylation/functional/G4s/'$group'.'$class'.txt
      intersectBed -a <(grep -f helpfiles/'$group'.txt functional/annotation/G4s/'$class'.onlyG4s.oppositestrand.bed) -b methylation/'$group'.methylation.bed -wao |awk -v OFS="\t" -v cl='$class' '"'"'($NF==1){print cl,"template",$10,$4}'"'"' >>methylation/functional/G4s/'$group'.'$class'.txt
    ' | sbatch -J $class.$group --ntasks=1 --cpus-per-task=1 --time=2:00:00 --out slurm/job.methylation.$group.$class.%j.out
  done
  # Run intergenic separately since it is not strand specific
  for class in "intergenic"
  do
      echo '#!/bin/bash
      module load bedtools/2.31.0
      intersectBed -a <(grep -f helpfiles/'$group'.txt functional/annotation/G4s/'$class'.exclG4s.bed) -b methylation/'$group'.methylation.bed -wao |awk -v OFS="\t" -v cl='$class' '"'"'($NF==1){print cl,"background",$(NF-1)}'"'"' >methylation/functional/G4s/'$group'.'$class'.txt
      intersectBed -a <(grep -f helpfiles/'$group'.txt functional/annotation/G4s/'$class'.onlyG4s.bed) -b methylation/'$group'.methylation.bed -wao |awk -v OFS="\t" -v cl='$class' '"'"'($NF==1){print cl,"ignorant",$(NF-1)}'"'"' >>methylation/functional/G4s/'$group'.'$class'.txt
      ' | sbatch -J $class.$group --ntasks=1 --cpus-per-task=1 --mem-per-cpu=4G --time=1:00:00 --out slurm/methylation.$group.$class.%j.out
  done
done

# Merge everything except intergenic
echo "Group Class Type Score Trx" |sed 's/ /\t/g' >methylation/functional/Group.allCpG.merged.txt
for group in "microdot" "macro" "micro"
do
  for class in "intronic" "CDS" "UTR5" "UTR3" "promoter" "nonprotcoding"
  do
    awk -v g=$group -v OFS="\t" '{print g,$0}' methylation/functional/G4s/$group.$class.txt |sed 's/;//g' |sed 's/"//g' >>methylation/functional/Group.allCpG.merged.txt
  done
done
# Merge intergenic
echo "Group Class Type Score" |sed 's/ /\t/g' >methylation/functional/Group.allCpG.intergenic.txt
for group in "microdot" "macro" "micro"
do
  awk -v g=$group -v OFS="\t" '{print g,$0}' methylation/functional/G4s/$group.intergenic.txt |sed 's/;//g' |sed 's/"//g' >>methylation/functional/Group.allCpG.intergenic.txt
done

# Figure out the proportion of genes that has/doesn't have methylation
prefix="bTaeGut7v0.4_MT_rDNA"
for group in "microdot" "macro" "micro"
do
  for class in "intronic" "nonprotcoding" "CDS" "UTR5" "UTR3" "promoter"
  do
      echo '#!/bin/bash
      module load bedtools/2.31.0
      intersectBed -a <(grep -f helpfiles/'$group'.txt functional/annotation/G4s/'$class'.exclG4s.bed) -b methylation/'$group'.methylation.bed -wao |awk -v OFS="\t" -v cl='$class' -v g='$group' '"'"'{if(NR==1){trx=$4; sum=$NF}else{if(trx==$4){sum+=$NF}else{print gr,cl,"background",$4,sum; trx=$4; sum=$NF}}}END{print gr,cl,"background",$4,sum;}'"'"' >methylation/functional/G4s/'$group'.'$class'.summary.txt
      intersectBed -a <(grep -f helpfiles/'$group'.txt functional/annotation/G4s/'$class'.onlyG4s.samestrand.bed) -b methylation/'$group'.methylation.bed -wao |awk -v OFS="\t" -v cl='$class' -v g='$group' '"'"'{if(NR==1){trx=$4; sum=$NF}else{if(trx==$4){sum+=$NF}else{print gr,cl,"coding",$4,sum; trx=$4; sum=$NF}}}END{print gr,cl,"coding",$4,sum;}'"'"' >>methylation/functional/G4s/'$group'.'$class'.summary.txt
      intersectBed -a <(grep -f helpfiles/'$group'.txt functional/annotation/G4s/'$class'.onlyG4s.oppositestrand.bed) -b methylation/'$group'.methylation.bed -wao |awk -v OFS="\t" -v cl='$class' -v g='$group' '"'"'{if(NR==1){trx=$4; sum=$NF}else{if(trx==$4){sum+=$NF}else{print gr,cl,"template",$4,sum; trx=$4; sum=$NF}}}END{print gr,cl,"template",$4,sum;}'"'"'  >>methylation/functional/G4s/'$group'.'$class'.summary.txt
    ' | sbatch -J $class.$group --ntasks=1 --cpus-per-task=1 --time=2:00:00 --partition=open --out slurm/methylation.$group.$class.summary.%j.out
  done
done
# Merge summary
echo "Group Class Type Trx Number" |sed 's/ /\t/g' >methylation/functional/Group.summaryCpG.merged.txt
for group in "microdot" "macro" "micro"
do
  for class in "intronic" "CDS" "UTR5" "UTR3" "promoter" "nonprotcoding"
  do
    awk -v g=$group -v OFS="\t" '{print g,$0}' methylation/functional/G4s/$group.$class.summary.txt |sed 's/;//g' |sed 's/"//g' >>methylation/functional/Group.summaryCpG.merged.txt
  done
done

