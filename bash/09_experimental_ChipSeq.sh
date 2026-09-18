################################################################################
### USING PREVIOUSLY MAPPED CHIP-SEQ DATA TO VALIDATE G4 PREDICTIONS IN CHICKEN
### CODE ADAPTED FROM QIAN ET AL. 2024, WITH THE HELP OF CLAUDE SONNET 5. 

##### TABLE OF CONTENTS
# DOWNLOAD DATA
# TRIM DATA WITH TRIMGALORE
# MAP DATA WITH BOWTIE2
# READ ENRICHMENT WITH DEEPTOOLS
#    - Generate profiles and heatmaps from bigWig files
# CALL PEAKS WITH MACS
# CALCULATE PEAK DENSITY OVER THE GENOME 
#   - GENOME WIDE AND REGION DENSITIES
#   - DENSITY ENRICHMENT IN FUNCTIONAL REGIONS
#   - DENSITY ENRICHMENT IN CENTROMERES


################################# DOWNLOAD DATA ################################
# DATA FROM CHICKEN CHIP-SEQ EXPERIMENT (Zheng et al 2020, NAR)

cd DATA
echo '#!/bin/bash
~/software/sratoolkit.3.2.1-ubuntu64/bin/fasterq-dump --split-files SRR9603970
#~/software/sratoolkit.3.2.1-ubuntu64/bin/fasterq-dump --split-files SRR9603969
' |sbatch -J sra  --ntasks=1 --cpus-per-task=1 --mem-per-cpu=2G --time=1-00:00:00 

# Zip files 
for fastq in DATA/*.fastq
do
    ls $fastq 
    echo '#!/bin/bash
    gzip '$fastq'
    ' |sbatch -J zip  --ntasks=1 --cpus-per-task=1 --mem-per-cpu=2G --time=1-00:00:00 
done 


########################### TRIM DATA WITH TRIMGALORE ##########################
module load python anaconda/2023.09
conda create --name mybiotools python=3.13
conda activate mybiotools
conda install -c bioconda trim-galore 

for sra in "SRR9603970" # "SRR9603969" # 
do
    echo '#!/bin/bash
    trim_galore --paired --cores 4 --fastqc DATA/'${sra}'_1.fastq.gz DATA/'${sra}'_2.fastq.gz
    ' |sbatch -J trim  --ntasks=1 --cpus-per-task=4 --mem-per-cpu=4G --time=1-00:00:00 
done

#conda deactivate 


############################# MAP DATA WITH BOWTIE2 ############################
#
mkdir ref/bowtie_index
# Build index 
echo '#!/bin/bash
module load bowtie2/2.5.2
bowtie2-build ref/chicken.v23.fa ref/bowtie_index/chicken.v23
' |sbatch -J index  --ntasks=1 --cpus-per-task=1 --mem-per-cpu=4G --time=1-00:00:00

# Map data 
mkdir logs experimental
mkdir -p experimental/bam experimental/bigwig experimental/matrix experimental/plots experimental/qc
INDEX="ref/bowtie_index/chicken.v23"
THREADS=8

# THIS PIPELINE ALSO FILTER AND REMOVE DUPLICATES
for sra in "SRR9603969"  "SRR9603970"
do
   echo '#!/bin/bash

    module load bowtie2/2.5.4 samtools/1.19.2
    bowtie2 \
    --sensitive-local \
    --no-unal \
    --no-discordant \
    --no-mixed \
    -p '$THREADS' \
    -x '$INDEX' \
    -1 DATA/'${sra}'_1_val_1.fq.gz \
    -2 DATA/'${sra}'_2_val_2.fq.gz \
    2> "logs/'${sra}'.bowtie2.log" \
    | samtools view -b -q 20 - \
    | samtools sort -n -@ '$THREADS' -o experimental/bam/'${sra}'.q20.namesort.bam -
    samtools fixmate -m experimental/bam/'${sra}'.q20.namesort.bam \
      experimental/bam/'${sra}'.q20.fixmate.bam
    samtools sort -@ '$THREADS' -o experimental/bam/'${sra}'.q20.sorted.bam \
      experimental/bam/'${sra}'.q20.fixmate.bam
    samtools markdup -r experimental/bam/'${sra}'.q20.sorted.bam \
      experimental/bam/'${sra}'.q20.markdup.bam
    samtools index experimental/bam/'${sra}'.q20.markdup.bam
' |sbatch -J $sra.bowtie2  --ntasks=1 --cpus-per-task=$THREADS --mem-per-cpu=4G --time=1-00:00:00
done


######################### READ ENRICHMENT WITH DEEPTOOLS #######################
# From the activated conda environment above
conda install bioconda::deeptools

# Testing the plotEnrichment function in a different gene regions: 
SAMPLE=SRR9603969
CONTROL=SRR9603970
for type in "CDS" "UTR5" "intergenic" "introns" "lncrna" "UTR3"
do 
  BED="annotation/chicken.v23.$type.bed"
  echo '#!/bin/bash
  plotEnrichment \
    -b experimental/bam/'${SAMPLE}'.q20.markdup.bam experimental/bam/'${CONTROL}'.q20.markdup.bam \
    --labels "'$SAMPLE'" "'$CONTROL'" \
    --BED '$BED' \
    --plotFile "experimental/qc/'${SAMPLE}'_vs_'${CONTROL}'.plotEnrichment.'$type'.pdf" \
    --outRawCounts "experimental/qc/'${SAMPLE}'_vs_'${CONTROL}'.plotEnrichment.'$type'.tsv"
  '| sbatch -J plotEnrichment --ntasks=1 --cpus-per-task=1 --mem-per-cpu=6G --time=1-00:00:00
done 

# Make normalized bigWig coverage files with deepTools bamCompare 
# Try both ratio and subtract mode
THREADS=8
echo '#!/bin/bash
bamCompare \
  -b1 "experimental/bam/'${SAMPLE}'.q20.markdup.bam" \
  -b2 "experimental/bam/'${CONTROL}'.q20.markdup.bam" \
  --operation ratio \
  --normalizeUsing RPKM --scaleFactorsMethod None \
  -p "'$THREADS'" \
  -o "experimental/bigwig/'${SAMPLE}'_vs_'${CONTROL}'.ratio.RPKM.bw"
' |sbatch -J bamCompare  --ntasks=1 --cpus-per-task=$THREADS --mem-per-cpu=4G --time=1-00:00:00

echo '#!/bin/bash
bamCompare \
  -b1 "experimental/bam/'${SAMPLE}'.q20.markdup.bam" \
  -b2 "experimental/bam/'${CONTROL}'.q20.markdup.bam" \
  --operation subtract \
  --normalizeUsing RPKM  --scaleFactorsMethod None \
  -p "'$THREADS'" \
  -o "experimental/bigwig/'${SAMPLE}'_vs_'${CONTROL}'.subtract.RPKM.bw"
' |sbatch -J bamCompare  --ntasks=1 --cpus-per-task=$THREADS --mem-per-cpu=4G --time=1-00:00:00


# ~~~~~~~~~~~~~ Generate profiles and heatmaps from bigWig files ~~~~~~~~~~~~~~~

# Only used if there are duplicated lines in the bed file
#awk 'BEGIN{OFS="\t"} !seen[$1 FS $2 FS $3]++' "$REGIONS_BED" > "experimental/bed/$(basename "${REGIONS_BED%.bed}").dedup.bed"

# Do this around the transcription start sites, so first create a 
# new bed file for chicken genes 
awk -F'\t' '$3=="gene"' ref/chicken.v23.gff | \
awk 'BEGIN{OFS="\t"} {
  name="NA"
  split($9, attrs, ";")
  for (i in attrs) {
    if (attrs[i] ~ /^ID=/) { split(attrs[i], kv, "="); name=kv[2] }
  }
  print $1, $4-1, $5, name, ".", $7
}' > annotation/chicken.v23.genes.bed
# Remove any duplicates
awk 'BEGIN{OFS="\t"} !seen[$1 FS $2 FS $3]++' annotation/chicken.v23.genes.bed > annotation/chicken.v23.genes.dedup.bed

BED="annotation/chicken.v23.genes.dedup.bed"
THREADS=8
echo '#!/bin/bash
computeMatrix reference-point \
  --referencePoint TSS \
  -S "experimental/bigwig/'${SAMPLE}'_vs_'${CONTROL}'.ratio.RPKM.bw" \
  -R "'$BED'" \
  -a 3000 \
  -b 3000 \
  -p "'$THREADS'" \
  -o "experimental/matrix/'${SAMPLE}'_vs_'${CONTROL}'.TSS.matrix.gz"
' |sbatch -J computeMat  --ntasks=1 --cpus-per-task=$THREADS --mem-per-cpu=4G --time=1-00:00:00

plotProfile \
  -m "experimental/matrix/${SAMPLE}_vs_${CONTROL}.TSS.matrix.gz" \
  -out "experimental/plots/${SAMPLE}_vs_${CONTROL}.profile.png" \
  --yAxisLabel "Sample vs Control ratio" \
  --perGroup \
  --samplesLabel " "

plotHeatmap \
  -m "experimental/matrix/${SAMPLE}_vs_${CONTROL}.TSS.matrix.gz" \
  -out "experimental/plots/${SAMPLE}_vs_${CONTROL}.heatmap.png" \
  --yAxisLabel "Sample vs Control ratio"



############################ CALL PEAKS WITH MACS  #############################

# In the same conda environment as above 
echo '#!/bin/bash
macs3 callpeak \
  -t "experimental/bam/'${SAMPLE}'.q20.markdup.bam" \
  -c "experimental/bam/'${CONTROL}'.q20.markdup.bam" \
  -f BAMPE \
  -g hs \
  -n "'$SAMPLE'" \
  --outdir experimental/peaks \
  -q 0.001 \
  --keep-dup 1
' |sbatch -J macs3  --ntasks=1 --cpus-per-task=1 --mem-per-cpu=6G --time=1-00:00:00



################### CALCULATE PEAK DENSITY OVER THE GENOME  ####################
# Note that we will work with density here and not coverage, since the peaks are
# narrow.

prefix="chicken.v23"

# ~~~~~~~~~~~~~~~~~~~ GENOME WIDE AND REGION DENSITIES ~~~~~~~~~~~~~~~~~~~~~~~~~
mkdir densities
module load bedtools/2.31.0
echo -e "Region\tDensity" >densities/${prefix}.chipSeq.tsv
totlen=`cat ref/${prefix}.fa.fai | awk '{sum+=$2}END{print sum}'`
echo "Totlen: "$totlen
awk -v tot=$totlen -v OFS="\t" 'END{d=NR/tot; print "genome_wide",d}' experimental/peaks/SRR9603969_peaks.narrowPeak >>densities/${prefix}.chipSeq.tsv
for group in "macro" "micro" "dot"
do
  grlen=`grep $group helpfiles/$prefix.groups.txt |cut -f1 |awk 'NR==FNR{a[$1]++;next}a[$1]' - ref/chicken.v23.fa.fai | awk '{sum+=$2}END{print sum}'`
  echo "Group length $group: $grlen"
  grep $group helpfiles/$prefix.groups.txt |cut -f1 |awk 'NR==FNR{a[$1]++;next}a[$1]' - experimental/peaks/SRR9603969_peaks.narrowPeak |\
  awk -v tot=$grlen -v g=$group -v OFS="\t" 'END{d=NR/tot; print g,d}' >>densities/${prefix}.chipSeq.tsv
done

# ~~~~~~~~~~~~~~~ DENSITY ENRICHMENT IN FUNCTIONAL REGIONS ~~~~~~~~~~~~~~~~~~~~~
prefix="chicken.v23"
gw=`grep genome_wide densities/${prefix}.chipSeq.tsv |cut -f2`
echo "Genome wide density: $gw"
for group in "macro" "micro" "dot"
do
  grd=`grep $group densities/${prefix}.chipSeq.tsv |cut -f2`
  echo "Group density is $grd"
  for class in "intergenic" "introns" "promoter" "CDS" "UTR5" "UTR3" 
  do
    echo '#!/bin/bash
    module load bedtools/2.31.0
    len=`awk '"'"'{sum+=$3-$2}END{print sum}'"'"' annotation/group_wise/'$prefix'.'$group'.'$class'.merged.bed`
    echo "Length of '$class' is $len"
    rm -f tmp.chipseq.'$prefix'.'$group'.'$class'
    d=`cut -f1-4 experimental/peaks/SRR9603969_peaks.narrowPeak  | \
    intersectBed -a - -b annotation/group_wise/'$prefix'.'$group'.'$class'.merged.bed -nonamecheck |\
    awk -v l=$len -v dgw='$gw' -v dgrp='$grd'  '"'"'END{d=NR/l; frac_gw=d/dgw; frac_grp=d/dgrp; print d,frac_gw,frac_grp}'"'"'`
          echo '$group'" "'$class'" "$d >>tmp.chipseq.'$prefix'.'$group'.'$class'
      ' | sbatch -J $class --ntasks=1 --cpus-per-task=1 --time=15:00 --out slurm/job.functional.chipseq.$group.$class.%j.out
  done
done

# Merge the tmp files
echo "Group Class Density Enrichment_gw Enrichment_grp" |sed "s/ /\t/g" >experimental/peaks/functional.enrichment.tsv
for group in "macro" "micro" "dot"
do
  for class in  "promoter" "intergenic" "introns" "CDS" "UTR5" "UTR3"
  do
      cat tmp.chipseq.$prefix.$group.$class |sed "s/ /\t/g"  >>experimental/peaks/functional.enrichment.tsv
  done
done


# And across all categories together:
prefix="chicken.v23"
gw=`grep genome_wide densities/${prefix}.chipSeq.tsv |cut -f2`
echo "Genome wide density: $gw"
for class in "intergenic" "introns" "promoter" "CDS" "UTR5" "UTR3" 
do
    echo '#!/bin/bash
    module load bedtools/2.31.0
    len=`awk '"'"'{sum+=$3-$2}END{print sum}'"'"' annotation/'$prefix'.'$class'.merged.bed`
    echo "Length of '$class' is $len"
    rm -f tmp.chipseq.'$prefix'.'$class'
    d=`cut -f1-4 experimental/peaks/SRR9603969_peaks.narrowPeak  | \
    intersectBed -a - -b annotation/'$prefix'.'$class'.merged.bed -nonamecheck |\
    awk -v l=$len -v dgw='$gw'  '"'"'END{d=NR/l; frac_gw=d/dgw; print d,frac_gw,"NA"}'"'"'`
          echo "genome "'$class'" "$d >>tmp.chipseq.'$prefix'.'$class'
      ' | sbatch -J $class --ntasks=1 --cpus-per-task=1 --time=15:00 --out slurm/job.functional.chipseq.$class.%j.out
done
# Merge with above
for class in  "promoter" "intergenic" "introns" "CDS" "UTR5" "UTR3"
do
  cat tmp.chipseq.$prefix.$class |sed "s/ /\t/g"  >>experimental/peaks/functional.enrichment.tsv
done


# ~~~~~~~~~~~~~~~~~~~ DENSITY ENRICHMENT IN CENTROMERES ~~~~~~~~~~~~~~~~~~~~~~~~
# I checked the overlap, but there are almost no peaks in centromeres (and 
# almost no G4 motifs either, according to former fig. 5) so I'll skip this.
intersectBed -a experimental/peaks/SRR9603969_peaks.narrowPeak  -b ref/chicken.v23.CEN.bed -nonamecheck |less




