################################################################################
### DETAILED ANALYSIS OF THE MOST COMMON PREDICTED G4 SEQUENCES AND THEIR 
### OVERLAP TO FUNCTIONAL REGIONS AND REPEATS
### CODE WRITTEN BY LINNÉA SMEDS

##### TABLE OF CONTENTS
# EXTRACT G4 SEQUENCES IN ZEBRA FINCH
# GET SCORES FOR COMMON G4s
# FIND OVERLAP WITH FUNCTIONAL REGIONS
# CHECK FOR PREVIOUS EXPERIMENTAL VALIDATION
# CHECKING G4 SCORES IN DIFFERENT FUNCTIONAL REGIONS


####################### EXTRACT G4 SEQUENCES IN ZEBRA FINCH ####################

# Extract the sequences of the predicted G4 regions 
module load bedtools/2.31.0
prefix="bTaeGut7v0.4_MT_rDNA"
for file in $(ls g4discovery/$prefix/*.bed.gz) 
do 
    chr=$(basename ${file%.bed.gz})
    echo $chr
    zcat $file | bedtools getfasta -fi ref/$prefix.fa -bed - -fo g4discovery/$prefix/$chr.G4.fasta
done

# Reverse complement sequences on the rev strand and print a list with
# position and sequence on one line 
for file in $(ls g4discovery/$prefix/*.fasta) 
do 
    chr=$(basename ${file%.fasta})
    echo $chr
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
        h=h":rev"
        seq=rc
    }else{
        h=h":fwd"
    }
    print h "\t" seq
    }' $file  >g4discovery/$prefix/$chr.revcomp.txt
done


# For each chromosome, sort and count the sequences 
for file in $(ls g4discovery/$prefix/*.G4.revcomp.txt) 
do 
    chr=$(basename ${file%.revcomp.txt})
    echo $chr
    cut -f2 $file | sort | uniq -c | sort -k1,1nr |awk '($1>1){print}' >g4discovery/$prefix/$chr.reocc_counts.txt
done 

# Merge all species and count the ocurrences of G4s in the diploid genome
cat g4discovery/$prefix/*.G4.reocc_counts.txt |awk '{print $2"\t"$1}' |sort -k1,1 |\
awk '{if(NR==1){seq=$1; sum=$2}else{if($1==seq){sum+=$2}else{print seq"\t"sum; seq=$1; sum=$2}}}END{print seq"\t"sum}' |\
sort -k2,2nr >g4discovery/$prefix/common.G4s.txt



########################## GET SCORES FOR COMMON G4s ###########################

# For each common G4, find positions in the genome and chose one to extract the
# scores from the original g4Discovery file (will try if all common G4s are 
# present on the largest chr2)

#chr="chr36_mat"
chr="chr2_mat"
cat g4discovery/$prefix/common.G4s.txt |head -n30 |while read line
do 
    seq=$(echo $line |cut -f1 -d" ")
    #echo "Processing sequence: "$seq
    no=`grep $seq g4discovery/$prefix/$chr.G4.revcomp.txt |wc -l`
    start=`grep $seq g4discovery/$prefix/$chr.G4.revcomp.txt |awk -F':|-|\t' -v seq=$seq '($5==seq){print $2}' - |shuf |head -n1`
    info=`zcat g4discovery/$prefix/$chr.bed.gz | awk -v st=$start '($2==st){print $0}'`
    echo "  "$info"  "$no" "$seq
done 
# All were present, except 4 of the last 5 ones. Seems like they are more 
# common on the dot chromosome, take them from chr36
  chr2_mat 152511581 152511602 64.0 21 + 1.71  672 GGGTTAGGGTTAGGGTTAGGG
  chr2_mat 53696573 53696591 69.0 18 + 2.0  3040 GGGAAGGGAAGGGAAGGG
  chr2_mat 27604804 27604822 69.0 18 + 2.0  197 GGGATGGGATGGGATGGG
  chr2_mat 76420846 76420861 73.0 15 + 2.4  227 GGGAGGGAGGGAGGG
  chr2_mat 97425459 97425478 113.0 19 + 3.37  360 GGGGAGGGGAGGGGAGGGG
  chr2_mat 151256632 151256649 70.0 17 + 2.06  23 GGGAGGGACTGGGAGGG
  chr2_mat 1379655 1379676 64.0 21 + 1.71  141 GGGAAAGGGAAAGGGAAAGGG
  chr2_mat 96472486 96472509 44.0 23 + 1.96  852 GGTGGGGGGCAGGAGGGAAGAGG
  chr2_mat 132738685 132738700 73.0 15 + 4.0  256 GGGGGGGGGGGGGGG
  chr2_mat 674264 674285 64.0 21 + 1.57  41 GGGACAGGGACAGGGACAGGG
  chr2_mat 126331423 126331445 109.0 22 + 2.91  93 GGGGAAGGGGAAGGGGAAGGGG
  chr2_mat 130867195 130867213 69.0 18 + 1.83  30 GGGCAGGGCAGGGCAGGG
  chr2_mat 1135955 1135980 59.0 25 + 1.8  28 GGGACTGGGGACAAGGGATGGAGGG
  chr2_mat 344095 344113 69.0 18 + 1.83  17 GGGCTGGGCTGGGCTGGG
  chr2_mat 152237824 152237847 61.0 23 + 1.87  22 GGGCTGGGGGAGATGGGCCTGGG
  chr2_mat 152513957 152513978 64.0 21 + 1.71  3 GGGTTAGGGTTTGGGTTAGGG
  chr2_mat 152510043 152510064 64.0 21 + 1.71  3 GGGTTAGGGTTAGGGTTTGGG
  chr2_mat 152513717 152513738 64.0 21 + 1.71  2 GGGTTTGGGTTAGGGTTAGGG
  chr2_mat 21847401 21847419 69.0 18 + 2.56  51 GGGCTGGGTGGGGGAGGG
  chr2_mat 7771786 7771807 64.0 21 + 1.71  52 GGGAATGGGAATGGGAATGGG
  chr2_mat 92695562 92695583 64.0 21 + 1.86  54 GGGAGAGGGAGAGGGAGAGGG
  chr2_mat 151385745 151385763 69.0 18 + 2.0  9 GGGTTGGGTTGGGTTGGG
  chr2_mat 102660117 102660135 115.0 18 + 4.0  81 GGGGGGGGGGGGGGGGGG
  chr2_mat 151263552 151263583 96.0 31 + 1.97  127 GGGGACATTGGGGACATTGGGGACATTGGGG
  chr2_mat 2019647 2019672 59.0 25 + 1.56  8 GGGCAGGGATGGATGGGATATTGGG
  chr36_mat 3575464 3575483 49.0 19 + 1.53  56 GGGCTGAGGGCGCGGCGGG
  chr2_mat 587928 587950 109.0 22 + 2.91  9 GGGGATGGGGATGGGGATGGGG
  chr36_mat 2707067 2707089 44.0 22 + 1.68  9 GGGGCAGTGGGCAGAGGCAGGG
  chr36_mat 1203562 1203604 114.0 42 + 2.14  50 GGGGGCTCGGGGCTGGCTGCGGGCGGAGGGGGCGGCAGGGGG
  chr36_mat 1212179 1212206 84.0 27 + 2.3  51 GGGGAATTTGGGGAGGTGGGAGTGGGG



##################### FIND OVERLAP WITH FUNCTIONAL REGIONS #####################
prefix="bTaeGut7v0.4_MT_rDNA"
# First make files with all positions for each common G4s
mkdir -p g4discovery/$prefix/common_g4s/
cat g4discovery/$prefix/common.G4s.txt |head -n30 |while read line
do 
    seq=$(echo $line |cut -f1 -d" ")
    rm -f g4discovery/$prefix/common_g4s/$seq.bed
    for chr in $(cat ref/$prefix.fa.fai |cut -f1)
    do 
        grep $seq g4discovery/$prefix/$chr.G4.revcomp.txt |awk -F':|-|\t' '{if($4=="rev"){strand="-"}else{strand="+"}; print $1"\t"$2"\t"$3"\t"$5"\t.\t"strand}' >>g4discovery/$prefix/common_g4s/$seq.bed
    done 
done

# Overlap with repeat and gene annotations, but first separate the different
# repeat types
awk '($4=="rRNA" || $4=="Simple_repeat" || $4=="unknown"){print $0}' annotation/$prefix.EDTA2.v0.2.bed >annotation/$prefix.Other.bed
awk '($4!="rRNA" && $4!="Simple_repeat" && $4!="unknown" && $4!="Satellite/Satellite"){print $0}' annotation/$prefix.EDTA2.v0.2.bed >annotation/$prefix.TEs.bed
# And make merged versions for later 
for type in "Other" "TEs" "Satellites" "TRF_withMers"
do
    sort -k1,1 -k2,2n annotation/$prefix.$type.bed|mergeBed -i - >annotation/$prefix.$type.merged.bed
done

# Make table that can be copied into the supplementary material 
# (first merge the regions, otherwise the fraction can exceed 100%)
module load bedtools/2.31.0
rm -f g4discovery/$prefix/common_g4s/common_g4s.annotation.table.txt
cat g4discovery/$prefix/common.G4s.txt |head -n30 |while read line
do 
    seq=$(echo $line |cut -f1 -d" ")
    file=g4discovery/$prefix/common_g4s/$seq.bed
    len=`awk '{sum+=($3-$2)}END{print sum}' $file`
    line=""
    for type in "TEs" "Satellites" "TRF_withMers" "Other" "UTR5" "promoter" "CDS" "introns" "UTR3" "lncrna" "intergenic"
    do
        #summary=`intersectBed  -a $file -b annotation/$type.bed -wao | awk '{len+=$3-$2; overlap+=$NF}END{print len"\t"overlap"\t"(overlap/len)}'`
        summary=`intersectBed  -a $file -b annotation/$prefix.$type.merged.bed -wao -nonamecheck | awk -v len=$len '{overlap+=$NF}END{print (overlap/len)}'`
        line=$line" "$summary
    done
    echo  -e $seq" "$line >>g4discovery/$prefix/common_g4s/common_g4s.annotation.table.txt
done

# Check more specific Satellites and TEs
module load bedtools/2.31.0
cat g4discovery/$prefix/common.G4s.txt |head -n30 |while read -r seq no
do 
    file=g4discovery/$prefix/common_g4s/$seq.bed
    len=`awk '{sum+=($3-$2)}END{print sum}' $file`
    line=""
    intersectBed -a $file -b annotation/$prefix.TEs.bed -wao | awk '{print $10,$13}' |sort -k1,1 |\
    awk -v l=$len -v seq=$seq '{if(NR==1){rep=$1;replen=$2}else{if(rep==$1){replen+=$2}else{frac=replen/l; if(frac>0.01){print seq"\t"rep"\t"frac}; rep=$1; replen=$2}}}END{frac=replen/l; if(frac>0.01){print seq"\t"rep"\t"frac}}'
done >g4discovery/$prefix/common_g4s/common_g4s.detailed_TE.table.txt

# Satellites
cat g4discovery/$prefix/common.G4s.txt |head -n30 |while read -r seq no
do 
    file=g4discovery/$prefix/common_g4s/$seq.bed
    len=`awk '{sum+=($3-$2)}END{print sum}' $file`
    line=""
    intersectBed -a $file -b annotation/$prefix.Satellites.bed -wao | awk '{print $10,$11}' |sort -k1,1 |\
    awk -v l=$len -v seq=$seq '{if(NR==1){rep=$1;replen=$2}else{if(rep==$1){replen+=$2}else{frac=replen/l; if(frac>0.01){print seq"\t"rep"\t"frac}; rep=$1; replen=$2}}}END{frac=replen/l; if(frac>0.01){print seq"\t"rep"\t"frac}}'
done >g4discovery/$prefix/common_g4s/common_g4s.detailed_Satellite.table.txt

# TRs
cat g4discovery/$prefix/common.G4s.txt |head -n30 |while read -r seq no
do 
    file=g4discovery/$prefix/common_g4s/$seq.bed
    len=`awk '{sum+=($3-$2)}END{print sum}' $file`
    line=""
    intersectBed -a $file -b annotation/$prefix.TRF_withMers.bed -wao -nonamecheck | awk '{print $11,$12}' |sort -k1,1 |\
    awk -v l=$len -v seq=$seq '{if(NR==1){rep=$1;replen=$2}else{if(rep==$1){replen+=$2}else{frac=replen/l; if(frac>0.01){print seq"\t"rep"\t"frac}; rep=$1; replen=$2}}}END{frac=replen/l; if(frac>0.01){print seq"\t"rep"\t"frac}}'
done >g4discovery/$prefix/common_g4s/common_g4s.detailed_TRF_withMers.table.txt



################# CHECK FOR PREVIOUS EXPERIMENTAL VALIDATION ###################

# Check if the most common G4s have been validated experimentally
cat g4discovery/$prefix/common.G4s.txt |head -n30|while read -r seq no
do 
    echo "Searching for "$seq
    grep -i $seq helpfiles/experimentally_validated.txt
    echo "#######################################################"
done 
# Some matches

#And check in chicken eG4 set from EndoQuad
cat g4discovery/$prefix/common.G4s.txt |head -n30|while read -r seq no
do 
    echo "Searching for "$seq
    grep -i $seq helpfiles/Chicken_eG4.txt
    echo "#######################################################"
done 
# Many matches!

# Check negative experimental validation 
cat g4discovery/$prefix/common.G4s.txt |head -n30|while read -r seq no
do 
    echo "Searching for "$seq
    grep -i $seq helpfiles/experimentally_not_G4s.txt
    echo "#######################################################"
done
# No matches, good!! 



############## CHECKING G4 SCORES IN DIFFERENT FUNCTIONAL REGIONS ##############

mkdir score_analysis
module load bedtools/2.31.0
prefix="bTaeGut7v0.4_MT_rDNA"
echo -e "Region\tChromosome_type\tPQS\tG4Hunter" >score_analysis/$prefix.G4.tsv
for region in "UTR5" "promoter" "CDS" "introns" "UTR3" "lncrna" "intergenic"
do 
    for type in "macro" "micro" "dot"
    do 
        echo "Processing region: "$region" on "$type" chromosomes"
        intersectBed -a <(grep $type helpfiles/$prefix.groups.txt |cut -f1 |grep -f - final_nonB/$prefix.G4.bed) -b annotation/$prefix.$region.merged.bed | awk -v t=$type -v r=$region '{print r"\t"t"\t"$4"\t"$7}' >>score_analysis/$prefix.G4.tsv
    done 
done 







