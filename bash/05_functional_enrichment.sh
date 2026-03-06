################################################################################
### CALCULATING NON-B DNA MOTIF ENRICHMENT IN GENES, AND USING ZEBRA FINCH 
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

############  PROCESS THE ANNOTATION FILE INTO GENE CATEGORIES  ################
# Only chicken and bustard have UTRs annotated in the gff files, so for the 
# others we will extract UTRs as exon minus CDS. Promoters will be 1k upstream 
# of 5' UTRs.

mkdir -p annotation
module load bedtools/2.31.0

# FIRST, GET LONGEST TRANSCRIPT PER GENE ONLY 
# There are three different styles of the gffs, so slightly different approaches
# are necessary
cat helpfiles/species_list.txt |head -n1 |while read -r sp longname prefix
do
  echo "Prepare annotation file"
  ls ref/${prefix}.gff
  echo "Get CDS, exons, UTRs, promoters, introns, non-coding RNA and intergenic" 
  #for chicken and bustard, nomenclature is different. Also, use UTRs directly from gff, 
  if [[ $sp == "chicken" || $sp == "great_bustard" ]]; 
  then
    ### === LiftOff style  GFF  === ###
    awk -F'\t' 'BEGIN{OFS="\t"}
    $3=="mRNA"{
      match($9,/ID=([^;]+)/,mid);        # save mRNA id
      tx = mid[1];

      match($9,/Parent=([^;]+)/,gid);    # save gene id
      gene = gid[1];

      tx_gene[tx] = gene;
    }

    $3=="exon"{
      match($9,/Parent=([^;]+)/,mid);
      tx = mid[1];

      len = $5 - $4 + 1;
      tx_len[tx] += len;
    }

    END{
      for (tx in tx_len){
        gene = tx_gene[tx];
        if (!(gene in best_len) || tx_len[tx] > best_len[gene]) {
          best_len[gene] = tx_len[tx];
          best_tx[gene] = tx;
        }
      }
      for (gene in best_tx)
        print gene, best_tx[gene], best_len[gene];
    }
    ' ref/${prefix}.gff > new_anno/${prefix}.longest_tx_per_gene.tsv

  elif [[ $sp == "zebra_finch" ]];
  then 
    ### ==== EGAPx Zebra finch style === ###
    awk -F'\t' 'BEGIN{OFS="\t"}
    ($3=="transcript"){
        match($9,/transcript_id ([^;]+)/,mtx);
        tx = mtx[1];

        match($9,/gene_id ([^;]+)/,mgene);
            gene = mgene[1];
        tx_gene[tx] = gene;
    }

    # sum exon lengths per transcript
    $3=="exon"{
        match($9,/transcript_id ([^;]+)/,mtx);
        tx = mtx[1];

        len = $5 - $4 + 1;
        tx_len[tx] += len;
    }

    END{
        for (tx in tx_len){
            gene = tx_gene[tx];
            if (!(gene in best_len) || tx_len[tx] > best_len[gene]) {
                best_len[gene] = tx_len[tx];
                best_tx[gene] = tx;
            }
        }
        for (gene in best_tx)
            print gene, best_tx[gene], best_len[gene];
    }
    ' ref/${prefix}.gff > new_anno/${prefix}.longest_tx_per_gene.tsv
  else
    ### === RefSeq / Gnomon style GFF  === ###
   awk -F'\t' 'BEGIN{OFS="\t"}
    ($3=="mRNA" || $3=="lnc_RNA"){
        match($9,/transcript_id=([^;]+)/,mtx);
        tx = mtx[1];

        # try gene= first
        if (match($9,/gene=([^;]+)/,mgene)) {
            gene = mgene[1];
        }
        # fallback to ID=gene-XXXX if gene= is missing
        else if (match($9,/ID=gene-([^;]+)/,gid)) {
            gene = gid[1];
        }
        # final fallback: use transcript ID as gene label (rare edge case)
        else {
            gene = tx;
        }

        tx_gene[tx] = gene;
    }

    # sum exon lengths per transcript
    $3=="exon"{
        match($9,/transcript_id=([^;]+)/,mtx);
        tx = mtx[1];

        len = $5 - $4 + 1;
        tx_len[tx] += len;
    }

    END{
        for (tx in tx_len){
            gene = tx_gene[tx];
            if (!(gene in best_len) || tx_len[tx] > best_len[gene]) {
                best_len[gene] = tx_len[tx];
                best_tx[gene] = tx;
            }
        }
        for (gene in best_tx)
            print gene, best_tx[gene], best_len[gene];
    }
    ' ref/${prefix}.gff > new_anno/${prefix}.longest_tx_per_gene.tsv
  fi
done

# Make new gffs with only longest transcripts
cat helpfiles/species_list.txt |head -n1 |while read -r sp longname prefix
do
  grep -F -f <(cut -f2 new_anno/${prefix}.longest_tx_per_gene.tsv |grep -v '^$') ref/${prefix}.gff \
  > ref/${prefix}.longest_only.gff
done 


# THEN, EXTRACT ALL DIFFERENT CLASSES, FROM LONGEST TRANSCRIPT ONLY 
cat helpfiles/species_list.txt |head -n1 |while read -r sp longname prefix
do
  echo "Prepare annotation file"
  ls ref/${prefix}.longest_only.gff
  echo "Get CDS, exons, UTRs, promoters, introns, non-coding RNA and intergenic" 
  #for chicken and bustard, nomenclature is different. Also, use UTRs directly from gff, 
  if [[ $sp == "chicken" || $sp == "great_bustard" ]]; 
  then
    #CDS FIRST:
    cat ref/$prefix.longest_only.gff |\
    awk -F'\t' '($3=="CDS"){match($9, /Parent=[^;]+/, type);
    gsub(/Parent=/, "", type[0]); s=$4-1; print $1,s,$5,type[0],".",$7}' |sed 's/ /\t/g' |sort -k1,1 -k2,2n |uniq >annotation/$prefix.CDS.bed
    echo "UTRs are annotated in the gff file" 
    # UTR5
    cat ref/$prefix.longest_only.gff  |\
    awk -F'\t' '($3=="five_prime_UTR"){match($9, /Parent=[^;]+/, type); gsub(/Parent=/, "", type[0]); s=$4-1; print $1,s,$5,type[0],".",$7}' |sed 's/ /\t/g' |sort -k1,1 -k2,2n |uniq |subtractBed -a - -b annotation/$prefix.CDS.bed >annotation/$prefix.UTR5.raw.bed
    # UTR3
    cat ref/$prefix.longest_only.gff  |\
    awk -F'\t' '($3=="three_prime_UTR"){match($9, /Parent=[^;]+/, type);gsub(/Parent=/, "", type[0]); s=$4-1; print $1,s,$5,type[0],".",$7}' |sed 's/ /\t/g' |sort -k1,1 -k2,2n |uniq |subtractBed -a - -b annotation/$prefix.CDS.bed |subtractBed -a - -b annotation/$prefix.UTR5.raw.bed  >annotation/$prefix.UTR3.raw.bed
    cat annotation/$prefix.UTR5.raw.bed annotation/$prefix.UTR3.raw.bed |sort -k1,1 -k2,2n >annotation/$prefix.UTR_all.bed
    # All exons 
    awk -F'\t' -v OFS="\t" '($3=="exon"){match($9, /Parent=[^;]+/, type);gsub(/Parent=/, "", type[0]);print $1, $4-1, $5, type[0], ".", $7}' ref/$prefix.longest_only.gff  |\
    sort -k1,1 -k4,4 -k2,2n |uniq >annotation/$prefix.exons_all.bed
  else
    if [[ $sp == "zebra_finch" ]];
    then 
      # Zebra finch has a different annotation style
      #CDS FIRST:
      cat ref/$prefix.longest_only.gff  |\
      awk -F'\t' '($3=="CDS"){match($9, /gene_id [^;]+/, type);
      gsub(/gene_id /, "", type[0]); s=$4-1; print $1,s,$5,type[0],".",$7}' |sed 's/ /\t/g' |sort -k1,1 -k2,2n |uniq >annotation/$prefix.CDS.bed
      #EXONS FOR PROTEIN CODING GENES 
      cat ref/$prefix.longest_only.gff  |\
        awk -F'\t' '($3=="exon" && $9 ~/transcript_biotype "mRNA";/){match($9, /gene_id [^;]+/, type); gsub(/gene_id /, "", type[0]); s=$4-1; print $1,s,$5,type[0],".",$7}' |sed 's/ /\t/g' |sort -k1,1 -k2,2n |uniq >annotation/$prefix.exons.bed
      # All exons (no pseudogenes are annotated in the zebra finch)
      awk -v OFS="\t" -F'\t' '($3=="exon"){match($9, /gene_id [^;]+/, type);gsub(/gene_id /, "", type[0]);print $1, $4-1, $5, type[0], ".", $7}' ref/$prefix.longest_only.gff  |\
        sort -k1,1 -k4,4 -k2,2n  |uniq >annotation/$prefix.exons_all.bed
    else
      #CDS FIRST:
      cat ref/$prefix.longest_only.gff  |\
      awk -F'\t' '($3=="CDS"){match($9, /gene=[^;]+/, type);
      gsub(/gene=/, "", type[0]); s=$4-1; print $1,s,$5,type[0],".",$7}' |sed 's/ /\t/g' |sort -k1,1 -k2,2n |uniq >annotation/$prefix.CDS.bed
       #EXONS FOR PROTEIN CODING GENES 
      cat ref/$prefix.longest_only.gff  |\
        awk -F'\t' '($3=="exon" && $9 ~/gbkey=mRNA;/){match($9, /gene=[^;]+/, type); gsub(/gene=/, "", type[0]); s=$4-1; print $1,s,$5,type[0],".",$7}' |sed 's/ /\t/g' |sort -k1,1 -k2,2n |uniq >annotation/$prefix.exons.bed
      # All exons except pseudogenes 
      awk -v OFS="\t" -F'\t' '($3=="exon" && $9!~/pseudo=true/){match($9, /;gene=[^;]+/, type);gsub(/;gene=/, "", type[0]);print $1, $4-1, $5, type[0], ".", $7}' ref/$prefix.longest_only.gff  |\
      sort -k1,1 -k4,4 -k2,2n  |uniq >annotation/$prefix.exons_all.bed
    fi
    echo "UTRs are NOT annotated in the gff file, extracting from exons minus CDS"
    #Extract all UTRs:
    subtractBed -s -a annotation/$prefix.exons.bed -b annotation/$prefix.CDS.bed >annotation/$prefix.UTR_all.bed
    # Define CDS bounds per transcript 
    bedtools groupby -i annotation/$prefix.CDS.bed -g 1,6,4 -c 2,3 -o min,max > annotation/$prefix.CDS_bounds.tsv
    # Divide UTRs into 3' and 5'
    awk 'BEGIN{OFS="\t"}
    NR==FNR {cs[$3]=$4; ce[$3]=$5; st[$3]=$2; next}
    {
      tid=$4
      if (!(tid in cs)) next
      if (st[tid]=="+" && $3 <= cs[tid])
        print > "annotation/'$prefix'.UTR5.raw.bed"
      else if (st[tid]=="+" && $2 >= ce[tid])
        print > "annotation/'$prefix'.UTR3.raw.bed"
      else if (st[tid]=="-" && $2 >= ce[tid])
        print > "annotation/'$prefix'.UTR5.raw.bed"
      else if (st[tid]=="-" && $3 <= cs[tid])
        print > "annotation/'$prefix'.UTR3.raw.bed"
    }' annotation/$prefix.CDS_bounds.tsv annotation/$prefix.UTR_all.bed     
  fi
  # ALL INTRONS (take the space between exons per gene)
  awk 'BEGIN{OFS="\t"}{if (NR==1) {
      chr=$1; start=$2; end=$3; gene=$4; strand=$6
      next
    }
    if ($1==chr && $4==gene) {
      if (end < $2)
        print chr, end, $2, gene, ".", strand
      end=$3
    } else {
      chr=$1; start=$2; end=$3; gene=$4; strand=$6
    }
  }' annotation/$prefix.exons_all.bed \
  > annotation/$prefix.introns_all.raw.bed
  # INTRONS between coding regions only 
  sort -k1,1 -k4,4 -k2,2n annotation/$prefix.CDS.bed | awk 'BEGIN{OFS="\t"}{if (NR==1) {
      chr=$1; start=$2; end=$3; gene=$4; strand=$6
      next
    }
    if ($1==chr && $4==gene) {
      if (end < $2)
        print chr, end, $2, gene, ".", strand
      end=$3
    } else {
      chr=$1; start=$2; end=$3; gene=$4; strand=$6
    }
  }'  > annotation/$prefix.introns_coding.raw.bed
  # PROMOTERS, Take 1k upstream from 5'UTRs)
  cat annotation/$prefix.UTR5.raw.bed |awk -v OFS="\t" '{if(NR==1){chr=$1; tr=$4; str=$6; if($6=="+"){start=$2}else if($6=="-"){start=$3}else{start="error"}}else{if($4==tr){if($6=="-"){start=$3}}else{if(str=="+"){e=start; s=e-1000}else{s=start; e=s+1000}; if(s<0){s=0}; print chr,s,e,tr,".",str; chr=$1; tr=$4; str=$6; if($6=="+"){start=$2}else if($6=="-"){start=$3}else{start="error"}}}}END{if(str=="+"){e=start; s=e-1000}else{s=start; e=s+1000}; if(s<0){s=0}; print chr,s,e,tr,".",str;}' >annotation/$prefix.promoter.raw.bed
  # INTERGENIC 
  awk -v OFS='\t' '{if($3=="gene"){s=$4-1; print $1,s,$5}}' ref/$prefix.gff  |sort -k1,1 -k2,2n |mergeBed -i -|awk -v OFS='\t' '{if(NR==1){chr=$1; end=$3}else{if($1==chr){print chr,end,$2; end=$3}else{chr=$1; end=$3}}}' >annotation/$prefix.intergenic.raw.bed
  # LONG NON-CODING RNAs - ZF doesn't have them annotated in column 3
  if [[ $sp == "zebra_finch" ]];
  then 
    awk -F'\t' -v OFS="\t" '($3=="exon" && $9 ~ /transcript_biotype "lnc_RNA"/){match($9, /gene_id [^;]+/, type);gsub(/gene_id /, "", type[0]);print $1, $4-1, $5, type[0], ".", $7}' ref/$prefix.longest_only.gff  >annotation/$prefix.lncrna.raw.bed
  else
    awk -F'\t' -v OFS="\t" '($3=="lnc_RNA"){match($9, /gene=[^;]+/, type);gsub(/gene=/, "", type[0]);print $1, $4-1, $5, type[0], ".", $7}' ref/$prefix.longest_only.gff  >annotation/$prefix.lncrna.raw.bed
  fi
  # ====================== CLEAN UP OVERLAPPING REGIONS ==========================
  echo "clean up and remove overlap"
  # UTRs (remove CDS)
  subtractBed -a annotation/$prefix.UTR5.raw.bed -b annotation/$prefix.CDS.bed -nonamecheck >annotation/$prefix.UTR5.bed
  subtractBed -a annotation/$prefix.UTR3.raw.bed -b annotation/$prefix.CDS.bed -nonamecheck|subtractBed -a - -b annotation/$prefix.UTR5.bed -nonamecheck >annotation/$prefix.UTR3.bed
  # Promoters (remove CDS and UTRs)
  subtractBed -a annotation/$prefix.promoter.raw.bed -b annotation/$prefix.CDS.bed -nonamecheck |subtractBed -a - -b annotation/$prefix.UTR_all.bed -nonamecheck >annotation/$prefix.promoter.bed
  # Introns (remove all above)
  subtractBed -a annotation/$prefix.introns_all.raw.bed -b annotation/$prefix.CDS.bed -nonamecheck |subtractBed -a - -b annotation/$prefix.UTR_all.bed -nonamecheck |subtractBed -a - -b annotation/$prefix.promoter.bed -nonamecheck >annotation/$prefix.introns.bed
  # introns coding 
   subtractBed -a annotation/$prefix.introns_coding.raw.bed -b annotation/$prefix.CDS.bed -nonamecheck |subtractBed -a - -b annotation/$prefix.UTR_all.bed -nonamecheck |subtractBed -a - -b annotation/$prefix.promoter.bed -nonamecheck >annotation/$prefix.introns_coding.bed
   # Introns noncoding 
   subtractBed -a annotation/$prefix.introns.bed -b annotation/$prefix.introns_coding.bed -nonamecheck >annotation/$prefix.introns_noncoding.bed
  # Long non-coding RNAs (remove CDS, UTRs, promoters)
  subtractBed -a annotation/$prefix.lncrna.raw.bed -b annotation/$prefix.CDS.bed -nonamecheck |subtractBed -a - -b annotation/$prefix.UTR_all.bed -nonamecheck |subtractBed -a - -b annotation/$prefix.promoter.bed -nonamecheck >annotation/$prefix.lncrna.bed
  # Intrergenic 
  subtractBed -a annotation/$prefix.intergenic.raw.bed -b annotation/$prefix.CDS.bed -nonamecheck |subtractBed -a - -b annotation/$prefix.UTR_all.bed -nonamecheck |subtractBed -a - -b annotation/$prefix.promoter.bed -nonamecheck |subtractBed -a - -b annotation/$prefix.introns_all.bed -nonamecheck >annotation/$prefix.intergenic.bed
done 


############### CREATE MERGED FILES FOR EACH FUNCTIONAL CATEGORY ###############

cat helpfiles/species_list.txt |head -n1 |while read -r sp longname prefix
do
  for class in "lncrna" # "introns" "promoter" "intergenic"  "CDS" "UTR5" "UTR3"  "introns_noncoding" "introns_coding"
  do
      sort -k1,1 -k2,2n annotation/$prefix.$class.bed | mergeBed -i - >annotation/$prefix.$class.merged.bed
  done
done 


######################## DIVIDE ANNOTATIONS PER GROUP ###########################

mkdir -p annotation/group_wise/
cat helpfiles/species_list.txt | head -n1 |while read -r sp longname prefix
do
  for group in "macro" "micro" "dot"
  do
    for class in "lncrna"  "introns" "promoter" "intergenic"  "CDS" "UTR5" "UTR3" "introns_noncoding" "introns_coding"
    do
      grep $group helpfiles/$prefix.groups.txt |awk 'NR==FNR{a[$1];next} ($1 in a){print}' - annotation/$prefix.$class.merged.bed >annotation/group_wise/$prefix.$group.$class.merged.bed
      grep $group helpfiles/$prefix.groups.txt |awk 'NR==FNR{a[$1];next} ($1 in a){print}' - annotation/$prefix.$class.bed >annotation/group_wise/$prefix.$group.$class.bed
    done
  done 
done


################### CHECK ANNOTATED LENGTHS FOR EACH SPECIES ###################
# Combined length
cat helpfiles/species_list.txt |grep zebra |while read -r sp longname prefix
do
  echo "========================================"
  echo "$sp annotation, in bp:"
  for group in "macro" "micro" "dot"
  do
    # Calc number of genes in total 
    if [[ $sp == "zebra_finch" ]]; then
      prot_genes=`cat annotation/group_wise/$prefix.$group.CDS.bed |cut -f1,4 |grep -E '_mat|chrZ_pat' |cut -f2 |sort |uniq |wc -l`
      nc_genes=`cat annotation/group_wise/$prefix.$group.lncrna.bed |cut -f1,4 |grep -E '_mat|chrZ_pat' |cut -f2 |sort |uniq |wc -l`
    else 
      prot_genes=`cat annotation/group_wise/$prefix.$group.CDS.bed |cut -f4|sort |uniq |wc -l`
      nc_genes=`cat annotation/group_wise/$prefix.$group.lncrna.bed |cut -f4|sort |uniq |wc -l`
    fi
    tmp=$prot_genes" "$nc_genes
    for class in "promoter" "UTR5" "CDS" "UTR3" "introns" "lncrna" "intergenic"
    do
      if [[ $sp == "zebra_finch" ]]; then
        mtot=`mergeBed -i annotation/group_wise/$prefix.$group.$class.merged.bed |grep -E '_mat|chrZ_pat' |awk '{sum+=$3-$2}END{print sum}'`
      else 
        #ltot=`awk '{sum+=$3-$2}END{print sum}' annotation/group_wise/$prefix.$group.$class.bed`
        mtot=`mergeBed -i annotation/group_wise/$prefix.$group.$class.merged.bed |awk '{sum+=$3-$2}END{print sum}'`
        #echo "$group $class $ltot, merged $mtot"
      fi
      tmp=$tmp" "$mtot
    done
    echo $tmp
  done 
done  

# Mean length per gene 
cat helpfiles/species_list.txt |head -n1 |while read -r sp longname prefix
do
  echo "========================================"
  echo "$sp annotation, in bp:"
  for class in "lncrna"  "introns" "promoter" "UTR5"  "CDS" "UTR3" "intergenic"
  do
    for group in "macro" "micro" "dot"
    do
      mean=`awk '{gene=$4; len=$3-$2; if(NR==1){cur=gene; sum=len; next}; 
      if(cur==gene){sum+=len}
      else{print cur,sum; cur=gene; sum=len}
      }END{print cur,sum}' annotation/group_wise/$prefix.$group.$class.bed |awk '{sum+=$2}END{print sum/NR}'`
       ltot=`awk '{sum+=$3-$2}END{print sum}' annotation/group_wise/$prefix.$group.$class.bed`
       mtot=`awk '{sum+=$3-$2}END{print sum}' annotation/group_wise/$prefix.$group.$class.merged.bed`
      echo "$group $class Mean: $mean, tot length: $ltot, merged $mtot"
    done
  done 
done  

#
# Lengths of chromosome groups 
prefix="bTaeGut7v0.4_MT_rDNA"
rm -f stats/$prefix.group.lengths.txt
for group in "macro" "micro" "dot" 
do 
  grep $group helpfiles/$prefix.groups.txt  |awk -v g=$group '{sum+=$2}END{print g"\t"sum'} >>stats/$prefix.group.lengths.txt
done

# Number of genes 
rm -f stats/$prefix.group.geneNo.txt
for group in "macro" "micro" "dot" 
do 
  no=`grep $group helpfiles/$prefix.groups.txt |cut -f1 |grep -f -  <(awk '($3=="gene"){print $1,$10}'  ref/bTaeGut7v0.4_MT_rDNA.EGAPx.v0.1.gtf) |wc -l`
  echo -e "$group\t$no" >>stats/$prefix.group.geneNo.txt
done


# Avg length of genes, introns and exons  
rm -f stats/$prefix.group.avggenelen.txt
for group in "macro" "micro" "dot" 
do 
  leng=`grep $group helpfiles/$prefix.groups.txt |cut -f1 |grep -f -   <(awk '($3=="gene"){print $1,$4,$5}' ref/bTaeGut7v0.4_MT_rDNA.EGAPx.v0.1.gtf)\
   |awk '{sum+=$3-$2; n++}END{print sum/n}'`
  lencds=`grep $group helpfiles/$prefix.groups.txt |cut -f1 |grep -f -  annotation/bTaeGut7v0.4_MT_rDNA.CDS.bed\
   |awk '{sum+=$3-$2; n++}END{print sum/n}'`
  lenintron=`grep $group helpfiles/$prefix.groups.txt |cut -f1 |grep -f -  annotation/bTaeGut7v0.4_MT_rDNA.introns.bed\
   |awk '{sum+=$3-$2; n++}END{print sum/n}'`
   exonno=`grep $group helpfiles/$prefix.groups.txt |cut -f1 |grep -f -  annotation/bTaeGut7v0.4_MT_rDNA.CDS.bed |cut -f4 |sort |uniq -c |awk '{sum+=$1; n++}END{print sum/n}'`
  echo -e "$group\t$leng\t$lencds\t$lenintron\t$exonno" >>stats/$prefix.group.avggenelen.txt
done

# Sum of repeats 
rm -f stats/$prefix.group.replen.txt
for group in "macro" "micro" "dot" 
do 
  len=`grep $group helpfiles/$prefix.groups.txt |cut -f1 |grep -f -   annotation/all_repeats.bed \
   |awk '{sum+=$3-$2}END{print sum}'`
  echo -e "$group\t$len" >>stats/$prefix.group.replen.txt
done





############################# CALCULATE ENRICHMENT #############################

mkdir -p functional/
cat helpfiles/species_list.txt |head -n1 |while read -r sp longname prefix
do
  for class in  "introns" "lncrna"  "promoter" "intergenic" "CDS" "UTR5" "UTR3"  #"introns_noncoding" "introns_coding"
  do
      echo '#!/bin/bash
      module load bedtools/2.31.0
      len=`awk '"'"'{sum+=$3-$2}END{print sum}'"'"' annotation/'$prefix'.'$class'.merged.bed`
      echo "Length of '$class' is $len"
      rm -f tmp.genome.'$prefix'.'$class'
      cat coverage/'${prefix}'.per_genome.tsv | while read -r non_b tot dens;
      do
          echo "looking at $non_b"
          d=`intersectBed -a annotation/'$prefix'.'$class'.merged.bed -b final_nonB/'${prefix}'.${non_b}.merged.bed -wo |awk -v l=$len -v dtot=$dens '"'"'{sum+=$7}END{d=sum/l; frac=d/dtot; print d,frac}'"'"'`
          echo '$sp'" "'$class'" "$non_b" "$d >>tmp.genome.'$prefix'.'$class'
      done
      ' | sbatch -J $class --ntasks=1 --cpus-per-task=1 --time=15:00 --out slurm/job.functional.$prefix.$class.%j.out
  done
done 
--account=kdm16_sc_default --partition=sla-prio 

# FOR TESTING - REMOVE 
len=`awk '{sum+=$3-$2}END{print sum}' annotation/$prefix.$class.merged.bed`
intersectBed -a annotation/$prefix.$class.merged.bed -b final_nonB/${prefix}.${non_b}.merged.bed -wo |awk -v l=$len -v dtot=$dens '{sum+=$7}END{d=sum/l; frac=d/dtot; print sum, l, d,frac}'

# Merge the tmp files
echo -e "Species\tClass\tnonB\tCoverage\tEnrichment_gw" >functional/8sp.enrichment.tsv
cat helpfiles/species_list.txt |while read -r sp longname prefix
do
  for class in "promoter" "intergenic" "introns" "CDS" "UTR5" "UTR3"  "lncrna"
  do
     # Chicken and great bustard dont have annotated lncrna annotated
     if [[ ($sp == "chicken" || $sp == "great_bustard") && $class == "lncrna" ]]; 
      then
        echo "Skipping $sp $class"
      else
        cat tmp.genome.$prefix.$class |sed "s/ /\t/g" >>functional/8sp.enrichment.tsv
    fi
  done
done 

# Remove temp files 
rm tmp.genome.*

# ~~~~~~~~~~~~~~~~~ CALCULATE ENRICHMENT FOR CHROMOSOME TYPES ~~~~~~~~~~~~~~~~~~

cat helpfiles/species_list.txt |head -n1 |while read -r sp longname prefix
do
  for group in "macro" "micro" "dot"
  do
    for class in "introns"  "lncrna"  "promoter" "intergenic"  "CDS" "UTR5" "UTR3" "introns_noncoding" "introns_coding"
    do
        echo '#!/bin/bash
        module load bedtools/2.31.0
        len=`awk '"'"'{sum+=$3-$2}END{print sum}'"'"' annotation/group_wise/'$prefix'.'$group'.'$class'.merged.bed`
        echo "Length of '$class' in '$prefix' '$group' chromosomes is $len"
        rm -f tmp.group.'$prefix'.'$group'.'$class'
        cat coverage/'${prefix}'.per_genome.tsv | while read -r non_b tot dens;
        do
            echo "looking at $non_b"
            d=`intersectBed -a annotation/group_wise/'$prefix'.'$group'.'$class'.merged.bed \
            -b final_nonB/'${prefix}'.${non_b}.merged.bed -wo -nonamecheck |\
            awk -v l=$len -v dtot=$dens '"'"'{sum+=$7}END{d=sum/l; frac=d/dtot; print d,frac}'"'"'`
            echo '$sp'" "'$group'" "'$class'" "$non_b" "$d >>tmp.group.'$prefix'.'$group'.'$class'
        done
        ' | sbatch -J $class.$prefix.$group --ntasks=1 --cpus-per-task=1 --time=15:00 --out slurm/job.funcEnrichGroup.$prefix.$group.$class.%j.out
    done
  done
done 
--account=kdm16_sc_default --partition=sla-prio


######### TESTING, REMOVE 
# BEFORE:
 head tmp.group.bTaeGut7v0.4_MT_rDNA.dot.introns
zebra_finch dot introns APR 0.0200092 3.48963
zebra_finch dot introns DR 0.0919522 3.81054
zebra_finch dot introns STR 0.0138861 1.12899
zebra_finch dot introns IR 0.0325626 1.02102
zebra_finch dot introns TRI 0.00278236 1.14661
zebra_finch dot introns G4 0.0562111 4.12571
zebra_finch dot introns Z 0.000853623 0.865961
zebra_finch dot introns Any 0.168652 2.23097

head tmp.group.chicken.v23.dot.introns
chicken dot introns APR 0.0188764 3.07243
chicken dot introns DR 0.250148 10.2212
chicken dot introns STR 0.0299511 1.93813
chicken dot introns IR 0.0206135 0.75255
chicken dot introns TRI 0.00386724 1.19675
chicken dot introns G4 0.252303 14.9569
chicken dot introns Z 0.00134522 0.943836
chicken dot introns Any 0.42789 5.61098

# AFTER 
zebra_finch dot introns APR 0.0484317 8.44657
zebra_finch dot introns DR 0.234337 9.71102
zebra_finch dot introns STR 0.0189764 1.54284
zebra_finch dot introns IR 0.0335849 1.05308
zebra_finch dot introns TRI 0.00291905 1.20294
zebra_finch dot introns G4 0.143902 10.5619
zebra_finch dot introns Z 0.000575377 0.583693
zebra_finch dot introns Any 0.359475 4.75522

chicken dot introns APR 0.0188764 3.07243
chicken dot introns DR 0.250148 10.2212
chicken dot introns STR 0.0299511 1.93813
chicken dot introns IR 0.0206135 0.75255
chicken dot introns TRI 0.00386724 1.19675
chicken dot introns G4 0.252303 14.9569
chicken dot introns Z 0.00134522 0.943836
chicken dot introns Any 0.42789 5.61098

 len=`grep $group helpfiles/$prefix.groups.txt |awk 'NR==FNR{a[$1];next} ($1 in a){sum+=$3-$2}END{print sum}' - annotation/$prefix.$class.merged.bed`

 intersectBed -a <(grep $group helpfiles/$prefix.groups.txt|awk 'NR==FNR{a[$1];next} ($1 in a){print}' - annotation/$prefix.$class.merged.bed) -b final_nonB/${prefix}.${non_b}.merged.bed -wo -nonamecheck |awk -v l=$len -v dtot=$dens '{sum+=$7}END{d=sum/l; frac=d/dtot; print sum, l, d,frac}'


# Merge the tmp files
echo -e "Species\tGroup\tClass\tnonB\tCoverage\tEnrichment_gw" >functional/8sp.enrichment.groups.tsv
cat helpfiles/species_list.txt |while read -r sp longname prefix
do
  for group in "macro" "micro" "dot"
  do
    for class in "promoter" "UTR5" "CDS" "UTR3" "introns" "lncrna" "intergenic"
    do
    # Chicken and great bustard dont have annotated lncrna annotated
     if [[ ($sp == "chicken" || $sp == "great_bustard") && $class == "lncrna" ]]; 
      then
        echo "Skipping $sp $class"
      else
        cat tmp.group.$prefix.$group.$class |sed 's/ /\t/g'  >>functional/8sp.enrichment.groups.tsv
      fi
    done
  done
done 
#rm tmp.group.*.*



# ~~~~~~~~ CALCULATE ENRICHMENT USING A GROUP BASED INTERGENIC AVG ~~~~~~~~~~
# Tried first dividing by group average, this was not included in the 
# paper as it doesn't make much sense for the dotchromosomes, where almost 
# everything is coding. 

cat helpfiles/species_list.txt |head -n1 |while read -r sp longname prefix
do
  for group in "macro" "micro" "dot"
  do
    for class in "lncrna"  "intergenic" "introns" "introns_noncoding" "introns_coding" "promoter" "CDS" "UTR5" "UTR3" 
    do
        echo '#!/bin/bash
        module load bedtools/2.31.0
        len=`awk '"'"'{sum+=$3-$2}END{print sum}'"'"' annotation/group_wise/'$prefix'.'$group'.'$class'.merged.bed`
        echo "Length of '$class' is $len"
        rm -f tmp.groupintergenic.'$prefix'.'$group'.'$class'
        cat functional/8sp.enrichment.groups.tsv |grep "intergenic" |awk -v s='$sp' -v g='$group' '"'"'(s==$1 && g==$2){print}'"'"' | while read -r spec gr reg non_b dens enr;
        do
            echo "looking at $non_b"
            d=`intersectBed -a annotation/group_wise/'$prefix'.'$group'.'$class'.merged.bed -b final_nonB/'${prefix}'.${non_b}.merged.bed -wo -nonamecheck |awk -v l=$len -v dtot=$dens '"'"'{sum+=$7}END{d=sum/l; frac=d/dtot; print d,frac}'"'"'`
            echo '$sp'" "'$group'" "'$class'" "$non_b" "$d >>tmp.groupintergenic.'$prefix'.'$group'.'$class'
        done
        ' | sbatch -J $class --ntasks=1 --cpus-per-task=1 --time=15:00 --out slurm/job.functional.groupavg.$group.$class.%j.out
    done
  done
done 

# Merge the tmp files
echo "Species Group Class nonB Coverage Enrichment_gw" |sed "s/ /\t/g" >functional/8sp.enrichment_groups.groupavg.tsv
cat helpfiles/species_list.txt |head -n1 |while read -r sp longname prefix
do
  for group in "macro" "micro" "dot"
  do
    for class in  "promoter" "intergenic" "introns" "CDS" "UTR5" "UTR3"  "lncrna"
    do
      if [[ ($sp == "chicken" || $sp == "great_bustard") && $class == "lncrna" ]]; 
      then
        echo "Skipping $sp $class"
      else
        cut -f2- -d " " tmp.groupintergenic.$prefix.$group.$class |sed "s/ /\t/g" |awk -v s=$sp '{print s"\t"$0}' >>functional/8sp.enrichment_groups.groupavg.tsv
      fi
    done
  done
done
#rm tmp.groupintergenic.*


################################################################################
# #################### SUBSAMBLE REGIONS TO GET ERROR BARS #####################
# CHECK SIGNIFICANCE OF THE BARS IN FIGURE 2A

# Subsample the functional regions to get a distribution of the enrichment
# values. If distribution overlaps with 1 (genome average), we deem this
# comparison NOT to be significally different from the genome average.

# Subsample 50 percent of the regions 100 times
cat helpfiles/species_list.txt |head -n1 |while read -r sp longname prefix
do
  echo "Starting subsampling for $prefix"
  mkdir -p functional/resample/$prefix/
  perc=50
  for group in "macro" "micro" "dot"
  do
    for class in "introns" "intergenic" "promoter" "CDS" "UTR5" "UTR3"  "lncrna"
    do
      num=`awk -v p=$perc '{}END{num=int(NR*(p/100)); print num}' annotation/group_wise/$prefix.$group.$class.merged.bed`
      echo "we will extract $num rows from $class on $group"
      # Subsample the sequences 100 times
      echo '#!/bin/bash
      module load bedtools/2.31.0
      for i in {1..100}
      do
          echo "working on $i..."
          cat annotation/group_wise/'$prefix'.'$group'.'$class'.merged.bed |shuf -n '$num' |sort -k1,1 -k2,2n  >functional/resample/'$prefix'/Resamp.'$perc'perc.'$class'.'$group'.$i.bed
          l=`wc -l functional/resample/'$prefix'/Resamp.'$perc'perc.'$class'.'$group'.$i.bed`
          echo "...extracting $l lines"
      done
      ' |sbatch -J $prefix.$class.$group  --ntasks=1 --mem-per-cpu=4G --cpus-per-task=1 --time=1:00:00  -o slurm/job.resamp.${perc}perc.$prefix.$class.$group.%j.out
    done
  done
done

# Subsample 50 percent of the region SIZE a 100 times
cat helpfiles/species_list.txt |head -n1 |while read -r sp longname prefix
do
  echo "Starting subsampling for $prefix"
  mkdir -p functional/resample/$prefix/
  perc=50
  for group in "macro" "micro" "dot"
  do
    for class in "introns" "intergenic"  "promoter" "CDS" "UTR5" "UTR3"  "lncrna"
    do
      TARGET=$(awk -v p=$perc '{sum+=$3-$2}END{print sum*(p/100)}' annotation/group_wise/$prefix.$group.$class.merged.bed)
       echo "we will extract $TARGET bp from $class on $group"
      # Subsample the sequences 100 times
      echo '#!/bin/bash
      module load bedtools/2.31.0
      for i in {1..100}
      do
          echo "working on $i..."
          shuf annotation/group_wise/'$prefix'.'$group'.'$class'.merged.bed |\
          awk -v target='$TARGET' '"'"'{
            len = $3-$2;
            cum += len;
            print $0;
            if (cum >= target) exit;
          }'"'"' >functional/resample/'$prefix'/RegionSampling.'$perc'perc.'$class'.'$group'.$i.bed
          l=`wc -l functional/resample/'$prefix'/RegionSampling.'$perc'perc.'$class'.'$group'.$i.bed`
          echo "...extracted $l lines"
      done
      ' |sbatch -J $prefix.$class.$group  --ntasks=1 --mem-per-cpu=4G --cpus-per-task=1 --time=10:00  -o slurm/job.RegionSampling.${perc}perc.$prefix.$class.$group.%j.out
    done
  done
done

# Calculate the enrichment for each of the subsamples
cat helpfiles/species_list.txt |head -n1 |while read -r sp longname prefix
do
  echo "Calculating enrichment for resamples of $prefix"
  for subset in "50perc" 
  do
    for group in "macro" "micro" "dot"
    do
      for class in   "introns" "promoter" "intergenic" "CDS" "UTR5" "UTR3"  "lncrna"
      do
          echo '#!/bin/bash
          module load bedtools/2.31.0
          for i in {1..100}
          do
              rm -f functional/resample/'$prefix'/tmp.RegionSampling.'$subset'.'$class'.'$group'.$i.txt
              len=`awk '"'"'{sum+=$3-$2}END{print sum}'"'"' functional/resample/'$prefix'/RegionSampling.'$subset'.'$class'.'$group'.$i.bed`
              cat coverage/'${prefix}'.per_genome.tsv |grep -v "Any" | while read -r non_b tot dens;
              do
                  d=`intersectBed -a functional/resample/'$prefix'/RegionSampling.'$subset'.'$class'.'$group'.$i.bed -b final_nonB/'$prefix'.${non_b}.merged.bed -wo |awk -v l=$len -v dtot=$dens '"'"'{sum+=$7}END{d=sum/l; frac=d/dtot; print frac}'"'"'`
                  echo $i '$group' '$class' $non_b $d >>functional/resample/'$prefix'/tmp.RegionSampling.'$subset'.'$class'.'$group'.$i.txt
              done
          done
          ' |sbatch -J $prefix.$class.$group -o slurm/job.RegionSampling-enrich.$subset.$class.$group.%j.out --ntasks=1 --cpus-per-task=1 --mem-per-cpu=4G --time=1:00:00
      done
    done
  done
done 

# Merge the results
rep="100rep"
subset="50perc"
cat helpfiles/species_list.txt |head -n1 |while read -r sp longname prefix
do
  echo "Rep Group Class nonB Enrichment" |sed "s/ /\t/g" >functional/$prefix.RegionSampling.summary.$subset.$rep.txt
  cat functional/resample/$prefix/tmp.RegionSampling.$subset.*.*.{1..100}.txt |sed "s/ /\t/g" >>functional/$prefix.RegionSampling.summary.$subset.$rep.txt
done 

# We want to find the min and max for each category, but remove the ~5% most
# extreme (meaning we remove the top and bottom two values, saving a 96% 'CI')
rep="100rep"
subset="50perc" 
cat helpfiles/species_list.txt |head -n1 |while read -r sp longname prefix
do
  echo "Group Class nonB Min Max" |sed "s/ /\t/g" >functional/$prefix.RegionSampling.minmax.$subset.$rep.96CI.txt
  for group in "macro" "micro" "dot"
  do
    for class in "promoter" "intergenic" "introns" "CDS" "UTR5" "UTR3"  "lncrna"
    do
        if [[ ($sp == "chicken" || $sp == "great_bustard") && $class == "lncrna" ]]; 
        then
          echo "Skipping $sp $class"
        else
          for nonb in APR DR G4 IR TRI STR Z
          do
              min=`grep $class functional/$prefix.RegionSampling.summary.$subset.$rep.txt |awk -v nb=$nonb -v g=$group -v OFS="\t" '($4==nb && $2==g){print $5}' |sort -n |head -n3 |tail -n1`
              max=`grep $class functional/$prefix.RegionSampling.summary.$subset.$rep.txt |awk -v nb=$nonb -v g=$group -v OFS="\t" '($4==nb && $2==g){print $5}' |sort -n |tail -n3 |head -n1`
              echo $group $class $nonb $min $max |sed "s/ /\t/g" >>functional/$prefix.RegionSampling.minmax.$subset.$rep.96CI.txt
          done
        fi
    done
  done
done



################################################################################
####################### SEPARATE G4s BASED ON STRANDNESS #######################
# We want to check if G4s occur more often on the coding or template strand

cat helpfiles/species_list.txt |head -n1 |while read -r sp longname prefix
do
  for group in "macro" "micro" "dot"
  do
    for class in "lncrna" "introns" "promoter" "intergenic" "introns_noncoding" "introns_coding" "CDS" "UTR5" "UTR3"  
    do
        echo '#!/bin/bash
        module load bedtools/2.31.0
        len=`awk '"'"'{sum+=$3-$2}END{print sum}'"'"' annotation/group_wise/'$prefix'.'$group'.'$class'.merged.bed`
        echo "Length of '$class' is $len"
        rm -f tmp.strand.'$prefix'.'$group'.'$class'
        dsame=`intersectBed -a annotation/group_wise/'$prefix'.'$group'.'$class'.bed -b final_nonB/'${prefix}'.G4.bed -s |sort -k1,1 -k2,2n |mergeBed -i - |awk -v l=$len '"'"'{sum+=$3-$2}END{d=sum/l; print d}'"'"'`
        doppo=`intersectBed -a annotation/group_wise/'$prefix'.'$group'.'$class'.bed -b final_nonB/'${prefix}'.G4.bed -S |sort -k1,1 -k2,2n |mergeBed -i - |awk -v l=$len '"'"'{sum+=$3-$2}END{d=sum/l; print d}'"'"'`
        echo '$sp'" "'$group'" "'$class'" Coding "$dsame >>tmp.strand.'$prefix'.'$group'.'$class'
        echo '$sp'" "'$group'" "'$class'" Template "$doppo >>tmp.strand.'$prefix'.'$group'.'$class'
        ' | sbatch -J $class --ntasks=1 --cpus-per-task=1 --mem-per-cpu=4G --time=1:00:00 --out slurm/functional.$group.$class.%j.out
    done
  done
done

#Merge
echo -e "Species\tGroup\tClass\tStrand\tCoverage" >functional/8sp.G4_strand.groups.tsv
cat helpfiles/species_list.txt |while read -r sp longname prefix
do
  for group in "macro" "micro" "dot"
  do
    for class in "promoter" "intergenic" "introns" "CDS" "UTR5" "UTR3"  "lncrna"
    do
     if [[ ($sp == "chicken" || $sp == "great_bustard") && $class == "lncrna" ]]; 
      then
        echo "Skipping $sp $class"
      else
        cut -f2- -d " " tmp.strand.$prefix.$group.$class |sed "s/ /\t/g" |awk -v s=$sp '{print s"\t"$0}' >>functional/8sp.G4_strand.groups.tsv
      fi
    done
  done
done 
  #rm tmp.strand.*

# To simplify downstream code, I want annotation files per group, where  
# overlapping fragments from the same gene are merged into one, but not 
# overlaps from DIFFERENT genes 
mkdir -p annotation/gene_merge/
cat helpfiles/species_list.txt |head -n1 |while read -r sp longname prefix
do
  for group in "macro" "micro" "dot"
  do
    for class in  "introns" "promoter" "CDS" "UTR5" "UTR3"  "lncrna" #"intergenic"
    do
      echo '#!/bin/bash
      cat annotation/group_wise/'$prefix'.'$group'.'$class'.bed | \
      sort -k1,1 -k6,6 -k4,4 -k2,2n | \
      awk '"'"'{
      chr=$1; start=$2; end=$3; gene=$4; strand=$6;
      if(NR==1){
          cur_chr=chr; cur_gene=gene; cur_strand=strand;
          mstart=start; mend=end;
          next;
      }
      # same gene + same chr + same strand AND overlapping?
      if(chr==cur_chr && gene==cur_gene && strand==cur_strand && start <= mend){
          if(end > mend) mend = end;
      } else {
          print cur_chr, mstart, mend, cur_gene, ".", cur_strand;
          cur_chr=chr; cur_gene=gene; cur_strand=strand;
          mstart=start; mend=end;
        }
      }
      END{
        print cur_chr, mstart, mend, cur_gene, ".", cur_strand;
      }'"'"' OFS="\t" >annotation/gene_merge/'$prefix'.'$group'.'$class'.bed
      '| sbatch -J $class --ntasks=1 --cpus-per-task=1 --mem-per-cpu=2G --time=5:00 --out slurm/job.genemerge.annotation.$prefix.$group.$class.%j.out
    done
  done 
done 



#TEST OUTSIDE SLURM - TO BE DELETED
 grep $group helpfiles/$prefix.groups.txt | cut -f1 | \
      awk 'NR==FNR{a[$1];next} $1 in a' - annotation/$prefix.$class.bed | \
      sort -k1,1 -k6,6 -k4,4 -k2,2n | \
      awk '{
      chr=$1; start=$2; end=$3; gene=$4; strand=$6;
      if(NR==1){
          cur_chr=chr; cur_gene=gene; cur_strand=strand;
          mstart=start; mend=end;
          next;
      }
      # same gene + same chr + same strand AND overlapping?
      if(chr==cur_chr && gene==cur_gene && strand==cur_strand && start <= mend){
          if(end > mend) mend = end;
      } else {
          print cur_chr, mstart, mend, cur_gene, ".", cur_strand;
          cur_chr=chr; cur_gene=gene; cur_strand=strand;
          mstart=start; mend=end;
        }
      }
      END{
        print cur_chr, mstart, mend, cur_gene, ".", cur_strand;
      }' OFS="\t"


# To get distributions for Coding and Template mean bars, I would like to get a
# coverage value for each gene
cat helpfiles/species_list.txt |head -n1 |while read -r sp longname prefix
do
  for group in "macro" "micro" "dot"
  do
    for class in "introns" "lncrna"  "promoter" "CDS" "UTR5" "UTR3"  
    do
        echo '#!/bin/bash
        module load bedtools/2.31.0
        rm -f tmp.strand.perGene.'$prefix'.'$group'.'$class'
        # SAME STRAND 
        intersectBed -a annotation/gene_merge/'$prefix'.'$group'.'$class'.bed -b final_nonB/'${prefix}'.G4.bed -s -wao |\
        awk -v sp='$sp' -v g='$group' -v c='$class' -v OFS="\t" '"'"'{tr=$4; len=$3-$2; g4=$NF;
         if(NR==1){
          cur=tr; prev_start=$2; tot=len; sum=g4; next;
         }
        if(tr==cur){
          sum+=g4; 
          if($2!=prev_start){
            tot+=len
          }
          prev_start=$2
        }else{
          d=sum/tot; 
          print sp,g,c,"Coding",cur,tot,sum,d; tr=$4;
          cur=tr; tot=len; sum=g4; prev_start=$2}
        }
        END{
          d=sum/tot; 
          print sp,g,c,"Coding",cur,tot,sum,d; tr=$4;}'"'"' >>tmp.strand.perGene.'$prefix'.'$group'.'$class'
      # OPPOSITE STRAND 
        intersectBed -a annotation/gene_merge/'$prefix'.'$group'.'$class'.bed -b final_nonB/'${prefix}'.G4.bed -S -wao |\
        awk -v sp='$sp' -v g='$group' -v c='$class' -v OFS="\t" '"'"'{if(NR==1){tr=$4; s=$2; tot=$3-$2; g4=$NF}
      else{if(tr==$4){g4+=$NF; if(s!=$2){tot+=$3-$2}; s=$2}
      else{d=g4/tot; print sp,g,c,"Template",tr,tot,g4,d; tr=$4; s=$2; tot=$3-$2; g4=$NF}
      }}END{d=g4/tot; print sp,g,c,"Template",tr,tot,g4,d}'"'"' >>tmp.strand.perGene.'$prefix'.'$group'.'$class'
        ' | sbatch -J $class --ntasks=1 --cpus-per-task=1 --mem-per-cpu=4G --time=5:00:00 --out slurm/job.functional.perGene.$prefix.$group.$class.%j.out
    done
  done
done 

#TEST OUTSIDE SLURM - TO BE DELETED 
intersectBed -a annotation/gene_merge/$prefix.$group.$class.bed -b final_nonB/${prefix}.G4.bed -s -wao |sort -k1,1 -k2,2n | mergeBed -i - -s -c 4 -o distinct


# Merge
echo "Species Group Class Strand Gene GenLen G4Len Coverage" |sed "s/ /\t/g" >functional/8sp.G4_strand.perGene.tsv
cat helpfiles/species_list.txt |while read -r sp longname prefix
do
  for group in "macro" "micro" "dot"
  do
    for class in  "promoter" "introns" "CDS" "UTR5" "UTR3"  "lncrna"
    do
      if [[ ($sp == "chicken" || $sp == "great_bustard") && $class == "lncrna" ]]; 
      then
        echo "Skipping $sp $class"
      else
        cat tmp.strand.perGene.$prefix.$group.$class |sed 's/"//g' |sed 's/;//g' >>functional/8sp.G4_strand.perGene.tsv
      fi
    done
  done
done 
#rm tmp.strand.perGene.*



################################################################################
################################# METHYLATION ##################################
mkdir -p methylation/functional

# Convert methylation bigwig file to bedgraph (assume UCSC bigWigToWig tool is
# installed in ~/software/)
~/software/bigWigToWig ref/bTaeGut7v0.4_MT_rDNA.PBmethylation.v0.1.bw ref/bTaeGut7v0.4_MT_rDNA.PBmethylation.v0.1.wig
# The Wig file is actually in BedGraph format, only need to remove the # lines
cat ref/bTaeGut7v0.4_MT_rDNA.PBmethylation.v0.1.wig |awk '($1!~/^#/){print}' >annotation/bTaeGut7v0.4_MT_rDNA.PBmethylation.v0.1.bed
# The file contains 22M methylated sites.

# Because the size of the file, I divide it into 3 parts before overlapping:
prefix="bTaeGut7v0.4_MT_rDNA"
for group in "macro" "micro" "dot"
do
  grep $group helpfiles/$prefix.groups.txt |cut -f1 |grep -f - annotation/$prefix.PBmethylation.v0.1.bed >methylation/$group.methylation.bed
done

# +++++++++++ TESTING 
head methylation/$group.methylation.bed
chr16_mat       25      26      7.4
chr16_mat       91      92      8.3
chr16_mat       21120   21121   59.1
chr16_mat       21233   21234   61.8
chr16_mat       21240   21241   55.2


# First I want to extract all G4s within functional regions, and the functional
# regions with G4s removed. Here, I don't merge overlapping genes because I want 
# to know each separate gene that contains a CpG site 
mkdir -p annotation/G4s
module load bedtools/2.31.0
mkdir -p methylation/G4s/
prefix="bTaeGut7v0.4_MT_rDNA"
for class in  "lncrna" "introns"  "promoter"  "CDS" "UTR5" "UTR3" 
do
  subtractBed -a annotation/$prefix.$class.bed -b final_nonB/$prefix.G4.bed -nonamecheck >annotation/G4s/$class.exclG4s.bed
  intersectBed -a annotation/$prefix.$class.bed -b final_nonB/$prefix.G4.bed -nonamecheck >annotation/G4s/$class.onlyG4s.bed
  intersectBed -a annotation/$prefix.$class.bed -b final_nonB/$prefix.G4.bed -s -nonamecheck >annotation/G4s/$class.onlyG4s.samestrand.bed
  intersectBed -a annotation/$prefix.$class.bed -b final_nonB/$prefix.G4.bed -S -nonamecheck >annotation/G4s/$class.onlyG4s.oppositestrand.bed
done
for class in "intergenic" 
do  
  subtractBed -a annotation/$prefix.$class.bed -b final_nonB/$prefix.G4.bed -nonamecheck |sort -k1,1 -k2,2n |mergeBed -i - >annotation/G4s/$class.exclG4s.bed
  intersectBed -a annotation/$prefix.$class.bed -b final_nonB/$prefix.G4.bed -nonamecheck |sort -k1,1 -k2,2n |mergeBed -i - >annotation/G4s/$class.onlyG4s.bed
done 

head annotation/G4s/$class.exclG4s.bed
chr10_mat       236728  237728  "chr10_mat_egapxtmp_032076"     .       +
chr10_mat       286799  287204  "chr10_mat_egapxtmp_032078"     .       +
chr10_mat       379545  379700  "chr10_mat_egapxtmp_032074"     .       -
chr10_mat       379718  379752  "chr10_mat_egapxtmp_032074"     .       -
chr10_mat       379780  379806  "chr10_mat_egapxtmp_032074"     .       -
chr10_mat       379825  379992  "chr10_mat_egapxtmp_032074"     .       -
chr10_mat       380019  380419  "chr10_mat_egapxtmp_032074"     .       -

# Make a list with all methylated sites within those regions. 
# We still keep overlaps so we can get the mean and median per gene later 
mkdir -p methylation/G4s/
prefix="bTaeGut7v0.4_MT_rDNA"
for group in "macro" "micro" "dot" 
do
  for class in "lncrna" "introns" "CDS" "UTR5" "UTR3" "promoter"  # 
  do
      echo '#!/bin/bash
      module load bedtools/2.31.0
      intersectBed -a <(grep '$group' helpfiles/'$prefix'.groups.txt |cut -f1 | grep -f - annotation/G4s/'$class'.exclG4s.bed) -b methylation/'$group'.methylation.bed -wao |awk -v OFS="\t" -v cl='$class' '"'"'($NF==1){print cl,"background",$(NF-1),$4}'"'"' >methylation/G4s/'$group'.'$class'.txt
      intersectBed -a <(grep '$group' helpfiles/'$prefix'.groups.txt |cut -f1 | grep -f - annotation/G4s/'$class'.onlyG4s.samestrand.bed) -b methylation/'$group'.methylation.bed -wao |awk -v OFS="\t" -v cl='$class' '"'"'($NF==1){print cl,"coding",$(NF-1),$4}'"'"' >>methylation/G4s/'$group'.'$class'.txt
      intersectBed -a <(grep '$group' helpfiles/'$prefix'.groups.txt |cut -f1 | grep -f - annotation/G4s/'$class'.onlyG4s.oppositestrand.bed) -b methylation/'$group'.methylation.bed -wao |awk -v OFS="\t" -v cl='$class' '"'"'($NF==1){print cl,"template",$(NF-1),$4}'"'"' >>methylation/G4s/'$group'.'$class'.txt
    ' | sbatch -J $class.$group --ntasks=1 --cpus-per-task=1  --time=10:00 --out slurm/job.methylation.$group.$class.%j.out
  done
  # Run intergenic separately since it is not strand specific
  for class in "intergenic"
  do
      echo '#!/bin/bash
      module load bedtools/2.31.0
      intersectBed -a annotation/G4s/'$class'.exclG4s.bed -b methylation/'$group'.methylation.bed -wao |awk -v OFS="\t" -v cl='$class' '"'"'($NF==1){print cl,"background",$(NF-1)}'"'"' >methylation/G4s/'$group'.'$class'.txt
      intersectBed -a annotation/G4s/'$class'.onlyG4s.bed -b methylation/'$group'.methylation.bed -wao |awk -v OFS="\t" -v cl='$class' '"'"'($NF==1){print cl,"ignorant",$(NF-1)}'"'"' >>methylation/G4s/'$group'.'$class'.txt
      ' | sbatch -J $class.$group --ntasks=1 --cpus-per-task=1 --mem-per-cpu=4G --time=15:00 --out slurm/job.methylation.$group.$class.%j.out
  done
done

# Test
 intersectBed -a <(grep $group helpfiles/$prefix.groups.txt |cut -f1 | grep -f - annotation/G4s/$class.exclG4s.bed) -b methylation/$group.methylation.bed -wao |awk -v OFS="\t" -v cl=$class '($NF==1){print cl,"background",$(NF-1),$4}'


# Merge everything except intergenic
prefix="bTaeGut7v0.4_MT_rDNA"
echo "Group Class Type Score Trx" |sed 's/ /\t/g' >methylation/$prefix.allCpG.group.txt
for group in "macro" "micro" "dot" 
do
  for class in  "promoter" "introns" "CDS" "UTR5" "UTR3"  "lncrna"
  do
    awk -v g=$group -v OFS="\t" '{print g,$0}' methylation/G4s/$group.$class.txt |sed 's/;//g' |sed 's/"//g' >>methylation/$prefix.allCpG.group.txt
  done
done
# Merge intergenic
echo "Group Class Type Score" |sed 's/ /\t/g' >methylation/$prefix.allCpG.intergenic.group.txt
for group in  "macro" "micro" "dot"
do
  awk -v g=$group -v OFS="\t" '{print g,$0}' methylation/G4s/$group.intergenic.txt |sed 's/;//g' |sed 's/"//g' >>methylation/$prefix.allCpG.intergenic.group.txt
done

# I THINK THIS IS WRONG!!!  OR IS IT THE R CODE??
# Figure out the proportion of genes that has/doesn't have methylation
prefix="bTaeGut7v0.4_MT_rDNA"
for group in  "macro" "micro" "dot"
do
  for class in "lncrna"  "introns" "promoter"  "CDS" "UTR5" "UTR3"  
  do
      echo '#!/bin/bash
      module load bedtools/2.31.0
      intersectBed -a <(grep '$group' helpfiles/'$prefix'.groups.txt |cut -f1 | grep -f - annotation/G4s/'$class'.exclG4s.bed) -b methylation/'$group'.methylation.bed -wao |\
      sort -k4,4 |awk -v OFS="\t" -v cl='$class' -v g='$group' '"'"'{if(NR==1){trx=$4; sum=$NF}else{if(trx==$4){sum+=$NF}else{print g,cl,"background",$4,sum; trx=$4; sum=$NF}}}END{print g,cl,"background",$4,sum;}'"'"' >methylation/G4s/'$group'.'$class'.summary.txt
      intersectBed -a <(grep '$group' helpfiles/'$prefix'.groups.txt |cut -f1 | grep -f - annotation/G4s/'$class'.onlyG4s.samestrand.bed) -b methylation/'$group'.methylation.bed -wao |\
      sort -k4,4 |awk -v OFS="\t" -v cl='$class' -v g='$group' '"'"'{if(NR==1){trx=$4; sum=$NF}else{if(trx==$4){sum+=$NF}else{print g,cl,"coding",$4,sum; trx=$4; sum=$NF}}}END{print g,cl,"coding",$4,sum;}'"'"' >>methylation/G4s/'$group'.'$class'.summary.txt
      intersectBed -a <(grep '$group' helpfiles/'$prefix'.groups.txt |cut -f1 | grep -f - annotation/G4s/'$class'.onlyG4s.oppositestrand.bed) -b methylation/'$group'.methylation.bed -wao |\
      sort -k4,4 |awk -v OFS="\t" -v cl='$class' -v g='$group' '"'"'{if(NR==1){trx=$4; sum=$NF}else{if(trx==$4){sum+=$NF}else{print g,cl,"template",$4,sum; trx=$4; sum=$NF}}}END{print g,cl,"template",$4,sum;}'"'"'  >>methylation/G4s/'$group'.'$class'.summary.txt
    ' | sbatch -J $class.$group --ntasks=1 --cpus-per-task=1 --time=10:00 --out slurm/job.methylation.$group.$class.summary.%j.out
  done
done

# FOR TESTING 
intersectBed -a annotation/G4s/$class.exclG4s.bed -b methylation/$group.methylation.bed -wao |sort -k4,4 |awk -v OFS="\t" -v cl=$class -v g=$group '{if(NR==1){trx=$4; sum=$NF}else{if(trx==$4){sum+=$NF}else{print g,cl,"background",$4,sum; trx=$4; sum=$NF}}}END{print g,cl,"background",$4,sum;}' 


# Merge summary
echo "Group Class Type Trx Number" |sed 's/ /\t/g' >methylation/$prefix.summaryCpG.group.txt
for group in  "macro" "micro" "dot"
do
  for class in  "promoter" "introns" "CDS" "UTR5" "UTR3"  "lncrna"
  do
    cat methylation/G4s/$group.$class.summary.txt |sed 's/;//g' |sed 's/"//g' >>methylation/$prefix.summaryCpG.group.txt
  done
done


################################################################################
# ##################### COMBINATION OF INTRONS AND REPEATS #####################
# We want to see if it is the mini satellites in introns that are driving the
# high non-B density here.

# Overlap introns with the TR repeat file
module load bedtools/2.31.0
mkdir repeats/intronTRF
prefix="bTaeGut7v0.4_MT_rDNA"
intersectBed -a annotation/$prefix.TRF_withMers.bed -b annotation/$prefix.introns.merged.bed >repeats/intronTRF/$prefix.introns_TRF.bed
subtractBed -a annotation/$prefix.introns.merged.bed -b repeats/intronTRF/$prefix.introns_TRF.bed >repeats/intronTRF/$prefix.introns_noTRF.bed

# For each chromosome category, check enrichment in these two
for group in "macro" "micro" "dot"
do
    echo '#!/bin/bash
    module load bedtools/2.31.0
    rm -f repeats/intronTRF/'$prefix'.enrichment.'$group'.txt
    for i in "TRF" "noTRF"
    do
        len=`grep '$group' helpfiles/'$prefix'.groups.txt | \
      awk '"'"'NR==FNR{a[$1];next} ($1 in a){sum+=$3-$2}END{print sum}'"'"' - repeats/intronTRF/'$prefix'.introns_${i}.bed`
        echo "length of $i: $len"
        cat coverage/'${prefix}'.per_genome.tsv | while read -r non_b tot dens;
        do
            d=`intersectBed -a <(grep '$group' helpfiles/'$prefix'.groups.txt | \
      awk '"'"'NR==FNR{a[$1];next} ($1 in a){print}'"'"' - repeats/intronTRF/'$prefix'.introns_${i}.bed) \
      -b final_nonB/'$prefix'.${non_b}.merged.bed -wo |awk -v l=$len -v dtot=$dens '"'"'{sum+=$NF}END{d=sum/l; frac=d/dtot; print d, frac}'"'"'`
            echo '$group' $i $non_b $d >>repeats/intronTRF/'$prefix'.enrichment.'$group'.txt
        done
    done
    ' |sbatch -J $group -o slurm/job.IntronTR.$group.%j.out --requeue --ntasks=1 --cpus-per-task=1 --mem-per-cpu=4G --time=10:00:00 
done
# Merge
echo "Group Subset NonB Coverage Enrichment_gw" |sed "s/ /\t/g" >repeats/$prefix.intron_enrichment.group.tsv
for group in "macro" "micro" "dot"
do
  cat repeats/intronTRF/$prefix.enrichment.$group.txt |sed "s/ /\t/g"
done >>repeats/$prefix.intron_enrichment.group.tsv



# TEST 

intersectBed -a <(grep $group helpfiles/$prefix.groups.txt | \
      awk 'NR==FNR{a[$1];next} ($1 in a){print}' - repeats/intronTRF/$prefix.introns_${i}.bed) \
      -b final_nonB/$prefix.${non_b}.merged.bed -wo |awk -v l=$len -v dtot=$dens '{sum+=$NF}END{d=sum/l; frac=d/dtot; print d, frac}'


# ALSO MAKE FILES WITH THE NUMBER OF TRFs / noTRFs PER CATEGORY AND COMPARTMENT
echo "Group Subset Length" |sed "s/ /\t/g" >repeats/$prefix.introns_TR.lengthsummary.tsv
for group in "macro" "micro" "dot"
do
  grep "length of" slurm/job.IntronTR.$group.487*.out |tail -n2 |sed 's/://g' |awk -v OFS="\t" -v g=$group '{print g,$3,$4}' >>repeats/$prefix.introns_TR.lengthsummary.tsv
done

# AND A FILE WITH TRF CLASSES
echo "Group Repeat Number" |sed "s/ /\t/g" >repeats/$prefix.introns_TR.repeatsummary.tsv
for group in "macro" "micro" "dot"
do
  grep $group helpfiles/$prefix.groups.txt |cut -f1 |grep -f - repeats/intronTRF/$prefix.introns_TRF.bed |\
  cut -f4 |sort |uniq -c |awk -v OFS="\t" -v g=$group '{print g,$2,$1}' >>repeats/$prefix.introns_TR.repeatsummary.tsv
done

# And exact length of repeats
cut -f1,5 repeats/intronTRF/$prefix.introns_TRF.bed |sed 's/-mer//g' >repeats/$prefix.introns_TRF.lengths.tsv







######################### A LOOK INTO THE HOX-D CLUSTER ########################
# December 8, 2025

# It seems like the HOX cluster have a lot of G4s in other species, let's see
# how it looks here.  
prefix="bTaeGut7v0.4_MT_rDNA"

# HOXD cluster: chr7_mat:17695980-17776819 and chr7_pat:17653526-17734186
# Make a file with these, including 1kb flanks:
echo -e "chr7_mat\t17694980\t17777819\nchr7_pat\t17652526\t17735186" >annotation/HOXD_cluster.bed

# Intersect with g4discovery 
module load bedtools/2.31.0
intersectBed -wao -a annotation/HOXD_cluster.bed -b <(cut -f1,2,3 final_nonB/$prefix.G4.bed) |cut -f1,2,3,7 |awk -v OFS="\t" '{if(NR==0){chr=$1; s=$2; e=$3; sum=$4}else{if($1==chr && $2==s){sum+=$4}else{print chr,s,e,sum; chr=$1; s=$2; e=$3; sum=$4}}}END{print chr,s,e,sum}' | sed "/^\s*$/d" >test.density.hoxd.bed
# And with Quadron
intersectBed -wao -a annotation/HOXD_cluster.bed -b <(cut -f1,2,3 final_nonB/$prefix.G4quadron.bed) |cut -f1,2,3,7 |awk -v OFS="\t" '{if(NR==0){chr=$1; s=$2; e=$3; sum=$4}else{if($1==chr && $2==s){sum+=$4}else{print chr,s,e,sum; chr=$1; s=$2; e=$3; sum=$4}}}END{print chr,s,e,sum}' | sed "/^\s*$/d" >test.density.Quadron.hoxd.bed

# Print coverage per region 
awk -v OFS="\t" '{len=$3-$2; dens=$4/len; print $1,$2,$3,len,$4,dens}' test.density.hoxd.bed
awk -v OFS="\t" '{len=$3-$2; dens=$4/len; print $1,$2,$3,len,$4,dens}' test.density.Quadron.hoxd.bed



