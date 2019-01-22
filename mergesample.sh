#!/bin/bash

svtype="INV"
svtype_full="inversion"
outname="NB_${svtype}"

destpth="/scratch1/irri-bioinformatics/StructuralVariants/MergedSamples"
sourcepth="/scratch1/irri-bioinformatics/StructuralVariants/MergedCallers/${svtype}"

#Change run type
runtype=CLUSTER  #CLUSTER or MERGE

if [[ $runtype == "MERGE" ]]; then

mergefile=$destpth/$outname"_mergesam.bed"
if [ -f $mergefile ]; then
    rm $mergefile
fi
    
touch $mergefile

module load bedtools

#MERGE all samples into one file
for i in $(find $sourcepth -name "*inversion.txt" -printf "%f\n");
do
  sample=${i/.${svtype_full}.txt/}
  echo ${sample}
  awk -v sam="${sample}" '{ print $0"\t"sam;}' $sourcepth/$i >> $mergefile  
done
  
#sort the mergefile and cluster intervals using bedtools (initial cluster)
sort -k1,1 -k2,2n $mergefile | bedtools cluster -i - > $destpth/$outname"_mergesam_intersect.txt"

#count number of intervals per bedtools cluster; 
#wk '{cluster[$6]++; if(!($6 in name)){ name[$6]=$1"\t"$2"\t"$3;}} END{for(i=1;i<=$6;i++) print i"\t"name[i]"\t"cluster[i];}' $destpth/${outname}"_mergesam_intersect.txt" > $destpth/$outname"_mergesam_intersect_sup.txt"

#compute distribution of support for all bedtools cluster
#wk '{tmp=sprintf("%d",$4/10); bin[tmp]++;} END{for(i=0;i<303;i++) if(bin[i]>0) print i*10"\t"bin[i]; else print i*10"\t0";}' $destpth/${chrom}"_mergesam_intersect_sup.txt" > $destpth/${chrom}"_mergesam_sup_dist.txt"

#Prepare input for hierarchical clustering
awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$6;}' $destpth/${outname}"_mergesam_intersect.txt" | uniq -c - | awk -v OFS='\t' '{print $2,$3,$4,$5,$4-$3}' - > $destpth/$outname"_mergesam_intersect_uniq.txt"

#count number of unique intervals per bedtools cluster; useful on estiamting reclustering runtime in R especially for clusters with many unique intervals
awk '{cluster[$4]++; if(!($4 in name)){ name[$4]=$1"\t"$2"\t"$3;}} END{for(i=1;i<=$4;i++) print i"\t"name[i]"\t"cluster[i];}' $destpth/${outname}"_mergesam_intersect_uniq.txt" > $destpth/$outname"_mergesam_intersect_sup.txt"

fi 

if [[ $runtype == "CLUSTER" ]]; then
#Adjust cluster ID after Hclustering
awk 'BEGIN{num=0;id="x";} $1~/^c/{tmp=$5"_"substr($6,0,length($6)-1); if(id!=tmp){id=tmp;num++;} print $1"\t"$2"\t"$3"\t"$4"\t"tmp"\t"num;}' $destpth/NB_${svtype}_cluster.txt | sort -k6,6n -k1,1 -k2,2n > $destpth/$outname"_mergesam_hclust.txt"  #-1 in length remove ^M char at the end

awk 'NR==FNR{tmp=$1"\t"$2"\t"$3; map[tmp]=$6;} NR!=FNR{x=$1"\t"$2"\t"$3; if(x in map){print x"\t"$4"\t"$5"\t"map[x];} }'  $destpth/${outname}"_mergesam_hclust.txt"  $destpth/${outname}"_mergesam_intersect.txt" >  $destpth/${outname}"_mergesam_clustered.txt"

#stat_brkpt.sh and stat_spread.sh for stat

#get length distribution
#awk '{tmp=sprintf("%d",$3-$2); bin[tmp]++;} END{for(i=0;i<=1000;i++) if(bin[i]==0) print i"\t0"; else print i"\t"bin[i];}' NB_mergesam_clustered.txt > length_dist.txt

fi
