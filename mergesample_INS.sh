#!/bin/bash

svtype="INS"
svtype_full="insertion"
outname="NB_${svtype}"

destpth=/scratch1/irri-bioinformatics/StructuralVariants/MergedSamples
sourcepth="/scratch1/irri-bioinformatics/StructuralVariants/MergedCallers/${svtype}"

mergefile=$destpth/$outname"_mergesam.bed"
if [ -f $mergefile ]; then
    rm $mergefile
fi

touch $mergefile

module load bedtools

for i in $(less /scratch1/irri-bioinformatics/StructuralVariants/3k_subset.txt );
do
  sample=${i}.insertion.txt
  echo $i

   awk -F"[\t;]" -v sam="${i}" 'BEGIN{pos="";} $5>=5{
	if(pos==""){ pos=$1"\t"$2"\t"$3; seq="";} 
	cur=$1"\t"$2"\t"$3; 
	if(cur==pos){ if($6!="" && seq!="") seq=seq","$6; else if($6!="") seq=$6; } 
	else{ print pos"\t"seq"\t"sam; pos=cur; seq=$6;}
   }'  $sourcepth/$sample >> $mergefile

  #f [ "$i" == "CX579" ]; then
	#break 
  #i
done

#sort the mergefile and cluster intervals using bedtools (initial cluster)
sort -k1,1 -k2,2n $mergefile | bedtools cluster -d 10 -i - > $destpth/$outname"_mergesam_clustered.txt"

#get list of insertion sites
awk '{print $1"\t"$2"\t"$6;}' $destpth/$outname"_mergesam_clustered.txt" | uniq > $destpth/$outname"_mergesam_INSsites.txt"

#awk -F"[\t;]" '{cur=$1"\t"$2; if(cur!=pos){pos=$1"\t"$2; bin[$5]++;}} END{for(i=1;i<100;i++){ if(bin[1]>0) print i"\t"bin[i]; else  print i"\t0";}}' ../MergedCallers/* > INS_length_dist.txt
