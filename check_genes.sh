#!/bin/bash

#THIS script checks whether a gene is totally deleted or partially(>50%) deleted
#or translocated to other genome region

dirpath=/scratch1/irri-bioinformatics/StructuralVariants/Analysis/Genes
delpath=/scratch1/irri-bioinformatics/StructuralVariants/MergedCallers/DEL
trapath=/scratch1/irri-bioinformatics/StructuralVariants/MergedCallers/TRA
genelist=/scratch1/irri-bioinformatics/StructuralVariants/Analysis/Genes/genes.bed
geneinfo=/scratch1/irri-bioinformatics/StructuralVariants/LOCUS_Info_MSU.txt

module load blast/2.4.0
module load bedtools/2.26.0

miscfolder=$dirpath/DELorTRA_Genes
if [ -d "$miscfolder" ]; then
  rm -r $miscfolder
fi
mkdir $miscfolder

statfile=DEL_TRA_gene_stat.txt
rm $statfile
touch $statfile

awk -v OFS="\t" 'NR>1{print $2,$3,$8,$9;}' $geneinfo | sort -k1,1 -k2,2 | awk '{if(tmp!=$1){print $0; tmp=$1;}}' - > ${miscfolder}/tmp4

for sample in $(less /scratch1/irri-bioinformatics/StructuralVariants/3k_sample.txt );
do
  file1=$delpath/${sample}.deletion.txt
  file2=$trapath/${sample}.translocation.txt
  #iscfile=$miscfolder/${sample}
  #kdir $miscfile
  echo $sample

  awk '{print $1"\t"$2"\t"$3;}' $file1 | bedtools intersect -f 1 -u -nonamecheck -a $genelist -b - > ${miscfolder}/tmp1

  if [ -f $file2 ]; then
 	awk '{print $1"\t"$2"\t"$3;}' $file2 | bedtools intersect -f 1 -u -nonamecheck -a $genelist -b - > ${miscfolder}/tmp2
        #exclude DEL that are already identified as TRA
  	awk 'NR==FNR{gene[$4];} 
	     NR!=FNR{if($4 in gene){ print $0"\tTranslocated";}
		else{ print $0"\tDeleted";}
  	     }' ${miscfolder}/tmp2 ${miscfolder}/tmp1 > ${miscfolder}/tmp3
  else
	#just write last row to sample without INS dataset or those excluded of the 562 samples
	awk 'END{print $0"\tDeleted";}' ${miscfolder}/tmp1 > ${miscfolder}/tmp3
  fi

  awk 'NR==FNR{TE[$1]=$3; Expr[$1]=$4;} NR!=FNR{print $0"\t"TE[$4]"\t"Expr[$4];}' ${miscfolder}/tmp4  ${miscfolder}/tmp3 > ${miscfolder}/$sample.txt  
  awk -v sam="${sample}" -v OFS="\t" 'BEGIN{count1=count2=count3=count4=0;}{if($5=="Deleted") count1++; if($5=="Translocated") count2++; if($6=="Y") count3++; if($7=="Y") count4++;} END{print sam,count1+count2,count2,count3,count4;}' ${miscfolder}/$sample.txt >> $statfile
done

rm ${miscfolder}/tmp1
rm ${miscfolder}/tmp2
rm ${miscfolder}/tmp3

subgrp=/scratch1/irri-bioinformatics/StructuralVariants/SUBGROUP_9_3k.txt
mergedfile=${miscfolder}/DELETEDgenes_merged.txt
rm $mergedfile
touch $mergedfile

for sample in $(less /scratch1/irri-bioinformatics/StructuralVariants/3k_sample.txt );
do
        awk -v sam="${sample}" 'NR==FNR{grp[$2]=$3;} NR!=FNR{print $0"\t"sam"\t"grp[sam];}' $subgrp ${miscfolder}/$sample.txt >> $mergedfile
done

#FIND deleted genes exclusive to certain subgroup
awk -F"\t" 'NR==FNR{loc1[$1]=0; loc2[$1]=0; loc3[$1]=0; loc4[$1]=0; loclist[NR]=$1;}
     NR!=FNR && $6=="N"{
                if($9=="japonica") loc1[$4]=1;
                if($9=="indica") loc2[$4]=1;
                if($9=="aromatic") loc3[$4]=1;
                if($9=="aus") loc4[$4]=1;
        } END{for(i=1;i<length(loclist);i++){
                if(loc1[loclist[i]]>0 || loc2[loclist[i]]>0 || loc3[loclist[i]]>0 || loc4[loclist[i]]>0)
                        print loclist[i],loc1[loclist[i]],loc2[loclist[i]],loc3[loclist[i]],loc4[loclist[i]];}
        }' ${miscfolder}/tmp4 $mergedfile >  ${miscfolder}/DELETEDgenes_subgrp.txt

rm ${miscfolder}/tmp4
