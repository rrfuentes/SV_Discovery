#!/bin/bash

svtype_full="insertion"
dirpath=/scratch1/irri-bioinformatics/StructuralVariants/Discovery/MindTheGap
FILE1=${dirpath}/qual.txt

if [ -f $FILE1 ]; then
    rm $FILE1
    touch $FILE1
else
    touch $FILE1
fi

for sample in $(less /scratch1/irri-bioinformatics/StructuralVariants/3k_subset.txt );
#or sample in $(less /scratch1/irri-bioinformatics/StructuralVariants/Analysis/INDELsubset/samples.txt );
do
  file=${dirpath}/${sample}/${sample}.insertions.fasta
  echo $sample
  awk -F'_' -v OFS="\t" -v ORS="" -v path="${dirpath}" 'BEGIN{ignore=0;} 
 	{
		if(NR%2>0){
			qual = sprintf("%d",$11);
			bin[qual]++;
			if( $11>0){ignore=0;} #ignore predictions with qual<=10
			else{ignore=1; next;}
		}else{
			if(ignore==1) next;
		} 

		if($1~/^>/)print $2,$4,$4+1,"INS;"$9";"; 
		#if($1~/^>/) print $2,$4,$4+1,$9"\t";
		else print $1"\n";
	}

	END{file=path"/qual.txt"; print bin[5]"\t"bin[10]"\t"bin[15]"\t"bin[25]"\t"bin[50]"\n" >> file;}' $file | sort -k1,1 -k2,2n > ${dirpath}/persample/${sample}.${svtype_full}.txt
 
#Compare GATK & MTG
#awk -v OFS="\t" 'NR==FNR && $4>15 && $4<21{print $0,"GATK"} NR!=FNR && $4>15 && $4<21{print $0,"MTG"}' /scratch1/irri-bioinformatics/StructuralVariants/Analysis/INDELsubset/GATK/${sample}.txt /scratch1/irri-bioinformatics/StructuralVariants/Analysis/INDELsubset/MindTheGap/filtered/${sample}.${svtype_full}.txt | sort -k1,1 -k2,2n -k4,4n -k6,6 | bedtools cluster -d 10 > /scratch1/irri-bioinformatics/StructuralVariants/Analysis/INDELsubset/Result/filtered/bin4/${sample}.txt
 
done

