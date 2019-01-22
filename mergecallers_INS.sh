#!/bin/bash

dirpath=/scratch1/irri-bioinformatics/StructuralVariants/MergedCallers/INS
metapath=/scratch1/irri-bioinformatics/StructuralVariants/Discovery/MetaSV
mindpath=/scratch1/irri-bioinformatics/StructuralVariants/Discovery/MindTheGap

for sample in $(less /scratch1/irri-bioinformatics/StructuralVariants/3k_subset.txt );
do
  file=persample/${sample}.insertion.txt
  echo $sample

  awk -F"[\t;]" -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}' $metapath/$file $mindpath/$file | sort -k1,1 -k2,2n -k5,5n | awk -v OFS="\t" '{if($6!="")print $1,$2,$3,$4";"$5";"$6; else print $1,$2,$3,$4";"$5;}' > $dirpath/${sample}.insertion.txt

done

