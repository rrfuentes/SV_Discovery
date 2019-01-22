#!/bin/bash

svtype_full="insertion"

dirpath=/scratch1/irri-bioinformatics/StructuralVariants/Discovery/MetaSV
for sample in $(less /scratch1/irri-bioinformatics/StructuralVariants/3k_subset.txt );
do
  file=${dirpath}/${sample}/variants.vcf.gz
  echo $sample
  less $file | \
        awk '$1~/^c/{ 
		if($7=="PASS" && $10!="0/0"){
		   split($8,x,";");
                   for(i in x){
                   	n=split(x[i],y,"=");
                        if(n==2)INFO[y[1]]= y[2];
                   }
		   seq="INSERTION_SEQUENCE";
		   if(INFO["SVTYPE"]=="INS"){ #exclude non-insertion
			if(seq in INFO) printf("%s\t%d\t%d\t%s;%s;%s\n", $1,$2,$2+1,INFO["SVTYPE"],INFO["SVLEN"],INFO[seq]);
			else printf("%s\t%d\t%d\t%s;%s\n", $1,$2,$2+1,INFO["SVTYPE"],INFO["SVLEN"]);
		   }
		   delete INFO;
		}
	     }' - > $dirpath/persample/${sample}.${svtype_full}.txt 
done

