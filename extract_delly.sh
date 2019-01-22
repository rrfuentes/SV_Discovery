#!/bin/bash

dirpath="/scratch1/irri-bioinformatics/StructuralVariants/Discovery/Delly"

idx=0;

for i in $(find /scratch1/irri-bioinformatics/StructuralVariants/Discovery/Delly/ -name "*.translocation.vcf.gz" -printf "%f\n");
do
  samname=${i/realigned.delly./};
  samname=${samname/.vcf.gz/.txt};
  less $dirpath/$i | \
	awk '$1~/^c/{ 
	    if($7=="PASS" && $10!~/^0\/0/ && $10!~/^\.\/\./){
		split($8,x,";"); 
		for(i in x){ 
		    n=split(x[i],y,"="); 
		    if(n==2)INFO[y[1]]= y[2];
		} 
		if(INFO["SVTYPE"]=="INS") printf("%s\t%d\t%d\t%s;%s\n", $1,$2,INFO["END"],INFO["SVTYPE"],INFO["INSLEN"]);
		else if(INFO["SVTYPE"]=="TRA") printf("%s\t%d\t%s\t%d\t%s\n", $1,$2,INFO["CHR2"],INFO["END"],INFO["SVTYPE"]);
		else  printf("%s\t%d\t%d\t%s;%s\n", $1,$2,INFO["END"],INFO["SVTYPE"],INFO["END"]-$2);
		delete INFO;
	     }
	}' - > $dirpath/pertype/$samname
   echo $i
done

