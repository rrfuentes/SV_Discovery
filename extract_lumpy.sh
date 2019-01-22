#!/bin/bash

svtype="DEL"
svtype_full="deletion"

dirpath="/scratch1/irri-bioinformatics/StructuralVariants/Discovery/Lumpy"

for i in $(find /scratch1/irri-bioinformatics/StructuralVariants/Lumpy/ -name "*.vcf" -printf "%f\n");
do
  samname=${i/.vcf/};
  echo $samname
  awk -v typ="${svtype}" '$1~/^c/{
                
		split($8,x,";");
		for(i in x){ 
		    	n=split(x[i],y,"="); 
		    	if(n==2)INFO[y[1]]= y[2];
		} 
		if(INFO["SVTYPE"]!=typ) next;
			
		if(INFO["SVTYPE"]=="DEL") printf("%s\t%d\t%d\t%s;%s\n", $1,$2,INFO["END"],INFO["SVTYPE"],-(INFO["SVLEN"]));
		else  printf("%s\t%d\t%d\t%s;%s\n", $1,$2,INFO["END"],INFO["SVTYPE"],INFO["SVLEN"]);
		delete INFO;
		
   }' /scratch1/irri-bioinformatics/StructuralVariants/Lumpy/${samname}/$i > $dirpath/pertype/${samname}.${svtype_full}.txt
   break
done
