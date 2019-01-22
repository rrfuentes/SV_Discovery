#!/bin/bash

#filename = find /shared/data/SV_rfuentes/Pindel_NB/ -name "*.gz"
svtype="DUP"
svtype_full="duplication"

for i in $(find /shared/data/SV_rfuentes/Pindel_NB/ -name "*.gz" -printf "%f\n");
do
  echo ${i:0:${#i}-24};
  less /shared/data/SV_rfuentes/Pindel_NB/$i | \
	awk -v typ="${svtype}" '$1~/^c/{
	    if($7=="PASS"){
                if($10!~/^0\/0/ && $10!~/\.\/\./){ 
			split($8,x,";");
			for(i in x){ 
		    		n=split(x[i],y,"="); 
		    		if(n==2)INFO[y[1]]= y[2];
			} 
			if(INFO["SVTYPE"]!~/^DUP/) next;
			
			if(INFO["SVTYPE"]=="INS") printf("%s\t%d\t%d\t%s;%s\n", $1,$2,$2+1,INFO["SVTYPE"],INFO["SVLEN"]);
			else if(INFO["SVTYPE"]=="DEL") printf("%s\t%d\t%d\t%s;%s\n", $1,$2,INFO["END"],INFO["SVTYPE"],-(INFO["SVLEN"]));
			else if(INFO["SVTYPE"]=="RPL") printf("%s\t%d\t%d\t%s;%s:%s\n", $1,$2,INFO["END"],INFO["SVTYPE"],-(INFO["SVLEN"]),INFO["NTLEN"]);
			else  printf("%s\t%d\t%d\t%s;%s\n", $1,$2,INFO["END"],INFO["SVTYPE"],INFO["SVLEN"]);
			delete INFO;
		}
	     }
	}' - > /home/rfuentes/Parsed_SV/Pindel/${i:0:${#i}-24}.${svtype_full}.txt
  echo $i 
done
