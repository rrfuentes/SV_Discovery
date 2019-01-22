#!/bin/bash

#filename = find /shared/data/SV_rfuentes/Pindel_NB/ -name "*.gz"
svtype="TRA"
svtype_full="translocation"

path="/shared/data/SV_rfuentes/GROM_NB/persample"
for i in $(find /shared/data/SV_rfuentes/GROM_NB/persample -name "*GROM_SV.txt" -printf "%f\n");
do
  samname=${i/.GROM_SV.txt/};
   awk -v typ="${svtype}" '$1~/^CTX/{
	   printf("%s\t%d\t%d\t%s;%d\n", $2,$3,$4,typ,$4-$3);
	}' $path/$i > /home/rfuentes/Parsed_SV/GROM/${samname}.${svtype_full}.txt
  echo $i 
done
