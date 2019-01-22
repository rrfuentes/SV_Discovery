#!/bin/bash

# Created by: Roven Rommel Fuentes
# Last Modified: Sept. 21, 2017
# For Merging INVERSION

svtype="INV"
svtype_full="inversion"

outpath="/scratch1/irri-bioinformatics/StructuralVariants/MergedCallers/${svtype}"
ro1="0.90"

#ORDER of breakpoint precision: Pindel,Delly,GROM,Lumpy
module load bedtools

for i in $(less /scratch1/irri-bioinformatics/StructuralVariants/3k_sample.txt);
do
   #if [ $i != "IRIS_313-8921" ]; 
   #then 
	#continue
   #fi

   samname=$i;
   dellyfile="/scratch1/irri-bioinformatics/StructuralVariants/Discovery/Delly/pertype/${samname}.$svtype_full.txt"
   pindelfile="/scratch1/irri-bioinformatics/StructuralVariants/Discovery/Pindel/pertype/${samname}.$svtype_full.txt"
   gromfile="/scratch1/irri-bioinformatics/StructuralVariants/Discovery/GROM/pertype/${samname}.$svtype_full.txt"
   lumpyfile="/scratch1/irri-bioinformatics/StructuralVariants/Discovery/Lumpy/pertype/${samname}.$svtype_full.txt"

   #SORTING breakpoints
   if [ ! -f $dellyfile ]; then 
	touch $dellyfile
   else
	sort -k1,1 -k2,2n $dellyfile -o $dellyfile
   fi
   if [ ! -f $pindelfile ]; then
	touch $pindelfile
   else
    	sort -k1,1 -k2,2n $pindelfile -o $pindelfile
   fi
   if [ ! -f $gromfile ]; then
	touch $gromfile
   else
   	sort -k1,1 -k2,2n $gromfile -o $gromfile
   fi
   if [ ! -f $lumpyfile ]; then
	touch $lumpyfile
   else
        sort -k1,1 -k2,2n $lumpyfile -o $lumpyfile
   fi

   if [ -f $outpath/${samname}_out000.txt ]; then
        rm $outpath/${samname}_out000.txt  #clear existing file
   fi
   touch $outpath/${samname}_out000.txt
   
   #MERGING callers
   if [ -f $pindelfile ] && [ -f $dellyfile ] && [ -f $gromfile ] && [ -f $lumpyfile ];
   then
        #Include all pindel and lumpy breakpoints
        awk -F"[;\t]" -v OFS='\t' '$5>=10 && $5<500000 {print}' $pindelfile >> $outpath/${samname}_out000.txt

        bedtools intersect -nonamecheck -sorted -v -f $ro1 -r  -a $lumpyfile -b $pindelfile \
           | awk -F"[;\t]" -v OFS='\t' '$5>=10 && $5<500000 {print}' >> $outpath/${samname}_out000.txt

        #Extracting breakpoints not intersecting with Pindel & Lumpy's
        #WARNING: Do not include -sorted since the list from Pindel and Lumpy is not sorted
        bedtools intersect -nonamecheck -v -f $ro1 -r  -a $dellyfile -b $outpath/${samname}_out000.txt \
           | awk -F"[;\t]" -v OFS='\t' '$5>=10 && $5<500000 {print}' > $outpath/${samname}_delly_000.txt

        bedtools intersect -nonamecheck -v -f $ro1 -r  -a $gromfile -b $outpath/${samname}_out000.txt \
           | awk -F"[;\t]" -v OFS='\t' '$5>=10 && $5<500000 {print}' > $outpath/${samname}_grom_000.txt

	#Get additional breakpoints
        bedtools intersect -nonamecheck -sorted -u -f $ro1 -r -a $outpath/${samname}_delly_000.txt -b $outpath/${samname}_grom_000.txt \
           >> $outpath/${samname}_out000.txt

   	sort -k1,1 -k2,2n -k3,3n $outpath/${samname}_out000.txt | uniq - >  $outpath/${samname}.${svtype_full}.txt
   	rm  $outpath/${samname}_*_000.txt
   	rm  $outpath/${samname}_out000.txt
   else
	echo "Missing File: ${samname}"
   fi
   echo $samname
   
done
