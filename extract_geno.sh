#!/bin/bash
svtype_s=DEL
input1=/scratch1/irri/irri-bioinformatics/StructuralVariants/Discovery/Pindel
input2=/scratch1/irri/irri-bioinformatics/StructuralVariants/Discovery/Delly
genofile=/scratch1/irri/irri-bioinformatics/StructuralVariants/Discovery/Genotypes
outdir=/scratch1/irri/irri-bioinformatics/StructuralVariants/MergedSamples/Genotyped
sam=/scratch1/irri/irri-bioinformatics/StructuralVariants/MergedCallers/$svtype_s

if [ $svtype_s == "DEL" ]
then
    svtype=deletion
elif [ $svtype_s == "INS" ]
then
   svtype=insertion
elif [ $svtype_s == "INV" ]
then
   svtype=inversion
elif [ $svtype_s == "DUP" ]
then
   svtype=duplication
else
   echo "ERROR: Incorrect SV type."
   exit 1
fi

for i in `less /scratch1/irri/irri-bioinformatics/StructuralVariants/3k_sample.txt`;
do
  less $input1/${i}.realigned.pindel.vcf.gz | awk 'BEGIN{tmp1="END"; tmp2="SVTYPE";} $1!~/^#/{n=split($8,a,";"); for(i=1;i<=n;i++){m=split(a[i],b,"="); if(m>1){INFO[b[1]]=b[2];}} split($10,c,":"); print $1"\t"$2"\t"INFO[tmp1]"\t"INFO[tmp2]"\t"c[1]; }' - > $genofile/${i}.Pindel

  less $input2/${i}.realigned.delly.*.vcf.gz | awk 'BEGIN{tmp1="END"; tmp2="SVTYPE";} $1!~/^#/{n=split($8,a,";"); for(i=1;i<=n;i++){m=split(a[i],b,"="); if(m>1){INFO[b[1]]=b[2];}} split($10,c,":"); print $1"\t"$2"\t"INFO[tmp1]"\t"INFO[tmp2]"\t"c[1]; }' - > $genofile/${i}.Delly

  cat $genofile/${i}.Pindel $genofile/${i}.Delly | awk 'NR==FNR{tmp=$1"\t"$2"\t"$3; geno[tmp]=$5; type[tmp]=$4;} NR!=FNR{tmp=$1"\t"$2"\t"$3; if(tmp in geno){print tmp"\t"$4"\t"geno[tmp];}else print tmp;}' - $sam/${i}.${svtype}.txt > $outdir/${i}.${svtype}.SVgeno.txt
  echo $i" "$sam/${i}.${svtype}.txt
done

