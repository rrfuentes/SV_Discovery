#!/bin/bash
#This tries to find matches between deletion and insertion sequences
#This do not fully address TRA detection and was only used to remove deletions that may be related to TEs or TRAs

svtype_full=translocation

dirpath=/scratch1/irri-bioinformatics/StructuralVariants/MergedCallers/TRA
delpath=/scratch1/irri-bioinformatics/StructuralVariants/MergedCallers/DEL

#Only use MetaSV insertions to avoid redundancies (size>100bp)
inspath=/scratch1/irri-bioinformatics/StructuralVariants/Discovery/MetaSV/persample

#Use TRA signals from Delly
delly=/scratch1/irri-bioinformatics/StructuralVariants/Discovery/Delly/pertype

reference=/scratch1/irri-bioinformatics/StructuralVariants/Benchmark/IRGSP-1.0_genome.fa

module load blast/2.4.0
module load bedtools/2.26.0

miscfolder=$dirpath/misc
if [ -d "$miscfolder" ]; then
  rm -r $miscfolder
fi
mkdir $miscfolder

for sample in $(less /scratch1/irri-bioinformatics/StructuralVariants/3k_subset.txt );
do
  file1=$delpath/${sample}.deletion.txt
  #Only use MetaSV insertions to avoid redundancies (size>100bp)
  file2=$inspath/${sample}.insertion.txt
  file3=$delly/${sample}.${svtype_full}.txt
  miscfile=$miscfolder/${sample}
  mkdir $miscfile
  echo $sample

  awk -F"[\t;]" -v OFS="\t" '$5>=100{print $1,$2,$3;}' $file1 | bedtools getfasta -fi $reference -bed - -fo $miscfile/${sample}.delseq.fa

  awk -F"[\t;]" -v OFS="\t" -v ORS="" '$5>=100 && NF==6{print ">"$1":"$2"-"$3"\n"$6"\n";}' $file2 > $miscfile/${sample}.insseq.fa

  makeblastdb -in $miscfile/${sample}.insseq.fa -dbtype nucl

  blastn -db $miscfile/${sample}.insseq.fa -query $miscfile/${sample}.delseq.fa -outfmt "6 qseqid sseqid qlen slen pident nident" | awk -F"\t" 'function max(a,b) {return a > b ? a : b} {if($6/max($3,$4)>0.90) print $0"\t"($6/max($3,$4))*100;}' - > $miscfile/${sample}.blastn_TRA.tsv

  awk -F"[\t:-]" -v OFS="\t" '{print $1,$2,$3,$7,$4,$5,$6,$8,$11;}' $miscfile/${sample}.blastn_TRA.tsv | sort -k1,1 -k2,2n -k5,5 -k6,6n | uniq > $dirpath/${sample}.${svtype_full}.txt
  
done
