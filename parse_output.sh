#!/bin/bash

dellypath=""
grompath=""
lumpypath=""
pindelpath=""
metapath=""
mindpath=""
outpath=""
svtype_s=""
helpmode=0

while getopts d:g:l:p:m:n:o:t:h option
do
  case "${option}"
  in
    d) dellypath=${OPTARG};;
    g) grompath=${OPTARG};;
    l) lumpypath=${OPTARG};;
    p) pindelpath=${OPTARG};;
    m) metapath=${OPTARG};;
    n) mindpath=${OPTARG};;
    o) outpath="${OPTARG}";;
    t) svtype_s="${OPTARG}";;
    h) helpmode=1;;
  esac
done

if [[ $svtype_s == "" || $outpath == "" || $helpmode == 1 || $dellypath == "" || $grompath == "" || $lumpypath == "" || $pindelpath == "" || $metapath == "" || $mindpath == "" ]]
then
    echo -e "\nThis script parses callers' output."
    echo -e "  -d\tDelly directory"
    echo -e "  -g\tGROM directory"
    echo -e "  -l\tLumpy directory"
    echo -e "  -p\tPindel directory"
    echo -e "  -m\tMetaSV directory (INS)"
    echo -e "  -n\tMindTheGap directory (INS)"
    echo -e "  -o\tOutput directory"
    echo -e "  -t\tSV Type (DEL/INS/DUP/INV)"
    echo -e "  -h\tHelp\n"
    exit 1
fi

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

if [ $dellypath != "" ]
then

  if [ ! -d "$outpath/Delly" ]; then 
   	mkdir $outpath/Delly 
  else
    	echo "Warning: Directory exists (${outpath}/Delly)."
  fi

  for i in $(find $dellypath -name "*.${svtype}.vcf.gz" -printf "%f\n");
  do
  samname=${i/.vcf.gz/.txt};
  less $dellypath/$i | \
	awk '$1!~/^#/{ 
	    if($7=="PASS" && $10!~/^0\/0/ && $10!~/^\.\/\./){
		split($8,x,";"); 
		for(i in x){ 
		    n=split(x[i],y,"="); 
		    if(n==2)INFO[y[1]]= y[2];
		} 
		if(INFO["SVTYPE"]=="INS") printf("%s\t%d\t%d\t%s;%s\n", $1,$2,INFO["END"],INFO["SVTYPE"],INFO["INSLEN"]);
		else if(INFO["SVTYPE"]=="TRA") printf("%s\t%d\t%s\t%d\t%s\n", $1,$2,INFO["CHR2"],INFO["END"],INFO["SVTYPE"]);
		else  printf("%s\t%d\t%d\t%s;%d\n", $1,$2,INFO["END"],INFO["SVTYPE"],INFO["END"]-$2);
		delete INFO;
	     }
	}' - > $outpath/Delly/$samname
   echo -e "Parsing Delly results for ${i}."
   done

fi

if [ $grompath != "" ]
then
   if [ ! -d "$outpath/GROM" ]; then
     	mkdir $outpath/GROM
   else
   	echo "Warning: Directory exists (${outpath}/GROM)."
   fi

   for i in $(find $grompath -name "*.vcf" -printf "%f\n");
   do
	samname=${i/.vcf/.${svtype}.txt};
	awk -v typ="$svtype_s" '$1!~/^#/ && $5~typ{
		gsub("END=","",$8); 
		gsub(/[<>\t]/,"",$5); 
		if($5~/INS/) printf("%s\t%d\t%d\t%s;\n", $1,$2,$8+1,$5);
		else  printf("%s\t%d\t%d\t%s;%d\n", $1,$2,$8,$5,$8-$2);
	}' $grompath/$i > $outpath/GROM/$samname
   echo -e "Parsing GROM results for ${i}."
   done

fi

if [ $lumpypath != "" ]
then
   if [ ! -d "$outpath/Lumpy" ]; then
        mkdir $outpath/Lumpy
   else
        echo "Warning: Directory exists (${outpath}/Lumpy)."
   fi

   for i in $(find $lumpypath -name "*.vcf" -printf "%f\n");
   do
	samname=${i/.vcf/.${svtype}.txt};
	awk -v typ="${svtype_s}" '$1!~/^#/{
                split($8,x,";");
                for(i in x){
                        n=split(x[i],y,"=");
                        if(n==2)INFO[y[1]]= y[2];
                }
                if(INFO["SVTYPE"]!=typ) next;

                if(INFO["SVTYPE"]=="DEL") printf("%s\t%d\t%d\t%s;%s\n", $1,$2,INFO["END"],INFO["SVTYPE"],-(INFO["SVLEN"]));
                else  printf("%s\t%d\t%d\t%s;%d\n", $1,$2,INFO["END"],INFO["SVTYPE"],INFO["SVLEN"]);
                delete INFO;

   	}' $lumpypath/$i > $outpath/Lumpy/$samname
   echo -e "Parsing Lumpy results for ${i}."
   done

fi

if [ $pindelpath != "" ]
then
   if [ ! -d "$outpath/Pindel" ]; then
        mkdir $outpath/Pindel
   else
        echo "Warning: Directory exists (${outpath}/Pindel)."
   fi

   for i in $(find $pindelpath -name "*.vcf.gz" -printf "%f\n");
   do
        samname=${i/.vcf.gz/.${svtype}.txt};
	less $pindelpath/$i | \
        awk -v typ="${svtype_s}" '$1!~/^#/{
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
        }' - > $outpath/Pindel/$samname 
   echo -e "Parsing Pindel results for ${i}."
   done

fi

if [[ $metapath != "" && $svtype_s == "INS" ]]
then
   if [ ! -d "$outpath/MetaSV" ]; then
        mkdir $outpath/MetaSV
   else
        echo "Warning: Directory exists (${outpath}/MetaSV)."
   fi

   for i in $(find $metapath -name "*.vcf.gz" -printf "%f\n");
   do
	samname=${i/.vcf.gz/.${svtype}.txt};
	less $metapath/$i | \
        awk '$1!~/^#/{
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
        }' - > $outpath/MetaSV/$samname
   echo -e "Parsing MetaSV results for ${i}."
   done

fi

if [[ $mindpath != "" && $svtype_s == "INS" ]]
then
   if [ ! -d "$outpath/MindTheGap" ]; then
        mkdir $outpath/MindTheGap
   else
        echo "Warning: Directory exists (${outpath}/MindTheGap)."
   fi

   for i in $(find $mindpath -name "*.insertions.fasta" -printf "%f\n");
   do
	samname=${i/.insertions.fasta/.${svtype}.txt};
	awk -F'_' -v OFS="\t" -v ORS="" 'BEGIN{ignore=0;}
        {
                if(NR%2>0){
                        if( $11>10){ignore=0;} #ignore predictions with qual<=10
                        else{ignore=1; next;}
                }else{
                        if(ignore==1) next;
                }

                if($1~/^>/)print $2,$4,$4+1,"INS;"$9";";
                else print $1"\n";
        }' $mindpath/$i | sort -k1,1 -k2,2n > $outpath/MindTheGap/$samname

   echo -e "Parsing MindTheGap results for ${i}."
   done

fi
