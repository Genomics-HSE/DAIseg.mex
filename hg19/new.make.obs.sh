#!/bin/bash

CHR=$1
panelfinal=$2

mex=$3
eu=$4
na=$5
af=$6
arch=$7
aa=$8
bed=$9

outtxt=${10}


bcftools query -f '%POS\t%REF\t%ALT\n' all.chr22.vcf.gz>positions.chr${CHR}.txt

bcftools query -S ${eu} -f '[%GT]\n' all.chr22.vcf.gz|awk '{for(i=1;i<=length($0);i++){a[substr($0,i,1)]=1} for(i in a){printf("%s",i)} print "";delete a}'|sed 's/|//g' |sed 's|/||g' |sed -e 's/\.//g'|sed 's/^$/\./' >al.chr${CHR}.${eu}
bcftools query -S ${na} -f '[%GT]\n' all.chr22.vcf.gz|awk '{for(i=1;i<=length($0);i++){a[substr($0,i,1)]=1} for(i in a){printf("%s",i)} print "";delete a}'|sed 's/|//g' |sed 's|/||g' |sed -e 's/\.//g'|sed 's/^$/\./' >al.chr${CHR}.${na}
bcftools query -S ${af} -f '[%GT]\n' all.chr22.vcf.gz|awk '{for(i=1;i<=length($0);i++){a[substr($0,i,1)]=1} for(i in a){printf("%s",i)} print "";delete a}'|sed 's/|//g' |sed 's|/||g' |sed -e 's/\.//g'|sed 's/^$/\./' >al.chr${CHR}.${af}

bcftools query -S ${arch} -f '[%GT]\n' all.chr22.vcf.gz|awk '{for(i=1;i<=length($0);i++){a[substr($0,i,1)]=1} for(i in a){printf("%s",i)} print "";delete a}'|sed 's|/||g' |sed 's/|//g'|sed -e 's/\.//g'|sed 's/^$/\./' >al.chr${CHR}.${arch}

bcftools query -S ${mex}  -f '[%GT ]\n'  ${panelfinal} |sed  's/|/ /g'|sed 's|/| |g' >  obs.chr${CHR}.${mex}

paste  positions.chr${CHR}.txt al.chr${CHR}.${eu} al.chr${CHR}.${na} al.chr${CHR}.${af} al.chr${CHR}.${arch} obs.chr${CHR}.${mex}> 3.chr${CHR}.txt


printf '#POSITIONS\t#REF\t#ALT\t#OUTGROUP\t#ARCHAIC\t#OBSERVATIONS\n'>header.txt
cat header.txt 3.chr${CHR}.txt > temp.allels.ref.and.obs.chr${CHR}.txt

rm header.txt
paste positions.chr${CHR}.txt obs.chr${CHR}.${mex} > obs.chr${CHR}.txt

rm al.*
rm positions.chr${CHR}.txt 
rm obs.chr${CHR}.ingroup.txt
rm 3.chr${CHR}.txt


python3 obs2.py ${CHR}  temp.allels.ref.and.obs.chr${CHR}.txt ${aa} ${bed} ${outtxt}

rm temp.allels.ref.and.obs.chr${CHR}.txt














