#!/bin/bash



CHR=$1
mex=$2
eu=$3
na=$4
af=$5
arch=$6
panelfinal=all.chr${CHR}.vcf.gz
for i in ${eu} ${na} ${af} ${arch} 
do
	bcftools query -S $i  -f '%POS\n'   ${panelfinal} > $i.positions.txt # print number of postions in separate text file
	bcftools query -S $i  -f '[%GT]\n'  ${panelfinal}  |sed 's/[^0]//g'  | awk '{ print length }'> $i.0.txt # считаем количество ноликов в строке
	bcftools query -S $i  -f '[%GT]\n'  ${panelfinal} |sed 's/[^1]//g'  | awk '{ print length }'> $i.1.txt # считаем количество единичек в строке
	paste  $i.0.txt $i.1.txt > $i.al.spec.txt #join 

	awk '{
		if ($1>"0"&& $2>"0"){print "1\t1"}
		else 
			if ($1>"0"&& $2=="0"){print "1\t-1"}
			else {print "-1\t1" }
		}' $i.al.spec.txt > $i.spec.txt	
done 

for i in ${eu} ${na} ${af} ${arch} 
do
paste $i.positions.txt $i.spec.txt  > chr${CHR}.$i.reference.txt
done

for i in ${eu} ${na} ${af} ${arch} 
do
rm $i.spec.*
rm $i.al.*
rm $i.0.*
rm $i.1.*
rm $i.positions.*
done


bcftools query -S ${mex}  -f '[%GT ]\n'  ${panelfinal} |sed  's/|/ /g' >  obs.chr${CHR}.ingroup.txt
rm samples.*

python3 obs.py ${CHR} chr${CHR}.${eu}.reference.txt chr${CHR}.${na}.reference.txt chr${CHR}.${af}.reference.txt chr${CHR}.${arch}.reference.txt obs.chr${CHR}.ingroup.txt 

for i in ${eu} ${na} ${af} ${arch} 
do
rm chr${CHR}.$i.reference.txt
done
rm obs.chr${CHR}.ingroup.txt














