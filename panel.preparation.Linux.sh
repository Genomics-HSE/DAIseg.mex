#!/bin/bash


#we download 1000GP and 3 ancient genomes to the folders 1000GP, neand/altai, neand/33.19, 
#full list of samples of ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz is in file all.samples.txt



#change the following names and directories
comp=scglab


CHR=$1
mex=$2
eu=$3
na=$4
af=$5
arch=$6

#change the following names and directories
NAME1000=/media/${comp}/T7/Work/Data/1000GP/${CHR}/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz #name of the 1000GP vcf files
n1=/media/${comp}/T7/Work/data/neand/33.19/chr${CHR}_mq25_mapab100.vcf.gz
n2=/media/${comp}/T7/Work/data/neand/altai/chr${CHR}_mq25_mapab100.vcf.gz
n3=/media/${comp}/T7/Work/data/neand/ChagyrskayaOkladnikov/split.${CHR}.vcf.gz



DIR=CHR${CHR}


panelfinal=all.chr${CHR}.vcf.gz
result0=1000GP.chr${CHR}final.vcf.gz
temporary=temp.panel.chr${CHR}.vcf.gz



cat ${mex} ${eu} ${na} ${af} > samples.for.hmm.txt




bcftools query -f '%POS\n' ${NAME1000}|sort|uniq -cd   > dublicated.snps.txt
sed -i 's/^ *//' dublicated.snps.txt
sed -i 's/.* //' dublicated.snps.txt 
sed -i -e 's/^/${CHR}\t/' dublicated.snps.txt 
bcftools view -v snps -T ^dublicated.snps.txt -S samples.for.hmm.txt  ${NAME1000} -Oz -o ${temporary}
tabix -p vcf ${temporary}
echo "removed dublicated snps and extract reference and observable samples"


bcftools query -f '%CHROM\t%POS\n' ${temporary} > positions.chr${CHR}.txt

rm dublicated.snps.txt







bcftools view -T positions.chr${CHR}.txt ${n1} -Oz -o filtered.snps.chr${CHR}.1.vcf.gz
bcftools view -T positions.chr${CHR}.txt ${n2} -Oz -o filtered.snps.chr${CHR}.2.vcf.gz
bcftools view -T positions.chr${CHR}.txt ${n3} -Oz -o filtered.snps.chr${CHR}.3.vcf.gz


tabix -p vcf filtered.snps.chr${CHR}.1.vcf.gz
tabix -p vcf filtered.snps.chr${CHR}.2.vcf.gz
tabix -p vcf filtered.snps.chr${CHR}.3.vcf.gz

bcftools merge filtered.snps.chr${CHR}.1.vcf.gz filtered.snps.chr${CHR}.2.vcf.gz filtered.snps.chr${CHR}.3.vcf.gz  -Oz -o merged.chr${CHR}.ancient.vcf.gz
tabix -p vcf merged.chr${CHR}.ancient.vcf.gz
rm filtered.snps.*

echo "we've glued ancient genomes"


bcftools query -f '%CHROM\t%POS\n' merged.chr${CHR}.ancient.vcf.gz > common.positions.txt
bcftools view -T common.positions.txt ${temporary} -Oz -o ${result0}
tabix -p vcf ${result0}

rm ${temporary}
rm common.positions.txt
 echo " merged 1000GP and ancient"


bcftools merge ${result0} merged.chr${CHR}.ancient.vcf.gz  -Oz -o merged.chr${CHR}.vcf.gz
tabix -p vcf merged.chr${CHR}.vcf.gz

bcftools view -m 2 -M 2 -v snps merged.chr${CHR}.vcf.gz -Oz -o ${panelfinal}
tabix -p vcf ${panelfinal}

echo "removed all bad snps"

rm merged.chr${CHR}.vcf.*
rm ${result0}
rm ${result0}.tbi
rm merged.chr${CHR}.ancient.vcf.*

rm ${temporary}.tbi










