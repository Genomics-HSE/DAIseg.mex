#!/bin/bash

dir=$1
CHR=$2
bed=$3
n1=$4
n2=$5
n3=$6
GP1000=$7

mex=$8
eu=$9
na=${10}
af=${11}
arch=${12}
outfilevcf=${13}
outtxt=${14}

aa=${15}


cd $1


./archaic.covering.sh ${CHR} ${bed} ${n1} ${n2} ${n3}


./new.panel.preparation.Linux.sh ${CHR} ${mex} ${eu} ${na} ${af} ${bed} ${GP1000} ${n1} ${n2} ${n3} ${outfilevcf} ${outtxt}


./new.make.obs.sh ${CHR} ${outfilevcf} ${mex} ${eu} ${na} ${af} ${arch}  ${15} ./regions/chr${CHR}.hg19.bed ${outtxt}

cd ../


