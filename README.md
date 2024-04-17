# DAIseg.mex
Highly accurate method for detecting archaic segments in the modern admixed genomes 


![Demography](https://github.com/Genomics-HSE/DAIseg.mex/blob/main/utilities/Mex.svg)

# DAIseg.mex
DAIseg.mex method is created to detect ancient introgressed segments using unadmixed outgroup population and several reference archaic genomes. 

__Input__: .vcf.gz{.tbi} file where merging all Neanderthal, Europeans, Native Americans, Africans and  ingroup observable samples together and five .txt files with ids of each group.

__Output__: .txt file where each two lines corresponds to the one admixed(mexican) sample and is  2d-array of tracts  by states: Modern Europeans ancestry, Neanderthal ancestry camed through Europeans, American ancestry, Neanderthal ancestry camed throug Native americans, African.


# Pipeline briefly
0. (optionally) Run __panel.preparation.mex.sh__ with samples' name files to merge 1000GP, neanderthal samples and obtain .vcf.gz file.
1. Using .vcf.gz{.tbi} and files with samples's names to run __./make.obs.sh__ to make observation files.
3. Run __daiseg.mex.py__ to obtain archaic tracts  with the posssibility of using EM algorithm.




# Files's summary
*  __outgroup.txt__(Africa), __archaic.txt__(Neanderthals)  and __obs.samples.txt__(European), are .txt files which consist of the samples' ids of reference Africans, Neanderthals and observable Europeans written in a column
```note
NA18484
NA18489
GM19129
```


*  __par.file.txt__
```note
29 # years per generation
1.25e-08    #mutation rate μ
1e-08    #recombination rate
1000    #window size
start end    #position of the first SNP in .vcf file
t_arch^c    #Coalescent time of AMH and Neanderthals
t_split^c    #Coalescent time out of Africa
t_intr^c    #coalescent time of archaic segments in modern genome with neanderthal samples
t_ea^c # coalescent time of Europeans and Asians
t_mex^c # modern coalescent time
t_intr #introgression time
t_mex # time of modern admixture
0.025    #admixture proportion of archaic introgression
0.45 # portion of European ancestry
0.45 # portion of American ancestry
0.1 # Portion of African ancestry
```

By default, the  time values are  550.000, 70.000, 55.000, 55.000 are used to make  initiall guess for the EM algorithm on Step 2. These values are good to find archqic segments but using EM algorithm allows to find short segments.


*  __all.chr22.vcf.gz{.tbi}__ files containing all reference genomes (Outgroup and Archaic) and observable samples with snps only (excluding indels, deletions etc.). The main reason of it is to avoid inconsistencies.
  
* __output.txt__ is a  file 
```note
[[[t_1,t_2], [t_3,t_4], [t_5,t_6]], [[t_2+1, t_3-1], [t_4+1, t_5-1]]]
[[[t'_1,t'_2], [t'_3,t'_4], [t'_5,t'_6]], [[t'_2+1, t'_3-1], [t'_4+1, t'_5-1]]]
...
...
```
where each two lines correspond to the one diploid sample from obs.samples.txt.





## Step 0. Merging 1000GP  and Archaic genomes
Download 1000GP panel 
>http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz 

and  archaic samples 
>http://cdna.eva.mpg.de/neandertal/Vindija/VCF/
>http://ftp.eva.mpg.de/neandertal/ChagyrskayaOkladnikov/

Make .txt files with samples' names  __eu.txt__, __na.txt__, __yri.txt__ and __mex.txt__ and  __archaic.txt__

Add full path to files  of 1000GP and three neanderthals to variables __$NAME1000__ and __$n1, $n2, $n3__ in  __panel.preparation.sh__ , change $CHR variable and run 
>./panel.preparation.sh
 
The resulting vcf.gz file is __all.chr22.vcf.gz{.tbi}__




## Step 1.  Make observations

You need  __all.chr22.vcf.gz{.tbi}__,  __outgroup.txt__, __observations.txt__, __archaic.txt__ to run  

>__./make.obs.sh__

and to make observation files __obs.neand.txt__, __obs.eu.txt__, __obs.na.txt__, __obs.yri.txt__ and the file with default parameters and start-end positions __par.file.txt__ (see the File's summary paragraph). 


## Step 2.0 Run DAI.seg without EM algorithm
>  python daiseg.mex.py  --EM no --HMM_par par.file.txt  --o output.tracts.txt

where file par.file.txt was created on the previous Step. 

## Step 3 (optional) Run DAI.seg using EM algorithm

par.file.txt obtained on the Step 1 could be used as the initial guess for EM algorithm.

There are two possible options to estimate parameters: 

> python daiseg.mex.py --EM yes --HMM_par par.file.txt --o output.txt

to obtain estimations only for coalescent times 

