# DAIseg.mex
Highly accurate method for detecting archaic segments in the modern admixed genomes 


![Demography](https://github.com/Genomics-HSE/DAIseg.mex/blob/main/utilities/Mex.svg)

# DAIseg.mex
DAIseg.mex method is created to detect ancient introgressed segments using unadmixed outgroup population and several reference archaic genomes. 



# Pipeline briefly
Run panel.preparation.Linux.sh with samples' name files to merge 1000GP, neanderthal samples and obtain .vcf.gz file (optionally). Using .vcf.gz{.tbi} and files with samples's names to run ./make.obs.sh to make observation files. Run daiseg.mex.py to obtain archaic tracts  with the posssibility of using EM algorithm.




# Files's summary
*  __eu.txt__(European), __na.txt__(American), __af.txt__(African),  __archaic.txt__(Neanderthals)  and __mex.txt__(European), are .txt files which consist of the samples' ids of reference Africans, Neanderthals and observable Europeans written in a column
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

*  __position.txt__ first-last positions of desired region
    ```note
    start_chr end_chr
    ```

By default, the  time values are  550.000, 70.000, 55.000, 55.000 are used to make  initiall guess for the EM algorithm on Step 2. These values are good to find archqic segments but using EM algorithm allows to find short segments.


*  __all.chr22.vcf.gz{.tbi}__ files containing all reference genomes (Outgroup and Archaic) and observable samples with snps only (excluding indels, deletions etc.). The main reason of it is to avoid inconsistencies.

  
* __out.eu.txt__ and __out.eu.txt__ are  files 
    ```note
    [[t_1,t_2], [t_3,t_4], [t_5,t_6]]
    [[t'_1,t'_2], [t'_3,t'_4]]
    ...
    ...
    ```
where each two lines correspond to the one diploid sample from mex.txt with European and Native American ancestries in each file respecting the order of admixed samples.



* __ancestral.alleles.file.txt__  file with information about ancestral allels ("-1"=="no information").
  ```note
  position1 -1
  position2 A
  position3 C
  ...
  ```
  The link on the ancestral alles files based on hg19 is https://drive.google.com/file/d/1Vw-QEG9uu1trkbGHpDVXhMlbGt-RQhbN/view?usp=sharing

* __gaps.by.pos.chr.22.txt__ is file with list of gaps
  ```note
  [[a,b], [c,d]. [e,f]]
  ```

## Step 0. Merging 1000GP  and Archaic genomes
Download [1000GP panel][1] and  archaic samples  [Link1][2] and [Link2][3]. Make .txt files with samples' names  __eu.txt__, __na.txt__, __yri.txt__ and __mex.txt__ and  __archaic.txt__

Add full path to files  of 1000GP,  Altai.Neanderthal, Vindija33.19, Chagyrskaya/Okladnikova to variables NAME1000 and n1, n2, n3 in  panel.preparation.*.sh and run 

```bash
./panel.preparation.*.sh 22 mex.txt eu.txt na.txt yri.txt archaic.txt 
```
 
The resulting vcf.gz file is all.chr22.vcf.gz{.tbi}




## Step 1.  Make observations

You need  __all.chr22.vcf.gz{.tbi}__,  __outgroup.txt__, __observations.txt__, __archaic.txt__ to run  

```bash
./make.obs.sh 22 mex.txt eu.txt na.txt yri.txt archaic.txt ancestral.alleles.file.txt  
```

and to make observation files 

obs.neand.chr22.txt, obs.eu.chr22.txt, obs.na.chr22txt, obs.yri.chr22.txt

and the file with default parameters and start-end positions __par.file.txt__ (see the File's summary paragraph) and file with start-end position __pos.chr22.txt. 


## Step 2.0 Run DAI.seg without EM algorithm

```bash
python3 daiseg.mex.empty.py --gaps ./GAPS.hg19/gaps.by.pos.chr.22.txt --location pos.chr22.txt --obs_eu obs.eu.chr22.txt --obs_na obs.na.chr22.txt --obs_af obs.yri.chr22.txt --obs_archaic obs.neand.chr22.txt --EM no  --HMM_par par.file.txt --o_eu out.eu.chr22.txt --o_na out.na.chr22.txt
```

par.file.txt with basic parameters is in directory.

--gaps ./GAPS.hg19/gaps.by.pos.chr.22.txt is optional (if there are some genome gaps in your desired region)

## Step 3 (optional) Run DAI.seg using EM algorithm

par.file.txt  could be used as the initial guess for EM algorithm.

There are two possible options to estimate parameters: 


```bash
python3 daiseg.mex.empty.py  --gaps ./GAPS.hg19/gaps.by.pos.chr.22.txt --location pos.chr22.txt --obs_eu obs.eu.chr22.txt --obs_na obs.na.chr22.txt --obs_af obs.yri.chr22.txt --obs_archaic obs.neand.chr22.txt --EM yes --EM_est coal --EM_steps 5  --HMM_par par.file.txt --o_eu out.eu.chr22.5.txt --o_na out.na.chr22.5.txt
```


[1]: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz 
[2]: http://cdna.eva.mpg.de/neandertal/Vindija/VCF/
[3]: http://ftp.eva.mpg.de/neandertal/ChagyrskayaOkladnikov/
to obtain estimations only for coalescent times 

