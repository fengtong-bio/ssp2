Getting Started
==========
```Bash
git clone https://github.com/fengtong2018/ssp2
cd ssp2 && chmod 755 *pl
#Quick start command
perl step0.pl ctg.list > ssp2.sh

#Step-by-step command
#Using kmer_count cut clean reads to K-mer
kmer_count -l 21 -f FQ -o FEMALE.out -i FEMALE_clean_R1.fq -i FEMALE_clean_R2.fq -G 10G
kmer_count -l 21 -f FQ -o MALE.out -i MALE_clean_R1.fq -i MALE_clean_R2.fq -G 10G

#Using filterX merge all K-mer data
filterx -k s -1 cnt>=1 FEMALE.out:grp=1 MALE.out:grp=1 > ALL_K-mer.out

#Extract specific reads from specific K-mer
perl step1.pl cfg.list ALL_K-mer.out

#Using SOAPdenovo2 to assembly special reads
SOAPdenovo127mer all -s FEMALE.cfg -o FEMALE_SOAPdenovo -p 10
SOAPdenovo127mer all -s MALE.cfg -o MALE_SOAPdenovo -p 10

#Using BWA and SamTools Computational coverage depth
bwa index FEMALE.fa
bwa mem -t 10 FEMALE.fa FEMALE_clean_R1.fq FEMALE_clean_R2.fq | samtools view -bS - > FEMALE_FEMALE.bam
bwa mem -t 10 FEMALE.fa MALE_clean_R1.fq MALE_clean_R2.fq | samtools view -bS - > MALE_FEMALE.bam
samtools sort -@ 10 FEMALE_FEMALE.bam > FEMALE_FEMALE_sort.bam
samtools sort -@ 10 MALE_FEMALE.bam > MALE_FEMALE_sort.bam
samtools index FEMALE_FEMALE_sort.bam
samtools index MALE_FEMALE_sort.bam
samtools depth -aa FEMALE_FEMALE_sort.bam MALE_FEMALE_sort.bam > FEMALE.depth
#For MALE is the same method

#Screening candidate specific sequences based on coverage depth
perl step2.pl FEMALE.fa MALE.fa 3 3 FEMALE.depth MALE.depth 0.9 20 100
```

Introduction
==========
Ssp2(Sex-specific-Perl-2.0) is the new version of `Sex-specific-Perl` which is a pipeline to rapid identification of sex-specific DNA fragments in Acipenseridae using comparative genomics with high-throughput sequencing. 

Installation
==========
Ssp2 is a program based on `Perl`, at the same time, it needs the support of several kinds of software.<br>
These are the software that needs to be used:`Perl` `kmer_count` `filterx` `SOAPdenovo` `BWA` `SamTools`.

Usage
==========
Step0.pl is recommended to get all command lines directly. Ssp2 needs to write the configuration file which named ctg.list, The details are as follows:<br>
```Bash
WORK_NAME=TEST	#Any name of the workflow
INSTALL_PATH=path_to/bin	#The installation path for ssp2
WORK_PATH=path_to/work	#The work path for ssp2
CFG_PATH=path_to/cfg.list	#The path for file cfg.list

KMER_COUNT_PATH=path_to/kmer_count	#The path for kmer_count
FILTERX_PATH=path_to/filterx	#The path for filterx
SOAPDENOVO_PATH=path_to/SOAPdenovo127mer	#The path for SOAPdenovo127mer
BWA_PATH=path_to/bwa	#The path for bwa
SAMTOOLS_PATH=path_to/samtools	#The path for samtools
THREAD=6	#Maximum number of threads used

G1_SAMPLE_NUMBER=3	#Number of samples in the first group/FEMALE
G2_SAMPLE_NUMBER=3	#Number of samples in the second group/MALE
GAP_LENGTH=20	#The base number of the largest gap in a specific sequence
MIX_LENGTH=100	#Minimum length of specific sequence
CONTIG_SPECIFIC_RATIO=0.9	#The minimum proportion of specific fragments in the same group of samples
SPECIFIC_KMER_RATIO=0.9	#The minimum proportion of specific K-mer in the same group of samples

G1_1_PATH=path_to/G1_0.R1.gz
G1_2_PATH=paht_to/G1_0.R2.gz
G2_1_PATH=path_to/G2_0.R1.gz
G2_2_PATH=paht_to/G2_0.R2.gz
#Provide all sample sequencing data paths in turn 
```
Getting Help
==========
Please use the [GitHub Issues page](https://github.com/fengtong2018/ssp2/issues) if you have questions.You can also contact Tong Feng by fengtong-bio@qq.com
