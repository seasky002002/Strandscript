## Introduction
&nbsp;&nbsp;&nbsp;&nbsp;One aspect of the Illumina array design that has always pestered researchers is the problem of strand consistency. The autosomal chromosomes are diploid and double stranded. In GWAS, the significant SNPs are always reported as risk allele (RA), occasionally referred to as effective allele (EA). The RA can be presented either on the forward or the reverse strand of the genome. When designing the probes to detect the two alleles of a SNP, the probes can be designed for either the forward or reverse strand. For the majority of the time, the reported RA for a SNP will be converted to the forward strand prior to publication. However, we cannot rely completely upon the authors on the subject of strand consistency.  Strand inconsistency will result in wrong interpretation of the direction of the association, and will unnecessarily convolute downstream genetics analyses based on these published SNPs, such as genetic risk scores or Mendelian randomization studies.  

&nbsp;&nbsp;&nbsp;&nbsp;In the case of a genotyping array, a simple definition of forward and reverse strand for each SNP the array would suffice. Instead, Illumina introduced a more unclear definition of strand, top and bottom, which have caused vast confusion with the forward and reverse strand. More disturbingly, when outputting the genotyping data from GenomeStudio, an option can be selected to convert all SNPs to the forward strand. This option should solve the problem of converting to the standard forward and reverse annotations. However, we computed that between 1 ~ 11% of SNPs are not presented in the forward strand after the GenomeStudio conversion on various Illumina genotyping arrays. The strand of a SNP can be checked by comparing the allele frequency in the dataset with those previously reported from an appropriate population database. However, when the allele frequency is near 50%, ambiguity can still arise. Additionally, strand issue can be resolved by comparing the alleles to a reference genome. However, when the two alleles of the SNPs are complementary (A/T or C/G), we are still left without the ability to determine the true strand. The only absolute solution for this problem is to determine the strand by comparing the probe sequences to the reference genome, provided that the probe sequences is correct.  For most users of microarray genotyping, accessing the probe sequences and comparing them to the reference sequence has been impractical, until now.  	

&nbsp;&nbsp;&nbsp;&nbsp;To provide a rapid and convenient solution to solve the Illumina Genotyping data strand inconsistency issue we introduce StrandScript, a toolset that can examine the accuracy of Illumina probe design, the accuracy of Illumina manifest file annotation and which will convert Illumina genotyping data consistently to the human genome reference forward strand. StrandScript works with all Illumina genotyping arrays. 


## Installation
* Requirements: Perl

  Install Perl (https://www.perl.org/) and add /bin directory to your executable path.
  

* Install Strandscript:

```
unzip Strandscript.zip       #Unzip the file
cd Strandscript/             #Change directories into the folder
chmod 755 bin/*.pl           #Change the mode of executable files

#Add Strandscript to Shell searching path ($PATH). This step is optional.
#If your NRSA is installed at /home/usrname/NRSA
export PATH=/home/usrname/NRSA/bin/:$PATH
```

* Download fasta files of reference genome

  Please download the fasta files for genome hg19, GRCh38, mm9, and mm10 from (http://hgdownload.cse.ucsc.edu/downloads.html). 
  Uncompress and save as hg19.fa/GRCh38.fa/mm9.fa/mm10.fa into folder /fasta under /Strandscript.


* Install required perl packages
```
#Check whether all the packages needed are installed by running test.modules,
./bin/test.modules

#the output of test.modules looks like,
 ok    Cwd 
 ok    File::Basename
 ok    File::Path
 fail   Getopt::Long
 
 
#If “fail” shows up before the package name, that means users don’t have this required package. In this case, Getopt::Long need to be installed. 
Users can install missing packages by running "./bin/install.modules packagename".
#e.g  
./bin/install.modules Getopt::Long
  
#Add the NRSA lib to your PERL5LIB environment variable,
#If your NRSA is installed at /home/usrname/NRSA
export PERL5LIB=$PERL5LIB:/home/usrname/NRSA/lib/
```

## Usage
**Step1: Illumina manifest file checking (input: manifest file)**
```
Usage: perl step1-mismatch.pl [options] -in manifest.csv
e.g: perl ./bin/step1-mismatch.pl -o test/ -n OncoChip -g hg19 -in OncoChip.csv
options:
-g [string]             define the genome: GRCh38, hg19, mm9, or mm10 (default: hg19)
-o [string]             output directory (default: current directory)
-in [string]            required, input manifest file in csv format
-n [string]             name for output
-h                      this help message
If you have added strandscript to Shell searching path ($PATH), please use the command "step1-mismatch.pl -o test/ -n OncoChip -g hg19 -in OncoChip.csv" instead.
```
Step1 outputs two files into the output directory. In the e.g., there would be two files (new_OncoChip.csv & outdated_OncoChip.csv) in test/ folder.  
&nbsp;&nbsp;&nbsp;&nbsp;new_OncoChip.csv: lists the basic information and mismathed number of base pairs for snps. This file would be provided as input for the step2.  
&nbsp;&nbsp;&nbsp;&nbsp;outdated_OncoChip.csv: lists the outdated snps.


**Step2: Strand flip for plink file (input: ped & map)**
```
Usage: perl step2-flip.pl [options] -in new_manifest.csv -map plink.map -ped plin.ped -o work_dir/
e.g: perl ./bin/step2-flip.pl -o test/ -in new_OncoChip.csv -map OncoChip.map -ped OncoChip.map -c 0.2
options:
-in [string]               required, input new manifest file in csv format from step1
-map [string]              required, plink map file for flipping
-ped [string]              required, plink ped file for flipping
-o [string]                output directory (default: current directory)
-c [real]                  cutoff for mismatch percentage of Probe sequence, SNPs would be flipped which mismatch<cutoff (default: 0.2)
-h                         this help message

If you have added Strandscript to Shell searching path ($PATH), please use the command "step2-flip.pl -o test/ -in new_OncoChip.csv -map OncoChip.map -ped OncoChip.map -c 0.2" instead.
```

Step2 outputs three files into the output directory. In the e.g., there would be three files (flipped_OncoChip.map,  flipped_OncoChip.ped & filtered-out-snp.txt) in test/ folder.  

&nbsp;&nbsp;&nbsp;&nbsp;flipped_OncoChip.map & flipped_OncoChip.ped: flipped plink files.  
&nbsp;&nbsp;&nbsp;&nbsp;filtered-out-snp.txt: lists the flitered out snps which mismatch > cutoff.
