# Introduction
aaaaaaaaaa

# Installation
* Requirements: Perl

```
Install Perl and add /bin directory to your executable path.
Perl: https://www.perl.org/
```


* Install Strandscript:

```
unzip Strandscript.zip       #Unzip the file
cd Strandscript/             #Change directories into the folder
chmod 755 bin/*.pl           #Change the mode of executable files

\#Add Strandscript to Shell searching path ($PATH). This step is optional.
\#If your NRSA is installed at /home/usrname/NRSA
export PATH=/home/usrname/NRSA/bin/:$PATH
```

* Download fasta files of reference genome
Please download the fasta files for genome hg19, GRCh38, mm9, and mm10 from http://hgdownload.cse.ucsc.edu/downloads.html. Uncompress and save as hg19.fa/GRCh38.fa/mm9.fa/mm10.fa into folder /fasta under Strandscript.

* Install required perl packages
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

Usage
Step1: Illumina manifest file checking (input: manifest file)
Usage: perl step1-mismatch.pl [options] -in manifest.csv
e.g: perl ./bin/step1-mismatch.pl -o test/ -n test1 -g hg19 -in OncoChip.csv
options:
-g [string]             define the genome: GRCh38, hg19, mm9, or mm10 (default: hg19)
-o [string]             output directory (default: current directory)
-in [string]            required, input manifest file in csv format
-n [string]             name for output
-h                      this help message
If you have added strandscript to Shell searching path ($PATH), please use the command "step1-mismatch.pl -o test/ -n test1 -g hg19 -in OncoChip.csv" instead.

Step2: Strand flip for plink file (input: ped & map)
Usage: perl step2-flip.pl [options] -in new_manifest.csv -map plink.map -ped plin.ped -o work_dir/
e.g: perl ./bin/step2-flip.pl -o test/ -in new_OncoChip.csv -map OncoChip.map -ped OncoChip.map -c 0.2
options:
-in [string]               required, input new manifest file in csv format from step1
-map [string]           required, plink map file for flipping
-ped [string]            required, plink ped file for flipping
-o [string]               output directory (default: current directory)
-c [real]                  cutoff for mismatch percentage of Probe sequence, SNPs would be flipped which
                               mismatch<cutoff (default: 0.2)
-h                           this help message

If you have added Strandscript to Shell searching path ($PATH), please use the command "step2-flip.pl -o test/ -in new_OncoChip.csv -map OncoChip.map -ped OncoChip.map -c 0.2" instead.

