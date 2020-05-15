cd ~
mkdir sra
cd sra

#6.3.5. Then enter the following command:

source /etc/profile.d/markcbm.sh
esearch -db sra -q 'slr-seq AND mus_musculus[organism]' | efetch -format runinfo > slr_seq_runinfo.csv

#This will create a file named ‘slr_seq_runinfo.csv’. Use the unix command such as ‘less’ to look at the information present in this file. 

#6.4
#6.4.1. Download the sra (SRR390728) run file for the cancer cell line RNA-Seq experiment SRX079566, using the below command;

prefetch -a "/usr/local/molbiocloud/aspera/cli/bin/ascp|/usr/local/molbiocloud/aspera/cli/etc/asperaweb_id_dsa.openssh" SRR390728

#-a option specifies the Aspera ascp command line utility and the required public key location.

#6.4.2. Run the sra-stat tools to look at the basic read statistics on the downloaded sra file

sra-stat --quick --xml SRR390728

#6.5. Extract data from sra file
#6.5.1. Run the fasterq-dump tool to download the fastq file, using the below command;

fasterq-dump SRR390728

#6.5.2. Dump SAM alignment from the sra file and write it in BAM format using the below command;

sam-dump SRR390728 | samtools19 view -bS - > SRR390728.bam

#if 30 didnt work, copy that bam file using the below command;
cp /home/manager/sra/SRR390728.bam .

#6.5.3. Visualize the dumped bam file using the below command;
#Sort the bam file and create an index
#look at the CPU usage when everyone is sorting

samtools19 sort --threads 16 SRR390728.bam -o SRR390728.sorted.bam
samtools19 index SRR390728.sorted.bam
#use igv to view 1:89051831 SRR390728.sorted.bam

#6.6 Align reads using Magic-BLAST
#From NCBI Magic-BLAST documentation page:
#Magic-BLAST is a tool for mapping large next-generation RNA or DNA sequencing runs against a whole genome or transcriptome. Each alignment optimizes a composite score, taking into account simultaneously the two reads of a pair, and in case of RNA-seq, locating the candidate introns and adding up the score of all exons. This is very different from other versions of BLAST, where each exon is scored as a separate hit and read-pairing is ignored.

#We will align the reads from one of the SRA runs from our first GEO example. 
#6.6.1 Go to NCBI GEO and search for ‘staphylococcus fosfomycin’. Click on the first result to bring up the GEO Accession view page. Find the link to the SRA project at the bottom of the page and navigate to the SRA results page.
#6.6.2 Which result has the fewest number of reads? Click on that hit to go to the SRA page. Copy the SRA run accession. We will align these reads. 
#6.6.3  In order to run magicblast, we will need a blastdb. Generate a blastdb for the reference genome as follows:

## create a directory in your home folder
cd
mkdir NCTC_8325
cd NCTC_8325

## fetch the reference genome sequence in FASTA format using Entrez Direct utility
efetch -db nuccore -id NC_007795.1 -format fasta > NC_007795.1.fa

## create a blastdb
makeblastdb -in NC_007795.1.fa -dbtype nucl -parse_seqids -out NCTC_8325 -title 'NCTC_8325'

## run magicblast 
## <skip -- long runtime>
## check precomputed results at ~/Desktop/unix/NCTC_8325
#magicblast -sra SRR3659236 -db NCTC_8325 -num_threads 4 -out SRR3659236.aligns.sam -paired

#magic blast for covid example
#download sars covid19 genome
efetch -db nuccore -id NC_045512.2 -format fasta > NC_045512.2.fasta

#make a blastdb for the downloaded covid19 reference genome
makeblastdb -in NC_045512.2.fasta -dbtype nucl -parse_seqids -out SARS_COV2 -title 'SARS_COV2'

#magic blast
magicblast -sra SRR11652915 -db SARS_COV2 -num_threads 8 -out SRR11652915.aligns.sam -paired

#sam to bam convert, sort
samtools19 view -u SRR11652915.aligns.sam | samtools19 sort - -o SRR11652915.aligns.bam

#index bam
samtools19 index SRR11652915.aligns.bam


