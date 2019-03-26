# Pipeline

## 01_FASTQ_download_ENA

### Fast download of FASTQ files and metadata from the European Nucleotide Archive (ENA)

1. How to perform bulk download of **fastq** files from the European Nucleitode Archive (ENA).

2. Typically people ask on how to get a certain SRA file from NCBI and how to convert it to fastq.

3. The common answer is prefetch followed by fastq-dump, but especially the latter is rather slow, so total file processing might take some time, especially if CPU (and disk) ressources are limited.

4. Luckily, most published (and unrestricted) sequencing data are mirrored at the ENA directly in fastq format, and there is a simple and efficient way to retrieve them.

5. In this tutorial, we will examplarily download an entire dataset of ChIP-seq and ATAC-seq data, requiring minimal preprocessing work.

6. We will use the Aspera client for download rates of several tens of Mb/s up to few hundred Mb/s (depending on the connection, I/O capacity and distance to the download location).

7. This example code should work on Linux and Mac.

```
$ awk 'FS="\t", OFS="\t" { gsub("ftp.sra.ebi.ac.uk", "era-fasp@fasp.sra.ebi.ac.uk:"); print }' accessions.txt | cut -f3 | awk -F ";" 'OFS="\n" {print $1, $2}' | awk NF | awk 'NR > 1, OFS="\n" {print "ascp -QT -l 300m -P33001 -i /Users/Akhil/Applications/Aspera_Connect.app/Contents/Resources/asperaweb_id_dsa.openssh" " " $1 " ."}' > download.txt
```

## 02_DNA_analysis_workflow

1. Download the file hg19.2bit file from https://hgdownload-test.gi.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit

2. Download the file twoBitToFa file from http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/twoBitToFa

3. Convert .2bit file to FASTA file. First make hg19.2bit file executable:

```sh
chmod 744 twoBitToFa
```

Then execute it like this:

```sh
./twoBitToFa hg19.2bit hg19.fasta
```

4. Use bwa to index the hg19.fasta file

```sh
bwa index hg19.fasta
```

5. Align fastq files using the reference sequence

```sh
bwa mem hf19.fasta ABC.1.fastq.gz ABC.2.fastq.gz > ABC.sam 
```

## 03_BWA_Alignment

#### Download Heng's BWA Aligner
So Heng's bwa aligner can perform the alignment. All BWA versions can be downloaded from here:

http://sourceforge.net/projects/bio-bwa/files/

The project page is here:

http://bio-bwa.sourceforge.net/

Below is how to get the latest version which recently was released:

```sh
wget http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.9a.tar.bz2/download
```
#### Unzip and Compile bwa
Next uncompress BWA and compile it:

```sh
tar xjvf bwa-0.7.9a.tar.bz2
cd bwa-0.7.9a
make
```
#### Data Download
Next download the data to be used for alignment - below are two example files to try out:

```sh
curl -O http://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/Escherichia_coli_O104_H4_2011C_3493_uid176127/NC_018658.fna


curl -O https://s3.amazonaws.com/public.ged.msu.edu/SRR390202.pe.qc.fq.gz
```
The first file is the reference file which has the Fasta format. The Fasta format is made of two lines, and there can be multiple two-line entries:

1. Description line with > at the front
2. Sequence

Below is an example with a portion of the sequence:

```sh
>gi|407479587|ref|NC_018658.1| Escherichia coli O104:H4 str. 2011C-3493 chromosome, complete genome
CATTATCGACTTTTGTTCGAGTGGAGTCCGCCGTGTCACTTTCGCTTTGGCAGCAGTGTCTTGCCCGATT
GCAGGATGAGTTACCAGCCACAGAATTCAGTATGTGGATACGCCCATTGCAGGCGGAACTGAGCGATAAC
ACGCTGGCCCTGTACGCGCCAAACCGTT
```
The second file is the paired-end reads file, which is described in the next section. To generate your custom file you will need to create both. The most important part for both the reads and reference files are the sequences and their respective descriptions/identifiers.

#### Process The Data
The reads to align are in Fastq format and below we extract a portion of them:

```sh
gunzip -c SRR390202.pe.qc.fq.gz | head -400000 > ecoli_pe.fq
```
The Fastq format are sequences which are listed in a 4-line format as follows, and there can be multiple four-line entries:

1. Sequence identifier
2. The sequence
3. Comments
4. Quality scores

##### Example of paired end (notice the /1 and /2):
```sh
@Example_read/1
CGCGTAACAAAAGTGTCTATAATCACGGCAGAAAA
+
HHHHHHHHHFHGGHHHHHHHHHHHHHHHHHHHHEH
@Example_read/2
TTTCCGGACCACATTAAATATACCAAAATTCCCCC
+
HHHHHHHHHFHGGHHHHHHHHHHHHHHHHHHHHEH
```
#### Create the Index of the reference for alignment
Next prepare the index to the reference for alignment:

```sh
bwa index -a bwtsw NC_018658.fna
```
#### Perform the Alignment
Next perform the alignment:

```sh
bwa mem -p NC_018658.fna ecoli_pe.fq > aln.x.ecoli_NC_018658.sam
```
#### Check for supplementary alignments
It's easiest to use samtools to determine which rows are supplementary alignments:

```sh
samtools view -f 0x800 -S aln.x.ecoli_NC_018658.sam
```
The alignment example does not contain any secondary alignments, but if it did you could retrieve them as such:

```sh
samtools view -f 0x100 -S aln.x.ecoli_NC_018658.sam
```
Below is a nice website where you can then type up the flag numbers, and it will checkbox for you which flags are set:

http://picard.sourceforge.net/explain-flags.html

To download samtools, below is the link:

```sh
http://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2/download
```
It will require compiling as with bwa, but it's fairly straightforward.