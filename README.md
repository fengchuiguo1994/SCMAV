## ChAIR library
![](img/tech.png)
Schematics of ChAIR procedure to generate ChAIR-RNA, -ATAC, and -PET datasets in single cells. In ChAIR, crosslinked cells are manipulated through in situ restriction digestion, proximity ligation, and Tn5 tagmentation. The processed cells are then individually encapsulated into droplets using the 10x Genomics Multiome system, which barcodes DNA and RNA molecules simultaneously. The corresponding DNA and RNA information is obtained through DNA sequencing. ChAIR data contain two library data (ChIATAC-DNA and cDNA) or three library data (ChIATAC-DNA, ATAC-DNA and cDNA). 

## scMAV
single-cell multi-omic analysis and visualization tools. scMAV contain a <strong>[ChAIR-PIPEv2](ChAIR-PIPEv2)</strong>, <strong>[ChIA-PIPEv2](ChIA-PIPEv2)</strong>, <strong>[scMAV](scMAV)</strong>, <strong>ChAIR-View ([shiny version](ChAIR-View-shiny) and [local version](ChAIR-View-local))</strong> and <strong>[basic browser](basic-browser)</strong>. <br/>

The <strong>core components</strong> of data processing have been uploaded to DockerHub to avoid tedious software installation. <br/>
**DockerHub Link: [https://hub.docker.com/repository/docker/fengchuiguo1994/scMAV:tagsname](https://hub.docker.com/repository/docker/fengchuiguo1994/scMAV)

## [System Requirements](https://sylabs.io/guides/3.8/admin-guide/installation.html)
#### Hardware
scMAV should be run on a Linux system that meets the following requirements:
* 8-core Intel or AMD processor (>20 cores recommended)
* 64GB RAM (>320GB recommended)
* 64-bit CentOS/RedHat 7.8 or Ubuntu 20.04

#### Software
* Singularity: a container platform, >=3.8

#### Quick installation of Singularity
```
## On Red Hat Enterprise Linux or CentOS install the following dependencies:
$ sudo yum update -y && sudo yum groupinstall -y 'Development Tools' && sudo yum install -y openssl-devel libuuid-devel libseccomp-devel wget squashfs-tools cryptsetup

## On Ubuntu or Debian install the following dependencies:
$ sudo apt-get update && sudo apt-get install -y build-essential uuid-dev libgpgme-dev squashfs-tools libseccomp-dev wget pkg-config git cryptsetup-bin

## Install Go
$ export VERSION=1.14.12 OS=linux ARCH=amd64 && wget https://dl.google.com/go/go$VERSION.$OS-$ARCH.tar.gz && sudo tar -C /usr/local -xzvf go$VERSION.$OS-$ARCH.tar.gz && rm go$VERSION.$OS-$ARCH.tar.gz

$ echo 'export GOPATH=${HOME}/go' >> ~/.bashrc && echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc && source ~/.bashrc

## Install singularity on CentOS without compile
$ yum install -y singularity
```

#### Quick download scMAV from DockerHub
You can download scMAV by running the following command:
```
singularity build scMAV.sif docker://fengchuiguo1994/scMAV:v1
```

## ChIATAC analysis
#### 1. Trim Nextera Transposase Sequence
![](img/fastqc.png)<br/>
<p align="center">generate by FastQC</p><br/>
ChIATAC/ChAIR-DNA/ATAC data's end always has Nextera Transposase Sequence. The Nextera Transposase Sequence is removed by using Trimmomatic.<br/>

```shell
trimmomatic LEADING:5 TRAILING:5 ILLUMINACLIP:${trimmomaticpath}/adapters/NexteraPE-PE.fa:2:30:7:8:true MINLEN:30
```

#### 2. Get truncated reads.
Linker calling by using cutadapt, and get the truncated reads.

```shell
cutadapt -b file:${linkerfile} -n 50
```

#### 3. Align all reads and get pairs
All reads from R1 (contain truncated reads and no-truncated reads) are aligned by bwa (bwa mem: >=70bp, bwa aln: <70bp). All reads from R2 (contain truncated reads and no-truncated reads) are aligned by bwa (bwa mem: >=70bp, bwa aln: <70bp). Low mapped quality reads are removed. If the reads has more than one loci on R1 and R2 aligned to the genome (R1 and R2 aligned), take the most terminal (5'end) loci of R1 and R2 to form a pair. If more than one loci from R1 or R2 reads are aligned to the genome (R1 or R2 aligned), take the most terminal (5'end and 3'end) loci to form a pair.

```shell
bwa mem -5
bwa aln
```

#### 4. Get ATAC signal and PET signal
According insert size and strand to classify. 
- ATAC: distance less than 1kb and FR type pairs. 
- PET: distance great than 1kb and inter-chromosome pairs. 
- TTANDEM: distance less than 1kb and non-FR type pairs.

#### 5. Peak calling
The ATAC signal can be treated as a standard ATAC-Seq data. Peak calling by using macs2 with "BAMPE" model.

```shell
macs2 -f BAMPE
```

#### 6. Loop calling
There are three strategies for loop calling by using PET signal. 
- Both ends of pairs overlap with peaks. Only interactions between peaks can be got. [recommend for ChIATAC or active histone modification or TF protein ChIAPET]
- Either ends of pairs overlap with peaks. Interactions with peak participation can be got. [recommend for less active TF protein ChIAPET]
- All pairs. [recommend for inactive histone modification ChIAPET]

## ChAIR data analysis
#### ChAIR-DNA analysis
ChAIR-DNA is a single-cell ChIATAC data. The data analysis process is the same as the standard ChIATAC analysis process, except that the STARsolo is used to identify the cell barcode of each reads and add it to the reads' information.
#### 1. Add cell barcode

```shell
STAR --soloCBstart $barcodeloci --soloCBlen $barcodelen --soloUMIstart 1 --soloUMIlen 1 --soloType CB_samTagOut --soloCBwhitelist $DNAbarcodefile --outSAMattributes CB --outSAMtype BAM Unsorted --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloCBmatchWLtype 1MM --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --outSAMunmapped Within
```

#### 2. The same data analysis like ChIATAC step1-step6

#### 3. Get the scATAC signal fragment and peak intensity matrix
scATAC signal fragment is generated with "+4 -5". Peak intensity matrix was generated by using fragments and peaks. Finally, cell calling by using STAR

```shell
STAR --runMode soloCellFiltering --soloCellFilter EmptyDrops_CR
```

#### 4. Get the scPET signal
Building a contact number-cell matrix to do the cell calling by using STAR

```shell
STAR --runMode soloCellFiltering --soloCellFilter EmptyDrops_CR
```

scbedpe and sc4dn are generated by using scPET signal. We choose the BandNorm to generating single-cell PET matrix.

## Hi-C data analysis
There are many tools for Hi-C data processing. Earlier tools have strictly filtered pairs based on the positional relationship between the sequence of the restriction enzyme cutting site and the position on the genome where such suquence is aligned (HiCUP, HiC-Pro ...). Now, some tools rely on the genomic distance between aligned loci to filter pairs (pairtools, 4dn pipeline ...). This aproach is employed in ChIATAC/ChIA-PET data processing.

#### 1. Trim illumina Universal Adapter

```shell
trimmomatic LEADING:5 TRAILING:5 ILLUMINACLIP:${trimmomaticpath}/adapters/TruSeq3-PE-2.fa:2:30:7:8:true MINLEN:30
```

#### 2. Align and get pairs

```shell
bwa mem -5SP
```

#### 3. Get PET signal
According insert size and strand to classify. 
- PET: distance great than 1kb and inter-chromosome pairs. 
- genome sequence (WGS) : distance less than 1kb and FR type pairs. 
- TTANDEM: distance less than 1kb and non-FR type pairs.
Convert the PET signal to hic file.

## scHi-C data analysis
The data analysis process is the same as the standard Hi-C analysis process, except that the STARsolo is used to identify the cell barcode of each reads and add it to the reads' information.
#### 1. Add cell barcode

```shell
STAR --soloCBstart $barcodeloci --soloCBlen $barcodelen --soloUMIstart 1 --soloUMIlen 1 --soloType CB_samTagOut --soloCBwhitelist $DNAbarcodefile --outSAMattributes CB --outSAMtype BAM Unsorted --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloCBmatchWLtype 1MM --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --outSAMunmapped Within
```

#### 2. Trim Nextera Transposase Sequence
<strong>CAUTION</strong>
scHi-C and Hi-C use the different sequence adapt.<br/>

```shell
trimmomatic LEADING:5 TRAILING:5 ILLUMINACLIP:${trimmomaticpath}/adapters/NexteraPE-PE.fa:2:30:7:8:true MINLEN:30
```

#### 3. The same data analysis like Hi-C step2-step3

#### 4. Get the scPET signal
Building a contact number-cell matrix to do the cell calling by using STAR

```shell
STAR --runMode soloCellFiltering --soloCellFilter EmptyDrops_CR
```

scbedpe and sc4dn are generated by using scPET signal. We choose the BandNorm to generating single-cell PET matrix.

## scATAC data analysis
There are many software can process the scATAC data, e.g. Cellranger-atac, chromap ... scMAV can call chromap or cellranger-atac to process scATAC data. scMAV also provide also provides the ability to analyze scATAC
#### 1. Add cell barcode

```shell
STAR --soloCBstart $barcodeloci --soloCBlen $barcodelen --soloUMIstart 1 --soloUMIlen 1 --soloType CB_samTagOut --soloCBwhitelist $DNAbarcodefile --outSAMattributes CB --outSAMtype BAM Unsorted --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloCBmatchWLtype 1MM --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --outSAMunmapped Within
```

#### 2. The same data analysis like ChIATAC step1-step6

#### 3. Get the scATAC signal fragment and peak intensity matrix
scATAC signal fragment is generated with "+4 -5". Peak intensity matrix was generated by using fragments and peaks. Finally, cell calling by using STAR

```shell
STAR --runMode soloCellFiltering --soloCellFilter EmptyDrops_CR
```

#### 4. Get the scPET signal
Building a contact number-cell matrix to do the cell calling by using STAR

```shell
STAR --runMode soloCellFiltering --soloCellFilter EmptyDrops_CR
```

scbedpe and sc4dn are generated by using scPET signal. We choose the BandNorm to generating single-cell PET matrix.










```
# samtools view ../5M/8kmousecortex-5Mchromaplocalbugbak/align/8kmousecortex-5M.bam | cut -f 1 | uniq -c | awk '$1==1' > name.list

# python extractFq.py name.list /data/home/ruanlab/huangxingyu/scChAIR/publicscATAC/8Kmousecortex/8k_mouse_cortex_ATACv2_nextgem_Chromium_X_fastqs/tmp/8kmousecortex-5M_S1_L001_R2_001.fastq.gz | gzip > 8kmousecortex-5M_S1_L001_R2_001.fastq.gz
# python extractFq.py name.list /data/home/ruanlab/huangxingyu/scChAIR/publicscATAC/8Kmousecortex/8k_mouse_cortex_ATACv2_nextgem_Chromium_X_fastqs/tmp/8kmousecortex-5M_S1_L001_R1_001.fastq.gz | gzip > 8kmousecortex-5M_S1_L001_R1_001.fastq.gz
# python extractFq.py name.list /data/home/ruanlab/huangxingyu/scChAIR/publicscATAC/8Kmousecortex/8k_mouse_cortex_ATACv2_nextgem_Chromium_X_fastqs/tmp/8kmousecortex-5M_S1_L001_R3_001.fastq.gz | gzip > 8kmousecortex-5M_S1_L001_R3_001.fastq.gz

# /data/home/ruanlab/huangxingyu/Tools/chromap/chromap --preset atac --remove-pcr-duplicates-at-cell-level --Tn5-shift -t 10 -r /data/home/ruanlab/huangxingyu/Tools/cellranger-atac-2.1.0/genome/refdata-cellranger-arc-mm10-2020-A-2.0.0/fasta/genome.fa -x /data/home/ruanlab/huangxingyu/Tools/cellranger-atac-2.1.0/genome/chromapindex/mm10chromapindex -1 8kmousecortex-5M_S1_L001_R1_001.fastq.gz -2 8kmousecortex-5M_S1_L001_R3_001.fastq.gz -b 8kmousecortex-5M_S1_L001_R2_001.fastq.gz --barcode-whitelist ../5M/8kmousecortex-5Mchromaplocal/barcodelist -o test.sam --SAM
# /data/home/ruanlab/huangxingyu/Tools/chromap/chromap --preset atac --remove-pcr-duplicates-at-cell-level --Tn5-shift -t 10 -r /data/home/ruanlab/huangxingyu/Tools/cellranger-atac-2.1.0/genome/refdata-cellranger-arc-mm10-2020-A-2.0.0/fasta/genome.fa -x /data/home/ruanlab/huangxingyu/Tools/cellranger-atac-2.1.0/genome/chromapindex/mm10chromapindex -1 8kmousecortex-5M_S1_L001_R1_001.fastq.gz -2 8kmousecortex-5M_S1_L001_R3_001.fastq.gz -b 8kmousecortex-5M_S1_L001_R2_001.fastq.gz --barcode-whitelist ../5M/8kmousecortex-5Mchromaplocal/barcodelist -o test.bed


# samtools view ../10M/8kmousecortex-10Mchromaplocalbugbak/align/8kmousecortex-10M.bam | cut -f 1 | uniq -c | awk '$1==1' > name.10M.list
# cat name.list name.10M.list | sort | uniq -c | awk '$1==2'
# grep A00836:1001:H5M52DSX3:3:1172:2935:6073 name.list 
# grep A00836:1001:H5M52DSX3:3:1172:2935:6073 name.10M.list
# samtools view ../10M/8kmousecortex-10Mchromaplocalbugbak/align/8kmousecortex-10M.bam|grep A00836:1001:H5M52DSX3:3:1172:2935:6073
```


旧版1-linker：
```
>linkerA
ACGCGATATCTTATCTGACT
>linkerA_
AGTCAGATAAGATATCGCGT
```

旧版2-linker（linker1）：
```
>linkerA
ACGTGATATTTCACGACTCT
>linkerA_
AGAGTCGTGAAATATCACGT
```

新版-linker（linker76）：
```
>linkerA
ATAGGCGCTGGTTGCGTTGT
>linkerA_
ACAACGCAACCAGCGCCTAT
```