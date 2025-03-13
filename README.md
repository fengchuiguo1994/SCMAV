## ChAIR library
![](img/tech.png)
Schematics of ChAIR procedure to generate ChAIR-RNA, -ATAC, and -PET datasets in single cells. In ChAIR, crosslinked cells are manipulated through in situ restriction digestion, proximity ligation, and Tn5 tagmentation. The processed cells are then individually encapsulated into droplets using the 10x Genomics Multiome system, which barcodes DNA and RNA molecules simultaneously. The corresponding DNA and RNA information is obtained through DNA sequencing. ChAIR data contain two library data (DNA and cDNA). 

## SCMAV
single-cell multi-omic analysis and visualization tools. SCMAV contain a <strong>[ChAIR-PIPEv2](ChAIR-PIPEv2)</strong>, <strong>[ChIA-PIPEv2](ChIA-PIPEv2)</strong>, <strong>ChAIR-View ([shiny version](ChAIR-View-shiny) and [local version](ChAIR-View-local))</strong> and <strong>[basic browser](basic-browser)</strong>. <br/>

The <strong>core components</strong> of data processing have been uploaded to DockerHub to avoid tedious software installation. <br/>
**DockerHub Link: [https://hub.docker.com/repository/docker/fengchuiguo1994/scmav:tagsname](https://hub.docker.com/repository/docker/fengchuiguo1994/scmav)

## [System Requirements](https://sylabs.io/guides/3.8/admin-guide/installation.html)
#### Hardware
SCMAV should be run on a Linux system that meets the following requirements:
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

#### Quick download SCMAV from DockerHub
You can download SCMAV by running the following command:
```
singularity build scmav.sif docker://fengchuiguo1994/scmav:v1
```




GSE253407
GSM8020244	MC_Droplet_Hi-C_10xATAC_v1.1 SRX23255282: GSM8020244: MC_Droplet_Hi-C_10xATAC_v1.1; Mus musculus; OTHER 2 ILLUMINA (Illumina NovaSeq 6000) runs: 3.1G spots, 1,007.3G bases, 315.8Gb downloads
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/081/SRR27586281/SRR27586281_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/081/SRR27586281/SRR27586281_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/080/SRR27586280/SRR27586280_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/080/SRR27586280/SRR27586280_2.fastq.gz

GSM8020245	MC_Droplet_Hi-C_10xATAC_v2 SRX23255283: GSM8020245: MC_Droplet_Hi-C_10xATAC_v2; Mus musculus; OTHER 1 ILLUMINA (Illumina NovaSeq 6000) run: 2.5G spots, 830.1G bases, 263.9Gb downloads
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/079/SRR27586279/SRR27586279_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/079/SRR27586279/SRR27586279_2.fastq.gz

GSM8020246	Cell line mixing_Droplet_Hi-C_10xATAC_v2 SRX23255284: GSM8020246: Cell line mixing_Droplet_Hi-C_10xATAC_v2; Homo sapiens; OTHER 1 ILLUMINA (Illumina NovaSeq 6000) run: 493.8M spots, 161G bases, 51.2Gb downloads
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/078/SRR27586278/SRR27586278_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/078/SRR27586278/SRR27586278_2.fastq.gz

GSM8020247	COLO320DM_Droplet_Hi-C_10xATAC_v2 SRX23255285: GSM8020247: COLO320DM_Droplet_Hi-C_10xATAC_v2; Homo sapiens; OTHER 2 ILLUMINA (Illumina NovaSeq 6000) runs: 684.3M spots, 146.5G bases, 46.9Gb downloads
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/069/SRR27586269/SRR27586269_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/069/SRR27586269/SRR27586269_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/077/SRR27586277/SRR27586277_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/077/SRR27586277/SRR27586277_2.fastq.gz

GSM8020248	COLO320HSR_Droplet_Hi-C_10xATAC_v2 SRX23255290: GSM8020248: COLO320HSR_Droplet_Hi-C_10xATAC_v2; Homo sapiens; OTHER 2 ILLUMINA (Illumina NovaSeq 6000) runs: 649.1M spots, 138.7G bases, 44.5Gb downloads
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/071/SRR27586271/SRR27586271_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/071/SRR27586271/SRR27586271_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/068/SRR27586268/SRR27586268_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/068/SRR27586268/SRR27586268_2.fastq.gz

GSM8020249	GBM39_Droplet_Hi-C_10xATAC_v1.1 SRX23255289: GSM8020249: GBM39_Droplet_Hi-C_10xATAC_v1.1; Homo sapiens; OTHER 2 ILLUMINA (Illumina NovaSeq 6000) runs: 1G spots, 312.9G bases, 98.6Gb downloads
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/070/SRR27586270/SRR27586270_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/070/SRR27586270/SRR27586270_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/073/SRR27586273/SRR27586273_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/073/SRR27586273/SRR27586273_2.fastq.gz

GSM8020250	GBM39-ER_Droplet_Hi-C_10xATAC_v1.1 SRX23255288: GSM8020250: GBM39-ER_Droplet_Hi-C_10xATAC_v1.1; Homo sapiens; OTHER 2 ILLUMINA (Illumina NovaSeq 6000) runs: 1.3G spots, 386.1G bases, 120.8Gb downloads
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/072/SRR27586272/SRR27586272_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/072/SRR27586272/SRR27586272_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/075/SRR27586275/SRR27586275_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/075/SRR27586275/SRR27586275_2.fastq.gz

GSM8020251	MC_Paired_Hi-C_rep1_DNA SRX23255287: GSM8020251: MC_Paired_Hi-C_rep1_DNA; Mus musculus; OTHER 1 ILLUMINA (Illumina NovaSeq 6000) run: 501.9M spots, 113.4G bases, 32.2Gb downloads
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/074/SRR27586274/SRR27586274_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/074/SRR27586274/SRR27586274_2.fastq.gz

GSM8020252	MC_Paired_Hi-C_rep1_RNA SRX23255286: GSM8020252: MC_Paired_Hi-C_rep1_RNA; Mus musculus; OTHER 1 ILLUMINA (NextSeq 2000) run: 121.2M spots, 12.4G bases, 3.7Gb downloads
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/076/SRR27586276/SRR27586276_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/076/SRR27586276/SRR27586276_2.fastq.gz

GSM8020253	MC_Paired_Hi-C_rep2_DNA SRX23255295: GSM8020253: MC_Paired_Hi-C_rep2_DNA; Mus musculus; OTHER 1 ILLUMINA (Illumina NovaSeq 6000) run: 934.1M spots, 304.5G bases, 95.3Gb downloads
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/063/SRR27586263/SRR27586263_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/063/SRR27586263/SRR27586263_2.fastq.gz

GSM8020254	MC_Paired_Hi-C_rep2_RNA SRX23255294: GSM8020254: MC_Paired_Hi-C_rep2_RNA; Mus musculus; OTHER 1 ILLUMINA (NextSeq 2000) run: 132.3M spots, 13.5G bases, 4Gb downloads
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/064/SRR27586264/SRR27586264_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/064/SRR27586264/SRR27586264_2.fastq.gz

GSM8020255	MC_Paired_Hi-C_rep3_DNA SRX23255293: GSM8020255: MC_Paired_Hi-C_rep3_DNA; Mus musculus; OTHER 1 ILLUMINA (Illumina NovaSeq 6000) run: 2.1G spots, 483.8G bases, 141.8Gb downloads
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/065/SRR27586265/SRR27586265_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/065/SRR27586265/SRR27586265_2.fastq.gz

GSM8020256	MC_Paired_Hi-C_rep3_RNA SRX23255292: GSM8020256: MC_Paired_Hi-C_rep3_RNA; Mus musculus; OTHER 1 ILLUMINA (NextSeq 2000) run: 290.3M spots, 29.6G bases, 8.9Gb downloads
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/066/SRR27586266/SRR27586266_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/066/SRR27586266/SRR27586266_2.fastq.gz

GSM8020257	GBM39_Paired_Hi-C_DNA SRX23255291: GSM8020257: GBM39_Paired_Hi-C_DNA; Homo sapiens; OTHER 1 ILLUMINA (Illumina NovaSeq 6000) run: 1.1G spots, 197.7G bases, 67.3Gb downloads
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/067/SRR27586267/SRR27586267_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/067/SRR27586267/SRR27586267_2.fastq.gz

GSM8020258	GBM39_Paired_Hi-C_RNA SRX23255296: GSM8020258: GBM39_Paired_Hi-C_RNA; Homo sapiens; OTHER 1 ILLUMINA (NextSeq 2000) run: 123.9M spots, 12.6G bases, 3.8Gb downloads
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/062/SRR27586262/SRR27586262_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/062/SRR27586262/SRR27586262_2.fastq.gz

GSM8020259	GBM39-ER_Paired_Hi-C_DNA SRX23255297: GSM8020259: GBM39-ER_Paired_Hi-C_DNA; Homo sapiens; OTHER 1 ILLUMINA (Illumina NovaSeq 6000) run: 1.9G spots, 434.6G bases, 130.3Gb downloads
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/061/SRR27586261/SRR27586261_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/061/SRR27586261/SRR27586261_2.fastq.gz

GSM8020260	GBM39-ER_Paired_Hi-C_RNA SRX23255298: GSM8020260: GBM39-ER_Paired_Hi-C_RNA; Homo sapiens; OTHER 1 ILLUMINA (NextSeq 2000) run: 273.6M spots, 27.9G bases, 8.4Gb downloads
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/060/SRR27586260/SRR27586260_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/060/SRR27586260/SRR27586260_2.fastq.gz

GSM8405621	Normal_like_cell line mixing_Droplet_Hi-C_10xATAC_v2 SRX25379911: GSM8405621: Normal_like_cell line mixing_Droplet_Hi-C_10xATAC_v2; Homo sapiens; OTHER 2 ILLUMINA (Illumina NovaSeq 6000) runs: 203.8M spots, 46.1G bases, 12.9Gb downloads
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR298/092/SRR29883592/SRR29883592_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR298/092/SRR29883592/SRR29883592_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR298/093/SRR29883593/SRR29883593_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR298/093/SRR29883593/SRR29883593_2.fastq.gz

GSM8476142	PBMC_Paired_Hi-C_DNA SRX25791833: GSM8476142: PBMC_Paired_Hi-C_DNA; Homo sapiens; OTHER 1 ILLUMINA (Illumina NovaSeq 6000) run: 291M spots, 66.9G bases, 20.4Gb downloads
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR303/039/SRR30331639/SRR30331639_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR303/039/SRR30331639/SRR30331639_2.fastq.gz

GSM8476144	PBMC_Paired_Hi-C_RNA SRX25791834: GSM8476144: PBMC_Paired_Hi-C_RNA; Homo sapiens; OTHER 1 ILLUMINA (Illumina NovaSeq 6000) run: 176.9M spots, 17.7G bases, 5.3Gb downloads
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR303/038/SRR30331638/SRR30331638_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR303/038/SRR30331638/SRR30331638_2.fastq.gz

GSM8476146	PBMC_Droplet_Hi-C_10xATAC_v2 SRX25791835: GSM8476146: PBMC_Droplet_Hi-C_10xATAC_v2; Homo sapiens; ATAC-seq 1 ILLUMINA (Illumina NovaSeq 6000) run: 52.9M spots, 12G bases, 3.7Gb downloads
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR303/037/SRR30331637/SRR30331637_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR303/037/SRR30331637/SRR30331637_2.fastq.gz

GSM8521393	BMMC_before_Paired_Hi-C_DNA
GSM8521394	BMMC_after_Paired_Hi-C_DNA
GSM8521395	GBM_P4CC3_Droplet_Hi-C

k562 hic
https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/1f922dbc-637b-4c43-b816-f604f0e79ad1/4DNFI8Y9SRP2.hic
https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/b03516d4-62f4-47ec-9d9b-6baf95623376/4DNFI244AS29.mcool
https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/5448b5f6-248c-4a4b-aee7-3ce97b64db58/4DNFIIX5BNC9.pairs.gz

GM12878 hic
https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/893f17d2-3a08-4687-b56e-5bd8a8ea68be/4DNFIIS73OJN.pairs.gz
https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/7f13cf57-f867-4a3a-89d4-f797f7936c1b/4DNFI3O82QVV.pairs.gz

HeLa-S3 hic
https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/c1a0bf52-2e65-4deb-871f-ede736babd41/4DNFIPZBEXCP.hic
https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/c4e41bbc-8cab-4475-b29f-b29452a7f97a/4DNFIT7QBK39.mcool
https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/8bf59da6-cb12-4cf1-a2b1-080d8dd2f2dd/4DNFIW458FJD.pairs.gz









10X scATAC 产品对应的是cellranger-atac-2.1.0/lib/python/atac/barcodes/737K-cratac-v1.txt.gz
10X multi ATAC+GEX 产品对应的是 cellranger-arc-2.0.2/lib/python/atac/barcodes/737K-arc-v1.txt.gz + cellranger-arc-2.0.2/lib/python/cellranger/barcodes/737K-arc-v1.txt.gz # https://kb.10xgenomics.com/hc/en-us/articles/4412343032205-Where-can-I-find-the-barcode-inclusion-list-formerly-barcode-whitelist-for-Single-Cell-Multiome-ATAC-GEX-product
10X multi ATAC+GEX 的barcode映射关系 # https://kb.10xgenomics.com/hc/en-us/articles/360049105612-Barcode-translation-in-Cell-Ranger-ARC
10X multi ATAC+GEX 早期测序的影响 # https://scg-lib-structs.readthedocs.io/en/latest/multi/10x_multiome.html

ls -lh cellranger-6.0.2/lib/python/cellranger/barcodes/737K-arc-v1.txt.gz
ls -lh cellranger-7.1.0/lib/python/cellranger/barcodes/737K-arc-v1.txt.gz
ls -lh cellranger-atac-2.1.0/lib/python/atac/barcodes/737K-cratac-v1.txt.gz
ls -lh cellranger-atac-2.1.0/lib/python/atac/barcodes/737K-arc-v1.txt.gz
ls -lh cellranger-arc-2.0.2/lib/python/atac/barcodes/737K-arc-v1.txt.gz
ls -lh cellranger-arc-2.0.2/lib/python/cellranger/barcodes/737K-arc-v1.txt.gz
