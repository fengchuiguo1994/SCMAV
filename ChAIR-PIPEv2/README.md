## ChAIR-PIPEv2

####  Introduction
ChAIR-PIPEv2 processes the sequencing data (ChAIR) to generate single-cell contact, ATAC fragment, peak count matrix, gene expression matrix ......, bulk-cell contact matrix, chromatin loops, coverage signal, peak ...... . ChAIR-PIPEv2 includes module: linker calling, mapping, duplicate removing, loop calling and more.

#### Usage
A ChAIR demo data from K562 was provided.<br/>
1. Prepare ChAIR data <br/>
```sh
wget https://download.cncb.ac.cn/gsa-human/HRA007068/HRR1686121/HRR1686121.tar.gz
wget https://download.cncb.ac.cn/gsa-human/HRA007068/HRR1686122/HRR1686122.tar.gz
tar zxf HRR1686121.tar.gz
tar zxf HRR1686122.tar.gz
mv K562_Rep1.DNA.R1_001.fastq.gz SCG0074-15872_SI-NA-D6_S5_L003_R1_001.fastq.gz
mv K562_Rep1.DNA.R2_001.fastq.gz SCG0074-15872_SI-NA-D6_S5_L003_R2_001.fastq.gz
mv K562_Rep1.DNA.R3_001.fastq.gz SCG0074-15872_SI-NA-D6_S5_L003_R3_001.fastq.gz
mv K562_Rep1.RNA.R1_001.fastq.gz SCG0081-15878_GTAACATGCG-AGGTAACACT_S1_L001_R1_001.fastq.gz
mv K562_Rep1.RNA.R2_001.fastq.gz SCG0081-15878_GTAACATGCG-AGGTAACACT_S1_L001_R2_001.fastq.gz

echo -e "$abspath/K562/data/SCG0074-15872_SI-NA-D6\t$abspath/K562/data/SCG0081-15878_GTAACATGCG-AGGTAACACT" > k562.input.list
```
2. Prepare genome <br/>
```sh
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.p13.genome.fa.gz
gunzip GRCh38.p13.genome.fa.gz
singularity run chiav2.sif bwa index GRCh38.p13.genome.fa
```
3. Modify the information (submit_cpu10x_for_scChAIR.hg38.pro.sh) <br/>
```sh
addLinkerSeq='ATAGGCGCTGGTTGCGTTGT' # linker sequence
genome=hg38 # assembly name
ARCref=~/Tools/cellranger-arc-2.0.2/genome/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 # cellranger-arc ref
GEXref=~/Tools/cellranger-7.1.0/genome/refdata-gex-GRCh38-2020-A # cellranger ref
ATACref=~/Tools/cellranger-atac-2.1.0/genome/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 # cellranger-atac ref
bwaref=~/refs/hg38/bwa/hg38.fa # bwa ref
blacklist=~/hg38.blacklist.bed # or null # download from encode
genebed=~/Tools/cellranger-7.1.0/genome/refdata-gex-GRCh38-2020-A/genes/hg38.info # gene info # bed format

siffile="~/scChAIR/ChAIR.sif"
codepath="~/scChAIR/codes"
```
4. Execute program. <br/>
```sh
bash submit_cpu10x_for_scChAIR.hg38.pro.sh k562.input.list # abspath
```

#### Example result
1. data stat
![](img/dataStat.png)

2. align insert size
![](img/NoLinker.png)
![](img/2tag.png)

3. profile
![](img/profile.png)

4. PET count distribution for loops
![](img/loop.PETcount.distribution.png)