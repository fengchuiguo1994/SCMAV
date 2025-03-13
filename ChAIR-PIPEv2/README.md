## ChAIR-PIPEv2

####  Introduction
ChAIR-PIPEv2 processes the sequencing data (ChAIR) to generate single-cell contact, ATAC fragment, peak count matrix, gene expression matrix ......, bulk-cell contact matrix, chromatin loops, coverage signal, peak ...... . ChAIR-PIPEv2 includes linker calling, mapping, duplicate removing, loop calling and more.


### Usage
A ChAIR demo data from mice was provided.
```
wget https://download.cncb.ac.cn/gsa-human/HRA007068/HRR1686121/HRR1686121.tar.gz
wget https://download.cncb.ac.cn/gsa-human/HRA007068/HRR1686122/HRR1686122.tar.gz
tar zxf HRR1686121.tar.gz
tar zxf HRR1686122.tar.gz
mv K562_Rep1.DNA.R1_001.fastq.gz SCG0074-15872_SI-NA-D6_S5_L003_R1_001.fastq.gz
mv K562_Rep1.DNA.R2_001.fastq.gz SCG0074-15872_SI-NA-D6_S5_L003_R2_001.fastq.gz
mv K562_Rep1.DNA.R3_001.fastq.gz SCG0074-15872_SI-NA-D6_S5_L003_R3_001.fastq.gz
mv K562_Rep1.RNA.R1_001.fastq.gz SCG0081-15878_GTAACATGCG-AGGTAACACT_S1_L001_R1_001.fastq.gz
mv K562_Rep1.RNA.R2_001.fastq.gz SCG0081-15878_GTAACATGCG-AGGTAACACT_S1_L001_R2_001.fastq.gz
    
# 1. Prepare genome
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M35/GRCm39.genome.fa.gz
gunzip GRCm39.genome.fa.gz
singularity run chiav2.sif bwa index GRCm39.genome.fa

# 2. Modify the information in the config file.

# 3. Execute program.
bash run.chiaV2.sh chiaV2.config.sh ${PWD}/demoData/ChIATACdata/SCG0192 result # abspath
chiaV2.config.sh: config file
${PWD}/demoData/ChIATACdata/SCG0192: the prefix of the input fastq data, must be absolute path
result: output dir
```

###   Example result
| file name | description |
| ----------- | ----------- |
| *.narrowPeak | DNA binding region (like ChIPSeq) |
| *.treat_pileup.NDP.bw | coverage signal (like ChIPSeq) |
| *.hic | contact matrix (Visualization in [juicer-box](https://aidenlab.org/juicebox/))|
| *.BE2.sigf.interactions | chromatin loops |
| *.final_stats.tsv | result report |
