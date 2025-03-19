# chiav2

ChiaV2 is a packaged version of the [ChIA-PIPE workflow](https://github.com/TheJacksonLaboratory/ChIA-PIPE), which eliminates the hassle of software installation. <br/>
**DockerHub Link: https://hub.docker.com/r/fengchuiguo1994/chiav2/tags**

##  Introduction
chiav2 processes the sequencing data (ChIATAC/ChIA-PET/Bridge-linker Hi-C) to generate chromatin loops, coverage signal and contact matrix. chiaV2 includes linker calling, mapping, duplicate removing, loop calling. It also contain the software for Hi-C analysis.

####   Quick download chiaV2 from DockerHub
You can download chiav2 by running the following command:
```
singularity build chiav2.sif docker://fengchuiguo1994/chiav2:1.0
singularity build chiav2.sif docker://docker.awsl9527.cn/fengchuiguo1994/chiav2:1.0
```

### Usage
A ChIATAC demo data from mice was provided.
```
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


# Hi-C pipeline
Hi-C pipeline is a packaged workflow that works for chromatin interaction data, which eliminates the hassle of software installation. <br/>
**DockerHub Link: https://hub.docker.com/r/fengchuiguo1994/chiav2/tags**

##  Introduction
Hi-C pipeline processes the sequencing data of chromation interactions to generate chromatin loops and contact matrix. Hi-C pipeline includes mapping, duplicate removing, pairing. chiaV2.sif also contain the software for Hi-C analysis.


### Usage
A Hi-C demo data from mice was provided.
```
# 1. Prepare genome
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M35/GRCm39.genome.fa.gz
gunzip GRCm39.genome.fa.gz
singularity run chiav2.sif bwa index GRCm39.genome.fa

# 2. Modify the information in the config file.

# 3. Execute program.
bash hic_pipeline.sh ${PWD}/demoData_hic/Bcell_hic ${PWD}/result # abspath
${PWD}/demoData_hic/Bcell_hic: the prefix of the input fastq data, must be absolute path
${PWD}/result: output dir, must be absolute path
```

###   Example result (in result folder)
| file name | description |
| ----------- | ----------- |
| *.hic | contact matrix (Visualization in [juicer-box](https://aidenlab.org/juicebox/))|
| *.flt1000.hic | contact matrix filtering cis contacts < 1kb (Visualization in [juicer-box](https://aidenlab.org/juicebox/))|
| *.mapped.pairs.cut.gz | contact txt|
| *.flt1000.mapped.pairs.cut.gz | contact txt filtering cis contacts < 1kb |
| *.out.qc | result report |
