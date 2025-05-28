## scRNA
In the software parameter, you can choose whether to call cellranger or starsolo. [cellranger](https://www.10xgenomics.com/support/software/cell-ranger/latest): slow. [STARsolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md): fast.

## scATAC
In the software parameter, you can choose whether to call cellranger-atac, chromap or custom. [cellranger-atac](https://www.10xgenomics.com/support/software/cell-ranger-atac/latest): slow. [chromap](https://github.com/haowenz/chromap): fast. custom: medium. <br>
#### Custom Design
- 1. Splitting large files into smaller chunk, trimming [Nextera Transposase Sequence](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) by using [Trimmomatic](https://github.com/usadellab/Trimmomatic).
- 2. Detecting barcode by using STARsolo, aligning by using [BWA-MEM](https://github.com/lh3/bwa).
- 3. Removing PCR reads by using [Picard](https://broadinstitute.github.io/picard/).
- 4. Peak calling by using [MACS2](https://github.com/macs3-project/MACS).
- 5. Generating the single cell fragment (+4, -5) result and the peak count matrix.

## scHi-C
In the software parameter, you can only choose custom. ChAIR-DNA data can also be processed with this flow.
#### Custom Design
- 1. Splitting large files into smaller chunk, trimming [Nextera Transposase Sequence](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) by using [Trimmomatic](https://github.com/usadellab/Trimmomatic).
- 2. Detecting barcode by using STARsolo, aligning by using [BWA-MEM](https://github.com/lh3/bwa), selecting uu|UU|Uu|uU|UR|RU type reads in [pairtools](https://github.com/open2c/pairtools) result.
- 3. Removing PCR reads by using [Picard](https://broadinstitute.github.io/picard/).
- 4. Generating the single cell contact and hic/cool file.

## ChAIR-DNA
In the software parameter, you can only choose custom.
#### Custom Design
- 1. Splitting large files into smaller chunk, trimming [Nextera Transposase Sequence](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) by using [Trimmomatic](https://github.com/usadellab/Trimmomatic), linker calling and removing by using [cutadapt](https://github.com/marcelm/cutadapt).
- 2. Detecting barcode by using STARsolo, aligning by using [BWA-MEM and BWA-aln](https://github.com/lh3/bwa).
- 3. Removing PCR reads by using [Picard](https://broadinstitute.github.io/picard/).
- 4. Extract the ATAC signal reads, Peak calling by using [MACS2](https://github.com/macs3-project/MACS), Generating the single cell fragment (+4, -5) result and the peak count matrix.
- 5. Extract the interaction signal reads, Generating the single cell contact and hic/cool file.

## scRNA+scATAC
In the software parameter, you can choose whether to call cellranger-arc or other combinations parameter.  [cellranger-arc](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest): slow. <br/>
Combinations: STARsolo+chromap: fast; STARsolo+custom: medium; ...

## scRNA+scHi-C
In the software parameter, you can only choose combinations parameter. ChAIR-DNA data can also be processed with this flow.<br/>
combinations: STARsolo+custom: fast; cellranger+custom: slow.

## ChAIR
In the software parameter, you can only choose combinations parameter. <br/>
combinations: STARsolo+custom: fast; cellranger+custom: slow.


## Tips
#### 1. Determine barcode location and orientation (eg. ChAIR).
```
zcat SCG0074-15872_SI-NA-D6_S5_L003_R2_001.fastq.gz | awk 'NR%4==2' | head -1000000 | cut -c 1-16 | sort | uniq | cat - ATAC.737K-arc-v1.txt | sort | uniq -d | wc -l
zcat SCG0074-15872_SI-NA-D6_S5_L003_R2_001.fastq.gz | awk 'NR%4==2' | head -1000000 | cut -c 9-24 | sort | uniq | cat - ATAC.737K-arc-v1.txt | sort | uniq -d | wc -l 
zcat SCG0074-15872_SI-NA-D6_S5_L003_R2_001.fastq.gz | awk 'NR%4==2' | head -1000000 | cut -c 1-16 | tr ACGT TGCA | rev | sort | uniq | cat - ATAC.737K-arc-v1.txt | sort | uniq -d | wc -l 
zcat SCG0074-15872_SI-NA-D6_S5_L003_R2_001.fastq.gz | awk 'NR%4==2' | head -1000000 | cut -c 9-24 | tr ACGT TGCA | rev | sort | uniq | cat - ATAC.737K-arc-v1.txt | sort | uniq -d | wc -l
```
