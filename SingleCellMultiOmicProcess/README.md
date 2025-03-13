## Transcriptome
#### STARsolo
[ref1](https://cumulus.readthedocs.io/en/latest/starsolo.html), [ref2](https://github.com/detrout/woldlab-rna-seq/blob/1.2.5/workflow/starsolo.snakefile), [ref3](https://broadinstitute.github.io/warp/docs/Pipelines/Multiome_Pipeline/README)
```
CB=16b，UMI=10b # 10X Chromium V2
CB=16b，UMI=12b # 10X Chromium V3 (ChAIR data)

--readFilesIn cDNAfragmentSequence.fastq.gz CellBarcodeUMIsequence.fastq.gz
--clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR

singularity exec STAR.sif STAR --genomeDir stargenome/hg38_index --outFileNamePrefix SCG0081_STARsoloCellRanger4 --readFilesCommand zcat --readFilesIn SCG0081-15878_GTAACATGCG-AGGTAACACT_S1_L001_R2_001.fastq.gz SCG0081-15878_GTAACATGCG-AGGTAACACT_S1_L001_R1_001.fastq.gz --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 --soloType CB_UMI_Simple --soloCBwhitelist 737K-arc-v1.txt --soloCellFilter EmptyDrops_CR --runThreadN 10 --outSAMattributes NH HI nM AS CR UR GX GN sS sQ sM MD CB UB --outSAMtype BAM SortedByCoordinate --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --outSAMunmapped Within --soloMultiMappers Unique EM
```

## ATAC
#### chromap

## PET (Hi-C like or ChIA-PET like)
#### SCMAV


## Tips
#### 1. Determine barcode location and orientation (eg. ChAIR).
```
zcat SCG0074-15872_SI-NA-D6_S5_L003_R2_001.fastq.gz | awk 'NR%4==2' | head -1000000 | cut -c 1-16 | sort | uniq | cat - ATAC.737K-arc-v1.txt | sort | uniq -d | wc -l
zcat SCG0074-15872_SI-NA-D6_S5_L003_R2_001.fastq.gz | awk 'NR%4==2' | head -1000000 | cut -c 9-24 | sort | uniq | cat - ATAC.737K-arc-v1.txt | sort | uniq -d | wc -l 
zcat SCG0074-15872_SI-NA-D6_S5_L003_R2_001.fastq.gz | awk 'NR%4==2' | head -1000000 | cut -c 1-16 | tr ACGT TGCA | rev | sort | uniq | cat - ATAC.737K-arc-v1.txt | sort | uniq -d | wc -l 
zcat SCG0074-15872_SI-NA-D6_S5_L003_R2_001.fastq.gz | awk 'NR%4==2' | head -1000000 | cut -c 9-24 | tr ACGT TGCA | rev | sort | uniq | cat - ATAC.737K-arc-v1.txt | sort | uniq -d | wc -l
```