# STARsolo

singleCellTK
SCTK-QC
https://www.nature.com/articles/s41467-022-29212-9
在 STAR 的 2.7 版本中 (2.7.6a) 出现了 STARsolo，可以进行单细胞数据的比对
需要注意 ReadFilesIn 先读入测序数据，再读入 barcode+UMI 文件，即先读入 R2 再读入 R1

#### code
[ref1](https://cloud.tencent.com/developer/article/2366729), [ref2](https://zhuanlan.zhihu.com/p/362727395)
```
## index
STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /data/users/Alignment/index/mm10 --genomeFastaFiles /data/users/Alignment/mm10/mm10.fa --sjdbGTFfile /data/users/Alignment/mm10/mm10.gtf
```

10X barcode
![](../img/scRNAbarcode.png)

```
## align
STAR --runThreadN 10 --genomeDir /data/users/Alignment/index/mm10 --readFilesIn $R1 $R2 
--runDirPerm All_RWX $GZIP $BAM --soloBarcodeMate 1 --clip5pNbases 39 0 --soloType CB_UMI_Simple --soloCBwhitelist $BC --soloCBstart 1 --soloCBlen $CBLEN --soloUMIstart $((CBLEN+1)) --soloUMIlen $UMILEN --soloStrand Forward --soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloCellFilter EmptyDrops_CR 
--outFilterScoreMin 30 
--soloFeatures Gene Velocyto 
--soloOutFileNames output/ features.tsv barcodes.tsv matrix.mtx --soloMultiMappers EM --outReadsUnmapped Fastx



STAR --runThreadN 10 --soloType CB_UMI_Simple --soloFeatures Gene Velocyto --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 --soloBarcodeReadLength 0 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS nM CR CY UR UY --soloOutFileNames Solo.out/ genes.tsv barcodes.tsv matrix.mtx --outFileNamePrefix ./out --genomeDir index/mm10 --soloCBwhitelist /software/STAR/BC/3M-february-2018.txt --readFilesIn CRRXXX_r2.fastq.gz CRRXXX_f1.fastq.gz

STAR  --genomeDir ghg38/ --readFilesCommand zcat --readFilesIn SRR7722939/SRR7722939_S1_L001_R2_001.fastq.gz SRR7722939/SRR7722939_S1_L001_R1_001.fastq.gz --soloType CB_UMI_Simple --soloCBwhitelist whitelist --runThreadN 8




我的是ICELL8平台测序，所用的命令如下。。
STAR --genomeDir /data/users/minmingw/Alignment/index --readFilesIn SRR11050949_2.fastq SRR11050949_1.fastq --soloCBwhitelist None --runThreadN 8 --soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 11 --soloUMIstart 12 --soloUMIlen 10 --soloBarcodeReadLength 26 --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM --outSAMtype BAM SortedByCoordinate --soloCellFilter EmptyDrops_CR --outFileNamePrefix SRR11050949
多个Lane: 注意fastq之间的空格
STAR --genomeDir /data/users/minmingw/Alignment/index --readFilesIn Chip91497_lane1_2.fastq,Chip91497_lane2_2.fastq,Chip91497_lane3_2.fastq Chip91497_lane1_1.fastq,Chip91497_lane2_1.fastq,Chip91497_lane3_1.fastq --soloCBwhitelist None --runThreadN 8 --soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 11 --soloUMIstart 12 --soloUMIlen 10 --soloBarcodeReadLength 26 --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM --outSAMtype BAM SortedByCoordinate --soloCellFilter EmptyDrops_CR --outFileNamePrefix Chip91497
```
