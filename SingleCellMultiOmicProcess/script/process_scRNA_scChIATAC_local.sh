Rbc=$1
Rrd=$2
Dbc=$3
Drd1=$4
Drd2=$5
outdir=$6
outprefix=$7
config=$8
source $config

SCRIPT_DIR=$(dirname "$0")
SCRIPT_DIR=$(realpath $SCRIPT_DIR)
source $SCRIPT_DIR/function.sh

RbcBName=$(basename $(realpath $Rbc))
ln -s -f $Rbc $outdir/$RbcBName
RrdBName=$(basename $(realpath $Rrd))
ln -s -f $Rrd $outdir/$RrdBName
DbcBName=$(basename $(realpath $Dbc))
ln -s $Dbc $outdir/$DbcBName
Drd1BName=$(basename $(realpath $Drd1))
ln -s $Drd1 $outdir/$Drd1BName
Drd2BName=$(basename $(realpath $Drd2))
ln -s $Drd2 $outdir/$Drd2BName

cat << ENDHERE > ${outdir}/${outprefix}.scRNA.header
#!/usr/bin/bash
ENDHERE

RNAsoftware=$(echo "$software" | awk -F"+" '{print $1}')
ATACsoftware=$(echo "$software" | awk -F"+" '{print $2}')
rnasample=$(echo ${RrdBName} | awk -F'_' '{for(i=1;i<=NF-4;i++) printf "%s%s", $i, (i==NF-4?"":"_")}')
case $RNAsoftware in
10X)
    jobname=${outprefix}.scRNA
    cat << ENDHERE > ${outdir}/${outprefix}.scRNA.taskrun

set -e
cd $outdir
start=\`date '+%s.%N'\`

echo \$(date) cellranger align >> jobdir/$jobname.log
/usr/bin/time -f 'real %E\nuser %U\nsys  %S\nMmem %M KB\ncpuP %P\nfsi %I\nfso %O' -o jobdir/$jobname.cellranger.time.log $cellranger count --id RNAalign --transcriptome $GEXref --fastqs=$outdir --include-introns=true --localcores=$NTHREAD --localmem=$GB --sample $rnasample --chemistry=ARC-v1 > jobdir/$jobname.cellranger.log 2> jobdir/$jobname.cellranger.log.err
mkdir -p RNA/raw_matrix
mkdir -p RNA/filtered_matrix
cp RNAalign/outs/raw_feature_bc_matrix/* RNA/raw_matrix
cp RNAalign/outs/filtered_feature_bc_matrix/* RNA/filtered_matrix
$Rscript $statRNA RNAalign/outs/metrics_summary.csv 10X > RNA/${outprefix}.stat

$bamCoverage -b RNAalign/outs/possorted_genome_bam.bam --minMappingQuality $mapq -o RNA/${outprefix}.gex.bw -bs 5 -p ${NTHREAD}
$bamCoverage -b RNAalign/outs/possorted_genome_bam.bam --minMappingQuality $mapq -o RNA/${outprefix}.gex.RPKM.bw -bs 5 -p ${NTHREAD} --normalizeUsing RPKM
$computeMatrix scale-regions --quiet -S RNA/${outprefix}.gex.RPKM.bw -R ${genebed} -a 3000 -b 3000 -m 3000 --skipZeros --missingDataAsZero --outFileName RNA/${outprefix}.gex.RPKM.genebody.gz -p ${NTHREAD}
$plotHeatmap -m RNA/${outprefix}.gex.RPKM.genebody.gz -out RNA/${outprefix}.gex.RPKM.genebody.Heatmap.pdf &
$plotProfile -m RNA/${outprefix}.gex.RPKM.genebody.gz -out RNA/${outprefix}.gex.RPKM.genebody.Profile.pdf
wait

echo \$(date) cellranger align finish >> jobdir/$jobname.log
dur=\`echo "\$(date +%s.%N) - \$start" | bc\`
printf "Execution time for running the scRNA pipline : %.2f seconds\n" \$dur >> jobdir/$jobname.log
ENDHERE
        
    jobnamefile=${outdir}/jobdir/${outprefix}.scRNA.sh
    cat ${outdir}/${outprefix}.scRNA.header ${outdir}/${outprefix}.scRNA.taskrun > $jobnamefile
    rm ${outdir}/${outprefix}.scRNA.taskrun
    runid=$(submitjobOptionalShell $jobnamefile)
    echo "mapGEX10x $runid $jobnamefile"
;;

STARsolo)
    var1=$(echo "$RNA_b_u" | awk '{print $1}') # barcode loci
    var2=$(echo "$RNA_b_u" | awk '{print $2}') # barcode length
    var3=$(echo "$RNA_b_u" | awk '{print $3}') # umi loci
    var4=$(echo "$RNA_b_u" | awk '{print $4}') # umi length
    var5=$(echo "$RNA_b_u" | awk '{print $5}') # Forward or Reverse
    jobname=${outprefix}.scRNA
    cat << ENDHERE > ${outdir}/${outprefix}.scRNA.taskrun

set -e
cd $outdir
start=\`date '+%s.%N'\`

echo \$(date) STARsolo align >> jobdir/$jobname.log
mkdir -p RNAalign
/usr/bin/time -f 'real %E\nuser %U\nsys  %S\nMmem %M KB\ncpuP %P\nfsi %I\nfso %O' -o jobdir/$jobname.STARsolo.time.log $STARsolo --genomeDir $STARref --outFileNamePrefix RNAalign/STAR --readFilesCommand zcat --readFilesIn $outdir/$RrdBName $outdir/$RbcBName --soloCBstart $var1 --soloCBlen $var2 --soloUMIstart $var3 --soloUMIlen $var4 --soloType CB_UMI_Simple --soloCBwhitelist $RNAbarcodefile --soloCellFilter EmptyDrops_CR --runThreadN $NTHREAD --outSAMattributes NH HI AS nM NM MD jM jI XS ch cN CR CY UR UY GX GN sM CB UB --outSAMtype BAM SortedByCoordinate --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --outSAMunmapped Within --soloStrand $var5 --chimOutType Junctions WithinBAM --chimSegmentMin 12 --chimJunctionOverhangMin 8 --chimOutJunctionFormat 1 --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --alignIntronMax 100000 --alignSJstitchMismatchNmax 5 -1 5 5 --chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG -4 --chimMultimapNmax 20 --chimNonchimScoreDropMin 10 --peOverlapNbasesMin 12 --peOverlapMMp 0.1 --alignInsertionFlush Right --alignSplicedMateMapLminOverLmate 0 --alignSplicedMateMapLmin 30 --soloFeatures Gene GeneFull SJ Velocyto > jobdir/$jobname.STARsolo.log 2> jobdir/$jobname.STARsolo.log.err

mkdir -p RNA/raw_matrix
mkdir -p RNA/filtered_matrix
awk '{print \$0"-1"}' RNAalign/STARSolo.out/Gene/raw/barcodes.tsv | $pigz -p $NTHREAD > RNA/raw_matrix/barcodes.tsv.gz
$pigz -p $NTHREAD -c RNAalign/STARSolo.out/Gene/raw/features.tsv > RNA/raw_matrix/features.tsv.gz
$pigz -p $NTHREAD -c RNAalign/STARSolo.out/Gene/raw/matrix.mtx > RNA/raw_matrix/matrix.mtx.gz
awk '{print \$0"-1"}' RNAalign/STARSolo.out/Gene/filtered/barcodes.tsv | $pigz -p $NTHREAD > RNA/filtered_matrix/barcodes.tsv.gz
$pigz -p $NTHREAD -c RNAalign/STARSolo.out/Gene/filtered/features.tsv > RNA/filtered_matrix/features.tsv.gz
$pigz -p $NTHREAD -c RNAalign/STARSolo.out/Gene/filtered/matrix.mtx > RNA/filtered_matrix/matrix.mtx.gz

$Rscript $statRNA RNAalign/STARSolo.out/Gene/Summary.csv STARsolo > RNA/${outprefix}.stat

$bamCoverage -b RNAalign/STARAligned.sortedByCoord.out.bam --minMappingQuality $mapq -o RNA/${outprefix}.gex.bw -bs 5 -p ${NTHREAD}
$bamCoverage -b RNAalign/STARAligned.sortedByCoord.out.bam --minMappingQuality $mapq -o RNA/${outprefix}.gex.RPKM.bw -bs 5 -p ${NTHREAD} --normalizeUsing RPKM
$computeMatrix scale-regions --quiet -S RNA/${outprefix}.gex.RPKM.bw -R ${genebed} -a 3000 -b 3000 -m 3000 --skipZeros --missingDataAsZero --outFileName RNA/${outprefix}.gex.RPKM.genebody.gz -p ${NTHREAD}
$plotHeatmap -m RNA/${outprefix}.gex.RPKM.genebody.gz -out RNA/${outprefix}.gex.RPKM.genebody.Heatmap.pdf &
$plotProfile -m RNA/${outprefix}.gex.RPKM.genebody.gz -out RNA/${outprefix}.gex.RPKM.genebody.Profile.pdf
wait

echo \$(date) STARsolo align finish >> jobdir/$jobname.log
dur=\`echo "\$(date +%s.%N) - \$start" | bc\`
printf "Execution time for running the scRNA pipline : %.2f seconds\n" \$dur >> jobdir/$jobname.log
ENDHERE
        
    jobnamefile=${outdir}/jobdir/${outprefix}.scRNA.sh
    cat ${outdir}/${outprefix}.scRNA.header ${outdir}/${outprefix}.scRNA.taskrun > $jobnamefile
    rm ${outdir}/${outprefix}.scRNA.taskrun
    runid=$(submitjobOptionalShell $jobnamefile)
    echo "mapGEXSTARsolo $runid $jobnamefile"
    ;;
    *)
    echo "Invalid parameter software: $software"
    ;;
esac


case $ATACsoftware in
custom)
    jobname=${outprefix}.scChIATAC
    var1=$(echo "$DNA_b" | awk '{print $1}')
    var2=$(echo "$DNA_b" | awk '{print $2}')
    var3=$(echo "$DNA_b" | awk '{print $3}')
    cat << ENDHERE > ${outdir}/${outprefix}.scChIATAC.taskrun
    
set -e
cd $outdir
start=\`date '+%s.%N'\`

begin=\`date '+%s.%N'\`
echo \$(date) prepare and split fastq >> jobdir/$jobname.log
mkdir -p fastq
/usr/bin/time -f 'real %E\nuser %U\nsys  %S\nMmem %M KB\ncpuP %P\nfsi %I\nfso %O' -o $outdir/fastq/$outprefix.splitFq.time.log $python $splitFq $Dbc $Drd1 $Drd2 $Nchunk $outdir/fastq $outprefix $DNA_b
nfile=\$(ls $outdir/fastq/$outprefix.barcode.*.fastq| wc -l)
dur=\`echo "\$(date +%s.%N) - \$begin" | bc\`
printf "\$(date) prepare and split fastq finish: %.2f seconds\n" \$dur >> jobdir/$jobname.log

begin=\`date '+%s.%N'\`
echo \$(date) trim, linker calling and align >> jobdir/$jobname.log
mkdir -p DNAalign
mkdir -p trim
mkdir -p linker
for i in \$(seq \$nfile); do 
    /usr/bin/time -f 'real %E\nuser %U\nsys  %S\nMmem %M KB\ncpuP %P\nfsi %I\nfso %O' -o $outdir/trim/$outprefix.\${i}.trim.time.log java -jar ${trimmomaticpath}/trimmomatic-0.39.jar PE -threads $NTHREAD -phred33 $outdir/fastq/$outprefix.R1.\${i}.fastq $outdir/fastq/$outprefix.R2.\${i}.fastq $outdir/trim/$outprefix.R1.\${i}.fastq $outdir/trim/$outprefix.R1.\${i}.unpair.fastq $outdir/trim/$outprefix.R2.\${i}.fastq $outdir/trim/$outprefix.R2.\${i}.unpair.fastq LEADING:5 TRAILING:5 ILLUMINACLIP:${trimmomaticpath}/adapters/NexteraPE-PE.fa:2:30:7:8:true MINLEN:30 > $outdir/trim/$outprefix.\${i}.trim.log 2> $outdir/trim/$outprefix.\${i}.trim.log.err 
    grep Surviving $outdir/trim/$outprefix.\${i}.trim.log.err | perl -lane 'if(/Both Surviving: (\d+) /){\$a=\$1;} if(/Forward Only Surviving: (\d+) /){\$b=\$1;} if(/Reverse Only Surviving: (\d+) /){\$c=\$1;} print "Trim_Reads\t".(\$a+\$b+\$c)' > $outdir/trim/$outprefix.\${i}.trim.readsnum

    $cutadapt -j $NTHREAD -b file:${linkerfile} -n 50 --info-file $outdir/linker/cutadapt.R1.\${i}.pair.Linker_info --discard -o $outdir/linker/cutadapt.R1.\${i}.pair.tmp.fastq -O 8 $outdir/trim/$outprefix.R1.\${i}.fastq > $outdir/linker/cutadapt.R1.\${i}.pair.cutadapt.stat 
    $cutadapt -j $NTHREAD -b file:${linkerfile} -n 50 --info-file $outdir/linker/cutadapt.R1.\${i}.unpair.Linker_info --discard -o $outdir/linker/cutadapt.R1.\${i}.unpair.tmp.fastq -O 8 $outdir/trim/$outprefix.R1.\${i}.unpair.fastq > $outdir/linker/cutadapt.R1.\${i}.unpair.cutadapt.stat 
    $cutadapt -j $NTHREAD -b file:${linkerfile} -n 50 --info-file $outdir/linker/cutadapt.R2.\${i}.pair.Linker_info --discard -o $outdir/linker/cutadapt.R2.\${i}.pair.tmp.fastq -O 8 $outdir/trim/$outprefix.R2.\${i}.fastq > $outdir/linker/cutadapt.R2.\${i}.pair.cutadapt.stat 
    $cutadapt -j $NTHREAD -b file:${linkerfile} -n 50 --info-file $outdir/linker/cutadapt.R2.\${i}.unpair.Linker_info --discard -o $outdir/linker/cutadapt.R2.\${i}.unpair.tmp.fastq -O 8 $outdir/trim/$outprefix.R2.\${i}.unpair.fastq > $outdir/linker/cutadapt.R2.\${i}.unpair.cutadapt.stat 
    $python $getReadPro $outdir/linker/cutadapt.R1.\${i}.pair.Linker_info $outdir/linker/cutadapt.R1.\${i}.unpair.Linker_info $outdir/linker/cutadapt.R2.\${i}.pair.Linker_info $outdir/linker/cutadapt.R2.\${i}.unpair.Linker_info $outdir/linker/cutadapt.\${i} 18 T > $outdir/linker/cutadapt.\${i}.splitReads.log
    cat $outdir/linker/cutadapt.R1.\${i}.pair.Linker_info  $outdir/linker/cutadapt.R1.\${i}.unpair.Linker_info  $outdir/linker/cutadapt.R2.\${i}.pair.Linker_info  $outdir/linker/cutadapt.R2.\${i}.unpair.Linker_info | awk -F"\t" '\$2!=-1' | awk '{aa[\$1]}END{print "Linker_Reads\t"length(aa)}' > $outdir/linker/cutadapt.\${i}.linker.readsnum
    $pigz -p $NTHREAD -f $outdir/linker/cutadapt.\${i}.R1.pair.long.fastq $outdir/linker/cutadapt.\${i}.R1.unpair.long.fastq $outdir/linker/cutadapt.\${i}.R2.pair.long.fastq $outdir/linker/cutadapt.\${i}.R2.unpair.long.fastq $outdir/linker/cutadapt.\${i}.R1.pair.short.fastq $outdir/linker/cutadapt.\${i}.R1.unpair.short.fastq $outdir/linker/cutadapt.\${i}.R2.pair.short.fastq $outdir/linker/cutadapt.\${i}.R2.unpair.short.fastq

    /usr/bin/time -f 'real %E\nuser %U\nsys  %S\nMmem %M KB\ncpuP %P\nfsi %I\nfso %O' -o DNAalign/STAR\${i}.time.log $STARsolo --genomeDir $barcoderef --outFileNamePrefix DNAalign/STAR\${i} --readFilesIn $outdir/fastq/$outprefix.R1.\${i}.tmp.fastq $outdir/fastq/$outprefix.barcode.\${i}.fastq --soloCBstart 1 --soloCBlen $var2 --soloUMIstart 1 --soloUMIlen 1 --soloType CB_samTagOut --soloCBwhitelist $DNAbarcodefile --runThreadN $NTHREAD --outSAMattributes CB --outSAMtype BAM Unsorted --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloCBmatchWLtype 1MM --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --outSAMunmapped Within --soloStrand Forward > DNAalign/STAR\${i}.log 2> DNAalign/STAR\${i}.log.err
    grep yes DNAalign/STAR\${i}Solo.out/Barcodes.stats | awk '{sum+=\$2}END{print "barcodeReads\t"sum}' > DNAalign/STAR\${i}.barcodereadsnum

    $bwa mem -5 -t $NTHREAD $bwaref $outdir/linker/cutadapt.\${i}.R1.pair.long.fastq.gz | $python $addMark R1 mem | $samtools sort -n -@ $NTHREAD -o $outdir/DNAalign/bwa\${i}.R1.pair.long.bam -
    $bwa mem -5 -t $NTHREAD $bwaref $outdir/linker/cutadapt.\${i}.R1.unpair.long.fastq.gz | $python $addMark R1 mem | $samtools sort -n -@ $NTHREAD -o $outdir/DNAalign/bwa\${i}.R1.unpair.long.bam -
    $bwa mem -5 -t $NTHREAD $bwaref $outdir/linker/cutadapt.\${i}.R2.pair.long.fastq.gz | $python $addMark R2 mem | $samtools sort -n -@ $NTHREAD -o $outdir/DNAalign/bwa\${i}.R2.pair.long.bam -
    $bwa mem -5 -t $NTHREAD $bwaref $outdir/linker/cutadapt.\${i}.R2.unpair.long.fastq.gz | $python $addMark R2 mem | $samtools sort -n -@ $NTHREAD -o $outdir/DNAalign/bwa\${i}.R2.unpair.long.bam -
    $bwa aln -t $NTHREAD -f $outdir/DNAalign/bwa\${i}.R1.pair.short.sai $bwaref $outdir/linker/cutadapt.\${i}.R1.pair.short.fastq.gz && $bwa samse $bwaref $outdir/DNAalign/bwa\${i}.R1.pair.short.sai $outdir/linker/cutadapt.\${i}.R1.pair.short.fastq.gz | $python $addMark R1 aln | $samtools sort -n -@ $NTHREAD -o $outdir/DNAalign/bwa\${i}.R1.pair.short.bam -
    $bwa aln -t $NTHREAD -f $outdir/DNAalign/bwa\${i}.R1.unpair.short.sai $bwaref $outdir/linker/cutadapt.\${i}.R1.unpair.short.fastq.gz && $bwa samse $bwaref $outdir/DNAalign/bwa\${i}.R1.unpair.short.sai $outdir/linker/cutadapt.\${i}.R1.unpair.short.fastq.gz | $python $addMark R1 aln | $samtools sort -n -@ $NTHREAD -o $outdir/DNAalign/bwa\${i}.R1.unpair.short.bam -
    $bwa aln -t $NTHREAD -f $outdir/DNAalign/bwa\${i}.R2.pair.short.sai $bwaref $outdir/linker/cutadapt.\${i}.R2.pair.short.fastq.gz && $bwa samse $bwaref $outdir/DNAalign/bwa\${i}.R2.pair.short.sai $outdir/linker/cutadapt.\${i}.R2.pair.short.fastq.gz | $python $addMark R2 aln | $samtools sort -n -@ $NTHREAD -o $outdir/DNAalign/bwa\${i}.R2.pair.short.bam - 
    $bwa aln -t $NTHREAD -f $outdir/DNAalign/bwa\${i}.R2.unpair.short.sai $bwaref $outdir/linker/cutadapt.\${i}.R2.unpair.short.fastq.gz && $bwa samse $bwaref $outdir/DNAalign/bwa\${i}.R2.unpair.short.sai $outdir/linker/cutadapt.\${i}.R2.unpair.short.fastq.gz | $python $addMark R2 aln | $samtools sort -n -@ $NTHREAD -o $outdir/DNAalign/bwa\${i}.R2.unpair.short.bam -
    
    if [ ${blacklist} != "null" ]; then
        $samtools merge -@ $NTHREAD -n -f -o - $outdir/DNAalign/bwa\${i}.R1.unpair.long.bam $outdir/DNAalign/bwa\${i}.R1.unpair.short.bam $outdir/DNAalign/bwa\${i}.R2.unpair.long.bam $outdir/DNAalign/bwa\${i}.R2.unpair.short.bam $outdir/DNAalign/bwa\${i}.R1.pair.long.bam $outdir/DNAalign/bwa\${i}.R1.pair.short.bam $outdir/DNAalign/bwa\${i}.R2.pair.long.bam $outdir/DNAalign/bwa\${i}.R2.pair.short.bam | $bedtools intersect -v -nonamecheck -a - -b ${blacklist} > $outdir/DNAalign/bwa\${i}.all.bam
    else
        $samtools merge -@ $NTHREAD -n -f -o $outdir/DNAalign/bwa\${i}.all.bam $outdir/DNAalign/bwa\${i}.R1.unpair.long.bam $outdir/DNAalign/bwa\${i}.R1.unpair.short.bam $outdir/DNAalign/bwa\${i}.R2.unpair.long.bam $outdir/DNAalign/bwa\${i}.R2.unpair.short.bam $outdir/DNAalign/bwa\${i}.R1.pair.long.bam $outdir/DNAalign/bwa\${i}.R1.pair.short.bam $outdir/DNAalign/bwa\${i}.R2.pair.long.bam $outdir/DNAalign/bwa\${i}.R2.pair.short.bam
    fi
    
    $python $get2pairFinal $outdir/DNAalign/bwa\${i}.all.bam $outdir/DNAalign/bwa\${i}.get2pair $mapq > $outdir/DNAalign/bwa\${i}.get2pair.log

    $samtools fixmate -c -m $outdir/DNAalign/bwa\${i}.get2pair.pair.withlinker.bam - | $samtools sort -@ $NTHREAD -o $outdir/DNAalign/bwa\${i}.withlinker.bam -
    $samtools fixmate -c -m $outdir/DNAalign/bwa\${i}.get2pair.pair.withoutlinker.bam - | $samtools sort -@ $NTHREAD -o $outdir/DNAalign/bwa\${i}.withoutlinker.bam -
    $samtools view -@ $NTHREAD $outdir/DNAalign/bwa\${i}.get2pair.pair.withlinker.bam | wc -l | awk '{print "unique_mapping_reads_with_linker\t"\$1/2}' > $outdir/DNAalign/bwa\${i}.withlinker.readsnum
    $samtools view -@ $NTHREAD $outdir/DNAalign/bwa\${i}.get2pair.pair.withoutlinker.bam | wc -l | awk '{print "unique_mapping_reads_without_linker\t"\$1/2}' > $outdir/DNAalign/bwa\${i}.withoutlinker.readsnum

    if [ ${barcodetranslate} == none ];then
        $python $giveBC -sb DNAalign/STAR\${i}Aligned.out.bam -bb DNAalign/bwa\${i}.withlinker.bam -o DNAalign/bwa\${i}.withlinker.addBC.bam -s DNAalign/bwa\${i}.withlinker.addBC.barcodeuniquereadsnum && rm DNAalign/bwa\${i}.withlinker.addBC.barcodeuniquereadsnum
        $python $giveBC -sb DNAalign/STAR\${i}Aligned.out.bam -bb DNAalign/bwa\${i}.withoutlinker.bam -o DNAalign/bwa\${i}.withoutlinker.addBC.bam -s DNAalign/bwa\${i}.withoutlinker.addBC.barcodeuniquereadsnum && rm DNAalign/bwa\${i}.withoutlinker.addBC.barcodeuniquereadsnum
    else
        $python $giveBC -sb DNAalign/STAR\${i}Aligned.out.bam -bb DNAalign/bwa\${i}.withlinker.bam -o DNAalign/bwa\${i}.withlinker.addBC.bam -s DNAalign/bwa\${i}.withlinker.addBC.barcodeuniquereadsnum -t $barcodetranslate && rm DNAalign/bwa\${i}.withlinker.addBC.barcodeuniquereadsnum
        $python $giveBC -sb DNAalign/STAR\${i}Aligned.out.bam -bb DNAalign/bwa\${i}.withoutlinker.bam -o DNAalign/bwa\${i}.withoutlinker.addBC.bam -s DNAalign/bwa\${i}.withoutlinker.addBC.barcodeuniquereadsnum -t $barcodetranslate && rm DNAalign/bwa\${i}.withoutlinker.addBC.barcodeuniquereadsnum
    fi
    
    rm $outdir/fastq/$outprefix.R1.\${i}.fastq $outdir/fastq/$outprefix.R2.\${i}.fastq
    rm $outdir/linker/cutadapt.R1.\${i}.pair.tmp.fastq $outdir/trim/$outprefix.R1.\${i}.fastq
    rm $outdir/linker/cutadapt.R1.\${i}.unpair.tmp.fastq $outdir/trim/$outprefix.R1.\${i}.unpair.fastq
    rm $outdir/linker/cutadapt.R2.\${i}.pair.tmp.fastq $outdir/trim/$outprefix.R2.\${i}.fastq
    rm $outdir/linker/cutadapt.R2.\${i}.unpair.tmp.fastq $outdir/trim/$outprefix.R2.\${i}.unpair.fastq
    rm $outdir/linker/cutadapt.R1.\${i}.pair.Linker_info $outdir/linker/cutadapt.R1.\${i}.unpair.Linker_info $outdir/linker/cutadapt.R2.\${i}.pair.Linker_info $outdir/linker/cutadapt.R2.\${i}.unpair.Linker_info
    rm $outdir/fastq/$outprefix.R1.\${i}.tmp.fastq $outdir/fastq/$outprefix.barcode.\${i}.fastq $outdir/DNAalign/STAR\${i}SJ.out.tab
    rm $outdir/DNAalign/bwa\${i}.R1.pair.short.sai $outdir/DNAalign/bwa\${i}.R1.unpair.short.sai $outdir/DNAalign/bwa\${i}.R2.pair.short.sai $outdir/DNAalign/bwa\${i}.R2.unpair.short.sai
    rm $outdir/DNAalign/bwa\${i}.R1.unpair.long.bam $outdir/DNAalign/bwa\${i}.R1.unpair.short.bam $outdir/DNAalign/bwa\${i}.R2.unpair.long.bam $outdir/DNAalign/bwa\${i}.R2.unpair.short.bam $outdir/DNAalign/bwa\${i}.R1.pair.long.bam $outdir/DNAalign/bwa\${i}.R1.pair.short.bam $outdir/DNAalign/bwa\${i}.R2.pair.long.bam $outdir/DNAalign/bwa\${i}.R2.pair.short.bam
    rm $outdir/DNAalign/bwa\${i}.all.bam
    rm $outdir/DNAalign/bwa\${i}.get2pair.pair.withlinker.bam $outdir/DNAalign/bwa\${i}.get2pair.pair.withoutlinker.bam $outdir/DNAalign/bwa\${i}.get2pair.single.withlinker.bam $outdir/DNAalign/bwa\${i}.get2pair.single.withoutlinker.bam
    rm $outdir/DNAalign/bwa\${i}.withlinker.bam $outdir/DNAalign/bwa\${i}.withoutlinker.bam
    rm -rf $outdir/DNAalign/STAR\${i}Solo.out
    rm $outdir/DNAalign/STAR\${i}Log.progress.out $outdir/DNAalign/STAR\${i}Log.out $outdir/DNAalign/STAR\${i}Log.final.out
    rm $outdir/linker/cutadapt.R1.\${i}.pair.cutadapt.stat $outdir/linker/cutadapt.R1.\${i}.unpair.cutadapt.stat $outdir/linker/cutadapt.R2.\${i}.pair.cutadapt.stat $outdir/linker/cutadapt.R2.\${i}.unpair.cutadapt.stat
    rm $outdir/linker/cutadapt.\${i}.R1.pair.long.fastq.gz $outdir/linker/cutadapt.\${i}.R1.unpair.long.fastq.gz $outdir/linker/cutadapt.\${i}.R2.pair.long.fastq.gz $outdir/linker/cutadapt.\${i}.R2.unpair.long.fastq.gz $outdir/linker/cutadapt.\${i}.R1.pair.short.fastq.gz $outdir/linker/cutadapt.\${i}.R1.unpair.short.fastq.gz $outdir/linker/cutadapt.\${i}.R2.pair.short.fastq.gz $outdir/linker/cutadapt.\${i}.R2.unpair.short.fastq.gz

done
wait
dur=\`echo "\$(date +%s.%N) - \$begin" | bc\`
printf "\$(date) trim, linker calling and align finish: %.2f seconds\n" \$dur >> jobdir/$jobname.log

begin=\`date '+%s.%N'\`
echo \$(date) merge align result and remove PCR >> jobdir/$jobname.log
$samtools merge -@ $NTHREAD -f -o DNAalign/$outprefix.withoutlinker.addBC.bam \$(for i in \$(seq \$nfile); do echo -e "$outdir/DNAalign/bwa\${i}.withoutlinker.addBC.bam \c"; done)
$samtools merge -@ $NTHREAD -f -o DNAalign/$outprefix.withlinker.addBC.bam \$(for i in \$(seq \$nfile); do echo -e "$outdir/DNAalign/bwa\${i}.withlinker.addBC.bam \c"; done)

cat DNAalign/STAR*.barcodereadsnum | awk '{sum+=\$2}END{print "Barcode_Reads\t"sum}' > DNAalign/STAR.barcodereadsnum
cat trim/$outprefix.*.trim.readsnum | awk '{sum+=\$2}END{print "Trim_Reads\t"sum}' > trim/$outprefix.trim.readsnum
cat linker/cutadapt.*.linker.readsnum | awk '{sum+=\$2}END{print "Linker_Reads\t"sum}' > linker/cutadapt.linker.readsnum

cat DNAalign/bwa*.withlinker.readsnum | awk '{sum+=\$2}END{print "unique_mapping_reads_with_linker\t"sum}' > DNAalign/bwa.withlinker.readsnum
cat DNAalign/bwa*.withoutlinker.readsnum | awk '{sum+=\$2}END{print "unique_mapping_reads_without_linker\t"sum}' > DNAalign/bwa.withoutlinker.readsnum

cat <(paste linker/cutadapt.*.splitReads.log | cut -f 1 | paste \$(for i in \$(seq \$(paste linker/cutadapt.*.splitReads.log | wc -l)); do echo -e "- \c"; done)) <(awk '{print \$2}' linker/cutadapt.*.splitReads.log | paste \$(for i in \$(seq \$(paste linker/cutadapt.*.splitReads.log | wc -l)); do echo -e "- \c"; done) ) | awk '{for(i=1;i<=NF;i++){a[FNR,i]=\$i}}END{for(i=1;i<=NF;i++){for(j=1;j<=FNR;j++){printf a[j,i]"\t"}print ""}}' - | awk '{for(i=2;i<=NF;i++){sum+=\$i; print \$1"\t"sum}}' > linker/cutadapt.splitReads.log
cat <(paste DNAalign/bwa*.get2pair.log | cut -f 1 | paste \$(for i in \$(seq \$(paste DNAalign/bwa*.get2pair.log | wc -l)); do echo -e "- \c"; done)) <(awk '{print \$2}' DNAalign/bwa*.get2pair.log | paste \$(for i in \$(seq \$(paste DNAalign/bwa*.get2pair.log | wc -l)); do echo -e "- \c"; done) ) | awk '{for(i=1;i<=NF;i++){a[FNR,i]=\$i}}END{for(i=1;i<=NF;i++){for(j=1;j<=FNR;j++){printf a[j,i]"\t"}print ""}}' - | awk '{for(i=2;i<=NF;i++){sum+=\$i; print \$1"\t"sum}}' > DNAalign/bwa.get2pair.log

$gatk MarkDuplicates --BARCODE_TAG=ZB --REMOVE_DUPLICATES=true --ADD_PG_TAG_TO_READS=false -M DNAalign/$outprefix.withlinker.addBC.dedup.MarkDuplicates.txt -O DNAalign/$outprefix.withlinker.addBC.dedup.bam -I DNAalign/$outprefix.withlinker.addBC.bam --ASSUME_SORT_ORDER=coordinate > DNAalign/$outprefix.withlinker.addBC.dedup.log 2>&1
$gatk MarkDuplicates --BARCODE_TAG=ZB --REMOVE_DUPLICATES=true --ADD_PG_TAG_TO_READS=false -M DNAalign/$outprefix.withoutlinker.addBC.dedup.MarkDuplicates.txt -O DNAalign/$outprefix.withoutlinker.addBC.dedup.bam -I DNAalign/$outprefix.withoutlinker.addBC.bam --ASSUME_SORT_ORDER=coordinate > DNAalign/$outprefix.withoutlinker.addBC.dedup.log 2>&1
awk '{if(/Unknown Library/){print "barcode_unique_mapping_reads_without_linker\t"\$4}}' DNAalign/$outprefix.withoutlinker.addBC.dedup.MarkDuplicates.txt > DNAalign/bwa.withoutlinker.addBC.barcodeuniquereadsnum
awk '{if(/Unknown Library/){print "barcode_unique_mapping_rmPCR_reads_without_linker\t"\$8}}' DNAalign/$outprefix.withoutlinker.addBC.dedup.MarkDuplicates.txt > DNAalign/bwa.withoutlinker.addBC.barcodeuniquermPCRreadsnum
awk '{if(/Unknown Library/){print "barcode_unique_mapping_reads_with_linker\t"\$4}}' DNAalign/$outprefix.withlinker.addBC.dedup.MarkDuplicates.txt > DNAalign/bwa.withlinker.addBC.barcodeuniquereadsnum
awk '{if(/Unknown Library/){print "barcode_unique_mapping_rmPCR_reads_with_linker\t"\$8}}' DNAalign/$outprefix.withlinker.addBC.dedup.MarkDuplicates.txt > DNAalign/bwa.withlinker.addBC.barcodeuniquermPCRreadsnum
$samtools sort -n -@ $NTHREAD -m 1000M -o DNAalign/$outprefix.withlinker.bam DNAalign/$outprefix.withlinker.addBC.dedup.bam
$samtools sort -n -@ $NTHREAD -m 1000M -o DNAalign/$outprefix.withoutlinker.bam DNAalign/$outprefix.withoutlinker.addBC.dedup.bam
$python $classPET DNAalign/$outprefix.withlinker.bam $contactMinDis DNAalign/$outprefix.withlinker.FRlt.bam DNAalign/$outprefix.withlinker.TANDEM.bam DNAalign/$outprefix.withlinker.scPET.bam DNAalign/$outprefix.withlinker.classPET.log > DNAalign/$outprefix.withlinker.classPET.num.log
$python ${classPETplot} DNAalign/$outprefix.withlinker.classPET.num.log &
$python $classPET DNAalign/$outprefix.withoutlinker.bam $contactMinDis DNAalign/$outprefix.withoutlinker.FRlt.bam DNAalign/$outprefix.withoutlinker.TANDEM.bam DNAalign/$outprefix.withoutlinker.scPET.bam DNAalign/$outprefix.withoutlinker.classPET.log > DNAalign/$outprefix.withoutlinker.classPET.num.log
$python ${classPETplot} DNAalign/$outprefix.withoutlinker.classPET.num.log &
$samtools merge -f -@ $NTHREAD -n -o DNAalign/$outprefix.FRlt.bam DNAalign/$outprefix.withoutlinker.FRlt.bam DNAalign/$outprefix.withlinker.FRlt.bam
$samtools merge -f -@ $NTHREAD -n -o DNAalign/$outprefix.scPET.bam DNAalign/$outprefix.withoutlinker.scPET.bam DNAalign/$outprefix.withlinker.scPET.bam

for i in \$(seq \$nfile); do rm $outdir/DNAalign/bwa\${i}.withoutlinker.addBC.bam $outdir/DNAalign/bwa\${i}.withlinker.addBC.bam $outdir/DNAalign/STAR\${i}Aligned.out.bam $outdir/DNAalign/STAR\${i}.barcodereadsnum $outdir/DNAalign/bwa\${i}.withlinker.readsnum $outdir/DNAalign/bwa\${i}.withoutlinker.readsnum $outdir/DNAalign/bwa\${i}.get2pair.log $outdir/linker/cutadapt.\${i}.linker.readsnum $outdir/linker/cutadapt.\${i}.splitReads.log $outdir/trim/$outprefix.\${i}.trim.readsnum; done

dur=\`echo "\$(date +%s.%N) - \$begin" | bc\`
printf "\$(date) merge align result and remove PCR finish: %.2f seconds\n" \$dur >> jobdir/$jobname.log

begin=\`date '+%s.%N'\`
echo \$(date) generate contact >> jobdir/$jobname.log
mkdir -p PET/raw_matrix
mkdir -p PET/filtered_matrix
/usr/bin/time -f 'real %E\nuser %U\nsys  %S\nMmem %M KB\ncpuP %P\nfsi %I\nfso %O' -o PET/${outprefix}.bam2bedpe.time.log $python $bam2bedpe -b DNAalign/$outprefix.scPET.bam -o PET/$outprefix.scPET.alldata.bedpe -q $mapq -d PET/raw_matrix > PET/${outprefix}.bam2bedpe.log
$python $CountBarcodePET PET/$outprefix.scPET.alldata.bedpe PET/$outprefix.scPET.alldata.bedpe.stat
$python ${forBandNorm} PET/$outprefix.scPET.alldata.bedpe PET/$outprefix.scPET.alldata.bedpe.mat 1000 $outprefix
$Rscript ${petmat} PET/$outprefix.scPET.alldata.bedpe.mat PET/$outprefix.scPET.alldata.bedpe.mat.svd.rds

mkdir -p tmpdir
bash $sortHuge4dn PET/$outprefix.scPET.alldata.bedpe.4dn PET/$outprefix.scPET.alldata.bedpe.sorted.4dn tmpdir $NTHREAD
$juicer_tools pre -r 2500000,2000000,1000000,500000,250000,200000,100000,50000,25000,20000,10000,5000,2500,2000,1000 PET/$outprefix.scPET.alldata.bedpe.sorted.4dn PET/$outprefix.scPET.alldata.bedpe.sorted.4dn.hic ${chrom_sizes} > PET/$outprefix.scPET.alldata.bedpe.sorted.4dn.hic.log 2>&1 &

$STARsolo --runMode soloCellFiltering PET/raw_matrix/ PET/filtered_matrix/tmp --soloCellFilter EmptyDrops_CR
awk '{if(NR==FNR){aa[\$1]}else{if(\$11 in aa){print \$0}}}' PET/filtered_matrix/tmpbarcodes.tsv PET/$outprefix.scPET.alldata.bedpe > PET/$outprefix.scPET.bedpe
$python $bedpe24dn PET/$outprefix.scPET.bedpe | bash $sortHuge4dn - PET/$outprefix.scPET.bedpe.4dn tmpdir $NTHREAD
$juicer_tools pre -r 2500000,2000000,1000000,500000,250000,200000,100000,50000,25000,20000,10000,5000,2500,2000,1000 PET/$outprefix.scPET.bedpe.4dn PET/$outprefix.scPET.bedpe.4dn.hic ${chrom_sizes} > PET/$outprefix.scPET.bedpe.4dn.hic.log 2>&1 &

$pigz -p $NTHREAD PET/raw_matrix/features.tsv PET/raw_matrix/matrix.mtx PET/raw_matrix/barcodes.tsv
mv PET/filtered_matrix/tmpfeatures.tsv PET/filtered_matrix/features.tsv
mv PET/filtered_matrix/tmpmatrix.mtx PET/filtered_matrix/matrix.mtx
mv PET/filtered_matrix/tmpbarcodes.tsv PET/filtered_matrix/barcodes.tsv
$pigz -p $NTHREAD PET/filtered_matrix/features.tsv PET/filtered_matrix/matrix.mtx PET/filtered_matrix/barcodes.tsv
$python $statChAIRPET fastq/$outprefix.readsnum DNAalign/STAR.barcodereadsnum trim/$outprefix.trim.readsnum linker/cutadapt.linker.readsnum DNAalign/bwa.withlinker.readsnum DNAalign/bwa.withoutlinker.readsnum DNAalign/bwa.withlinker.addBC.barcodeuniquereadsnum DNAalign/bwa.withoutlinker.addBC.barcodeuniquereadsnum DNAalign/bwa.withlinker.addBC.barcodeuniquermPCRreadsnum DNAalign/bwa.withoutlinker.addBC.barcodeuniquermPCRreadsnum PET/${outprefix}.bam2bedpe.log PET/filtered_matrix/matrix.mtx.gz custom > PET/${outprefix}.stat
wait
dur=\`echo "\$(date +%s.%N) - \$begin" | bc\`
printf "\$(date) generate contact finish: %.2f seconds\n" \$dur >> jobdir/$jobname.log

begin=\`date '+%s.%N'\`
echo \$(date) generate ATAC signal fragment and peak intensity matrix >> jobdir/$jobname.log
mkdir -p fragment
/usr/bin/time -f 'real %E\nuser %U\nsys  %S\nMmem %M KB\ncpuP %P\nfsi %I\nfso %O' -o fragment/${outprefix}.bam2fragment.time.log $python $bam2fragment -p DNAalign/$outprefix.FRlt.bam -s DNAalign/$outprefix.scPET.bam  -o fragment/${outprefix}_fragments.unsorted.tsv > fragment/${outprefix}_fragments.unsorted.tsv.log
mkdir -p tmpdir
wait
bash $sortHugeBed fragment/${outprefix}_fragments.unsorted.tsv fragment/${outprefix}_fragments.tsv tmpdir $NTHREAD
$bgzip -f fragment/${outprefix}_fragments.tsv
$tabix -f -p bed fragment/${outprefix}_fragments.tsv.gz

mkdir -p peak
$samtools merge -f -@ $NTHREAD -n -o DNAalign/$outprefix.scATAC.bam DNAalign/$outprefix.FRlt.bam DNAalign/$outprefix.scPET.bam
$samtools sort -m 1000M -@ $NTHREAD -o DNAalign/$outprefix.scATAC.sorted.bam DNAalign/$outprefix.scATAC.bam
samtools index DNAalign/$outprefix.scATAC.sorted.bam &
/usr/bin/time -f 'real %E\nuser %U\nsys  %S\nMmem %M KB\ncpuP %P\nfsi %I\nfso %O' -o peak/$outprefix.macs2.time.log $macs2 callpeak -f BAM --keep-dup all --nomodel --shift $shiftsize --extsize $peakext -B --SPMR -t DNAalign/$outprefix.scATAC.sorted.bam -g $glen -n peak/$outprefix --qvalue $macsq > peak/$outprefix.macs2.log 2> peak/$outprefix.macs2.log.err
awk '\$7+0.0 > 1.5' peak/${outprefix}_peaks.narrowPeak | grep -v -i "chrM" | cut -f 1-3 > peak/${outprefix}.peak

mkdir -p ATAC/raw_matrix
mkdir -p ATAC/filtered_matrix
/usr/bin/time -f 'real %E\nuser %U\nsys  %S\nMmem %M KB\ncpuP %P\nfsi %I\nfso %O' -o fragment/${outprefix}.getMat.time.log $python $getMat peak/${outprefix}.peak fragment/${outprefix}_fragments.tsv.gz ATAC/raw_matrix/features.tsv ATAC/raw_matrix/barcodes.tsv ATAC/raw_matrix/matrix.mtx
$STARsolo --runMode soloCellFiltering ATAC/raw_matrix/ ATAC/filtered_matrix/tmp --soloCellFilter EmptyDrops_CR

$pigz -p $NTHREAD ATAC/raw_matrix/features.tsv ATAC/raw_matrix/matrix.mtx ATAC/raw_matrix/barcodes.tsv
mv ATAC/filtered_matrix/tmpfeatures.tsv ATAC/filtered_matrix/features.tsv
mv ATAC/filtered_matrix/tmpmatrix.mtx ATAC/filtered_matrix/matrix.mtx
mv ATAC/filtered_matrix/tmpbarcodes.tsv ATAC/filtered_matrix/barcodes.tsv
$pigz -p $NTHREAD ATAC/filtered_matrix/features.tsv ATAC/filtered_matrix/matrix.mtx ATAC/filtered_matrix/barcodes.tsv

$python $statChAIRATAC fastq/$outprefix.readsnum DNAalign/STAR.barcodereadsnum trim/$outprefix.trim.readsnum linker/cutadapt.linker.readsnum DNAalign/bwa.withlinker.readsnum DNAalign/bwa.withoutlinker.readsnum DNAalign/bwa.withlinker.addBC.barcodeuniquereadsnum DNAalign/bwa.withoutlinker.addBC.barcodeuniquereadsnum DNAalign/bwa.withlinker.addBC.barcodeuniquermPCRreadsnum DNAalign/bwa.withoutlinker.addBC.barcodeuniquermPCRreadsnum fragment/${outprefix}_fragments.unsorted.tsv.log peak/${outprefix}.peak ATAC/filtered_matrix/matrix.mtx.gz custom > ATAC/${outprefix}.stat

wait
$bamCoverage -b DNAalign/$outprefix.scATAC.sorted.bam --minMappingQuality $mapq -o ATAC/${outprefix}.ATAC.bw -bs 5 -p ${NTHREAD}
$bamCoverage -b DNAalign/$outprefix.scATAC.sorted.bam --minMappingQuality $mapq -o ATAC/${outprefix}.ATAC.RPKM.bw -bs 5 -p ${NTHREAD} --normalizeUsing RPKM
$computeMatrix scale-regions --quiet -S ATAC/${outprefix}.ATAC.RPKM.bw -R ${genebed} -a 3000 -b 3000 -m 3000 --skipZeros --missingDataAsZero --outFileName ATAC/${outprefix}.ATAC.RPKM.genebody.gz -p ${NTHREAD}
$computeMatrix reference-point --quiet --referencePoint TSS -S ATAC/${outprefix}.ATAC.RPKM.bw -R ${genebed} -a 3000 -b 3000 --skipZeros --missingDataAsZero --outFileName ATAC/${outprefix}.ATAC.RPKM.geneTSS.gz -p ${NTHREAD}
$plotHeatmap -m ATAC/${outprefix}.ATAC.RPKM.genebody.gz -out ATAC/${outprefix}.ATAC.RPKM.genebody.Heatmap.pdf &
$plotProfile -m ATAC/${outprefix}.ATAC.RPKM.genebody.gz -out ATAC/${outprefix}.ATAC.RPKM.genebody.Profile.pdf &
$plotHeatmap -m ATAC/${outprefix}.ATAC.RPKM.geneTSS.gz -out ATAC/${outprefix}.ATAC.RPKM.geneTSS.Heatmap.pdf &
$plotProfile -m ATAC/${outprefix}.ATAC.RPKM.geneTSS.gz -out ATAC/${outprefix}.ATAC.RPKM.geneTSS.Profile.pdf
wait

dur=\`echo "\$(date +%s.%N) - \$begin" | bc\`
printf "\$(date) generate ATAC signal fragment and peak intensity matrix finish: %.2f seconds\n" \$dur >> jobdir/$jobname.log

echo \$(date) finish >> jobdir/$jobname.log
dur=\`echo "\$(date +%s.%N) - \$start" | bc\`
printf "Execution time for running the pipline : %.2f seconds\n" \$dur >> jobdir/$jobname.log
ENDHERE

    jobnamefile=${outdir}/jobdir/${outprefix}.scChIATAC.sh
    cat ${outdir}/${outprefix}.scRNA.header ${outdir}/${outprefix}.scChIATAC.taskrun > $jobnamefile
    rm ${outdir}/${outprefix}.scRNA.header ${outdir}/${outprefix}.scChIATAC.taskrun
    runid=$(submitjobOptionalShell $jobnamefile)
    echo "mapChIATACcustom $runid $jobnamefile"
;;

*)
    echo "Invalid parameter software: $software"
    ;;
esac
