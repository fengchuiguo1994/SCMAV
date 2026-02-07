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
    jobname=${outprefix}.scHi-C
    var1=$(echo "$DNA_b" | awk '{print $1}')
    var2=$(echo "$DNA_b" | awk '{print $2}')
    var3=$(echo "$DNA_b" | awk '{print $3}')
    cat << ENDHERE > ${outdir}/${outprefix}.scHi-C.taskrun
    
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
for i in \$(seq \$nfile); do 
    /usr/bin/time -f 'real %E\nuser %U\nsys  %S\nMmem %M KB\ncpuP %P\nfsi %I\nfso %O' -o $outdir/trim/$outprefix.\${i}.trim.time.log java -jar ${trimmomaticpath}/trimmomatic-0.39.jar PE -threads $NTHREAD -phred33 $outdir/fastq/$outprefix.R1.\${i}.fastq $outdir/fastq/$outprefix.R2.\${i}.fastq $outdir/trim/$outprefix.R1.\${i}.fastq $outdir/trim/$outprefix.R1.\${i}.unpair.fastq $outdir/trim/$outprefix.R2.\${i}.fastq $outdir/trim/$outprefix.R2.\${i}.unpair.fastq LEADING:5 TRAILING:5 ILLUMINACLIP:${trimmomaticpath}/adapters/NexteraPE-PE.fa:2:30:7:8:true MINLEN:30 > $outdir/trim/$outprefix.\${i}.trim.log 2> $outdir/trim/$outprefix.\${i}.trim.log.err 
    grep Surviving $outdir/trim/$outprefix.\${i}.trim.log.err | perl -lane 'if(/Both Surviving: (\d+) /){\$a=\$1;} if(/Forward Only Surviving: (\d+) /){\$b=\$1;} if(/Reverse Only Surviving: (\d+) /){\$c=\$1;} print "Trim_Reads\t".(\$a+\$b+\$c)' > $outdir/trim/$outprefix.\${i}.trim.readsnum

    /usr/bin/time -f 'real %E\nuser %U\nsys  %S\nMmem %M KB\ncpuP %P\nfsi %I\nfso %O' -o DNAalign/STAR\${i}.time.log $STARsolo --genomeDir $barcoderef --outFileNamePrefix DNAalign/STAR\${i} --readFilesIn $outdir/fastq/$outprefix.R1.\${i}.tmp.fastq $outdir/fastq/$outprefix.barcode.\${i}.fastq --soloCBstart 1 --soloCBlen $var2 --soloUMIstart 1 --soloUMIlen 1 --soloType CB_samTagOut --soloCBwhitelist $DNAbarcodefile --runThreadN $NTHREAD --outSAMattributes CB --outSAMtype BAM Unsorted --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloCBmatchWLtype 1MM --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --outSAMunmapped Within --soloStrand Forward > DNAalign/STAR\${i}.log 2> DNAalign/STAR\${i}.log.err
    grep yes DNAalign/STAR\${i}Solo.out/Barcodes.stats | awk '{sum+=\$2}END{print "barcodeReads\t"sum}' > DNAalign/STAR\${i}.barcodereadsnum
    
    if [ ${blacklist} != "null" ]; then
        $bwa mem -5SP -t $NTHREAD $bwaref $outdir/trim/$outprefix.R1.\${i}.fastq $outdir/trim/$outprefix.R2.\${i}.fastq | $samtools view -@ $NTHREAD -Sb -F 2048 - | $python $BWAmemFilter 1 $mapq DNAalign/bwa\${i}.uniquereadsnum | $bedtools pairtobed -abam - -b ${blacklist} -type neither | $samtools sort -@ $NTHREAD -o DNAalign/bwa\${i}.bam -
    else
        $bwa mem -5SP -t $NTHREAD $bwaref $outdir/trim/$outprefix.R1.\${i}.fastq $outdir/trim/$outprefix.R2.\${i}.fastq | $samtools view -@ $NTHREAD -Sb -F 2048 - | $python $BWAmemFilter 1 $mapq DNAalign/bwa\${i}.uniquereadsnum | $samtools sort -@ $NTHREAD -o DNAalign/bwa\${i}.bam -
    fi
    
    if [ ${barcodetranslate} == none ];then
        $python $giveBC -sb DNAalign/STAR\${i}Aligned.out.bam -bb DNAalign/bwa\${i}.bam -o DNAalign/bwa\${i}.addBC.bam -s DNAalign/bwa\${i}.addBC.barcodeuniquereadsnum
    else
        $python $giveBC -sb DNAalign/STAR\${i}Aligned.out.bam -bb DNAalign/bwa\${i}.bam -o DNAalign/bwa\${i}.addBC.bam -s DNAalign/bwa\${i}.addBC.barcodeuniquereadsnum -t $barcodetranslate
    fi
    
    rm $outdir/fastq/$outprefix.R1.\${i}.fastq $outdir/fastq/$outprefix.R2.\${i}.fastq
    rm $outdir/fastq/$outprefix.R1.\${i}.tmp.fastq $outdir/fastq/$outprefix.barcode.\${i}.fastq $outdir/DNAalign/STAR\${i}SJ.out.tab
    rm $outdir/DNAalign/bwa\${i}.bam 
    rm -rf $outdir/DNAalign/STAR\${i}Solo.out
    rm $outdir/DNAalign/STAR\${i}Log.progress.out $outdir/DNAalign/STAR\${i}Log.out $outdir/DNAalign/STAR\${i}Log.final.out
    rm $outdir/trim/$outprefix.R1.\${i}.fastq $outdir/trim/$outprefix.R1.\${i}.unpair.fastq $outdir/trim/$outprefix.R2.\${i}.fastq $outdir/trim/$outprefix.R2.\${i}.unpair.fastq

done
wait
dur=\`echo "\$(date +%s.%N) - \$begin" | bc\`
printf "\$(date) trim, linker calling and align finish: %.2f seconds\n" \$dur >> jobdir/$jobname.log

begin=\`date '+%s.%N'\`
echo \$(date) merge align result and remove PCR >> jobdir/$jobname.log

$samtools merge -@ $NTHREAD -f -o DNAalign/$outprefix.addBC.bam \$(for i in \$(seq \$nfile); do echo -e "$outdir/DNAalign/bwa\${i}.addBC.bam \c"; done)

cat DNAalign/STAR*.barcodereadsnum | awk '{sum+=\$2}END{print "Barcode_Reads\t"sum}' > DNAalign/STAR.barcodereadsnum
cat trim/$outprefix.*.trim.readsnum | awk '{sum+=\$2}END{print "Trim_Reads\t"sum}' > trim/$outprefix.trim.readsnum

cat DNAalign/bwa*.uniquereadsnum | awk '{sum+=\$2}END{print "unique_mapping_reads\t"sum}' > DNAalign/bwa.uniquereadsnum
cat DNAalign/bwa*.addBC.barcodeuniquereadsnum | awk '{sum+=\$2}END{print "barcode_unique_mapping_reads\t"sum}' > DNAalign/bwa.addBC.barcodeuniquereadsnum

/usr/bin/time -f 'real %E\nuser %U\nsys  %S\nMmem %M KB\ncpuP %P\nfsi %I\nfso %O' -o DNAalign/$jobname.MarkDuplicates.time.log $gatk MarkDuplicates --BARCODE_TAG=ZB --REMOVE_DUPLICATES=true --ADD_PG_TAG_TO_READS=false -M DNAalign/$outprefix.addBC.dedup.MarkDuplicates.txt -O DNAalign/$outprefix.addBC.dedup.bam -I DNAalign/$outprefix.addBC.bam --ASSUME_SORT_ORDER=coordinate > DNAalign/$outprefix.addBC.dedup.log 2>&1

awk '{if(/Unknown Library/){print "barcode_unique_mapping_reads\t"\$4}}' DNAalign/$outprefix.addBC.dedup.MarkDuplicates.txt > DNAalign/bwa.addBC.barcodeuniquereadsnum
awk '{if(/Unknown Library/){print "barcode_unique_mapping_rmPCR_reads\t"\$8}}' DNAalign/$outprefix.addBC.dedup.MarkDuplicates.txt > DNAalign/bwa.addBC.barcodeuniquermPCRreadsnum

$samtools sort -n -@ $NTHREAD -m 1000M -o DNAalign/$outprefix.bam DNAalign/$outprefix.addBC.dedup.bam
$python $classPET DNAalign/$outprefix.bam $contactMinDis DNAalign/$outprefix.FRlt.bam DNAalign/$outprefix.TANDEM.bam DNAalign/$outprefix.scPET.bam DNAalign/$outprefix.classPET.log > DNAalign/$outprefix.classPET.num.log
$python ${classPETplot} DNAalign/$outprefix.classPET.num.log &

for i in \$(seq \$nfile); do rm $outdir/DNAalign/bwa\${i}.addBC.bam $outdir/DNAalign/STAR\${i}Aligned.out.bam $outdir/DNAalign/STAR\${i}.barcodereadsnum $outdir/trim/$outprefix.\${i}.trim.readsnum; done
wait

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
$python $statHiC fastq/$outprefix.readsnum DNAalign/STAR.barcodereadsnum DNAalign/bwa.uniquereadsnum DNAalign/bwa.addBC.barcodeuniquereadsnum DNAalign/bwa.addBC.barcodeuniquermPCRreadsnum PET/${outprefix}.bam2bedpe.log PET/filtered_matrix/matrix.mtx.gz custom > PET/${outprefix}.stat

wait
dur=\`echo "\$(date +%s.%N) - \$begin" | bc\`
printf "\$(date) generate contact finish: %.2f seconds\n" \$dur >> jobdir/$jobname.log

echo \$(date) finish >> jobdir/$jobname.log
dur=\`echo "\$(date +%s.%N) - \$start" | bc\`
printf "Execution time for running the pipline : %.2f seconds\n" \$dur >> jobdir/$jobname.log
ENDHERE

    jobnamefile=${outdir}/jobdir/${outprefix}.scHi-C.sh
    cat ${outdir}/${outprefix}.scRNA.header ${outdir}/${outprefix}.scHi-C.taskrun > $jobnamefile
    rm ${outdir}/${outprefix}.scRNA.header ${outdir}/${outprefix}.scHi-C.taskrun
    runid=$(submitjobOptionalShell $jobnamefile)
    echo "mapHi-Ccustom $runid $jobnamefile"
;;

*)
    echo "Invalid parameter software: $software"
    ;;
esac
