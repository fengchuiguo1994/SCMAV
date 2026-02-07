Rbc=$1
Rrd=$2
outdir=$3
outprefix=$4
config=$5
source $config

SCRIPT_DIR=$(dirname "$0")
SCRIPT_DIR=$(realpath $SCRIPT_DIR)
source $SCRIPT_DIR/function.sh

RbcBName=$(basename $(realpath $Rbc))
ln -s -f $Rbc $outdir/$RbcBName
RrdBName=$(basename $(realpath $Rrd))
ln -s -f $Rrd $outdir/$RrdBName

RNAsoftware=$(echo "$software" | awk -F"+" '{print $1}')
rnasample=$(echo ${RrdBName} | awk -F'_' '{for(i=1;i<=NF-4;i++) printf "%s%s", $i, (i==NF-4?"":"_")}')
case $RNAsoftware in
10X)
    jobname=${outprefix}.scRNA
    cat << ENDHERE > ${outdir}/${outprefix}.scRNA.taskrun

set -e
cd $outdir
start=\`date '+%s.%N'\`

echo \$(date) cellranger align >> jobdir/$jobname.log
/usr/bin/time -f 'real %E\nuser %U\nsys  %S\nMmem %M KB\ncpuP %P\nfsi %I\nfso %O' -o jobdir/$jobname.cellranger.time.log $cellranger count --id RNAalign --transcriptome $GEXref --fastqs=$outdir --include-introns=true --localcores=$NTHREAD --localmem=$GB --sample $rnasample > jobdir/$jobname.cellranger.log 2> jobdir/$jobname.cellranger.log.err
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

    jobnamefile=${outdir}/jobdir/${outprefix}.scRNA.slurm
    cat << ENDHERE > ${outdir}/${outprefix}.scRNA.header
#!/usr/bin/bash
#SBATCH --job-name $jobname
#SBATCH -N 1
#SBATCH -n ${NTHREAD} 
#SBATCH --mem ${GB}G
#SBATCH --partition=${QUEUE}
#SBATCH --output=$jobnamefile.out
#SBATCH --error=$jobnamefile.err
ENDHERE
    cat ${outdir}/${outprefix}.scRNA.header ${outdir}/${outprefix}.scRNA.taskrun > $jobnamefile
    rm ${outdir}/${outprefix}.scRNA.header ${outdir}/${outprefix}.scRNA.taskrun
    runid=$(submitjobOptionalSLURM $jobnamefile)
    echo "mapGEX10x $runid $jobnamefile"
;;

STARsolo)
    var1=$(echo "$RNA_b_u" | awk '{print $1}') # barcode
    var2=$(echo "$RNA_b_u" | awk '{print $2}')
    var3=$(echo "$RNA_b_u" | awk '{print $3}') # umi
    var4=$(echo "$RNA_b_u" | awk '{print $4}')
    var5=$(echo "$RNA_b_u" | awk '{print $5}')
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

    jobnamefile=${outdir}/jobdir/${outprefix}.scRNA.slurm
    cat << ENDHERE > ${outdir}/${outprefix}.scRNA.header
#!/usr/bin/bash
#SBATCH --job-name $jobname
#SBATCH -N 1
#SBATCH -n ${NTHREAD} 
#SBATCH --mem ${GB}G
#SBATCH --partition=${QUEUE}
#SBATCH --output=$jobnamefile.out
#SBATCH --error=$jobnamefile.err
ENDHERE
    cat ${outdir}/${outprefix}.scRNA.header ${outdir}/${outprefix}.scRNA.taskrun > $jobnamefile
    rm ${outdir}/${outprefix}.scRNA.header ${outdir}/${outprefix}.scRNA.taskrun
    runid=$(submitjobOptionalSLURM $jobnamefile)
    echo "mapGEXSTARsolo $runid $jobnamefile"
;;
*)
    echo "Invalid parameter software: $software"
    ;;
esac
