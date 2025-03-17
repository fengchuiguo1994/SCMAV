#!/usr/bin/bash

set -e
#################################### must
# addLinkerSeq='ACGTGATATTTCACGACTCT' # add linker sequence for cpu stag '-A ACGTGATATTTCACGACTCT' # Important: if using default linker just leave it empty .ie. '' !!!
addLinkerSeq='ATAGGCGCTGGTTGCGTTGT' # linker76

genome=hg38 # uniq ID
ARCref=/data/home/ruanlab/huangxingyu/Tools/cellranger-arc-2.0.2/genome/refdata-cellranger-arc-GRCh38-2020-A-2.0.0
GEXref=/data/home/ruanlab/huangxingyu/Tools/cellranger-7.1.0/genome/refdata-gex-GRCh38-2020-A
ATACref=/data/home/ruanlab/huangxingyu/Tools/cellranger-atac-2.1.0/genome/refdata-cellranger-arc-GRCh38-2020-A-2.0.0
bwaref=/data/home/ruanlab/huangxingyu/Haoxi20230315/20230205_pipeline/refs/hg38/bwa/hg38.fa
blacklist=/data/home/ruanlab/huangxingyu/Haoxi20230315/20230205_pipeline/pipelineTools/hg38.blacklist.bed # or null
genebed=/data/home/ruanlab/huangxingyu/Tools/cellranger-7.1.0/genome/refdata-gex-GRCh38-2020-A/genes/hg38.info

siffile="/data/home/ruanlab/huangxingyu/scChAIR/ChAIR.sif"
codepath="/data/home/ruanlab/huangxingyu/scChAIR/codes"

####################### for slurm system
NTHREADmap=48  # number of cpu for reads align (for cellranger/cellranger-arc/cellranger-atac/bwa)
NTHREAD=20     # number of cpu for other task
highGB=300 # high memory task, advice: 400
GB=200 # regular jobs memory
smallGB=100 # small memory task,
highQUEUE=cpu
QUEUE=cpu
smallQUEUE=cpu

####################### for macs2
glen='hs' #MACS2: hs=human, mm=mouse, dm=fruitfly
shiftsize=-75
macsq=0.001
peakext=150
alldata=T # or F. T: use all bam, F: bam from self-lagation reads with none linker # no used

####################### for cpu (ChIATAC)
minMapScore=30 # CPU mappping output filter
mapq=20 # mapping quality filter by samtools
selfbp=1 # min PET distance in single cell
bulkSelfbp=8000 # CPU PET interactions min distance
minPETsignif=3 # for ChiaSig
#define names
pairlabel='singlelinker.paired'  #PET or 2tag file names
singlabel='singlelinker.single'  #1tag file names
nonelabel='none'                 #0tag file names
dumpfolder=dump_intermediate     #where intermediate files to be dumped

####################### tools
cellranger="singularity run ${siffile} cellranger"
cellrangerarc="singularity run ${siffile} cellranger-arc"
cellrangeratac="singularity run ${siffile} cellranger-atac"
Rscript="singularity exec ${siffile} Rscript"
perl="singularity exec ${siffile} perl"
python="singularity run ${siffile} python"
pigz="singularity run ${siffile} pigz -f"
gatk="singularity run ${siffile} gatk"
bedtools="singularity run ${siffile} bedtools"
samtools="singularity run ${siffile} samtools"
juicer_tools="singularity run ${siffile} juicer_tools"
cpu="singularity run ${siffile} cpu"
bedClip="singularity run ${siffile} bedClip"
bedGraphToBigWig="singularity run ${siffile} bedGraphToBigWig"
bgzip="singularity run ${siffile} bgzip"
sortBed="singularity run ${siffile} sortBed"
ChiaSig="singularity run ${siffile} ChiaSig"
tabix="singularity run ${siffile} tabix"
macs2="singularity run ${siffile} macs2"
bamCoverage="singularity run ${siffile} bamCoverage"
computeMatrix="singularity run ${siffile} computeMatrix"
plotHeatmap="singularity run ${siffile} plotHeatmap"
plotProfile="singularity run ${siffile} plotProfile"


statpy_noTag=${codepath}/statBarcode_givenCBlist.v2.py
addCBpy=${codepath}/addCBZB_givenCBlist.v2.py
injectCBreadID=${codepath}/useCBtag_injectReadID.pl
pickUpy=${codepath}/pickingUxxU.xyh.py
statpy_bamTag=${codepath}/statCB_PET.py
buildCBpy=${codepath}/generateCBlist_for_cpuRead.v2.py
mergeTable=${codepath}/mergeCBstats.v3.py
bamtofrag=${codepath}/bam2fragment.py # https://github.com/wbaopaul/scATAC-pro # xyh
sortfrag=${codepath}/sort_frags.R
sumsh=${codepath}/extract_summary_cpu10x_1genome.v3.sh
filterinsertsizeless1000=${codepath}/classPET.py
CountBarcodePET=${codepath}/CountBarcodePET.py # xyh
getfrom10xG=${codepath}/getfrom10xG.R
getfrom10xA=${codepath}/getfrom10xA.py
getGEXbarcode=${codepath}/getGEXbarcode.R
barcodePlot=${codepath}/barcodePlot.R
plotPETnumber=${codepath}/plotPETnumber.R
plotFragmentlen=${codepath}/plotFragmentlen.R
plotStat=${codepath}/plotStat.R
gex2fragment=${codepath}/gex2fragment.py
bedpe2mat=${codepath}/bedpe2mat.forBandNorm.v2.py
petpca=${codepath}/petpca.pro.v2.R
catBam=${codepath}/catBam.py

maindir=$PWD

#------------------------------------
# Auxilary functions
#------------------------------------
die () {
    echo >&2 "$@"
    exit 1
}
[ $# -ge 1 ] || die "Please give 1 arg: fqlist file!"

function odir {
    local outputDir=$1
    if [ ! -d $outputDir ]; then
        echo "Create $outputDir"
        mkdir $outputDir
    fi
}

function submitjob { # one depend job submit
    local slurmfile=$1 # jobs
    local jid=$2 # depend
    [[ ! -z "$jid" ]] && jid1=$(sbatch --dependency=afterok:$jid  $slurmfile) || jid1=$(sbatch $slurmfile) # -z is empty or not, empty return true
    echo $jid1 | awk '{print $NF}'
}

function submitjobOptional { # more then one depend job submit
    local slurmfile=$1 # jobs
    local jid1=$2
    local jid2optional=$3

    if [ -z "$jid1" ]; then
        jidnow=$(sbatch $slurmfile)
    elif [ -z "$jid2optional" ]; then
        jidnow=$(sbatch --dependency=afterok:$jid1  $slurmfile)
    else
        jidnow=$(sbatch --dependency=afterok:$jid1:$jid2optional  $slurmfile)
    fi
    echo $jidnow | awk '{print $NF}'
}

function submitjobOptional8 { # more
    local slurmfile=$1
    local jid1optional=$2
    local jid2optional=$3
    local jid3optional=$4
    local jid4optional=$5
    local jid5optional=$6
    local jid6optional=$7
    local jid7optional=$8
    local jid8optional=$9

    if [ -z "$jid1optional" ]; then
        jidnow=$(sbatch $slurmfile)
    elif [ -z "$jid2optional" ]; then
        jidnow=$(sbatch --dependency=afterany:$jid1optional  $slurmfile)
    elif [ -z "$jid3optional" ]; then
        jidnow=$(sbatch --dependency=afterany:$jid1optional:$jid2optional  $slurmfile)
    elif [ -z "$jid4optional" ]; then
        jidnow=$(sbatch --dependency=afterany:$jid1optional:$jid2optional:$jid3optional  $slurmfile)
    elif [ -z "$jid5optional" ]; then
        jidnow=$(sbatch --dependency=afterany:$jid1optional:$jid2optional:$jid3optional:$jid4optional  $slurmfile)
    elif [ -z "$jid6optional" ]; then
        jidnow=$(sbatch --dependency=afterany:$jid1optional:$jid2optional:$jid3optional:$jid4optional:$jid5optional  $slurmfile)
    elif [ -z "$jid7optional" ]; then
        jidnow=$(sbatch --dependency=afterany:$jid1optional:$jid2optional:$jid3optional:$jid4optional:$jid5optional:$jid6optional  $slurmfile)
    elif [ -z "$jid8optional" ]; then
        jidnow=$(sbatch --dependency=afterany:$jid1optional:$jid2optional:$jid3optional:$jid4optional:$jid5optional:$jid6optional:$jid7optional  $slurmfile)
    else
        jidnow=$(sbatch --dependency=afterany:$jid1optional:$jid2optional:$jid3optional:$jid4optional:$jid5optional:$jid6optional:$jid7optional:$jid8optional  $slurmfile)
    fi
    echo $jidnow | awk '{print $NF}'
}

############################################### just for get barcode
function mapARC10x {
    # run 10x pipeline for multiome, just for get the barcode list
    local fqDNApref=$1
    local fqRNApref=$2
    dnaRUN=$(basename $fqDNApref)
    rnaRUN=$(basename $fqRNApref)
    dnaDir=$(dirname $fqDNApref)
    rnaDir=$(dirname $fqRNApref)
    jobname=${dnaRUN}.arc
    pbs=${jobname}.slurm
    wdir=$maindir/$dnaRUN
    libcsv=$wdir/$jobname.csv
    echo 'fastqs,sample,library_type' > $libcsv
    echo "${dnaDir},${dnaRUN},Chromatin Accessibility" >> $libcsv
    echo "${rnaDir},${rnaRUN},Gene Expression" >> $libcsv

    cat << ENDHERE > $pbs
#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -n ${NTHREADmap}
#SBATCH --mem ${GB}G 
#SBATCH --partition=${QUEUE}
#SBATCH --job-name $jobname
#SBATCH --output="$jobname.slurm.o"
#SBATCH --error="$jobname.slurm.e"

set -e 
cd $wdir
export SINGULARITY_BIND=$bwarefabs,$blacklistabs,$sifpathabs,$ARCrefabs,$GEXrefabs,$ATACrefabs,$codepathabs,$maindir,$wdir,$genebedabs
date

/usr/bin/time -f 'real %E\nuser %U\nsys  %S\nmem  %M\n' $cellrangerarc count --reference $ARCref --id 10xARC --libraries=$libcsv --localcores=$NTHREADmap --jobmode=local && echo "success" > slurm/$jobname.log

date
ENDHERE

    jidarc=$(submitjob $pbs)
    echo "mapARC10x $jidarc $pbs" >&2
    buildCBlist $fqDNApref $jidarc
    
}


function buildCBlist {
    local fqDNApref=$1
    dnaRUN=$(basename $fqDNApref)
    local prevJobid=$2
    jobname=${dnaRUN}.buildCBpy
    pbs=${jobname}.slurm
    wdir=$maindir/$dnaRUN
    
    cat << ENDHERE > $pbs
#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem ${highGB}G
#SBATCH --partition=${highQUEUE}
#SBATCH --job-name $jobname
#SBATCH --output="$jobname.slurm.o"
#SBATCH --error="$jobname.slurm.e"

set -e 
cd $wdir
export SINGULARITY_BIND=$bwarefabs,$blacklistabs,$sifpathabs,$ARCrefabs,$GEXrefabs,$ATACrefabs,$codepathabs,$maindir,$wdir,$genebedabs
date 

/usr/bin/time -f 'real %E\nuser %U\nsys  %S\nmem  %M\n' $python $buildCBpy $maindir/$dnaRUN/10xARC/outs/atac_possorted_bam.bam $maindir/$dnaRUN/bulkChiatac/${dnaRUN}_R1.fastq.gz $dnaRUN && echo "success" > slurm/$jobname.log

date 
ENDHERE

    jidpkl=$(submitjob $pbs $prevJobid)
    echo "buildCBlist $jidpkl $pbs wait $prevJobid" >&2
}


############################################### deal scRNA data
function mapGEX10x {
    #run 10x pipeline for scRNA
    local fqDNApref=$1
    local fqRNApref=$2
    dnaRUN=$(basename $fqDNApref)
    rnaRUN=$(basename $fqRNApref)
    dnaDir=$(dirname $fqDNApref)
    rnaDir=$(dirname $fqRNApref)
    jobname=${rnaRUN}.gex
    pbs=${jobname}.slurm
    wdir=$maindir/$dnaRUN

    cat << ENDHERE > $pbs
#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -n ${NTHREADmap} 
#SBATCH --mem ${GB}G
#SBATCH --partition=${QUEUE}
#SBATCH --job-name $jobname
#SBATCH --output="$jobname.slurm.o"
#SBATCH --error="$jobname.slurm.e"

set -e 
cd $wdir
export SINGULARITY_BIND=$bwarefabs,$blacklistabs,$sifpathabs,$ARCrefabs,$GEXrefabs,$ATACrefabs,$codepathabs,$maindir,$wdir,$genebedabs
date

/usr/bin/time -f 'real %E\nuser %U\nsys  %S\nmem  %M\n' $cellranger count --transcriptome $GEXref --id 10xGEX --sample ${rnaRUN} --fastqs=$rnaDir --localcores=$NTHREADmap --localmem=$GB --include-introns --chemistry=ARC-v1 && echo "success" > slurm/$jobname.log # # v6

date
ENDHERE

    runid=$(submitjob $pbs)
    echo "mapGEX10x $runid $pbs" >&2
}


############################################### treat as ChIATAC
function cpuLinker {
    # linker calling and split data
    local fqDNApref=$1
    RUN=$(basename $fqDNApref)
    fastq1=${RUN}_R1.fastq.gz
    fastq2=${RUN}_R3.fastq.gz
    nf=$(/bin/ls ${fqDNApref}_*L*_R1*.fastq.gz |wc -l)
    r1=$(/bin/ls ${fqDNApref}_*L*_R1*.fastq.gz)
    r3=$(/bin/ls ${fqDNApref}_*L*_R3*.fastq.gz)
    cpudir=$maindir/$RUN/bulkChiatac
    jobname=${RUN}.linker
    pbs=${jobname}.slurm

    cat << ENDHERE > $pbs
#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -n $NTHREAD
#SBATCH --mem ${smallGB}G
#SBATCH --partition=${smallQUEUE}
#SBATCH --job-name $jobname
#SBATCH --output="$jobname.slurm.o"
#SBATCH --error="$jobname.slurm.e"

set -e 
cd $cpudir
export SINGULARITY_BIND=$bwarefabs,$blacklistabs,$sifpathabs,$ARCrefabs,$GEXrefabs,$ATACrefabs,$codepathabs,$maindir,$wdir,$genebedabs
date

if [ "$nf" -gt  1 ]; then
   ls ${fqDNApref}_*L*_R1*.fastq.gz 
   cat ${fqDNApref}_*L*_R1*.fastq.gz > $fastq1
   ls ${fqDNApref}_*L*_R3*.fastq.gz
   cat ${fqDNApref}_*L*_R3*.fastq.gz > $fastq2
else
    ln -s -f $r1 $fastq1
    ln -s -f $r3 $fastq2
fi

/usr/bin/time -f 'real %E\nuser %U\nsys  %S\nmem  %M\n' $cpu stag -W -T 18 -t ${NTHREAD} -A $addLinkerSeq -O ${RUN} $fastq1 $fastq2 && \
$cpu stat -s -p -T 18 -t 1 ${RUN}.cpu 1> ${RUN}.stat && \
$pigz -p ${NTHREAD} ${RUN}.singlelinker.paired.fastq && \
$pigz -p ${NTHREAD} ${RUN}.singlelinker.single.fastq && \
$pigz -p ${NTHREAD} ${RUN}.none.fastq && \
rm ${RUN}.conflict.fastq && \
rm ${RUN}.tied.fastq && echo "success" > ../slurm/$jobname.log

date
ENDHERE

    jidLinker=$(submitjob $pbs)
    echo "cpuLinker $jidLinker $pbs" >&2
}


function cpuMapping {
    # mapping
    local fqDNApref=$1
    local LABEL=$2 #either pairlabel, singlabel, or nonelabel
    local prevJobid=$3 #sbatch dependency
    pkl=$maindir/$RUN/${RUN}.CB_forCPU.pkl
    RUN=$(basename $fqDNApref)
    jobname=${RUN}.map.$LABEL
    pbs=${jobname}.slurm
    cpudir=$maindir/$RUN/bulkChiatac
    halfN=$(echo $NTHREADmap/2 | bc )

    cat << ENDHERE > $pbs
#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -n $NTHREADmap
#SBATCH --mem ${GB}G
#SBATCH --partition=${QUEUE}
#SBATCH --job-name $jobname
#SBATCH --output="$jobname.slurm.o"
#SBATCH --error="$jobname.slurm.e"

set -e 
cd $cpudir
export SINGULARITY_BIND=$bwarefabs,$blacklistabs,$sifpathabs,$ARCrefabs,$GEXrefabs,$ATACrefabs,$codepathabs,$maindir,$wdir,$genebedabs
date

/usr/bin/time -f 'real %E\nuser %U\nsys  %S\nmem  %M\n' $cpu memaln -T $minMapScore -t ${NTHREADmap} $bwaref ${RUN}.${LABEL}.fastq.gz | $python $addCBpy $pkl > ${RUN}.${LABEL}.sam && $pigz -p ${NTHREADmap} ${RUN}.${LABEL}.sam
$python $statpy_noTag ${RUN}.${LABEL}.sam.gz | gzip > ${RUN}.${LABEL}.statCB.tsv.gz &
$cpu pair -s $selfbp -S -t 2 -q $mapq ${RUN}.${LABEL}.sam.gz 1> ${RUN}.${LABEL}.stat.xls
$python $statpy_noTag ${RUN}.${LABEL}.UU.bam | gzip > ${RUN}.${LABEL}.UU.statCB.tsv.gz &
$python $statpy_noTag ${RUN}.${LABEL}.UxxU.bam | gzip > ${RUN}.${LABEL}.UxxU.statCB.tsv.gz &

if [ ${LABEL} == ${pairlabel} ];then
    if [ ${blacklist} != "null" ];then
        /usr/bin/time -f 'real %E\nuser %U\nsys  %S\nmem  %M\n' $bedtools intersect -nonamecheck -a ${RUN}.${LABEL}.UxxU.bam -b ${blacklist} -v | $samtools view -@ ${halfN} -F 2048 -F 4 -q ${mapq} -h - | $python $pickUpy ${mapq} | $samtools view -@ ${halfN} -Sb -o ${RUN}.${LABEL}.UxxU.hq.bam - &
        /usr/bin/time -f 'real %E\nuser %U\nsys  %S\nmem  %M\n' $samtools view -@ ${halfN} -F 2048 -Sb ${RUN}.${LABEL}.UU.bam | $bedtools pairtobed -abam - -b ${blacklist} -type neither > ${RUN}.${LABEL}.UU.hq.bam
    else
        ln -s -f ${RUN}.${LABEL}.UU.bam ${RUN}.${LABEL}.UU.hq.bam
        /usr/bin/time -f 'real %E\nuser %U\nsys  %S\nmem  %M\n' $samtools view -@ ${halfN} -F 2048 -F 4 -q ${mapq} -h ${RUN}.${LABEL}.UxxU.bam | $python $pickUpy ${mapq} | $samtools view -@ ${halfN} -Sb -o ${RUN}.${LABEL}.UxxU.hq.bam -
    fi
    wait
    $python $statpy_noTag ${RUN}.${LABEL}.UU.hq.bam | gzip > ${RUN}.${LABEL}.UU.hq.statCB.tsv.gz &
    $python $statpy_noTag ${RUN}.${LABEL}.UxxU.hq.bam | gzip > ${RUN}.${LABEL}.UxxU.hq.statCB.tsv.gz &
    $samtools view ${RUN}.${LABEL}.UU.hq.bam | $python $statpy_bamTag | gzip > ${RUN}.${LABEL}.UU.hq.CB.statPET.tsv.gz &
    $samtools view -@ ${NTHREADmap} -h ${RUN}.${LABEL}.UU.hq.bam | grep -v noCB | $samtools sort -@ $NTHREADmap -o ${RUN}.${LABEL}.UU.hq.CB.csorted.bam -
    $samtools view -@ ${NTHREADmap} -h ${RUN}.${LABEL}.UxxU.hq.bam | grep -v noCB | $samtools sort -@ $NTHREADmap -o ${RUN}.${LABEL}.UxxU.hq.CB.csorted.bam -
    wait
    $gatk MarkDuplicates --BARCODE_TAG=ZB --REMOVE_DUPLICATES=true --ADD_PG_TAG_TO_READS=false -M ${RUN}.${LABEL}.UU.hq.CB.nr.MarkDuplicates.txt -O ${RUN}.${LABEL}.UU.hq.CB.csorted.dedup.bam -I ${RUN}.${LABEL}.UU.hq.CB.csorted.bam && $samtools sort -n -@ $NTHREADmap -o ${RUN}.${LABEL}.UU.hq.CB.nr.bam ${RUN}.${LABEL}.UU.hq.CB.csorted.dedup.bam &
    $gatk MarkDuplicates --BARCODE_TAG=ZB --REMOVE_DUPLICATES=true --ADD_PG_TAG_TO_READS=false -M ${RUN}.${LABEL}.UxxU.hq.CB.nr.MarkDuplicates.txt -O ${RUN}.${LABEL}.UxxU.hq.CB.nr.bam -I ${RUN}.${LABEL}.UxxU.hq.CB.csorted.bam
    wait
    $gatk CollectInsertSizeMetrics -I ${RUN}.${LABEL}.UU.hq.CB.csorted.dedup.bam -O ${RUN}.${LABEL}.UU.hq.CB.csorted.dedup.insert_size_metrics.txt -H ${RUN}.${LABEL}.UU.hq.CB.csorted.dedup.insert_size_histogram.pdf -M 0.05 &
    $python $statpy_noTag ${RUN}.${LABEL}.UU.hq.CB.nr.bam | gzip > ${RUN}.${LABEL}.UU.hq.CB.nr.statCB.tsv.gz &
    $python $statpy_noTag ${RUN}.${LABEL}.UxxU.hq.CB.nr.bam | gzip > ${RUN}.${LABEL}.UxxU.hq.CB.nr.statCB.tsv.gz &
    $samtools view ${RUN}.${LABEL}.UU.hq.CB.nr.bam | $python $statpy_bamTag | gzip > ${RUN}.${LABEL}.UU.hq.CB.nr.statPET.tsv.gz &
    $python ${filterinsertsizeless1000} ${RUN}.${LABEL}.UU.hq.CB.nr.bam ${RUN}.${LABEL}.UU.hq.CB.nr &
    $samtools view -h ${RUN}.${LABEL}.UU.hq.CB.nr.bam | $perl $injectCBreadID | $bedtools bamtobed -bedpe | sort -k7,7 > ${RUN}.${LABEL}.UU.hq.CB.nr.contain.all.bedpe &
    $samtools index ${RUN}.${LABEL}.UU.hq.CB.csorted.dedup.bam
    $bamCoverage -b ${RUN}.${LABEL}.UU.hq.CB.csorted.dedup.bam -o ${RUN}.${LABEL}.UU.hq.CB.nr.bw -p $NTHREADmap -bs 5 && $bamCoverage -b ${RUN}.${LABEL}.UU.hq.CB.csorted.dedup.bam -o ${RUN}.${LABEL}.UU.hq.CB.nr.RPKM.bw -p $NTHREADmap -bs 5 --normalizeUsing RPKM && $bamCoverage -b ${RUN}.${LABEL}.UU.hq.CB.csorted.dedup.bam -o ${RUN}.${LABEL}.UU.hq.CB.nr.bedgraph -p $NTHREADmap -bs 5 -of bedgraph && $pigz -p ${NTHREADmap} ${RUN}.${LABEL}.UU.hq.CB.nr.bedgraph && $bamCoverage -b ${RUN}.${LABEL}.UU.hq.CB.csorted.dedup.bam -o ${RUN}.${LABEL}.UU.hq.CB.nr.RPKM.bedgraph -p $NTHREADmap -bs 5 -of bedgraph --normalizeUsing RPKM && $pigz -p ${NTHREADmap} ${RUN}.${LABEL}.UU.hq.CB.nr.RPKM.bedgraph
    $samtools index ${RUN}.${LABEL}.UxxU.hq.CB.nr.bam 
    $bamCoverage -b ${RUN}.${LABEL}.UxxU.hq.CB.nr.bam -o ${RUN}.${LABEL}.UxxU.hq.CB.nr.bw -p $NTHREADmap -bs 5 && $bamCoverage -b ${RUN}.${LABEL}.UxxU.hq.CB.nr.bam -o ${RUN}.${LABEL}.UxxU.hq.CB.nr.RPKM.bw -p $NTHREADmap -bs 5 --normalizeUsing RPKM && $bamCoverage -b ${RUN}.${LABEL}.UxxU.hq.CB.nr.bam -o ${RUN}.${LABEL}.UxxU.hq.CB.nr.bedgraph -p $NTHREADmap -bs 5 -of bedgraph && $pigz -p ${NTHREADmap} ${RUN}.${LABEL}.UxxU.hq.CB.nr.bedgraph && $bamCoverage -b ${RUN}.${LABEL}.UxxU.hq.CB.nr.bam -o ${RUN}.${LABEL}.UxxU.hq.CB.nr.RPKM.bedgraph -p $NTHREADmap -bs 5 -of bedgraph --normalizeUsing RPKM && $pigz -p ${NTHREADmap} ${RUN}.${LABEL}.UxxU.hq.CB.nr.RPKM.bedgraph
    wait
    $python $statpy_noTag ${RUN}.${LABEL}.UU.hq.CB.nr.pet.bam | gzip > ${RUN}.${LABEL}.UU.hq.CB.nr.pet.statCB.tsv.gz &
    $python $statpy_noTag ${RUN}.${LABEL}.UU.hq.CB.nr.other.bam | gzip > ${RUN}.${LABEL}.UU.hq.CB.nr.other.statCB.tsv.gz &
    $python $statpy_noTag ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.bam | gzip > ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.statCB.tsv.gz &
    $samtools view ${RUN}.${LABEL}.UU.hq.CB.nr.pet.bam | $python $statpy_bamTag | gzip > ${RUN}.${LABEL}.UU.hq.CB.nr.pet.statPET.tsv.gz &
    $samtools view ${RUN}.${LABEL}.UU.hq.CB.nr.other.bam | $python $statpy_bamTag | gzip > ${RUN}.${LABEL}.UU.hq.CB.nr.other.statPET.tsv.gz &
    $samtools view ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.bam | $python $statpy_bamTag | gzip > ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.statPET.tsv.gz &
    $samtools view -h ${RUN}.${LABEL}.UU.hq.CB.nr.pet.bam | $perl $injectCBreadID | $bedtools bamtobed -bedpe | sort -k7,7 > ${RUN}.${LABEL}.UU.hq.CB.nr.bedpe &
    $samtools sort -@ $NTHREADmap -o ${RUN}.${LABEL}.UU.hq.CB.nr.pet.csorted.bam ${RUN}.${LABEL}.UU.hq.CB.nr.pet.bam
    $samtools sort -@ $NTHREADmap -o ${RUN}.${LABEL}.UU.hq.CB.nr.other.csorted.bam ${RUN}.${LABEL}.UU.hq.CB.nr.other.bam
    $samtools sort -@ $NTHREADmap -o ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.csorted.bam ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.bam
    $gatk CollectInsertSizeMetrics -I ${RUN}.${LABEL}.UU.hq.CB.nr.pet.csorted.bam -O ${RUN}.${LABEL}.UU.hq.CB.nr.pet.csorted.insert_size_metrics.txt -H ${RUN}.${LABEL}.UU.hq.CB.nr.pet.csorted.insert_size_histogram.pdf -M 0.05 &
    $gatk CollectInsertSizeMetrics -I ${RUN}.${LABEL}.UU.hq.CB.nr.other.csorted.bam -O ${RUN}.${LABEL}.UU.hq.CB.nr.other.csorted.insert_size_metrics.txt -H ${RUN}.${LABEL}.UU.hq.CB.nr.other.csorted.insert_size_histogram.pdf -M 0.05 &
    $gatk CollectInsertSizeMetrics -I ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.csorted.bam -O ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.csorted.insert_size_metrics.txt -H ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.csorted.insert_size_histogram.pdf -M 0.05 &
    $samtools index ${RUN}.${LABEL}.UU.hq.CB.nr.pet.csorted.bam
    $samtools index ${RUN}.${LABEL}.UU.hq.CB.nr.other.csorted.bam
    $samtools index ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.csorted.bam
    $bamCoverage -b ${RUN}.${LABEL}.UU.hq.CB.nr.pet.csorted.bam -o ${RUN}.${LABEL}.UU.hq.CB.nr.pet.csorted.bw -p $NTHREADmap -bs 5 && $bamCoverage -b ${RUN}.${LABEL}.UU.hq.CB.nr.pet.csorted.bam -o ${RUN}.${LABEL}.UU.hq.CB.nr.pet.csorted.bedgraph -p $NTHREADmap -bs 5 -of bedgraph && $pigz -p ${NTHREADmap} ${RUN}.${LABEL}.UU.hq.CB.nr.pet.csorted.bedgraph && $bamCoverage -b ${RUN}.${LABEL}.UU.hq.CB.nr.pet.csorted.bam -o ${RUN}.${LABEL}.UU.hq.CB.nr.pet.csorted.RPKM.bw -p $NTHREADmap -bs 5 --normalizeUsing RPKM && $bamCoverage -b ${RUN}.${LABEL}.UU.hq.CB.nr.pet.csorted.bam -o ${RUN}.${LABEL}.UU.hq.CB.nr.pet.csorted.RPKM.bedgraph -p $NTHREADmap -bs 5 -of bedgraph --normalizeUsing RPKM && $pigz -p ${NTHREADmap} ${RUN}.${LABEL}.UU.hq.CB.nr.pet.csorted.RPKM.bedgraph
    $bamCoverage -b ${RUN}.${LABEL}.UU.hq.CB.nr.other.csorted.bam -o ${RUN}.${LABEL}.UU.hq.CB.nr.other.csorted.bw -p $NTHREADmap -bs 5 && $bamCoverage -b ${RUN}.${LABEL}.UU.hq.CB.nr.other.csorted.bam -o ${RUN}.${LABEL}.UU.hq.CB.nr.other.csorted.bedgraph -p $NTHREADmap -bs 5 -of bedgraph && $pigz -p ${NTHREADmap} ${RUN}.${LABEL}.UU.hq.CB.nr.other.csorted.bedgraph && $bamCoverage -b ${RUN}.${LABEL}.UU.hq.CB.nr.other.csorted.bam -o ${RUN}.${LABEL}.UU.hq.CB.nr.other.csorted.RPKM.bw -p $NTHREADmap -bs 5 --normalizeUsing RPKM && $bamCoverage -b ${RUN}.${LABEL}.UU.hq.CB.nr.other.csorted.bam -o ${RUN}.${LABEL}.UU.hq.CB.nr.other.csorted.RPKM.bedgraph -p $NTHREADmap -bs 5 -of bedgraph --normalizeUsing RPKM && $pigz -p ${NTHREADmap} ${RUN}.${LABEL}.UU.hq.CB.nr.other.csorted.RPKM.bedgraph
    $bamCoverage -b ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.csorted.bam -o ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.csorted.bw -p $NTHREADmap -bs 5 && $bamCoverage -b ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.csorted.bam -o ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.csorted.bedgraph -p $NTHREADmap -bs 5 -of bedgraph && $pigz -p ${NTHREADmap} ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.csorted.bedgraph && $bamCoverage -b ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.csorted.bam -o ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.csorted.RPKM.bw -p $NTHREADmap -bs 5 --normalizeUsing RPKM && $bamCoverage -b ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.csorted.bam -o ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.csorted.RPKM.bedgraph -p $NTHREADmap -bs 5 -of bedgraph --normalizeUsing RPKM && $pigz -p ${NTHREADmap} ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.csorted.RPKM.bedgraph
    wait

elif [ ${LABEL} == ${nonelabel} ];then
    if [ ${blacklist} != "null" ];then
        /usr/bin/time -f 'real %E\nuser %U\nsys  %S\nmem  %M\n' $bedtools intersect -nonamecheck -a ${RUN}.${LABEL}.UxxU.bam -b ${blacklist} -v | $samtools view -@ ${halfN} -F 2048 -F 4 -q ${mapq} -h - | $python $pickUpy ${mapq} | $samtools view -@ ${halfN} -Sb -o ${RUN}.${LABEL}.UxxU.hq.bam - &
        /usr/bin/time -f 'real %E\nuser %U\nsys  %S\nmem  %M\n' $samtools view -@ ${halfN} -F 2048 -Sb ${RUN}.${LABEL}.UU.bam | $bedtools pairtobed -abam - -b ${blacklist} -type neither > ${RUN}.${LABEL}.UU.hq.bam
    else
        ln -s -f ${RUN}.${LABEL}.UU.bam ${RUN}.${LABEL}.UU.hq.bam
        /usr/bin/time -f 'real %E\nuser %U\nsys  %S\nmem  %M\n' $samtools view -@ ${halfN} -F 2048 -F 4 -q ${mapq} -h ${RUN}.${LABEL}.UxxU.bam | $python $pickUpy ${mapq} | $samtools view -@ ${halfN} -Sb -o ${RUN}.${LABEL}.UxxU.hq.bam -
    fi
    wait
    $python $statpy_noTag ${RUN}.${LABEL}.UU.hq.bam | gzip > ${RUN}.${LABEL}.UU.hq.statCB.tsv.gz &
    $python $statpy_noTag ${RUN}.${LABEL}.UxxU.hq.bam | gzip > ${RUN}.${LABEL}.UxxU.hq.statCB.tsv.gz &
    $samtools view ${RUN}.${LABEL}.UU.hq.bam | $python $statpy_bamTag | gzip > ${RUN}.${LABEL}.UU.hq.CB.statPET.tsv.gz &
    $samtools view -@ ${NTHREADmap} -h ${RUN}.${LABEL}.UU.hq.bam | grep -v noCB | $samtools sort -@ $NTHREADmap -o ${RUN}.${LABEL}.UU.hq.CB.csorted.bam -
    $samtools view -@ ${NTHREADmap} -h ${RUN}.${LABEL}.UxxU.hq.bam | grep -v noCB | $samtools sort -@ $NTHREADmap -o ${RUN}.${LABEL}.UxxU.hq.CB.csorted.bam -
    wait
    $gatk MarkDuplicates --BARCODE_TAG=ZB --REMOVE_DUPLICATES=true --ADD_PG_TAG_TO_READS=false -M ${RUN}.${LABEL}.UU.hq.CB.nr.MarkDuplicates.txt -O ${RUN}.${LABEL}.UU.hq.CB.csorted.dedup.bam -I ${RUN}.${LABEL}.UU.hq.CB.csorted.bam && $samtools sort -n -@ $NTHREADmap -o ${RUN}.${LABEL}.UU.hq.CB.nr.bam ${RUN}.${LABEL}.UU.hq.CB.csorted.dedup.bam &
    $gatk MarkDuplicates --BARCODE_TAG=ZB --REMOVE_DUPLICATES=true --ADD_PG_TAG_TO_READS=false -M ${RUN}.${LABEL}.UxxU.hq.CB.nr.MarkDuplicates.txt -O ${RUN}.${LABEL}.UxxU.hq.CB.nr.bam -I ${RUN}.${LABEL}.UxxU.hq.CB.csorted.bam
    wait
    $gatk CollectInsertSizeMetrics -I ${RUN}.${LABEL}.UU.hq.CB.csorted.dedup.bam -O ${RUN}.${LABEL}.UU.hq.CB.csorted.dedup.insert_size_metrics.txt -H ${RUN}.${LABEL}.UU.hq.CB.csorted.dedup.insert_size_histogram.pdf -M 0.05 &
    $python $statpy_noTag ${RUN}.${LABEL}.UU.hq.CB.nr.bam | gzip > ${RUN}.${LABEL}.UU.hq.CB.nr.statCB.tsv.gz &
    $python $statpy_noTag ${RUN}.${LABEL}.UxxU.hq.CB.nr.bam | gzip > ${RUN}.${LABEL}.UxxU.hq.CB.nr.statCB.tsv.gz
    $samtools view ${RUN}.${LABEL}.UU.hq.CB.nr.bam | $python $statpy_bamTag | gzip > ${RUN}.${LABEL}.UU.hq.CB.nr.statPET.tsv.gz &
    $python ${filterinsertsizeless1000} ${RUN}.${LABEL}.UU.hq.CB.nr.bam ${RUN}.${LABEL}.UU.hq.CB.nr &
    $samtools index ${RUN}.${LABEL}.UU.hq.CB.csorted.dedup.bam
    $bamCoverage -b ${RUN}.${LABEL}.UU.hq.CB.csorted.dedup.bam -o ${RUN}.${LABEL}.UU.hq.CB.nr.bw -p $NTHREADmap -bs 5 && $bamCoverage -b ${RUN}.${LABEL}.UU.hq.CB.csorted.dedup.bam -o ${RUN}.${LABEL}.UU.hq.CB.nr.bedgraph -p $NTHREADmap -bs 5 -of bedgraph && $pigz -p ${NTHREADmap} ${RUN}.${LABEL}.UU.hq.CB.nr.bedgraph && $bamCoverage -b ${RUN}.${LABEL}.UU.hq.CB.csorted.dedup.bam -o ${RUN}.${LABEL}.UU.hq.CB.nr.RPKM.bw -p $NTHREADmap -bs 5 --normalizeUsing RPKM && $bamCoverage -b ${RUN}.${LABEL}.UU.hq.CB.csorted.dedup.bam -o ${RUN}.${LABEL}.UU.hq.CB.nr.RPKM.bedgraph -p $NTHREADmap -bs 5 -of bedgraph --normalizeUsing RPKM && $pigz -p ${NTHREADmap} ${RUN}.${LABEL}.UU.hq.CB.nr.RPKM.bedgraph
    $samtools index ${RUN}.${LABEL}.UxxU.hq.CB.nr.bam 
    $bamCoverage -b ${RUN}.${LABEL}.UxxU.hq.CB.nr.bam -o ${RUN}.${LABEL}.UxxU.hq.CB.nr.bw -p $NTHREADmap -bs 5 && $bamCoverage -b ${RUN}.${LABEL}.UxxU.hq.CB.nr.bam -o ${RUN}.${LABEL}.UxxU.hq.CB.nr.bedgraph -p $NTHREADmap -bs 5 -of bedgraph && $pigz -p ${NTHREADmap} ${RUN}.${LABEL}.UxxU.hq.CB.nr.bedgraph && $bamCoverage -b ${RUN}.${LABEL}.UxxU.hq.CB.nr.bam -o ${RUN}.${LABEL}.UxxU.hq.CB.nr.RPKM.bw -p $NTHREADmap -bs 5 --normalizeUsing RPKM && $bamCoverage -b ${RUN}.${LABEL}.UxxU.hq.CB.nr.bam -o ${RUN}.${LABEL}.UxxU.hq.CB.nr.RPKM.bedgraph -p $NTHREADmap -bs 5 -of bedgraph --normalizeUsing RPKM && $pigz -p ${NTHREADmap} ${RUN}.${LABEL}.UxxU.hq.CB.nr.RPKM.bedgraph
    wait
    $python $statpy_noTag ${RUN}.${LABEL}.UU.hq.CB.nr.pet.bam | gzip > ${RUN}.${LABEL}.UU.hq.CB.nr.pet.statCB.tsv.gz &
    $python $statpy_noTag ${RUN}.${LABEL}.UU.hq.CB.nr.other.bam | gzip > ${RUN}.${LABEL}.UU.hq.CB.nr.other.statCB.tsv.gz &
    $python $statpy_noTag ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.bam | gzip > ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.statCB.tsv.gz &
    $samtools view ${RUN}.${LABEL}.UU.hq.CB.nr.pet.bam | $python $statpy_bamTag | gzip > ${RUN}.${LABEL}.UU.hq.CB.nr.pet.statPET.tsv.gz &
    $samtools view ${RUN}.${LABEL}.UU.hq.CB.nr.other.bam | $python $statpy_bamTag | gzip > ${RUN}.${LABEL}.UU.hq.CB.nr.other.statPET.tsv.gz &
    $samtools view ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.bam | $python $statpy_bamTag | gzip > ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.statPET.tsv.gz &
    $samtools view -h ${RUN}.${LABEL}.UU.hq.CB.nr.pet.bam | $perl $injectCBreadID | $bedtools bamtobed -bedpe | sort -k7,7 > ${RUN}.${LABEL}.UU.hq.CB.nr.bedpe &
    $samtools sort -@ $NTHREADmap -o ${RUN}.${LABEL}.UU.hq.CB.nr.pet.csorted.bam ${RUN}.${LABEL}.UU.hq.CB.nr.pet.bam
    $samtools sort -@ $NTHREADmap -o ${RUN}.${LABEL}.UU.hq.CB.nr.other.csorted.bam ${RUN}.${LABEL}.UU.hq.CB.nr.other.bam
    $samtools sort -@ $NTHREADmap -o ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.csorted.bam ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.bam
    $gatk CollectInsertSizeMetrics -I ${RUN}.${LABEL}.UU.hq.CB.nr.pet.csorted.bam -O ${RUN}.${LABEL}.UU.hq.CB.nr.pet.csorted.insert_size_metrics.txt -H ${RUN}.${LABEL}.UU.hq.CB.nr.pet.csorted.insert_size_histogram.pdf -M 0.05 &
    $gatk CollectInsertSizeMetrics -I ${RUN}.${LABEL}.UU.hq.CB.nr.other.csorted.bam -O ${RUN}.${LABEL}.UU.hq.CB.nr.other.csorted.insert_size_metrics.txt -H ${RUN}.${LABEL}.UU.hq.CB.nr.other.csorted.insert_size_histogram.pdf -M 0.05 &
    $gatk CollectInsertSizeMetrics -I ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.csorted.bam -O ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.csorted.insert_size_metrics.txt -H ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.csorted.insert_size_histogram.pdf -M 0.05 &
    $samtools index ${RUN}.${LABEL}.UU.hq.CB.nr.pet.csorted.bam
    $samtools index ${RUN}.${LABEL}.UU.hq.CB.nr.other.csorted.bam
    $samtools index ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.csorted.bam
    $bamCoverage -b ${RUN}.${LABEL}.UU.hq.CB.nr.pet.csorted.bam -o ${RUN}.${LABEL}.UU.hq.CB.nr.pet.csorted.bw -p $NTHREADmap -bs 5 && $bamCoverage -b ${RUN}.${LABEL}.UU.hq.CB.nr.pet.csorted.bam -o ${RUN}.${LABEL}.UU.hq.CB.nr.pet.csorted.bedgraph -p $NTHREADmap -bs 5 -of bedgraph && $pigz -p ${NTHREADmap} ${RUN}.${LABEL}.UU.hq.CB.nr.pet.csorted.bedgraph && $bamCoverage -b ${RUN}.${LABEL}.UU.hq.CB.nr.pet.csorted.bam -o ${RUN}.${LABEL}.UU.hq.CB.nr.pet.csorted.RPKM.bw -p $NTHREADmap -bs 5 --normalizeUsing RPKM && $bamCoverage -b ${RUN}.${LABEL}.UU.hq.CB.nr.pet.csorted.bam -o ${RUN}.${LABEL}.UU.hq.CB.nr.pet.csorted.RPKM.bedgraph -p $NTHREADmap -bs 5 -of bedgraph --normalizeUsing RPKM && $pigz -p ${NTHREADmap} ${RUN}.${LABEL}.UU.hq.CB.nr.pet.csorted.RPKM.bedgraph
    $bamCoverage -b ${RUN}.${LABEL}.UU.hq.CB.nr.other.csorted.bam -o ${RUN}.${LABEL}.UU.hq.CB.nr.other.csorted.bw -p $NTHREADmap -bs 5 && $bamCoverage -b ${RUN}.${LABEL}.UU.hq.CB.nr.other.csorted.bam -o ${RUN}.${LABEL}.UU.hq.CB.nr.other.csorted.bedgraph -p $NTHREADmap -bs 5 -of bedgraph && $pigz -p ${NTHREADmap} ${RUN}.${LABEL}.UU.hq.CB.nr.other.csorted.bedgraph && $bamCoverage -b ${RUN}.${LABEL}.UU.hq.CB.nr.other.csorted.bam -o ${RUN}.${LABEL}.UU.hq.CB.nr.other.csorted.RPKM.bw -p $NTHREADmap -bs 5 --normalizeUsing RPKM && $bamCoverage -b ${RUN}.${LABEL}.UU.hq.CB.nr.other.csorted.bam -o ${RUN}.${LABEL}.UU.hq.CB.nr.other.csorted.RPKM.bedgraph -p $NTHREADmap -bs 5 -of bedgraph --normalizeUsing RPKM && $pigz -p ${NTHREADmap} ${RUN}.${LABEL}.UU.hq.CB.nr.other.csorted.RPKM.bedgraph
    $bamCoverage -b ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.csorted.bam -o ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.csorted.bw -p $NTHREADmap -bs 5 && $bamCoverage -b ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.csorted.bam -o ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.csorted.bedgraph -p $NTHREADmap -bs 5 -of bedgraph && $pigz -p ${NTHREADmap} ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.csorted.bedgraph && $bamCoverage -b ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.csorted.bam -o ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.csorted.RPKM.bw -p $NTHREADmap -bs 5 --normalizeUsing RPKM && $bamCoverage -b ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.csorted.bam -o ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.csorted.RPKM.bedgraph -p $NTHREADmap -bs 5 -of bedgraph --normalizeUsing RPKM && $pigz -p ${NTHREADmap} ${RUN}.${LABEL}.UU.hq.CB.nr.singleton.csorted.RPKM.bedgraph
    wait

else
    if [ ${blacklist} != "null" ];then
        /usr/bin/time -f 'real %E\nuser %U\nsys  %S\nmem  %M\n' $bedtools intersect -nonamecheck -a ${RUN}.${LABEL}.UxxU.bam -b ${blacklist} -v | $samtools view -@ ${halfN} -F 2048 -F 4 -q ${mapq} -h - | $python $pickUpy ${mapq} | $samtools view -@ ${halfN} -Sb -o ${RUN}.${LABEL}.UxxU.hq.bam -
    else
        /usr/bin/time -f 'real %E\nuser %U\nsys  %S\nmem  %M\n' $samtools view -@ ${halfN} -F 2048 -F 4 -q ${mapq} -h ${RUN}.${LABEL}.UxxU.bam | $python $pickUpy ${mapq} | $samtools view -@ ${halfN} -Sb -o ${RUN}.${LABEL}.UxxU.hq.bam -
    fi
    $python $statpy_noTag ${RUN}.${LABEL}.UxxU.hq.bam | gzip > ${RUN}.${LABEL}.UxxU.hq.statCB.tsv.gz &
    $samtools view -@ ${NTHREADmap} -h ${RUN}.${LABEL}.UxxU.hq.bam | grep -v noCB | $samtools sort -@ $NTHREADmap -o ${RUN}.${LABEL}.UxxU.hq.CB.csorted.bam -
    wait
    $gatk MarkDuplicates --BARCODE_TAG=ZB --REMOVE_DUPLICATES=true --ADD_PG_TAG_TO_READS=false -M ${RUN}.${LABEL}.UxxU.hq.CB.nr.MarkDuplicates.txt -O ${RUN}.${LABEL}.UxxU.hq.CB.nr.bam -I ${RUN}.${LABEL}.UxxU.hq.CB.csorted.bam
    $python $statpy_noTag ${RUN}.${LABEL}.UxxU.hq.CB.nr.bam | gzip > ${RUN}.${LABEL}.UxxU.hq.CB.nr.statCB.tsv.gz &
    $samtools index ${RUN}.${LABEL}.UxxU.hq.CB.nr.bam
    $bamCoverage -b ${RUN}.${LABEL}.UxxU.hq.CB.nr.bam -o ${RUN}.${LABEL}.UxxU.hq.CB.nr.bw -p $NTHREADmap -bs 5 && $bamCoverage -b ${RUN}.${LABEL}.UxxU.hq.CB.nr.bam -o ${RUN}.${LABEL}.UxxU.hq.CB.nr.bedgraph -p $NTHREADmap -bs 5 -of bedgraph && $pigz -p ${NTHREADmap} ${RUN}.${LABEL}.UxxU.hq.CB.nr.bedgraph && $bamCoverage -b ${RUN}.${LABEL}.UxxU.hq.CB.nr.bam -o ${RUN}.${LABEL}.UxxU.hq.CB.nr.RPKM.bw -p $NTHREADmap -bs 5 --normalizeUsing RPKM && $bamCoverage -b ${RUN}.${LABEL}.UxxU.hq.CB.nr.bam -o ${RUN}.${LABEL}.UxxU.hq.CB.nr.RPKM.bedgraph -p $NTHREADmap -bs 5 -of bedgraph --normalizeUsing RPKM && $pigz -p ${NTHREADmap} ${RUN}.${LABEL}.UxxU.hq.CB.nr.RPKM.bedgraph
    wait

fi

echo "success" > ../slurm/$jobname.log
date 
ENDHERE

    if [ ${LABEL} == ${pairlabel} ];then
        jidCBdedup1=$(submitjobOptional8 $pbs $prevJobid $jidpkl)
        echo "cpuMapping $LABEL $jidCBdedup1 $pbs wait $prevJobid $jidpkl" >&2
    elif [ ${LABEL} == ${nonelabel} ];then
        jidCBdedup2=$(submitjobOptional8 $pbs $prevJobid $jidpkl)
        echo "cpuMapping $LABEL $jidCBdedup2 $pbs wait $prevJobid $jidpkl" >&2
    else
        jidCBdedup3=$(submitjobOptional8 $pbs $prevJobid $jidpkl)
        echo "cpuMapping $LABEL $jidCBdedup3 $pbs wait $prevJobid $jidpkl" >&2
    fi
}


function cpuCluster {
    # cluster calling
    fqDNApref=$1
    RUN=$(basename $fqDNApref)
    local prevJobid1=$2
    local prevJobid2=$3
    cpudir=$maindir/$RUN/bulkChiatac
    local cis_prefix="${RUN}.e500.clusters.cis.chiasig"
    local juicefile=${RUN}.e500.juice
    jobname=${RUN}.cpuCluster
    pbs=${jobname}.slurm

    cat << ENDHERE > $pbs
#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -n $NTHREAD
#SBATCH --mem ${GB}G
#SBATCH --partition=${QUEUE}
#SBATCH --job-name $jobname
#SBATCH --output="$jobname.slurm.o"
#SBATCH --error="$jobname.slurm.e"

set -e 
cd $cpudir
export SINGULARITY_BIND=$bwarefabs,$blacklistabs,$sifpathabs,$ARCrefabs,$GEXrefabs,$ATACrefabs,$codepathabs,$maindir,$wdir,$genebedabs
date

$samtools view -H ${RUN}.${pairlabel}.UU.hq.CB.nr.pet.bam | grep '^@SQ' | cut -f 2-3 | sed s?SN:?? | sed s?LN:?? > ${RUN}.chrlength.txt &
$cpu span -s $bulkSelfbp -t 4 ${RUN}.${pairlabel}.UU.hq.CB.nr.pet.bam 1> ${RUN}.$pairlabel.UU.nr.span.xls &
$cpu span -s $bulkSelfbp -t 4 ${RUN}.${nonelabel}.UU.hq.CB.nr.pet.bam 1> ${RUN}.$nonelabel.UU.nr.span.xls &
$samtools merge -@ $NTHREAD -f -n -o ${RUN}.hq.CB.nr.bam ${RUN}.${pairlabel}.UU.hq.CB.nr.pet.bam ${RUN}.${nonelabel}.UU.hq.CB.nr.pet.bam
$cpu span -s $bulkSelfbp -t 4 ${RUN}.hq.CB.nr.bam 1>${RUN}.hq.CB.nr.span.xls &
cat ${RUN}.${pairlabel}.UU.hq.CB.nr.bedpe ${RUN}.${nonelabel}.UU.hq.CB.nr.bedpe > ${RUN}.hq.CB.nr.bedpe
$perl -lane '@tmp=split(/:/,\$F[6],2);print join("\t",@F[0..5])."\t".join("\t",@tmp)."\t".join("\t",@F[7..9])' ${RUN}.hq.CB.nr.bedpe > ${RUN}.hq.CB.nr.tmp.bedpe
$python ${CountBarcodePET} ${RUN}.hq.CB.nr.tmp.bedpe ${RUN}.hq.CB.nr.bedpe.stat
wait
tmpid=$(echo ${RUN} | awk -F"-" '{print $1}')
$python ${bedpe2mat} ${RUN}.hq.CB.nr.tmp.bedpe ${RUN}.1000k.bedpe.mat 1000 $tmpid && $Rscript ${petpca} ${RUN}.1000k.bedpe.mat ${RUN}.1000k.bedpe.matrix.rds &


cat ${RUN}.hq.CB.nr.bedpe | awk -v OFS="\t" '{print 0,\$1,int((\$2+\$3)/2),0,0,\$4,int((\$5+\$6)/2),1}' | grep -v -e chrM | LC_COLLATE=C sort --parallel=${NTHREAD} -k2,2 -k6,6 -k3,3n -k7,7n > ${RUN}.hq.CB.nr.sorted.4dn
/usr/bin/time -f 'real %E\nuser %U\nsys  %S\nmem  %M\n' $juicer_tools pre -r 2500000,2000000,1000000,500000,250000,200000,100000,50000,25000,20000,10000,5000,2500,2000,1000 ${RUN}.hq.CB.nr.sorted.4dn ${RUN}.hq.CB.nr.sorted.4dn.hic ${RUN}.chrlength.txt &

/usr/bin/time -f 'real %E\nuser %U\nsys  %S\nmem  %M\n' $cpu cluster -s $bulkSelfbp -M -B 1000 -5 5,-20 -3 5,480 -t ${NTHREAD} -O ${RUN}.e500 -j -x -v 1 -g ${RUN}.hq.CB.nr.bam
zcat ${RUN}.e500.clusters.cis.chiasig.gz | awk '{ if ( \$7 > 1 ) print }' > ${RUN}.e500.clusters.cis.chiasig.BE2
/usr/bin/time -f 'real %E\nuser %U\nsys  %S\nmem  %M\n' $ChiaSig -c $minPETsignif -t $NTHREAD -p -m $bulkSelfbp ${RUN}.e500.clusters.cis.chiasig.BE2 &
zcat ${RUN}.e500.juice.gz | grep -v -e chrM > ${RUN}.e500.juice
/usr/bin/time -f 'real %E\nuser %U\nsys  %S\nmem  %M\n' $juicer_tools pre -r 2500000,2000000,1000000,500000,250000,200000,100000,50000,25000,20000,10000,5000,2500,2000,1000 ${RUN}.e500.juice $RUN.disgt8k.hic ${RUN}.chrlength.txt &
wait

date
echo "success" > ../slurm/$jobname.log
ENDHERE

    jidChiasig=$(submitjobOptional8 $pbs $prevJobid1 $prevJobid2)
    echo "cpuCluster $jidChiasig $pbs wait $prevJobid1 $prevJobid2" >&2
}


############################################### generate scATAC matrix
function mergeFinalBAM {
    # merge all bam file to generate scATAC matrix
    local fqDNApref=$1
    RUN=$(basename $fqDNApref)
    local j1=$2
    local j2=$3
    local j3=$4
    wdir=$maindir/$RUN
    cpudir=$maindir/$RUN/bulkChiatac
    jobdir=$maindir/$RUN/slurm
    jobname=${RUN}.mergeFinalBAM
    pbs=${jobname}.slurm
    tsvfile=${RUN}.fragments.tsv
    halfN=$(echo $NTHREADmap/2 | bc )

    cat << ENDHERE > $pbs
#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -n $NTHREADmap
#SBATCH --mem ${GB}G
#SBATCH --partition=${QUEUE}
#SBATCH --job-name $jobname
#SBATCH --output="$jobname.slurm.o"
#SBATCH --error="$jobname.slurm.e"

set -e 
cd $cpudir
export SINGULARITY_BIND=$bwarefabs,$blacklistabs,$sifpathabs,$ARCrefabs,$GEXrefabs,$ATACrefabs,$codepathabs,$maindir,$wdir,$genebedabs
date

$samtools view -H ${RUN}.${nonelabel}.UU.hq.CB.nr.bam | grep '^@SQ' | cut -f 2-3 | sed s?SN:?? | sed s?LN:?? > ${genome}.chrlength.txt &
# for peak calling
if [ ${alldata} != "T" ];then # just use ${RUN}.${nonelabel}.UU.hq.CB.nr.singleton.bam/${RUN}.${nonelabel}.UU.hq.CB.nr.singleton.csorted.bam
    $samtools view -@ $halfN -h ${RUN}.${nonelabel}.UU.hq.CB.nr.singleton.csorted.bam | $perl $injectCBreadID | $samtools view -@ $halfN -o ${RUN}.forpeak.bam -
else # use all align reads
    $samtools merge -f -@ $NTHREADmap -o ${RUN}.forpeak.tmp.bam ${RUN}.${pairlabel}.UU.hq.CB.csorted.dedup.bam ${RUN}.${pairlabel}.UxxU.hq.CB.nr.bam ${RUN}.${singlabel}.UxxU.hq.CB.nr.bam ${RUN}.${nonelabel}.UU.hq.CB.csorted.dedup.bam ${RUN}.${nonelabel}.UxxU.hq.CB.nr.bam
    $samtools view -@ $halfN -h ${RUN}.forpeak.tmp.bam | $perl $injectCBreadID | $samtools view -@ $halfN -o ${RUN}.forpeak.bam -
    rm ${RUN}.forpeak.tmp.bam
    samtools index ${RUN}.forpeak.bam
fi
$macs2 callpeak --keep-dup all --nomodel --shift $shiftsize --extsize $peakext -B --SPMR -t ${RUN}.forpeak.bam -g $glen -n $RUN  --qvalue $macsq &
$bamCoverage -b ${RUN}.forpeak.bam -o ${RUN}.atac.allreads.bw -p $NTHREADmap -bs 5 && $bamCoverage -b ${RUN}.forpeak.bam -o ${RUN}.atac.allreads.bedgraph -p $NTHREADmap -bs 5 -of bedgraph && $pigz -p ${NTHREADmap} ${RUN}.atac.allreads.bedgraph && $bamCoverage -b ${RUN}.forpeak.bam -o ${RUN}.atac.allreads.RPKM.bw -p $NTHREADmap -bs 5 --normalizeUsing RPKM && $bamCoverage -b ${RUN}.forpeak.bam -o ${RUN}.atac.allreads.RPKM.bedgraph -p $NTHREADmap -bs 5 -of bedgraph --normalizeUsing RPKM && $pigz -p ${NTHREADmap} ${RUN}.atac.allreads.RPKM.bedgraph
ln -s -f ${RUN}.none.UU.hq.CB.nr.singleton.csorted.bw ${RUN}.atac.perfectreads.bw
ln -s -f ${RUN}.none.UU.hq.CB.nr.singleton.csorted.RPKM.bw ${RUN}.atac.perfectreads.RPKM.bw
ln -s -f ${RUN}.none.UU.hq.CB.nr.singleton.csorted.bedgraph.gz ${RUN}.atac.perfectreads.bedgraph.gz
ln -s -f ${RUN}.none.UU.hq.CB.nr.singleton.csorted.RPKM.bedgraph.gz ${RUN}.atac.perfectreads.RPKM.bedgraph.gz

# For generating fragment file, ${RUN}.${nonelabel}.UU.hq.CB.nr.singleton.bam is the perfect file to generate the single cell fragment. For singleton reads: paired-end, startpos+4, endpos-5. For other, pet and single end: single-end, +:+4; -:-5.
$python $bamtofrag -p ${RUN}.${pairlabel}.UU.hq.CB.nr.singleton.bam,${RUN}.${nonelabel}.UU.hq.CB.nr.singleton.bam -s ${RUN}.${pairlabel}.UU.hq.CB.nr.pet.bam,${RUN}.${pairlabel}.UU.hq.CB.nr.other.bam,${RUN}.${pairlabel}.UxxU.hq.CB.nr.bam,${RUN}.${nonelabel}.UU.hq.CB.nr.pet.bam,${RUN}.${nonelabel}.UU.hq.CB.nr.other.bam,${RUN}.${nonelabel}.UxxU.hq.CB.nr.bam,${RUN}.${singlabel}.UxxU.hq.CB.nr.bam -o ${tsvfile} &
$python $catBam ${RUN}.${pairlabel}.UU.hq.CB.nr.singleton.bam,${RUN}.${nonelabel}.UU.hq.CB.nr.singleton.bam ${RUN}.fragment.paired-end.bam
$python $catBam ${RUN}.${pairlabel}.UU.hq.CB.nr.pet.bam,${RUN}.${pairlabel}.UU.hq.CB.nr.other.bam,${RUN}.${pairlabel}.UxxU.hq.CB.nr.bam,${RUN}.${nonelabel}.UU.hq.CB.nr.pet.bam,${RUN}.${nonelabel}.UU.hq.CB.nr.other.bam,${RUN}.${nonelabel}.UxxU.hq.CB.nr.bam,${RUN}.${singlabel}.UxxU.hq.CB.nr.bam ${RUN}.fragment.single-end.bam
wait
$pigz -p ${NTHREADmap} ${tsvfile}.len
LC_COLLATE=C sort --parallel=${NTHREADmap} -k1,1 -k2,2n -k3,3n ${tsvfile} > ${tsvfile}.tmp && mv ${tsvfile}.tmp ${tsvfile}
$bgzip -f ${tsvfile}
$tabix -f -p bed ${tsvfile}.gz
zcat ${RUN}.fragments.tsv.len.gz | awk '{aa[\$1]+=1}END{for(kk in aa){print kk"\t"aa[kk]}}' | sort -k1,1n > ${RUN}.ATAC.fragment.length.distribution.txt
cp ${tsvfile}.gz atac_fragments.tsv.gz
cp ${tsvfile}.gz.tbi atac_fragments.tsv.gz.tbi

$bedClip ${RUN}_treat_pileup.bdg ${genome}.chrlength.txt ${RUN}_treat_pileup.clip.bdg
LC_COLLATE=C sort --parallel=${NTHREADmap} -k1,1 -k2,2n -k3,3n ${RUN}_treat_pileup.clip.bdg > ${RUN}_treat_pileup.clip.sorted.bdg
$bedGraphToBigWig ${RUN}_treat_pileup.clip.sorted.bdg ${genome}.chrlength.txt ${RUN}.treat_pileup.NDP.bw &
cut -f 1-3 ${RUN}_peaks.narrowPeak | awk -v OFS="\t" '{if(NR==FNR){if(\$1=="chrM"){aa = \$2}}else{if(\$1=="chrM"){if(\$3>aa){\$3=aa;}} print \$0}}' ${genome}.chrlength.txt - > ${RUN}.peak
wait

cat <($samtools view -H ../10xGEX/outs/possorted_genome_bam.bam) <($samtools view -@ ${halfN} ../10xGEX/outs/possorted_genome_bam.bam | grep xf:i:25) | $samtools view -@ ${halfN} -Sb -o ${RUN}.gex.bam - 
$samtools index ${RUN}.gex.bam
$python $gex2fragment ${RUN}.gex.bam ${RUN}.gex.fragments.unsorted.tsv 2> ${RUN}.gex.fragments.unsorted.tsv.log
LC_COLLATE=C sort --parallel=${NTHREADmap} -k1,1 -k2,2n -k3,3n ${RUN}.gex.fragments.unsorted.tsv | $pigz -p ${NTHREADmap} > ${RUN}.gex.fragments.tsv.gz
$bamCoverage -b ${RUN}.gex.bam -o ${RUN}.gex.bw -p $NTHREADmap -bs 5 && $bamCoverage -b ${RUN}.gex.bam -o ${RUN}.gex.bedgraph -p $NTHREADmap -bs 5 -of bedgraph && $pigz -p ${NTHREADmap} ${RUN}.gex.bedgraph && $bamCoverage -b ${RUN}.gex.bam -o ${RUN}.gex.RPKM.bw -p $NTHREADmap -bs 5 --normalizeUsing RPKM && $bamCoverage -b ${RUN}.gex.bam -o ${RUN}.gex.RPKM.bedgraph -p $NTHREADmap -bs 5 -of bedgraph --normalizeUsing RPKM && $pigz -p ${NTHREADmap} ${RUN}.gex.RPKM.bedgraph

cd $wdir
if [ -d "10xATAC" ]; then
    rm -rf 10xATAC
fi
/usr/bin/time -f 'real %E\nuser %U\nsys  %S\nmem  %M\n' ${cellrangeratac} reanalyze --id 10xATAC --peaks $cpudir/${RUN}.peak --fragments $cpudir/atac_fragments.tsv.gz --reference $ATACref --localcores ${NTHREADmap}
wait

date
cd $cpudir && echo "success" > ../slurm/$jobname.log
ENDHERE

    jm=$(submitjobOptional8  $pbs $j1 $j2 $j3)
    echo "mergeFinalBAM $jm $pbs wait $j1 $j2 $j3" >&2
}


############################################### summary
function summaryData {
    # summary
    local fqDNApref=$1
    local j1=$2
    RUN=$(basename $fqDNApref)
    wdir=$maindir/$RUN
    cpudir=$maindir/$RUN/bulkChiatac
    jobdir=$maindir/$RUN/slurm
    jobname=${RUN}.summaryData
    pbs=${jobname}.slurm
    halfN=$(echo $NTHREAD/2 | bc )

    cat << ENDHERE > $pbs
#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -n $NTHREAD
#SBATCH --mem ${smallGB}G
#SBATCH --partition=${smallQUEUE}
#SBATCH --job-name $jobname
#SBATCH --output="$jobname.slurm.o"
#SBATCH --error="$jobname.slurm.e"

set -e 
cd $cpudir
export SINGULARITY_BIND=$bwarefabs,$blacklistabs,$sifpathabs,$ARCrefabs,$GEXrefabs,$ATACrefabs,$codepathabs,$maindir,$wdir,$genebedabs
date

cd ${cpudir}

$python $mergeTable ${RUN}
bash $sumsh $RUN $genome ${siffile} ${getfrom10xG} ${getfrom10xA} ${mapq}
$Rscript $getGEXbarcode ../10xGEX/outs/raw_feature_bc_matrix ${RUN}.gex.barcode.txt
# awk '{aa[\$7]+=1}END{print "nPET\tbarcode"; for(k in aa){print k"\t"aa[k]"\t"k}}' ${RUN}.hq.CB.nr.tmp.bedpe > ${RUN}.pet.barcode.txt
awk '{if(\$1==\$4){aa[\$7]["cis"]+=1}else{aa[\$7]["trans"]+=1}}END{print "nPET\tbarcode\tcis\ttrans"; for(k in aa){aa[k]["cis"]+=0; aa[k]["trans"]+=0; print k"\t"aa[k]["cis"]+aa[k]["trans"]"\t"k"\t"aa[k]["cis"]"\t"aa[k]["trans"]}}' ${RUN}.hq.CB.nr.tmp.bedpe > ${RUN}.pet.barcode.txt
zcat ${RUN}.fragments.tsv.gz | awk '{aa[\$4]+=1}END{print "nFragment\tbarcode"; for(k in aa){print k"\t"aa[k]"\t"k}}' > ${RUN}.fragment.barcode.txt
$Rscript $barcodePlot ${RUN}.gex.barcode.txt ${RUN}.fragment.barcode.txt ${RUN}.pet.barcode.txt ${RUN}.barcodeRankPlot.pdf
$Rscript $plotPETnumber ${RUN}.${genome}.bulk_stats.tsv ${RUN}.loop.PETcount.distribution.pdf
$Rscript $plotFragmentlen ${RUN}.ATAC.fragment.length.distribution.txt ${RUN}.ATAC.fragment.length.distribution.pdf
$Rscript $plotStat ${RUN}.${genome}.bulk_stats.tsv ${RUN}.dataStat.pdf

$computeMatrix scale-regions -q -S ${RUN}.gex.RPKM.bw -R ${genebed} -a 3000 -b 3000 -m 3000 --skipZeros --missingDataAsZero --outFileName ${RUN}.gex.RPKM.genebody.gz -p ${NTHREAD} --maxThreshold 50000
$computeMatrix scale-regions -q -S ${RUN}.atac.perfectreads.RPKM.bw -R ${genebed} -a 3000 -b 3000 -m 3000 --skipZeros --missingDataAsZero --outFileName ${RUN}.atac.perfectreads.RPKM.genebody.gz -p ${NTHREAD} --maxThreshold 50000
$computeMatrix scale-regions -q -S ${RUN}.atac.allreads.RPKM.bw -R ${genebed} -a 3000 -b 3000 -m 3000 --skipZeros --missingDataAsZero --outFileName ${RUN}.atac.allreads.RPKM.genebody.gz -p ${NTHREAD} --maxThreshold 50000
$computeMatrix reference-point -q --referencePoint TSS -S ${RUN}.atac.perfectreads.RPKM.bw -R ${genebed} -a 3000 -b 3000 --skipZeros --missingDataAsZero --outFileName ${RUN}.atac.perfectreads.RPKM.geneTSS.gz -p ${NTHREAD} --maxThreshold 50000
$computeMatrix reference-point -q --referencePoint TSS -S ${RUN}.atac.allreads.RPKM.bw -R ${genebed} -a 3000 -b 3000 --skipZeros --missingDataAsZero --outFileName ${RUN}.atac.allreads.RPKM.geneTSS.gz -p ${NTHREAD} --maxThreshold 50000

plotHeatmap -m ${RUN}.gex.RPKM.genebody.gz -out ${RUN}.gex.RPKM.genebody.Heatmap.pdf && plotProfile -m ${RUN}.gex.RPKM.genebody.gz -out ${RUN}.gex.RPKM.genebody.Profile.pdf
plotHeatmap -m ${RUN}.atac.perfectreads.RPKM.genebody.gz -out ${RUN}.atac.perfectreads.RPKM.genebody.Heatmap.pdf && plotProfile -m ${RUN}.atac.perfectreads.RPKM.genebody.gz -out ${RUN}.atac.perfectreads.RPKM.genebody.Profile.pdf
plotHeatmap -m ${RUN}.atac.allreads.RPKM.genebody.gz -out ${RUN}.atac.allreads.RPKM.genebody.Heatmap.pdf && plotProfile -m ${RUN}.atac.allreads.RPKM.genebody.gz -out ${RUN}.atac.allreads.RPKM.genebody.Profile.pdf
plotHeatmap -m ${RUN}.atac.perfectreads.RPKM.geneTSS.gz -out ${RUN}.atac.perfectreads.RPKM.geneTSS.Heatmap.pdf && plotProfile -m ${RUN}.atac.perfectreads.RPKM.geneTSS.gz -out ${RUN}.atac.perfectreads.RPKM.geneTSS.Profile.pdf
plotHeatmap -m ${RUN}.atac.allreads.RPKM.geneTSS.gz -out ${RUN}.atac.allreads.RPKM.geneTSS.Heatmap.pdf && plotProfile -m ${RUN}.atac.allreads.RPKM.geneTSS.gz -out ${RUN}.atac.allreads.RPKM.geneTSS.Profile.pdf

wait

###### mkdir aligndir
###### mv *.none.fastq.gz *.singlelinker.paired.fastq.gz *.singlelinker.single.fastq.gz ${RUN}.stat ${RUN}.cpu aligndir
###### mv *.none.sam.gz *.singlelinker.paired.sam.gz *.singlelinker.single.sam.gz aligndir
###### mv *perfectreads*bedgraph.gz *perfectreads*bw *.hq.CB.nr*bw *.hq.CB.nr*bedgraph.gz *atac.perfectreads* aligndir
###### mv *.UxxU.bam *.UU.bam *.hq.bam *.hq.CB.bam aligndir
###### mv *.UxxU.statCB.tsv.gz *.UU.statCB.tsv.gz *.hq.statCB.tsv.gz *.hq.CB.statCB.tsv.gz aligndir
###### mv *.none.statCB.tsv.gz *.none.stat.xls aligndir

###### mkdir peakdir
###### mv *_summits.bed *_treat_pileup.bdg *_treat_pileup.clip.bdg *_treat_pileup.clip.sorted.bdg *_control_lambda.bdg peakdir

###### mv *.xx.bam *.nn.bam *.failed.bam *.discordant.bam *.chiasig *.e500.juice *.e500.juice.gz $dumpfolder

###### rm atac_fragments.tsv.gz atac_fragments.tsv.gz.tbi

date
echo "success" > ../slurm/$jobname.log
ENDHERE

    sm=$(submitjobOptional8  $pbs $j1)
    echo "summaryData $sm $pbs wait $j1" >&2
}


#------- loop & submit
fqlist=$1
bwarefabs=$(dirname $(realpath $bwaref))
blacklistabs=$(dirname $(realpath $blacklist))
sifpathabs=$(dirname $(realpath $siffile))
genebedabs=$(dirname $(realpath $genebed))
ARCrefabs=$(realpath $ARCref)
GEXrefabs=$(realpath $GEXref)
ATACrefabs=$(realpath $ATACref)
codepathabs=$(realpath $codepath)
while read dna rna ; do
   RUN=$(basename $dna)
   wdir=$maindir/$RUN
   cpudir=$wdir/bulkChiatac
   jobdir=$wdir/slurm
   export SINGULARITY_BIND=$bwarefabs,$blacklistabs,$sifpathabs,$ARCrefabs,$GEXrefabs,$ATACrefabs,$codepathabs,$maindir,$wdir,$genebedabs
   cd $maindir
   odir $RUN
   odir $jobdir
   odir $cpudir && odir $cpudir/$dumpfolder 
   cd $jobdir
   # cpuLinker $dna #will run from linker detection all the way to dedupCB
   # mapARC10x $dna $rna #run10x 
   # mapGEX10x $dna $rna #run10x scRNA if needed 
   

   # cpuMapping $dna $pairlabel $jidLinker
   # cpuMapping $dna $singlabel $jidLinker
   # cpuMapping $dna $nonelabel $jidLinker


   # cpuCluster $dna $jidCBdedup1 $jidCBdedup2


   mergeFinalBAM $dna $jidCBdedup1 $jidCBdedup2 $jidCBdedup3


   summaryData $dna $jm

   sleep 2
   echo "All jobs for $dna submitted"
done < $fqlist 
