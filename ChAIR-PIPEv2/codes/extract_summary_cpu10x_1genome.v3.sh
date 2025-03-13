#!/bin/bash
set -e

RUN=$1
GENOME1=$2
siffile=$3
Rscript="singularity run ${siffile} Rscript"
getfrom10xG=$4
getfrom10xA=$5
mapq=$6
pairlabel='singlelinker.paired'
singlelabel='singlelinker.single'
nonelabel='none'
extbp=500
grun=${RUN}.${GENOME1}
out_file="${grun}.bulk_stats.tsv"

#------------------------------------
echo Starting to sumarize ${grun}
echo "Output file (overwrite mode): $out_file" 

######################################################################### Get library ID
cat ${RUN}.stat | grep "Filename" | awk -F'[ \t.]' '{print "Library_ID\t"$2}' > ${out_file}
# Get reference genome
echo -e "Genome_origin\t"${GENOME1} >> ${out_file}

######################################################################### Get scATAC/PET count
n_pet=$( cat ${RUN}.stat | grep "Total pairs" | awk -F'[ \t]' '{print $3}' )
# Write PET count
n_pet=$( printf "%'.f" ${n_pet} )
echo -e "DNA_raw_reads\t"${n_pet} >> ${out_file}

## Get PET count with linker
nolinker_link=$( cat ${RUN}.stat | grep "No Linker" | awk -F '[ \t]' '{print $3}')
nolinker_linkf=$( printf "%'.f" ${nolinker_link} )
echo -e "DNA_raw_without_linker\t"${nolinker_linkf} >> ${out_file}
read_link=$( cat ${RUN}.stat | grep "Linker detected" | awk -F '[ \t]' '{print $3}')
read_linkf=$( printf "%'.f" ${read_link} )
echo -e "DNA_linker_detected_reads\t"${read_linkf} >> ${out_file}
pet_link=$( cat ${RUN}.stat | grep "Single Linker 2 tags" | awk -F '[ \t]' '{print $6}' )
pet_linkf=$( printf "%'.f" ${pet_link} )
echo -e "DNA_raw_2tag_with_linker\t"${pet_linkf} >> ${out_file}
single_link=$( cat ${RUN}.stat | grep "Single Linker 1 tag" | awk -F '[ \t]' '{print $6}' )
single_linkf=$( printf "%'.f" ${single_link} )
echo -e "DNA_raw_1tag_with_linker\t"${single_linkf} >> ${out_file}
other_link=$( echo $read_link - $pet_link - $single_link | bc )
other_linkf=$( printf "%'.f" ${other_link} )
echo -e "DNA_raw_others_with_linker\t"${other_linkf} >> ${out_file}


######################################################################### 2 tags
echo "" >> ${out_file}

pairlabeltotal=$(cat ${RUN}.$pairlabel.stat.xls | grep "TOTAL" | awk -F '[\t]' '{print $3}' | xargs printf "%'.f" )
echo -e "paired_reads/PET/2tag\t"${pairlabeltotal} >> ${out_file}

unique=$( cat ${RUN}.$pairlabel.stat.xls | grep "Uniquely Mapped" | awk -F '[\t]' '{print $3}' | xargs printf "%'.f" )
echo -e "PET_UU_all\t"${unique} >> ${out_file}

# pairUUmap=$(zcat ${RUN}.$pairlabel.UU.statCB.tsv.gz | awk '{sum+=$2}END{print sum}' | xargs printf "%'.f")
# echo -e "PET_UU_mapped\t"${pairUUmap} >> ${out_file}

pairUUuniq=$(zcat ${RUN}.$pairlabel.UU.hq.statCB.tsv.gz | awk '{sum+=$2}END{print sum}' | xargs printf "%'.f")
echo -e "PET_UU_mapped_q${mapq}\t"${pairUUuniq} >> ${out_file}

pairUUuniqCB=$(zcat ${RUN}.$pairlabel.UU.hq.statCB.tsv.gz | grep -v noCB | awk '{sum+=$2}END{print sum}' | xargs printf "%'.f")
echo -e "PET_UU_mapped_q${mapq}_validCB\t"${pairUUuniqCB} >> ${out_file}
pairUUuniqCB=$(zcat ${RUN}.$pairlabel.UU.hq.statCB.tsv.gz | grep -v noCB | awk '{sum+=$2}END{print sum}')

pairUUuniqnr=$(zcat ${RUN}.$pairlabel.UU.hq.CB.nr.statCB.tsv.gz | awk '{sum+=$2}END{print sum}' | xargs printf "%'.f")
echo -e "PET_UU_mapped_q${mapq}_validCB_dedup\t"${pairUUuniqnr} >> ${out_file}
pairUUuniqnr=$(zcat ${RUN}.$pairlabel.UU.hq.CB.nr.statCB.tsv.gz | awk '{sum+=$2}END{print sum}' )

uuPairDupRate=$(grep "Unknown Library"  ${RUN}.${pairlabel}.UU.hq.CB.nr.MarkDuplicates.txt | awk -F '[\t]' '{print $9*100}' |  xargs printf "%.1f%%" )
echo -e "PET_UU_Redundancy\t"${uuPairDupRate} >> ${out_file}

PETtrans=$(head -2 ${RUN}.full.CBsums.tsv | tail -1 | cut -f 59 | xargs printf "%'.f" )
PETcis=$(head -2 ${RUN}.full.CBsums.tsv | tail -1 | cut -f 60 | xargs printf "%'.f" )
PETcis1k=$(head -2 ${RUN}.full.CBsums.tsv | tail -1 | cut -f 61 | xargs printf "%'.f" )
PETcis8k=$(head -2 ${RUN}.full.CBsums.tsv | tail -1 | cut -f 62 | xargs printf "%'.f" )
echo -e "PET_trans\t"${PETtrans} >> ${out_file}
echo -e "PET_cis\t"${PETcis} >> ${out_file}
echo -e "PET_cis>1kb\t"${PETcis1k} >> ${out_file}
echo -e "PET_cis>8kb\t"${PETcis8k} >> ${out_file}
PETtrans=$(head -2 ${RUN}.full.CBsums.tsv | tail -1 | cut -f 59 )
PETcis=$(head -2 ${RUN}.full.CBsums.tsv | tail -1 | cut -f 60 )
PETcis1k=$(head -2 ${RUN}.full.CBsums.tsv | tail -1 | cut -f 61 )
PETcis8k=$(head -2 ${RUN}.full.CBsums.tsv | tail -1 | cut -f 62 )
ratio=$( echo "${PETcis} / ${PETtrans}" | bc -l )
ratio=$( printf %.2f ${ratio} )
echo -e "PET_ratio_cis/trans\t"${ratio} >> ${out_file}
ratio=$( echo "${PETcis1k} / ${PETtrans}" | bc -l )
ratio=$( printf %.2f ${ratio} )
echo -e "PET_ratio_cis>1kb/trans\t"${ratio} >> ${out_file}

uxxu=$( cat ${RUN}.$pairlabel.stat.xls | grep "Ux/xU" | awk -F '[\t]' '{print $3}' | xargs printf "%'.f" )
echo -e "PET_UxxU_all\t"${uxxu} >> ${out_file}

# pairUxxUmap=$(zcat ${RUN}.$pairlabel.UxxU.statCB.tsv.gz  | awk '{sum+=$2}END{print sum}' | xargs printf "%'.f")
# echo -e "PET_UxxU_mapped\t"${pairUxxUmap} >> ${out_file}

pairUxxUuniq=$(zcat ${RUN}.$pairlabel.UxxU.hq.statCB.tsv.gz  | awk '{sum+=$2}END{print sum}' | xargs printf "%'.f")
echo -e "PET_UxxU_mapped_q${mapq}\t"${pairUxxUuniq} >> ${out_file}

pairUxxUuniqCB=$(zcat ${RUN}.$pairlabel.UxxU.hq.statCB.tsv.gz  | grep -v noCB | awk '{sum+=$2}END{print sum}' | xargs printf "%'.f")
echo -e "PET_UxxU_mapped_q${mapq}_validCB\t"${pairUxxUuniqCB} >> ${out_file}
pairUxxUuniqCB=$(zcat ${RUN}.$pairlabel.UxxU.hq.statCB.tsv.gz  | grep -v noCB | awk '{sum+=$2}END{print sum}' )

pairUxxUuniqnr=$(zcat ${RUN}.$pairlabel.UxxU.hq.CB.nr.statCB.tsv.gz  | awk '{sum+=$2}END{print sum}' | xargs printf "%'.f")
echo -e "PET_UxxU_mapped_q${mapq}_validCB_dedup\t"${pairUxxUuniqnr} >> ${out_file}
pairUxxUuniqnr=$(zcat ${RUN}.$pairlabel.UxxU.hq.CB.nr.statCB.tsv.gz  | awk '{sum+=$2}END{print sum}' )

uxxuPairDupRate=$(grep "Unknown Library"  ${RUN}.${pairlabel}.UxxU.hq.CB.nr.MarkDuplicates.txt | awk -F '[\t]' '{print $9*100}' |  xargs printf "%.1f%%" )
echo -e "PET_UxxU_Redundancy\t"${uxxuPairDupRate} >> ${out_file} 


######################################################################### 1 tags
echo "" >> ${out_file}

singlelabeltotal=$(cat ${RUN}.$singlelabel.stat.xls | grep "TOTAL" | awk -F '[\t]' '{print $3}' | xargs printf "%'.f" )
# echo -e "single_reads (through filter)\t"${singlelabeltotal} >> ${out_file}
echo -e "single_reads\t"${singlelabeltotal} >> ${out_file}

uxxu=$( cat ${RUN}.$singlelabel.stat.xls | grep "Ux/xU" | awk -F '[\t]' '{print $3}' | xargs printf "%'.f" )
echo -e "single_UxxU_all\t"${uxxu} >> ${out_file}

# singleUxxUmap=$(zcat ${RUN}.$singlelabel.UxxU.statCB.tsv.gz  | awk '{sum+=$2}END{print sum}' | xargs printf "%'.f")
# echo -e "single_UxxU_mapped\t"${singleUxxUmap} >> ${out_file}

singleUxxUuniq=$(zcat ${RUN}.$singlelabel.UxxU.hq.statCB.tsv.gz  | awk '{sum+=$2}END{print sum}' | xargs printf "%'.f")
echo -e "single_UxxU_mapped_q${mapq}\t"${singleUxxUuniq} >> ${out_file}

singleUxxUuniqCB=$(zcat ${RUN}.$singlelabel.UxxU.hq.statCB.tsv.gz  | grep -v noCB | awk '{sum+=$2}END{print sum}' | xargs printf "%'.f")
echo -e "single_UxxU_mapped_q${mapq}_validCB\t"${singleUxxUuniqCB} >> ${out_file}
singleUxxUuniqCB=$(zcat ${RUN}.$singlelabel.UxxU.hq.statCB.tsv.gz  | grep -v noCB | awk '{sum+=$2}END{print sum}' )

singleUxxUuniqnr=$(zcat ${RUN}.$singlelabel.UxxU.hq.CB.nr.statCB.tsv.gz  | awk '{sum+=$2}END{print sum}' | xargs printf "%'.f")
echo -e "single_UxxU_mapped_q${mapq}_validCB_dedup\t"${singleUxxUuniqnr} >> ${out_file}
singleUxxUuniqnr=$(zcat ${RUN}.$singlelabel.UxxU.hq.CB.nr.statCB.tsv.gz  | awk '{sum+=$2}END{print sum}' )

uxxusingleDupRate=$(grep "Unknown Library"  ${RUN}.${singlelabel}.UxxU.hq.CB.nr.MarkDuplicates.txt | awk -F '[\t]' '{print $9*100}' |  xargs printf "%.1f%%" )
echo -e "single_UxxU_Redundancy\t"${uxxusingleDupRate} >> ${out_file}

######################################################################### 0 tags
echo "" >> ${out_file}

nonelabeltotal=$(cat ${RUN}.$nonelabel.stat.xls | grep "TOTAL" | awk -F '[\t]' '{print $3}' | xargs printf "%'.f" )
# echo -e "NoLinker_reads (through filter)\t"${nonelabeltotal} >> ${out_file}
echo -e "NoLinker_reads\t"${nonelabeltotal} >> ${out_file}

unique=$( cat ${RUN}.$nonelabel.stat.xls | grep "Uniquely Mapped" | awk -F '[\t]' '{print $3}' | xargs printf "%'.f" )
echo -e "NoLinker_UU_all\t"${unique} >> ${out_file}

# nlUUmap=$(zcat ${RUN}.$nonelabel.UU.statCB.tsv.gz | awk '{sum+=$2}END{print sum}' | xargs printf "%'.f")
# echo -e "NoLinker_UU_mapped\t"${nlUUmap} >> ${out_file}

nlUUuniq=$(zcat ${RUN}.$nonelabel.UU.hq.statCB.tsv.gz | awk '{sum+=$2}END{print sum}' | xargs printf "%'.f")
echo -e "NoLinker_UU_mapped_q${mapq}\t"${nlUUuniq} >> ${out_file}

nlUUuniqCB=$(zcat ${RUN}.$nonelabel.UU.hq.statCB.tsv.gz | grep -v noCB | awk '{sum+=$2}END{print sum}' | xargs printf "%'.f")
echo -e "NoLinker_UU_mapped_q${mapq}_validCB\t"${nlUUuniqCB} >> ${out_file}
nlUUuniqCB=$(zcat ${RUN}.$nonelabel.UU.hq.statCB.tsv.gz | grep -v noCB | awk '{sum+=$2}END{print sum}' )

nlUUuniqnr=$(zcat ${RUN}.$nonelabel.UU.hq.CB.nr.statCB.tsv.gz | awk '{sum+=$2}END{print sum}' | xargs printf "%'.f")
echo -e "NoLinker_UU_mapped_q${mapq}_validCB_dedup\t"${nlUUuniqnr} >> ${out_file}
nlUUuniqnr=$(zcat ${RUN}.$nonelabel.UU.hq.CB.nr.statCB.tsv.gz | awk '{sum+=$2}END{print sum}' )

uunlDupRate=$(grep "Unknown Library"  ${RUN}.${nonelabel}.UU.hq.CB.nr.MarkDuplicates.txt | awk -F '[\t]' '{print $9*100}' |  xargs printf "%.1f%%" )
echo -e "NoLinker_UU_Redundancy\t"${uunlDupRate} >> ${out_file}

nlUUpet=$(zcat ${RUN}.$nonelabel.UU.hq.CB.nr.pet.statCB.tsv.gz | awk '{sum+=$2}END{print sum}' | xargs printf "%'.f")
echo -e "NoLinker_UU_PET\t"${nlUUpet} >> ${out_file}
nlUUpet=$(zcat ${RUN}.$nonelabel.UU.hq.CB.nr.pet.statCB.tsv.gz | awk '{sum+=$2}END{print sum}' )

NoLinkertrans=$(head -2 ${RUN}.full.CBsums.tsv | tail -1 | cut -f 79 | xargs printf "%'.f" )
NoLinkercis=$(head -2 ${RUN}.full.CBsums.tsv | tail -1 | cut -f 80 | xargs printf "%'.f" )
NoLinkercis1k=$(head -2 ${RUN}.full.CBsums.tsv | tail -1 | cut -f 81 | xargs printf "%'.f" )
NoLinkercis8k=$(head -2 ${RUN}.full.CBsums.tsv | tail -1 | cut -f 82 | xargs printf "%'.f" )
echo -e "NoLinker_trans\t"${NoLinkertrans} >> ${out_file}
echo -e "NoLinker_cis\t"${NoLinkercis} >> ${out_file}
echo -e "NoLinker_cis>1kb\t"${NoLinkercis1k} >> ${out_file}
echo -e "NoLinker_cis>8kb\t"${NoLinkercis8k} >> ${out_file}
NoLinkertrans=$(head -2 ${RUN}.full.CBsums.tsv | tail -1 | cut -f 79 )
NoLinkercis=$(head -2 ${RUN}.full.CBsums.tsv | tail -1 | cut -f 80 )
NoLinkercis1k=$(head -2 ${RUN}.full.CBsums.tsv | tail -1 | cut -f 81 )
NoLinkercis8k=$(head -2 ${RUN}.full.CBsums.tsv | tail -1 | cut -f 82 )
ratio=$( echo "${NoLinkercis} / ${NoLinkertrans}" | bc -l )
ratio=$( printf %.2f ${ratio} )
echo -e "NoLinker_ratio_cis/trans\t"${ratio} >> ${out_file}
ratio=$( echo "${NoLinkercis1k} / ${NoLinkertrans}" | bc -l )
ratio=$( printf %.2f ${ratio} )
echo -e "NoLinker_ratio_cis>1kb/trans\t"${ratio} >> ${out_file}

uxxu=$( cat ${RUN}.$nonelabel.stat.xls | grep "Ux/xU" | awk -F '[\t]' '{print $3}' | xargs printf "%'.f" )
echo -e "NoLinker_UxxU_all\t"${uxxu} >> ${out_file}

# nlUxxUuniq=$(zcat ${RUN}.$nonelabel.UxxU.statCB.tsv.gz  | awk '{sum+=$2}END{print sum}' | xargs printf "%'.f")
# echo -e "NoLinker_UxxU_mapped\t"${nlUxxUuniq} >> ${out_file}

nlUxxUuniq=$(zcat ${RUN}.$nonelabel.UxxU.hq.statCB.tsv.gz  | awk '{sum+=$2}END{print sum}' | xargs printf "%'.f")
echo -e "NoLinker_UxxU_mapped_q${mapq}\t"${nlUxxUuniq} >> ${out_file}

nlUxxUuniqCB=$(zcat ${RUN}.$nonelabel.UxxU.hq.statCB.tsv.gz  | grep -v noCB | awk '{sum+=$2}END{print sum}' | xargs printf "%'.f")
echo -e "NoLinker_UxxU_mapped_q${mapq}_validCB\t"${nlUxxUuniqCB} >> ${out_file}
nlUxxUuniqCB=$(zcat ${RUN}.$nonelabel.UxxU.hq.statCB.tsv.gz  | grep -v noCB | awk '{sum+=$2}END{print sum}' )

nlUxxUuniqnr=$(zcat ${RUN}.$nonelabel.UxxU.hq.CB.nr.statCB.tsv.gz  | awk '{sum+=$2}END{print sum}' | xargs printf "%'.f")
echo -e "NoLinker_UxxU_mapped_q${mapq}_validCB_dedup\t"${nlUxxUuniqnr} >> ${out_file}
nlUxxUuniqnr=$(zcat ${RUN}.$nonelabel.UxxU.hq.CB.nr.statCB.tsv.gz  | awk '{sum+=$2}END{print sum}' )

uxxunlDupRate=$(grep "Unknown Library"  ${RUN}.${nonelabel}.UxxU.hq.CB.nr.MarkDuplicates.txt | awk -F '[\t]' '{print $9*100}' |  xargs printf "%.1f%%" )
echo -e "NoLinker_UxxU_Redundancy\t"${uxxunlDupRate} >> ${out_file} 

######################################################################### dup stat
echo "" >> ${out_file}
redun=$( echo "(1 - (${pairUUuniqnr}+${pairUxxUuniqnr}+${singleUxxUuniqnr}+${nlUUuniqnr}+${nlUxxUuniqnr}) / (${pairUUuniqCB}+${pairUxxUuniqCB}+${singleUxxUuniqCB}+${nlUUuniqCB}+${nlUxxUuniqCB})) * 100" | bc -l )
redun=$( printf %.2f ${redun} )
echo -e "Average_Redundancy\t"${redun}"%" >> ${out_file}

trans=$( echo "${PETtrans} + ${NoLinkertrans}" | bc )
cis=$( echo "${PETcis} + ${NoLinkercis}" | bc )
cis1k=$( echo "${PETcis1k} + ${NoLinkercis1k}" | bc )
cis8k=$( echo "${PETcis8k} + ${NoLinkercis8k}" | bc )
allpet=$( echo "${trans} + ${cis}" | bc )
ratio=$( echo "${cis} / ${trans}" | bc -l )
ratio1=$( printf %.2f ${ratio} )
ratio=$( echo "${cis1k} / ${trans}" | bc -l )
ratio2=$( printf %.2f ${ratio} )
allpet=$( printf "%'.f" ${allpet} )
trans=$( printf "%'.f" ${trans} )
cis=$( printf "%'.f" ${cis} )
cis1k=$( printf "%'.f" ${cis1k} )
cis8k=$( printf "%'.f" ${cis8k} )
echo -e "all_pet\t"${allpet} >> ${out_file}
echo -e "all_trans\t"${trans} >> ${out_file}
echo -e "all_cis\t"${cis} >> ${out_file}
echo -e "all_cis>1kb\t"${cis1k} >> ${out_file}
echo -e "all_cis>8kb\t"${cis8k} >> ${out_file}
echo -e "all_ratio_cis/trans\t"${ratio1} >> ${out_file}
echo -e "all_ratio_cis>1kb/trans\t"${ratio2} >> ${out_file}

######################################################################### cluster
echo "" >> ${out_file}
total_cluster_number=$(zcat ${RUN}*chiasig.gz | awk '$7 != 1{print}' | wc -l)
total_cluster_number=$( printf "%'.f" ${total_cluster_number} )
echo -e "PET_cluster\t"${total_cluster_number} >> ${out_file}
intra_cluster=$( zcat ${RUN}*cis.chiasig.gz | awk '$7 >=2 {print}' | wc -l )
inter_cluster=$( zcat ${RUN}*trans.chiasig.gz | awk '$7 >=2 {print}' | wc -l)
intra_cluster=$( printf "%'.f" ${intra_cluster} )
echo -e "Intra-chr_PET_cluster\t"${intra_cluster} >> ${out_file}
for i in $(seq 2 10)
do
    intra_pets_number=$(zcat ${RUN}*cis.chiasig.gz | awk -v cutoff=${i} '$7 == cutoff {print}' | wc -l | xargs printf "%'.f")
    echo -e "pets_number_"${i}"\t"${intra_pets_number} >> ${out_file}
done
echo -e "pets_number>10\t"$(zcat ${RUN}*cis.chiasig.gz | awk '$7 >10 {print}' | wc -l | xargs printf "%'.f") >> ${out_file}

inter_cluster=$( printf "%'.f" ${inter_cluster} )
echo -e "Inter-chr_PET_cluster\t"${inter_cluster} >> ${out_file}
for i in $(seq 2 10)
do
    inter_pets_number=$(zcat ${RUN}*trans.chiasig.gz | awk -v cutoff=${i} '$7 == cutoff {print}' | wc -l | xargs printf "%'.f")
    echo -e "pets_number_"${i}"\t"${inter_pets_number} >> ${out_file}
done
echo -e "pets_number>10\t"$(zcat ${RUN}*trans.chiasig.gz | awk '$7 >10 {print}' | wc -l | xargs printf "%'.f") >> ${out_file}

ns=$( grep -v : ${RUN}.e${extbp}.clusters.cis.chiasig.BE2.sigf.interactions | wc -l )
ns=$( printf "%'.f" ${ns} )
echo -e "Significant_interaction \t"$ns >> ${out_file}

######################################################################### Get scRNA stat
echo "" >> ${out_file}
$Rscript $getfrom10xG ../10xGEX/outs/metrics_summary.csv >> ${out_file} 2>&1

######################################################################### Get scATAC stat
echo "" >> ${out_file}
np=$( cat ${RUN}_peaks.narrowPeak | wc -l )
n_peak=$( printf "%'.f" ${np} )
echo -e "ATAC_Peak\t"$n_peak >> ${out_file}

nf=$( zcat atac_fragments.tsv.gz | wc -l )
n_fragment=$( printf "%'.f" ${nf} )
echo -e "ATAC_Fragment\t"$n_fragment >> ${out_file}

nb=$( cat ../10xATAC/outs/raw_peak_bc_matrix/barcodes.tsv | wc -l )
n_barcode=$( printf "%'.f" ${nb} )
echo -e "ATAC_Barcode\t"$n_barcode >> ${out_file}

python $getfrom10xA ../10xATAC/outs >> ${out_file}
