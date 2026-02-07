# zcat ~/20240614/MOSCCA0015-15878_GTAACATGCG-AGGTAACACT_S1_L001_R1_001.fastq.gz | head -12000000 | gzip > testinput/MOSCCA0015-15878_GTAACATGCG-AGGTAACACT_S1_L001_R1_001.fastq.gz
# zcat ~/20240614/MOSCCA0015-15878_GTAACATGCG-AGGTAACACT_S1_L001_R2_001.fastq.gz | head -12000000 | gzip > testinput/MOSCCA0015-15878_GTAACATGCG-AGGTAACACT_S1_L001_R2_001.fastq.gz
# zcat ~/20240614/MOSCCA0020-15872_SI-NA-D6_S5_L003_R3_001.fastq.gz | head -12000000 | gzip > testinput/MOSCCA0020-15872_SI-NA-D6_S5_L003_R3_001.fastq.gz
# zcat ~/20240614/MOSCCA0020-15872_SI-NA-D6_S5_L003_R2_001.fastq.gz | head -12000000 | gzip > testinput/MOSCCA0020-15872_SI-NA-D6_S5_L003_R2_001.fastq.gz
# zcat ~/20240614/MOSCCA0020-15872_SI-NA-D6_S5_L003_R1_001.fastq.gz | head -12000000 | gzip > testinput/MOSCCA0020-15872_SI-NA-D6_S5_L003_R1_001.fastq.gz


# cp /data/home/ruanlab/huangxingyu/scChAIR/dropletHiC/GSE253407/MC_Paired_Hi-C_rep3_DNAMusmusculus/tmp/MC-RNA-5M_S1_L001_R1_001.fastq.gz testinput/MC-RNA-5M_S1_L001_R1_001.fastq.gz
# cp /data/home/ruanlab/huangxingyu/scChAIR/dropletHiC/GSE253407/MC_Paired_Hi-C_rep3_DNAMusmusculus/tmp/MC-RNA-5M_S1_L001_R2_001.fastq.gz testinput/MC-RNA-5M_S1_L001_R2_001.fastq.gz
# cp /data/home/ruanlab/huangxingyu/scChAIR/dropletHiC/GSE253407/MC_Paired_Hi-C_rep3_DNAMusmusculus/tmp/MC-DNA-5M_S1_L001_R2_001.fastq.gz testinput/MC-DNA-5M_S1_L001_R2_001.fastq.gz
# cp /data/home/ruanlab/huangxingyu/scChAIR/dropletHiC/GSE253407/MC_Paired_Hi-C_rep3_DNAMusmusculus/tmp/MC-DNA-5M_S1_L001_R1_001.fastq.gz testinput/MC-DNA-5M_S1_L001_R1_001.fastq.gz
# cp /data/home/ruanlab/huangxingyu/scChAIR/dropletHiC/GSE253407/MC_Paired_Hi-C_rep3_DNAMusmusculus/tmp/MC-DNA-5M_S1_L001_R3_001.fastq.gz testinput/MC-DNA-5M_S1_L001_R3_001.fastq.gz
# zcat testinput/MC-DNA-5M_S1_L001_R2_001.fastq.gz | awk 'NR%4==2' | head -1000000 | cut -c 1-16 | sort | uniq | wc -l # cut -c 9-24 | tr ACGT TGCA
# for i in 737K-arc-v1.txt.gz 737K-cratac-v1.txt.gz
# do 
#     echo ${i}
#     cat <(zcat testinput/MC-DNA-5M_S1_L001_R2_001.fastq.gz | awk 'NR%4==2' | head -1000000 | cut -c 1-16 | sort | uniq) <(zcat /data/home/ruanlab/huangxingyu/Tools/cellranger-atac-2.1.0raw/lib/python/atac/barcodes/${i}) | sort | uniq -d | wc -l
#     cat <(zcat testinput/MC-DNA-5M_S1_L001_R2_001.fastq.gz | awk 'NR%4==2' | head -1000000 | cut -c 1-16 | tr ACGT TGCA | rev | sort | uniq) <(zcat /data/home/ruanlab/huangxingyu/Tools/cellranger-atac-2.1.0raw/lib/python/atac/barcodes/${i}) | sort | uniq -d | wc -l
# done
# for i in 737K-arc-v1.txt.gz 737K-cratac-v1.txt.gz
# do
#     echo ${i}
#     cat <(zcat testinput/MC-DNA-5M_S1_L001_R2_001.fastq.gz | awk 'NR%4==2' | head -1000000 | cut -c 9-24 | sort | uniq) <(zcat /data/home/ruanlab/huangxingyu/Tools/cellranger-atac-2.1.0raw/lib/python/atac/barcodes/${i}) | sort | uniq -d | wc -l
#     cat <(zcat testinput/MC-DNA-5M_S1_L001_R2_001.fastq.gz | awk 'NR%4==2' | head -1000000 | cut -c 9-24 | tr ACGT TGCA | rev | sort | uniq) <(zcat /data/home/ruanlab/huangxingyu/Tools/cellranger-atac-2.1.0raw/lib/python/atac/barcodes/${i}) | sort | uniq -d | wc -l
# done



# sbatch -p ruan_cpu2 -J sh.run -o sh.run.out -e sh.run.err -N 1 -n 10 --mem=100G --wrap=" bash /data/home/ruanlab/huangxingyu/scChAIR/testdir/mainPro.sh config/ChAIRconfig10Xcustomslurm.txt input/ChAIRdataslurm.input.txt "
# bash /data/home/ruanlab/huangxingyu/scChAIR/testdir/mainPro.sh config/ChAIRconfig10Xcustomlocal.txt input/ChAIRdatalocal.input.txt > ChAIRdatalocal.log 2>&1

# sbatch -p ruan_cpu2 -J PairedSeq.run -o PairedSeq.run.out -e PairedSeq.run.err -N 1 -n 10 --mem=100G --wrap=" bash /data/home/ruanlab/huangxingyu/scChAIR/testdir/mainPro.sh config/PairedSeqconfig10Xcustomslurm.txt input/PairedSeqdataslurm.input.txt "
# bash /data/home/ruanlab/huangxingyu/scChAIR/testdir/mainPro.sh config/PairedSeqconfig10Xcustomlocal.txt input/PairedSeqdatalocal.input.txt > PairedSeqdatalocal.log 2>&1




# cp /data/home/ruanlab/huangxingyu/scChAIR/ISSAAC-seq/K56210xMultiome/tmp/K562_10x_Multiome_RNA-5M_S1_L001_R1_001.fastq.gz /data/home/ruanlab/huangxingyu/scChAIR/ISSAAC-seq/K56210xMultiome/tmp/K562_10x_Multiome_RNA-5M_S1_L001_R2_001.fastq.gz /data/home/ruanlab/huangxingyu/scChAIR/ISSAAC-seq/K56210xMultiome/tmp/K562_10x_Multiome_ATAC-5M_S1_L001_I2_001.fastq.gz /data/home/ruanlab/huangxingyu/scChAIR/ISSAAC-seq/K56210xMultiome/tmp/K562_10x_Multiome_ATAC-5M_S1_L001_R1_001.fastq.gz /data/home/ruanlab/huangxingyu/scChAIR/ISSAAC-seq/K56210xMultiome/tmp/K562_10x_Multiome_ATAC-5M_S1_L001_R2_001.fastq.gz testinput
# sbatch -p ruan_cpu3 -J RNAATAC.run -o RNAATAC.run.out -e RNAATAC.run.err -N 1 -n 10 --mem=100G --wrap=" bash /data/home/ruanlab/huangxingyu/scChAIR/testdir/mainPro.sh config/RNAATACconfig10Xcustomslurm.txt input/RNAATACdataslurm.input.txt "
# bash /data/home/ruanlab/huangxingyu/scChAIR/testdir/mainPro.sh config/RNAATACconfig10Xcustomlocal.txt input/RNAATACdatalocal.input.txt > RNAATACdatalocal.log 2>&1



# zcat /data/home/ruanlab/huangxingyu/Tools/cellranger-atac-2.1.0raw/lib/python/atac/barcodes/737K-cratac-v1.txt.gz > barcodefile/DNA.737K-cratac-v1.txt
# cp /data/home/ruanlab/huangxingyu/scChAIR/dropletHiC/GSE253407/MC_Droplet_Hi-C_10xATAC_v1.1Musmusculus/tmp/MC-5M_S1_L001_R2_001.fastq.gz /data/home/ruanlab/huangxingyu/scChAIR/dropletHiC/GSE253407/MC_Droplet_Hi-C_10xATAC_v1.1Musmusculus/tmp/MC-5M_S1_L001_R1_001.fastq.gz /data/home/ruanlab/huangxingyu/scChAIR/dropletHiC/GSE253407/MC_Droplet_Hi-C_10xATAC_v1.1Musmusculus/tmp/MC-5M_S1_L001_R3_001.fastq.gz testinput
# sbatch -p ruan_cpu2 -J HiC.run -o HiC.run.out -e HiC.run.err -N 1 -n 10 --mem=100G --wrap=" bash /data/home/ruanlab/huangxingyu/scChAIR/testdir/mainPro.sh config/HiCconfig10Xcustomslurm.txt input/HiCdataslurm.input.txt "
# bash /data/home/ruanlab/huangxingyu/scChAIR/testdir/mainPro.sh config/HiCconfig10Xcustomlocal.txt input/HiCdatalocal.input.txt > HiCdatalocal.log 2>&1

# cp /data/home/ruanlab/huangxingyu/scChAIR/publicscATAC/8Kmousecortex/8k_mouse_cortex_ATACv2_nextgem_Chromium_X_fastqs/tmp/8kmousecortex-5M_S1_L001_R2_001.fastq.gz /data/home/ruanlab/huangxingyu/scChAIR/publicscATAC/8Kmousecortex/8k_mouse_cortex_ATACv2_nextgem_Chromium_X_fastqs/tmp/8kmousecortex-5M_S1_L001_R1_001.fastq.gz /data/home/ruanlab/huangxingyu/scChAIR/publicscATAC/8Kmousecortex/8k_mouse_cortex_ATACv2_nextgem_Chromium_X_fastqs/tmp/8kmousecortex-5M_S1_L001_R3_001.fastq.gz testinput
# sbatch -p ruan_cpu2 -J ATAC.run -o ATAC.run.out -e ATAC.run.err -N 1 -n 10 --mem=100G --wrap=" bash /data/home/ruanlab/huangxingyu/scChAIR/testdir/mainPro.sh config/ATACconfig10Xcustomslurm.txt input/ATACdataslurm.input.txt "
# bash /data/home/ruanlab/huangxingyu/scChAIR/testdir/mainPro.sh config/ATACconfig10Xcustomlocal.txt input/ATACdatalocal.input.txt > ATACdatalocal.log 2>&1





# cp /data/home/ruanlab/huangxingyu/scChAIR/publicscRNA/PRJNA515497/tmp/SRR8448140-5M_S1_L001_R1_001.fastq.gz /data/home/ruanlab/huangxingyu/scChAIR/publicscRNA/PRJNA515497/tmp/SRR8448140-5M_S1_L001_R2_001.fastq.gz testinput
# sbatch -p ruan_cpu2 -J RNA.run -o RNA.run.out -e RNA.run.err -N 1 -n 10 --mem=100G --wrap=" bash /data/home/ruanlab/huangxingyu/scChAIR/testdir/mainPro.sh config/RNAconfig10Xcustomslurm.txt input/RNAdataslurm.input.txt "
# bash /data/home/ruanlab/huangxingyu/scChAIR/testdir/mainPro.sh config/RNAconfig10Xcustomlocal.txt input/RNAdatalocal.input.txt > RNAdatalocal.log 2>&1






zcat e18_mouse_brain_DNA_S1_L001_R2_001.fastq.gz | awk 'NR%4==2' | head -1000000 | cut -c 1-16 | sort | uniq | wc -l
for i in 737K-arc-v1.txt.gz 737K-cratac-v1.txt.gz
do 
    echo ${i}
    cat <(zcat e18_mouse_brain_DNA_S1_L001_R2_001.fastq.gz | awk 'NR%4==2' | head -1000000 | cut -c 1-16 | sort | uniq) <(zcat /data/home/ruanlab/huangxingyu/Tools/cellranger-atac-2.1.0raw/lib/python/atac/barcodes/${i}) | sort | uniq -d | wc -l
    cat <(zcat e18_mouse_brain_DNA_S1_L001_R2_001.fastq.gz | awk 'NR%4==2' | head -1000000 | cut -c 1-16 | tr ACGT TGCA | rev | sort | uniq) <(zcat /data/home/ruanlab/huangxingyu/Tools/cellranger-atac-2.1.0raw/lib/python/atac/barcodes/${i}) | sort | uniq -d | wc -l
done
for i in 737K-arc-v1.txt.gz 737K-cratac-v1.txt.gz
do
    echo ${i}
    cat <(zcat e18_mouse_brain_DNA_S1_L001_R2_001.fastq.gz | awk 'NR%4==2' | head -1000000 | cut -c 9-24 | sort | uniq) <(zcat /data/home/ruanlab/huangxingyu/Tools/cellranger-atac-2.1.0raw/lib/python/atac/barcodes/${i}) | sort | uniq -d | wc -l
    cat <(zcat e18_mouse_brain_DNA_S1_L001_R2_001.fastq.gz | awk 'NR%4==2' | head -1000000 | cut -c 9-24 | tr ACGT TGCA | rev | sort | uniq) <(zcat /data/home/ruanlab/huangxingyu/Tools/cellranger-atac-2.1.0raw/lib/python/atac/barcodes/${i}) | sort | uniq -d | wc -l
done

