# zcat ~/20240614/MOSCCA0015-15878_GTAACATGCG-AGGTAACACT_S1_L001_R1_001.fastq.gz | head -12000000 | gzip > testinput/MOSCCA0015-15878_GTAACATGCG-AGGTAACACT_S1_L001_R1_001.fastq.gz
# zcat ~/20240614/MOSCCA0015-15878_GTAACATGCG-AGGTAACACT_S1_L001_R2_001.fastq.gz | head -12000000 | gzip > testinput/MOSCCA0015-15878_GTAACATGCG-AGGTAACACT_S1_L001_R2_001.fastq.gz
# zcat ~/20240614/MOSCCA0020-15872_SI-NA-D6_S5_L003_R3_001.fastq.gz | head -12000000 | gzip > testinput/MOSCCA0020-15872_SI-NA-D6_S5_L003_R3_001.fastq.gz
# zcat ~/20240614/MOSCCA0020-15872_SI-NA-D6_S5_L003_R2_001.fastq.gz | head -12000000 | gzip > testinput/MOSCCA0020-15872_SI-NA-D6_S5_L003_R2_001.fastq.gz
# zcat ~/20240614/MOSCCA0020-15872_SI-NA-D6_S5_L003_R1_001.fastq.gz | head -12000000 | gzip > testinput/MOSCCA0020-15872_SI-NA-D6_S5_L003_R1_001.fastq.gz



sbatch -p ruan_cpu2 -J sh.run -o sh.run.out -e sh.run.err -N 1 -n 10 --mem=100G --wrap=" bash /data/home/ruanlab/huangxingyu/scChAIR/testdir/mainPro.sh config/config10Xcustomslurm.txt input/ChAIRdataslurm.input.txt "
bash /data/home/ruanlab/huangxingyu/scChAIR/testdir/mainPro.sh config/config10Xcustomlocal.txt input/ChAIRdatalocal.input.txt > ChAIRdatalocal.log 2>&1
