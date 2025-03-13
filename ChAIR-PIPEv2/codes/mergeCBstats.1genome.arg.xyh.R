args <- commandArgs(trailingOnly = TRUE);
# traceback() 
RUN <- args[1]; 

message('Combining stats for: ', RUN)
singlabel='.singlelinker.single'
pairlabel='.singlelinker.paired'
nonelabel='.none'

readFile = function(infile) {
  if (length(grep(".gz$",infile)) || length(grep(".GZ$",infile)) || length(grep(".gzip$",infile)) || length(grep(".GZIP$",infile))) {
    gf = gzfile(infile,'rt')
    daf = read.delim(gf, check.names=FALSE, header=FALSE)
    close(gf)
    return(daf)
  } else {
    return(read.delim(infile, check.names=FALSE, header=FALSE))
  }
}

########################## total
none = readFile(paste0(RUN,nonelabel,'.statCB.tsv.gz')) # *.none.statCB.tsv.gz
stag = readFile(paste0(RUN,singlabel,'.statCB.tsv.gz')) # *.singlelinker.single.statCB.tsv.gz
tpet = readFile(paste0(RUN,pairlabel,'.statCB.tsv.gz')) # *.singlelinker.paired.statCB.tsv.gz
colnames(none) = c('CorrectBarcode', 'NoLinker_tot')
colnames(stag) = c('CorrectBarcode', '1tag_tot')
colnames(tpet) = c('CorrectBarcode', '2tag_tot')
alldata1=merge(merge(none, stag, all=TRUE), tpet, all=TRUE)
alldata1[is.na(alldata1)] <- 0
print(gc())

########################## split
# noneUU = readFile(paste0(RUN,nonelabel,'.UU.statCB.tsv.gz')) # *.none.statCB.tsv.gz
# noneUx = readFile(paste0(RUN,nonelabel,'.UxxU.statCB.tsv.gz')) # *.none.statCB.tsv.gz
# stagUx = readFile(paste0(RUN,singlabel,'.UxxU.statCB.tsv.gz')) # *.singlelinker.single.statCB.tsv.gz
# tpetUU = readFile(paste0(RUN,pairlabel,'.UU.statCB.tsv.gz')) # *.singlelinker.paired.statCB.tsv.gz
# tpetUx = readFile(paste0(RUN,pairlabel,'.UxxU.statCB.tsv.gz')) # *.singlelinker.paired.statCB.tsv.gz
# colnames(noneUU) = c('CorrectBarcode', 'NoLinker_UU')
# colnames(noneUx) = c('CorrectBarcode', 'NoLinker_UxxU')
# colnames(stagUx) = c('CorrectBarcode', '1tag_UxxU')
# colnames(tpetUU) = c('CorrectBarcode', '2tag_UU')
# colnames(tpetUx) = c('CorrectBarcode', '2tag_UxxU')
# alldata2=merge(merge(merge(merge(noneUU, noneUx, all=TRUE), stagUx, all=TRUE), tpetUU, all=TRUE), tpetUx, all=TRUE)
# alldata2[is.na(alldata2)] <- 0
# print(gc())

########################## high quality
noneUU_hq = readFile(paste0(RUN,nonelabel,'.UU.hq.statCB.tsv.gz')) # *.none.statCB.tsv.gz
noneUx_hq = readFile(paste0(RUN,nonelabel,'.UxxU.hq.statCB.tsv.gz')) # *.none.statCB.tsv.gz
stagUx_hq = readFile(paste0(RUN,singlabel,'.UxxU.hq.statCB.tsv.gz')) # *.singlelinker.single.statCB.tsv.gz
tpetUU_hq = readFile(paste0(RUN,pairlabel,'.UU.hq.statCB.tsv.gz')) # *.singlelinker.paired.statCB.tsv.gz
tpetUx_hq = readFile(paste0(RUN,pairlabel,'.UxxU.hq.statCB.tsv.gz')) # *.singlelinker.paired.statCB.tsv.gz
colnames(noneUU_hq) = c('CorrectBarcode', 'NoLinker_UU_hq')
colnames(noneUx_hq) = c('CorrectBarcode', 'NoLinker_UxxU_hq')
colnames(stagUx_hq) = c('CorrectBarcode', '1tag_UxxU_hq')
colnames(tpetUU_hq) = c('CorrectBarcode', '2tag_UU_hq')
colnames(tpetUx_hq) = c('CorrectBarcode', '2tag_UxxU_hq')
alldata3=merge(merge(merge(merge(noneUU_hq, noneUx_hq, all=TRUE), stagUx_hq, all=TRUE), tpetUU_hq, all=TRUE), tpetUx_hq, all=TRUE)
alldata3[is.na(alldata3)] <- 0
alldata3$NoLinker_hq = alldata3$NoLinker_UU_hq + alldata3$NoLinker_UxxU_hq
alldata3$`1tag_hq` = alldata3$`1tag_UxxU_hq`
alldata3$`2tag_hq` = alldata3$`2tag_UU_hq` + alldata3$`2tag_UxxU_hq`
alldata3 = alldata3[c("NoLinker_hq","1tag_hq","2tag_hq")]
atab = merge(alldata1, alldata3, all=TRUE)
rm(alldata1)
rm(alldata3)
gc()
print(gc())

########################## barcode
noneUU_CB = readFile(paste0(RUN,nonelabel,'.UU.hq.CB.statCB.tsv.gz')) # *.none.statCB.tsv.gz
noneUx_CB = readFile(paste0(RUN,nonelabel,'.UxxU.hq.CB.statCB.tsv.gz')) # *.none.statCB.tsv.gz
stagUx_CB = readFile(paste0(RUN,singlabel,'.UxxU.hq.CB.statCB.tsv.gz')) # *.singlelinker.single.statCB.tsv.gz
tpetUU_CB = readFile(paste0(RUN,pairlabel,'.UU.hq.CB.statCB.tsv.gz')) # *.singlelinker.paired.statCB.tsv.gz
tpetUx_CB = readFile(paste0(RUN,pairlabel,'.UxxU.hq.CB.statCB.tsv.gz')) # *.singlelinker.paired.statCB.tsv.gz
colnames(noneUU_CB) = c('CorrectBarcode', 'NoLinker_UU_CB')
colnames(noneUx_CB) = c('CorrectBarcode', 'NoLinker_UxxU_CB')
colnames(stagUx_CB) = c('CorrectBarcode', '1tag_UxxU_CB')
colnames(tpetUU_CB) = c('CorrectBarcode', '2tag_UU_CB')
colnames(tpetUx_CB) = c('CorrectBarcode', '2tag_UxxU_CB')
alldata4=merge(merge(merge(merge(noneUU_CB, noneUx_CB, all=TRUE), stagUx_CB, all=TRUE), tpetUU_CB, all=TRUE), tpetUx_CB, all=TRUE)
alldata4[is.na(alldata4)] <- 0
alldata4$NoLinker_CB = alldata4$NoLinker_UU_CB + alldata4$NoLinker_UxxU_CB
alldata4$`1tag_CB` = alldata4$`1tag_UxxU_CB`
alldata4$`2tag_CB` = alldata4$`2tag_UU_CB` + alldata4$`2tag_UxxU_CB`
alldata4 = alldata4[c("NoLinker_CB","1tag_CB","2tag_CB")]
atab = merge(atab, alldata4, all=TRUE)
rm(alldata4)
gc()
print(gc())

########################## nr
noneUU_nr = readFile(paste0(RUN,nonelabel,'.UU.hq.CB.nr.statCB.tsv.gz')) # *.none.statCB.tsv.gz
noneUx_nr = readFile(paste0(RUN,nonelabel,'.UxxU.hq.CB.nr.statCB.tsv.gz')) # *.none.statCB.tsv.gz
stagUx_nr = readFile(paste0(RUN,singlabel,'.UxxU.hq.CB.nr.statCB.tsv.gz')) # *.singlelinker.single.statCB.tsv.gz
tpetUU_nr = readFile(paste0(RUN,pairlabel,'.UU.hq.CB.nr.statCB.tsv.gz')) # *.singlelinker.paired.statCB.tsv.gz
tpetUx_nr = readFile(paste0(RUN,pairlabel,'.UxxU.hq.CB.nr.statCB.tsv.gz')) # *.singlelinker.paired.statCB.tsv.gz
colnames(noneUU_nr) = c('CorrectBarcode', 'NoLinker_UU_nr')
colnames(noneUx_nr) = c('CorrectBarcode', 'NoLinker_UxxU_nr')
colnames(stagUx_nr) = c('CorrectBarcode', '1tag_UxxU_nr')
colnames(tpetUU_nr) = c('CorrectBarcode', '2tag_UU_nr')
colnames(tpetUx_nr) = c('CorrectBarcode', '2tag_UxxU_nr')
alldata5=merge(merge(merge(merge(noneUU_nr, noneUx_nr, all=TRUE), stagUx_nr, all=TRUE), tpetUU_nr, all=TRUE), tpetUx_nr, all=TRUE)
alldata5[is.na(alldata5)] <- 0
alldata5$NoLinker_nr = alldata5$NoLinker_UU_nr + alldata5$NoLinker_UxxU_nr
alldata5$`1tag_nr` = alldata5$`1tag_UxxU_nr`
alldata5$`2tag_nr` = alldata5$`2tag_UU_nr` + alldata5$`2tag_UxxU_nr`
alldata5 = alldata5[c("NoLinker_nr","1tag_nr","2tag_nr")]
atab = merge(atab, alldata5, all=TRUE)
rm(alldata5)
gc()
print(gc())

########################## singleton other pet
noneUU_singleton = readFile(paste0(RUN,nonelabel,'.UU.hq.CB.nr.singleton.statCB.tsv.gz')) # *.none.statCB.tsv.gz
noneUU_other = readFile(paste0(RUN,nonelabel,'.UU.hq.CB.nr.other.statCB.tsv.gz')) # *.none.statCB.tsv.gz
noneUU_pet = readFile(paste0(RUN,nonelabel,'.UU.hq.CB.nr.pet.statCB.tsv.gz')) # *.none.statCB.tsv.gz
tpetUU_singleton = readFile(paste0(RUN,pairlabel,'.UU.hq.CB.nr.singleton.statCB.tsv.gz')) # *.none.statCB.tsv.gz
tpetUU_other = readFile(paste0(RUN,pairlabel,'.UU.hq.CB.nr.other.statCB.tsv.gz')) # *.none.statCB.tsv.gz
tpetUU_pet = readFile(paste0(RUN,pairlabel,'.UU.hq.CB.nr.pet.statCB.tsv.gz')) # *.none.statCB.tsv.gz

colnames(noneUU_singleton) = c('CorrectBarcode', 'NoLinker_UU_singleton')
colnames(noneUU_other) = c('CorrectBarcode', 'NoLinker_UU_other')
colnames(noneUU_pet) = c('CorrectBarcode', 'NoLinker_UU_pet')
colnames(tpetUU_singleton) = c('CorrectBarcode', '2tag_UU_singleton')
colnames(tpetUU_other) = c('CorrectBarcode', '2tag_UU_other')
colnames(tpetUU_pet) = c('CorrectBarcode', '2tag_UU_pet')
alldata6=merge(merge(merge(merge(merge(noneUU_singleton, noneUU_other, all=TRUE), noneUU_pet, all=TRUE), tpetUU_singleton, all=TRUE), tpetUU_other, all=TRUE), tpetUU_pet, all=TRUE)
alldata6[is.na(alldata6)] <- 0
atab = merge(atab, alldata6, all=TRUE)
rm(alldata6)
gc()
print(gc())

# CorrectBarcode	nPET	nTrans	nCis	>1kb	>8kb
readFile2 = function(infile) {
  if (length(grep(".gz$",infile)) || length(grep(".GZ$",infile)) || length(grep(".gzip$",infile)) || length(grep(".GZIP$",infile))) {
    gf = gzfile(infile,'rt')
    daf = read.delim(gf, check.names=FALSE, header=TRUE)
    close(gf)
    return(daf)
  } else {
    return(read.delim(infile, check.names=FALSE, header=TRUE))
  }
}

########################## all
noneUU_allpet=readFile2(paste0(RUN, nonelabel,'.UU.hq.CB.statPET.tsv.gz'))
tpetUU_allpet=readFile2(paste0(RUN, pairlabel,'.UU.hq.CB.statPET.tsv.gz'))
colnames(noneUU_allpet)[2:6] = paste0(colnames(noneUU_allpet)[2:6],'_NoLinker_all_pairs')
colnames(tpetUU_allpet)[2:6] = paste0(colnames(tpetUU_allpet)[2:6],'_2tag_all_pairs')
alldata7 = merge(noneUU_allpet, tpetUU_allpet, all=TRUE)
alldata7[is.na(alldata7)] <- 0
atab = merge(atab, alldata7, all=TRUE)
rm(alldata7)
gc()
print(gc())

########################## nr
noneUU_nrpet=readFile2(paste0(RUN, nonelabel,'.UU.hq.CB.nr.statPET.tsv.gz'))
tpetUU_nrpet=readFile2(paste0(RUN, pairlabel,'.UU.hq.CB.nr.statPET.tsv.gz'))
colnames(noneUU_nrpet)[2:6] = paste0(colnames(noneUU_nrpet)[2:6],'_NoLinker_nr_pairs')
colnames(tpetUU_nrpet)[2:6] = paste0(colnames(tpetUU_nrpet)[2:6],'_2tag_nr_pairs')
alldata8 = merge(noneUU_nrpet, tpetUU_nrpet, all=TRUE)
alldata8[is.na(alldata8)] <- 0
atab = merge(atab, alldata8, all=TRUE)
rm(alldata8)
gc()
print(gc())

########################## singleton other pet
noneUU_singletonpet = readFile2(paste0(RUN,nonelabel,'.UU.hq.CB.nr.singleton.statPET.tsv.gz')) # *.none.statCB.tsv.gz
noneUU_otherpet = readFile2(paste0(RUN,nonelabel,'.UU.hq.CB.nr.other.statPET.tsv.gz')) # *.none.statCB.tsv.gz
noneUU_petpet = readFile2(paste0(RUN,nonelabel,'.UU.hq.CB.nr.pet.statPET.tsv.gz')) # *.none.statCB.tsv.gz
tpetUU_singletonpet = readFile2(paste0(RUN,pairlabel,'.UU.hq.CB.nr.singleton.statPET.tsv.gz')) # *.none.statCB.tsv.gz
tpetUU_otherpet = readFile2(paste0(RUN,pairlabel,'.UU.hq.CB.nr.other.statPET.tsv.gz')) # *.none.statCB.tsv.gz
tpetUU_petpet = readFile2(paste0(RUN,pairlabel,'.UU.hq.CB.nr.pet.statPET.tsv.gz')) # *.none.statCB.tsv.gz
colnames(noneUU_singletonpet)[2:6] = paste0(colnames(noneUU_singletonpet)[2:6],'_NoLinker_singleton_pairs')
colnames(noneUU_otherpet)[2:6] = paste0(colnames(noneUU_otherpet)[2:6],'_NoLinker_other_pairs')
colnames(noneUU_petpet)[2:6] = paste0(colnames(noneUU_petpet)[2:6],'_NoLinker_pet_pairs')
colnames(tpetUU_singletonpet)[2:6] = paste0(colnames(tpetUU_singletonpet)[2:6],'_2tag_singleton_pairs')
colnames(tpetUU_otherpet)[2:6] = paste0(colnames(tpetUU_otherpet)[2:6],'_2tag_other_pairs')
colnames(tpetUU_petpet)[2:6] = paste0(colnames(tpetUU_petpet)[2:6],'_2tag_pet_pairs')
alldata9=merge(merge(merge(merge(merge(noneUU_singletonpet, noneUU_otherpet, all=TRUE), noneUU_petpet, all=TRUE), tpetUU_singletonpet, all=TRUE), tpetUU_otherpet, all=TRUE), tpetUU_petpet, all=TRUE)
alldata9[is.na(alldata9)] <- 0
atab = merge(atab, alldata9, all=TRUE)
rm(alldata9)
gc()
print(gc())

# atab = merge(merge(merge(merge(merge(merge(merge(merge(alldata1, alldata2, all=TRUE), alldata3, all=TRUE), alldata4, all=TRUE), alldata5, all=TRUE), alldata6, all=TRUE), alldata7, all=TRUE), alldata8, all=TRUE), alldata9, all=TRUE)
# atab = merge(merge(merge(merge(merge(alldata1, alldata5, all=TRUE), alldata6, all=TRUE), alldata7, all=TRUE), alldata8, all=TRUE), alldata9, all=TRUE)
# atab = merge(merge(merge(alldata1, alldata4, all=TRUE), alldata5, all=TRUE), alldata6, all=TRUE)
write.table(atab, file=paste0(RUN, '.CBcounts.tsv'), quote=F, row.names=F, col.names=T, sep="\t")

#sums
nctab = ncol(atab)
tot = colSums(atab[,2:nctab])
totnocb = colSums(subset(atab, CorrectBarcode!='noCB')[,2:nctab])
write.table(rbind(tot, totnocb), file=paste0(RUN, '.CBsums.tsv'), quote=F, row.names=F, col.names=T, sep="\t")
