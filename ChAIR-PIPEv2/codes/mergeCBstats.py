import sys
import pandas as pd

def readFile(infile, header):
    return pd.read_csv(infile, sep='\t', names=header, header=None,index_col=0)

def readFile2(infile):
    return pd.read_csv(infile, sep='\t', index_col=0)

import gzip        
def writeFile(outfile):
    if outfile.endswith((".gz","gzip")):
        fout = gzip.open(outfile,'wt')
    else:
        fout = open(outfile,'w')
    return fout

RUN = sys.argv[1]
print('Combining stats for: {0}'.format(RUN))
singlabel='.singlelinker.single'
pairlabel='.singlelinker.paired'
nonelabel='.none'


########################## total
none = readFile("{0}{1}.statCB.tsv.gz".format(RUN,nonelabel), ['CorrectBarcode', 'NoLinker_tot'])
stag = readFile("{0}{1}.statCB.tsv.gz".format(RUN,singlabel), ['CorrectBarcode', '1tag_tot'])
tpet = readFile("{0}{1}.statCB.tsv.gz".format(RUN,pairlabel), ['CorrectBarcode', '2tag_tot'])
alldata1 = pd.concat([none, stag, tpet], axis=1)
alldata1 = alldata1.fillna(0).astype(int)

########################## split
noneUU = readFile("{0}{1}.UU.statCB.tsv.gz".format(RUN,nonelabel),['CorrectBarcode', 'NoLinker_UU']) 
noneUx = readFile("{0}{1}.UxxU.statCB.tsv.gz".format(RUN,nonelabel),['CorrectBarcode', 'NoLinker_UxxU']) 
stagUx = readFile("{0}{1}.UxxU.statCB.tsv.gz".format(RUN,singlabel),['CorrectBarcode', '1tag_UxxU']) 
tpetUU = readFile("{0}{1}.UU.statCB.tsv.gz".format(RUN,pairlabel),['CorrectBarcode', '2tag_UU']) 
tpetUx = readFile("{0}{1}.UxxU.statCB.tsv.gz".format(RUN,pairlabel),['CorrectBarcode', '2tag_UxxU']) 
alldata2 = pd.concat([noneUU, noneUx, stagUx, tpetUU, tpetUx], axis=1)
alldata2 = alldata2.fillna(0).astype(int)
alldata2['NoLinker_mapped'] = alldata2['NoLinker_UU'] + alldata2['NoLinker_UxxU']
alldata2['1tag_mapped'] = alldata2['1tag_UxxU'] 
alldata2['2tag_mapped'] = alldata2['2tag_UU'] + alldata2['2tag_UxxU']
alldata22 = alldata2[["NoLinker_mapped","1tag_mapped","2tag_mapped"]]

########################## high quality
noneUU_hq = readFile("{0}{1}.UU.hq.statCB.tsv.gz".format(RUN,nonelabel), ['CorrectBarcode', 'NoLinker_UU_hq']) 
noneUx_hq = readFile("{0}{1}.UxxU.hq.statCB.tsv.gz".format(RUN,nonelabel), ['CorrectBarcode', 'NoLinker_UxxU_hq']) 
stagUx_hq = readFile("{0}{1}.UxxU.hq.statCB.tsv.gz".format(RUN,singlabel), ['CorrectBarcode', '1tag_UxxU_hq']) 
tpetUU_hq = readFile("{0}{1}.UU.hq.statCB.tsv.gz".format(RUN,pairlabel), ['CorrectBarcode', '2tag_UU_hq']) 
tpetUx_hq = readFile("{0}{1}.UxxU.hq.statCB.tsv.gz".format(RUN,pairlabel), ['CorrectBarcode', '2tag_UxxU_hq']) 
alldata3 = pd.concat([noneUU_hq, noneUx_hq, stagUx_hq, tpetUU_hq, tpetUx_hq], axis=1)
alldata3 = alldata3.fillna(0).astype(int)
alldata3['NoLinker_hq'] = alldata3['NoLinker_UU_hq'] + alldata3['NoLinker_UxxU_hq']
alldata3['1tag_hq'] = alldata3['1tag_UxxU_hq'] 
alldata3['2tag_hq'] = alldata3['2tag_UU_hq'] + alldata3['2tag_UxxU_hq']
alldata33 = alldata3[["NoLinker_hq","1tag_hq","2tag_hq"]]

########################## barcode
noneUU_CB = readFile("{0}{1}.UU.hq.CB.statCB.tsv.gz".format(RUN,nonelabel), ['CorrectBarcode', 'NoLinker_UU_CB']) 
noneUx_CB = readFile("{0}{1}.UxxU.hq.CB.statCB.tsv.gz".format(RUN,nonelabel), ['CorrectBarcode', 'NoLinker_UxxU_CB']) 
stagUx_CB = readFile("{0}{1}.UxxU.hq.CB.statCB.tsv.gz".format(RUN,singlabel), ['CorrectBarcode', '1tag_UxxU_CB']) 
tpetUU_CB = readFile("{0}{1}.UU.hq.CB.statCB.tsv.gz".format(RUN,pairlabel), ['CorrectBarcode', '2tag_UU_CB']) 
tpetUx_CB = readFile("{0}{1}.UxxU.hq.CB.statCB.tsv.gz".format(RUN,pairlabel), ['CorrectBarcode', '2tag_UxxU_CB']) 
alldata4 = pd.concat([noneUU_CB, noneUx_CB, stagUx_CB, tpetUU_CB, tpetUx_CB], axis=1)
alldata4 = alldata4.fillna(0).astype(int)
alldata4['NoLinker_CB'] = alldata4['NoLinker_UU_CB'] + alldata4['NoLinker_UxxU_CB']
alldata4['1tag_CB'] = alldata4['1tag_UxxU_CB'] 
alldata4['2tag_CB'] = alldata4['2tag_UU_CB'] + alldata4['2tag_UxxU_CB']
alldata44 = alldata4[["NoLinker_CB","1tag_CB","2tag_CB"]]

########################## nr
noneUU_nr = readFile("{0}{1}.UU.hq.CB.nr.statCB.tsv.gz".format(RUN,nonelabel), ['CorrectBarcode', 'NoLinker_UU_nr']) 
noneUx_nr = readFile("{0}{1}.UxxU.hq.CB.nr.statCB.tsv.gz".format(RUN,nonelabel), ['CorrectBarcode', 'NoLinker_UxxU_nr']) 
stagUx_nr = readFile("{0}{1}.UxxU.hq.CB.nr.statCB.tsv.gz".format(RUN,singlabel), ['CorrectBarcode', '1tag_UxxU_nr']) 
tpetUU_nr = readFile("{0}{1}.UU.hq.CB.nr.statCB.tsv.gz".format(RUN,pairlabel), ['CorrectBarcode', '2tag_UU_nr']) 
tpetUx_nr = readFile("{0}{1}.UxxU.hq.CB.nr.statCB.tsv.gz".format(RUN,pairlabel), ['CorrectBarcode', '2tag_UxxU_nr']) 
alldata5 = pd.concat([noneUU_nr, noneUx_nr, stagUx_nr, tpetUU_nr, tpetUx_nr], axis=1)
alldata5 = alldata5.fillna(0).astype(int)
alldata5['NoLinker_nr'] = alldata5['NoLinker_UU_nr'] + alldata5['NoLinker_UxxU_nr']
alldata5['1tag_nr'] = alldata5['1tag_UxxU_nr'] 
alldata5['2tag_nr'] = alldata5['2tag_UU_nr'] + alldata5['2tag_UxxU_nr']
alldata55 = alldata5[["NoLinker_nr","1tag_nr","2tag_nr"]]

########################## singleton other pet
noneUU_singleton = readFile("{0}{1}.UU.hq.CB.nr.singleton.statCB.tsv.gz".format(RUN,nonelabel), ['CorrectBarcode', 'NoLinker_UU_singleton']) 
noneUU_other = readFile("{0}{1}.UU.hq.CB.nr.other.statCB.tsv.gz".format(RUN,nonelabel), ['CorrectBarcode', 'NoLinker_UU_other']) 
noneUU_pet = readFile("{0}{1}.UU.hq.CB.nr.pet.statCB.tsv.gz".format(RUN,nonelabel), ['CorrectBarcode', 'NoLinker_UU_pet']) 
tpetUU_singleton = readFile("{0}{1}.UU.hq.CB.nr.singleton.statCB.tsv.gz".format(RUN,pairlabel), ['CorrectBarcode', '2tag_UU_singleton']) 
tpetUU_other = readFile("{0}{1}.UU.hq.CB.nr.other.statCB.tsv.gz".format(RUN,pairlabel), ['CorrectBarcode', '2tag_UU_other']) 
tpetUU_pet = readFile("{0}{1}.UU.hq.CB.nr.pet.statCB.tsv.gz".format(RUN,pairlabel), ['CorrectBarcode', '2tag_UU_pet']) 
alldata6 = pd.concat([noneUU_singleton, noneUU_other, noneUU_pet, tpetUU_singleton, tpetUU_other, tpetUU_pet], axis=1)
alldata6 = alldata6.fillna(0).astype(int)
alldata6['singleton'] = alldata6['NoLinker_UU_singleton'] + alldata6['2tag_UU_singleton']
alldata6['other'] = alldata6['NoLinker_UU_other'] + alldata6['2tag_UU_other'] 
alldata6['pet'] = alldata6['NoLinker_UU_pet'] + alldata6['2tag_UU_pet']
alldata66 = alldata6[["singleton","other","pet"]]




########################## all
noneUU_allpet=readFile2("{0}{1}.UU.hq.CB.statPET.tsv.gz".format(RUN, nonelabel))
tpetUU_allpet=readFile2("{0}{1}.UU.hq.CB.statPET.tsv.gz".format(RUN, pairlabel))
noneUU_allpet.columns = ["{0}_NoLinker_all_pairs".format(x) for x in noneUU_allpet.columns]
tpetUU_allpet.columns = ["{0}_2tag_all_pairs".format(x) for x in tpetUU_allpet.columns]
alldata7 = pd.concat([noneUU_allpet, tpetUU_allpet], axis=1)
alldata7 = alldata7.fillna(0).astype(int)
for kk in ["nPET","nTrans","nCis",">1kb",">8kb"]:
    alldata7["{0}_all_pairs".format(kk)] = alldata7["{0}_NoLinker_all_pairs".format(kk)] + alldata7["{0}_2tag_all_pairs".format(kk)]
alldata77 = alldata7[["nPET_all_pairs","nTrans_all_pairs","nCis_all_pairs",">1kb_all_pairs",">8kb_all_pairs"]]

########################## nr
noneUU_nrpet=readFile2("{0}{1}.UU.hq.CB.nr.statPET.tsv.gz".format(RUN, nonelabel))
tpetUU_nrpet=readFile2("{0}{1}.UU.hq.CB.nr.statPET.tsv.gz".format(RUN, pairlabel))
noneUU_nrpet.columns = ["{0}_NoLinker_nr_pairs".format(x) for x in noneUU_nrpet.columns]
tpetUU_nrpet.columns = ["{0}_2tag_nr_pairs".format(x) for x in tpetUU_nrpet.columns]
alldata8 = pd.concat([noneUU_nrpet, tpetUU_nrpet], axis=1)
alldata8 = alldata8.fillna(0).astype(int)
for kk in ["nPET","nTrans","nCis",">1kb",">8kb"]:
    alldata8["{0}_nr_pairs".format(kk)] = alldata8["{0}_NoLinker_nr_pairs".format(kk)] + alldata8["{0}_2tag_nr_pairs".format(kk)]
alldata88 = alldata8[["nPET_nr_pairs","nTrans_nr_pairs","nCis_nr_pairs",">1kb_nr_pairs",">8kb_nr_pairs"]]

########################## singleton other pet
noneUU_singletonpet = readFile2("{0}{1}.UU.hq.CB.nr.singleton.statPET.tsv.gz".format(RUN,nonelabel))
noneUU_otherpet = readFile2("{0}{1}.UU.hq.CB.nr.other.statPET.tsv.gz".format(RUN,nonelabel))
noneUU_petpet = readFile2("{0}{1}.UU.hq.CB.nr.pet.statPET.tsv.gz".format(RUN,nonelabel))
tpetUU_singletonpet = readFile2("{0}{1}.UU.hq.CB.nr.singleton.statPET.tsv.gz".format(RUN,pairlabel))
tpetUU_otherpet = readFile2("{0}{1}.UU.hq.CB.nr.other.statPET.tsv.gz".format(RUN,pairlabel))
tpetUU_petpet = readFile2("{0}{1}.UU.hq.CB.nr.pet.statPET.tsv.gz".format(RUN,pairlabel))
noneUU_singletonpet.columns = ["{0}_NoLinker_singleton_pairs".format(x) for x in noneUU_singletonpet.columns]
noneUU_otherpet.columns = ["{0}_NoLinker_other_pairs".format(x) for x in noneUU_otherpet.columns]
noneUU_petpet.columns = ["{0}_NoLinker_pet_pairs".format(x) for x in noneUU_petpet.columns]
tpetUU_singletonpet.columns = ["{0}_2tag_singleton_pairs".format(x) for x in tpetUU_singletonpet.columns]
tpetUU_otherpet.columns = ["{0}_2tag_other_pairs".format(x) for x in tpetUU_otherpet.columns]
tpetUU_petpet.columns = ["{0}_2tag_pet_pairs".format(x) for x in tpetUU_petpet.columns]
alldata9 = pd.concat([noneUU_singletonpet, noneUU_otherpet, noneUU_petpet, tpetUU_singletonpet, tpetUU_otherpet, tpetUU_petpet], axis=1)
alldata9 = alldata9.fillna(0).astype(int)
for kk in ["nPET","nTrans","nCis",">1kb",">8kb"]:
    alldata9['{0}_singleton_pairs'.format(kk)] = alldata9['{0}_NoLinker_singleton_pairs'.format(kk)] + alldata9['{0}_2tag_singleton_pairs'.format(kk)]
    alldata9['{0}_other_pairs'.format(kk)] = alldata9['{0}_NoLinker_other_pairs'.format(kk)] + alldata9['{0}_2tag_other_pairs'.format(kk)] 
    alldata9['{0}_pet_pairs'.format(kk)] = alldata9['{0}_NoLinker_pet_pairs'.format(kk)] + alldata9['{0}_2tag_pet_pairs'.format(kk)]
alldata99 = alldata9[["nPET_singleton_pairs","nTrans_singleton_pairs","nCis_singleton_pairs",">1kb_singleton_pairs",">8kb_singleton_pairs", "nPET_other_pairs","nTrans_other_pairs","nCis_other_pairs",">1kb_other_pairs",">8kb_other_pairs", "nPET_pet_pairs","nTrans_pet_pairs","nCis_pet_pairs",">1kb_pet_pairs",">8kb_pet_pairs"]]

atab = pd.concat([alldata1, alldata2, alldata3, alldata4, alldata5, alldata6, alldata7, alldata8, alldata9], axis=1)
atab = atab.fillna(0).astype(int)
atab.to_csv("{0}.full.CBcounts.tsv".format(RUN), sep='\t', index=True)
#sums
aa = atab.sum()
fout = writeFile("{0}.full.CBsums.tsv".format(RUN))
fout.write("type\t")
fout.write("{0}\n".format("\t".join(aa.index)))
fout.write("all\t")
fout.write("{0}\n".format("\t".join(map(str,aa))))
fout.write("noCB\t")
fout.write("{0}\n".format("\t".join(map(str, atab.loc["noCB"]))))
bb = atab[~atab.index.isin(["noCB"])].sum()
fout.write("CB\t")
fout.write("{0}\n".format("\t".join(map(str,bb))))
fout.close()


atab = pd.concat([alldata1, alldata22, alldata33, alldata44, alldata55, alldata66, alldata77, alldata88, alldata99], axis=1)
atab = atab.fillna(0).astype(int)
atab.to_csv("{0}.simple.CBcounts.tsv".format(RUN), sep='\t', index=True)
#sums
aa = atab.sum()
fout = writeFile("{0}.simple.CBsums.tsv".format(RUN))
fout.write("type\t")
fout.write("{0}\n".format("\t".join(aa.index)))
fout.write("all\t")
fout.write("{0}\n".format("\t".join(map(str,aa))))
fout.write("noCB\t")
fout.write("{0}\n".format("\t".join(map(str, atab.loc["noCB"]))))
bb = atab[~atab.index.isin(["noCB"])].sum()
fout.write("CB\t")
fout.write("{0}\n".format("\t".join(map(str,bb))))
fout.close()