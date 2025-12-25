import sys

import gzip
def readFile(infile):
    if infile.endswith((".gz","gzip")):
        fin = gzip.open(infile,'rt')
    else:
        fin = open(infile,'r')
    return fin
        
def writeFile(outfile):
    if outfile.endswith((".gz","gzip")):
        fout = gzip.open(outfile,'wt')
    else:
        fout = open(outfile,'w')
    return fout

class LinkerRead:
    def __init__(self,infile):
        if infile == None or infile == "-":
            self.fin = sys.stdin
        else:
            self.fin = readFile(infile)
    def __iter__(self):
        flag = None
        outstr = []
        for line in self.fin:
            tmp = line.strip('\r\n').split("\t")
            if flag != None and flag != tmp[0]:
                yield outstr
                outstr = []
            flag = tmp[0]
            outstr.append(tmp)
        self.fin.close()
        yield outstr

def getRead(linker):
    seq = []
    qual = []
    linkerset = set()
    rid = linker[0][0].split()[0]
    if linker[0][1] == "-1":
        seq.append(linker[0][2])
        qual.append(linker[0][3])
        return rid, seq, qual, 1, "", linkerset
    for i in range(len(linker)):
        if i == 0:
            seq = [linker[i][4], linker[i][6]]
            qual = [linker[i][8], linker[i][10]]
            if linker[i][7].endswith("_"):
                linkerset.add(linker[i][7][:-1])
            else:
                linkerset.add(linker[i][7])
        else: # 
            tmpstr = "".join(linker[i][4:7])
            for iii,jjj in enumerate(seq):
                if jjj == tmpstr:
                    break
            seq[iii] = linker[i][4]
            seq.insert(iii+1, linker[i][6])
            qual[iii] = linker[i][8]
            qual.insert(iii+1, linker[i][10])
            if linker[i][7].endswith("_"):
                linkerset.add(linker[i][7][:-1])
            else:
                linkerset.add(linker[i][7])
    # print(seq)
    # if len(linkerset) == 1:
    #     linkerstr = list(linkerset)[0]
    # else:
    #     linkerstr = "mix{0}".format(";".join(sorted(linkerset)))
    linkerstr = "{0}".format(";".join(sorted(linkerset)))
    return rid, seq, qual, len(seq), linkerstr, linkerset

def printOut(flag, rid, seq, qual, count, linkerstr, linkerset, fout, foutshort):
    # print(seq)
    if flag and len(linkerset) != 0:
        for i in range(count):
            if len(seq[i]) >= minlen:
                if len(seq[i]) >= 70:
                    fout.write("@{0}_{1}\n{2}\n+\n{3}\n".format(rid,"{0}_{1}_{2}".format(i+1, count, linkerstr), seq[i], qual[i]))
                else:
                    foutshort.write("@{0}_{1}\n{2}\n+\n{3}\n".format(rid,"{0}_{1}_{2}".format(i+1, count, linkerstr), seq[i], qual[i]))
    else:
        for i in range(count):
            if len(seq[i]) >= minlen:
                if len(seq[i]) >= 70:
                    fout.write("@{0}_{1}\n{2}\n+\n{3}\n".format(rid,"{0}_{1}".format(i+1, count), seq[i], qual[i]))
                else:
                    foutshort.write("@{0}_{1}\n{2}\n+\n{3}\n".format(rid,"{0}_{1}".format(i+1, count), seq[i], qual[i]))

fin1 = LinkerRead(sys.argv[1])
finunpair1 = LinkerRead(sys.argv[2])
fin2 = LinkerRead(sys.argv[3])
finunpair2 = LinkerRead(sys.argv[4])
outprefix = sys.argv[5]
minlen = int(sys.argv[6]) # 18
flag = False
if len(sys.argv) > 7:
    if sys.argv[7] == "T":
        flag = True

fout1 = writeFile("{0}.R1.pair.long.fastq".format(outprefix))
fout1short = writeFile("{0}.R1.pair.short.fastq".format(outprefix))
foutunpair1 = writeFile("{0}.R1.unpair.long.fastq".format(outprefix))
foutunpair1short = writeFile("{0}.R1.unpair.short.fastq".format(outprefix))
fout2 = writeFile("{0}.R2.pair.long.fastq".format(outprefix))
fout2short = writeFile("{0}.R2.pair.short.fastq".format(outprefix))
foutunpair2 = writeFile("{0}.R2.unpair.long.fastq".format(outprefix))
foutunpair2short = writeFile("{0}.R2.unpair.short.fastq".format(outprefix))
linkerstat = {"Linker":0, "Linker1":0, "Linker2":0, "NoLinker":0, "pairData":0, "R1Data":0, "R2Data":0, "R1Data1":0, "R2Data2":0, "NoData":0, "Linker-pairData":0, "Linker-R1Data":0, "Linker-R2Data":0, "Linker-NoData":0, "NoLinker-pairData":0, "NoLinker-R1Data":0, "NoLinker-R2Data":0, "NoLinker-NoData":0}

for linker1, linker2 in zip(fin1, fin2):
    if len(linker1) == 0 and len(linker2) == 0:
        continue
    rid1, seq1, qual1, count1, linkerstr1, linkerset1 = getRead(linker1)
    rid2, seq2, qual2, count2, linkerstr2, linkerset2 = getRead(linker2)
    if len(linkerset1) > 0 or len(linkerset2) > 0:
        linkerstat["Linker"] += 1
        linkflag = "Linker"
    else:
        linkerstat["NoLinker"] += 1
        linkflag = "NoLinker"
    flag1 = False
    for k in seq1:
        if len(k) >= minlen:
            flag1 = True
            break
    flag2 = False
    for k in seq2:
        if len(k) >= minlen:
            flag2 = True
            break
    if flag1 and flag2:
        printOut(flag, rid1, seq1, qual1, count1, linkerstr1, linkerset1, fout1, fout1short)
        printOut(flag, rid2, seq2, qual2, count2, linkerstr2, linkerset2, fout2, fout2short)
        linkerstat["pairData"] += 1
        if linkflag == "Linker":
            linkerstat["Linker-pairData"] += 1
        else:
            linkerstat["NoLinker-pairData"] += 1
    elif flag1:
        printOut(flag, rid1, seq1, qual1, count1, linkerstr1, linkerset1, foutunpair1, foutunpair1short)
        linkerstat["R1Data"] += 1
        if linkflag == "Linker":
            linkerstat["Linker-R1Data"] += 1
        else:
            linkerstat["NoLinker-R1Data"] += 1
    elif flag2:
        printOut(flag, rid2, seq2, qual2, count2, linkerstr2, linkerset2, foutunpair2, foutunpair2short)
        linkerstat["R2Data"] += 1
        if linkflag == "Linker":
            linkerstat["Linker-R2Data"] += 1
        else:
            linkerstat["NoLinker-R2Data"] += 1
    else:
        linkerstat["NoData"] += 1
        if linkflag == "Linker":
            linkerstat["Linker-NoData"] += 1
        else:
            linkerstat["NoLinker-NoData"] += 1

for linker1 in finunpair1:
    if len(linker1) == 0:
        continue
    rid1, seq1, qual1, count1, linkerstr1, linkerset1 = getRead(linker1)
    if len(linkerset1) > 0:
        linkerstat["Linker"] += 1
        linkflag = "Linker"
    else:
        linkerstat["NoLinker"] += 1
        linkflag = "NoLinker"
    flag1 = False
    for k in seq1:
        if len(k) >= minlen:
            flag1 = True
            break
    if flag1:
        printOut(flag, rid1, seq1, qual1, count1, linkerstr1, linkerset1, foutunpair1, foutunpair1short)
        linkerstat["R1Data"] += 1
        if linkflag == "Linker":
            linkerstat["Linker-R1Data"] += 1
        else:
            linkerstat["NoLinker-R1Data"] += 1
    else:
        linkerstat["NoData"] += 1
        if linkflag == "Linker":
            linkerstat["Linker-NoData"] += 1
        else:
            linkerstat["NoLinker-NoData"] += 1

for linker2 in finunpair2:
    if len(linker2) == 0:
        continue
    rid2, seq2, qual2, count2, linkerstr2, linkerset2 = getRead(linker2)
    if len(linkerset2) > 0:
        linkerstat["Linker"] += 1
        linkflag = "Linker"
    else:
        linkerstat["NoLinker"] += 1
        linkflag = "NoLinker"
    flag2 = False
    for k in seq2:
        if len(k) >= minlen:
            flag2 = True
            break
    if flag2:
        printOut(flag, rid2, seq2, qual2, count2, linkerstr2, linkerset2, foutunpair2, foutunpair2short)
        linkerstat["R2Data"] += 1
        if linkflag == "Linker":
            linkerstat["Linker-R2Data"] += 1
        else:
            linkerstat["NoLinker-R2Data"] += 1
    else:
        linkerstat["NoData"] += 1
        if linkflag == "Linker":
            linkerstat["Linker-NoData"] += 1
        else:
            linkerstat["NoLinker-NoData"] += 1

fout1.close()
fout1short.close()
foutunpair1.close()
foutunpair1short.close()
fout2.close()
fout2short.close()
foutunpair2.close()
foutunpair2short.close()

print("Linker\t{0}".format(linkerstat["Linker"]))
print("NoLinker\t{0}".format(linkerstat["NoLinker"]))
print("pairData\t{0}".format(linkerstat["pairData"]))
print("R1Data\t{0}".format(linkerstat["R1Data"]))
print("R2Data\t{0}".format(linkerstat["R2Data"]))
print("NoData\t{0}".format(linkerstat["NoData"]))
print("Linker-pairData\t{0}".format(linkerstat["Linker-pairData"]))
print("Linker-R1Data\t{0}".format(linkerstat["Linker-R1Data"]))
print("Linker-R2Data\t{0}".format(linkerstat["Linker-R2Data"]))
print("Linker-NoData\t{0}".format(linkerstat["Linker-NoData"]))
print("NoLinker-pairData\t{0}".format(linkerstat["NoLinker-pairData"]))
print("NoLinker-R1Data\t{0}".format(linkerstat["NoLinker-R1Data"]))
print("NoLinker-R2Data\t{0}".format(linkerstat["NoLinker-R2Data"]))
print("NoLinker-NoData\t{0}".format(linkerstat["NoLinker-NoData"]))