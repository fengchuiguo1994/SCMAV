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

fin = readFile(sys.argv[1])
fout = writeFile(sys.argv[2])
mydict = {}
for line in fin:
    tmp = line.strip().split()
    if tmp[10] not in mydict:
        mydict[tmp[10]] = [0,0,0,0,0,0,0,0] # total inter intra intra<=1000 intra>1000 intra>8000 intra>10000 intra>20000
    mydict[tmp[10]][0] += 1
    if tmp[0] == tmp[3]:
        mydict[tmp[10]][2] += 1
        if int(tmp[5]) - int(tmp[1]) >= 20000:
            mydict[tmp[10]][7] += 1
        if int(tmp[5]) - int(tmp[1]) >= 10000:
            mydict[tmp[10]][6] += 1
        if int(tmp[5]) - int(tmp[1]) >= 8000:
            mydict[tmp[10]][5] += 1
        if int(tmp[5]) - int(tmp[1]) >= 1000:
            mydict[tmp[10]][4] += 1
        else:
            mydict[tmp[10]][3] += 1
    else:
        mydict[tmp[10]][1] += 1

fout.write("PET_total\tPET_inter\tPET_intra\tPET_intrale1000\tPET_intragt1000\tPET_intragt8000\tPET_intragt10000\tPET_intragt20000\n")
for k in mydict.keys():
    fout.write("{0}\t{1}\n".format(k, '\t'.join(map(str, mydict[k]))))
fin.close()
fout.close()
