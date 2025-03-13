import sys
import collections

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

from bx.intervals.intersection import Intersecter, Interval
def regionTree(resFrag,chrom,start,end,value):
    if chrom in resFrag:
        resFrag[chrom].add_interval(Interval(start,end,value))
    else:
        resFrag[chrom] = Intersecter()
        resFrag[chrom].add_interval(Interval(start,end,value))

def regionFind(tree,start,end):
    return tree.find(start,end)

resFrag = {}
chromall = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","genome"]
bindict = collections.OrderedDict()
step = 0.125
start = 10
for i in chromall:
    regionTree(resFrag,i,int(0),int(2**10),"bin{0}".format(1))
    bindict["bin{0}".format(1)] = 1
binindex = 2
while start <= 27.750:
    for i in chromall:
        regionTree(resFrag,i,int(2**start),int(2**(start+1)),"bin{0}".format(binindex))
        bindict["bin{0}".format(binindex)] = 1
    start += 0.125
    binindex += 1
binlist = bindict.keys()
# print(binlist)

fin = readFile(sys.argv[1])
files = {}
resultdict = {}
for i in chromall:
    files[i] = writeFile("{0}.{1}.txt".format(sys.argv[2],i))
    resultdict[i] = {}
barcodedict = collections.OrderedDict()

for line in fin:
    tmp = line.strip().split()
    if tmp[0] != tmp[3]:
        continue
    if tmp[0] not in chromall:
        continue
    dis = (int(tmp[5]) + int(tmp[4]) - int(tmp[2]) - int(tmp[1]))/2
    if dis < 1:
        dis = 1 
    barcodedict[tmp[6]] = 1
    hits = resFrag[tmp[0]].find(dis,dis+1)
    if len(hits) == 0:
        print(line)
        print(dis)
        sys.exit(1)
    else:
        for r in hits: # dict -> chrom -> bin -> barcode
            if r.value not in resultdict[tmp[0]]:
                resultdict[tmp[0]][r.value] = {}
            if tmp[6] not in resultdict[tmp[0]][r.value]:
                resultdict[tmp[0]][r.value][tmp[6]] = 0
            resultdict[tmp[0]][r.value][tmp[6]] += 1
            if r.value not in resultdict["genome"]:
                resultdict["genome"][r.value] = {}
            if tmp[6] not in resultdict["genome"][r.value]:
                resultdict["genome"][r.value][tmp[6]] = 0
            resultdict["genome"][r.value][tmp[6]] += 1
fin.close()

barcodelist = barcodedict.keys()
for i in chromall:
    files[i].write("bin\t")
    files[i].write("\t".join(barcodelist))
    files[i].write("\n")
    for k in binlist:
        files[i].write("{0}".format(k))
        for v in barcodelist:
            if k not in resultdict[i]:
                resultdict[i][k] = {}
            if v not in resultdict[i][k]:
                files[i].write("\t0")
            else:
                files[i].write("\t{0}".format(resultdict[i][k][v]))
        files[i].write("\n")
    files[i].close()
