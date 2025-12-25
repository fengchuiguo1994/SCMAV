import sys
from collections import Counter, OrderedDict

# ~/software/cellranger-atac/mro/atac/stages/metrics/estimate_library_complexity
# ~/software/cellranger-atac/mro/atac/stages/processing/generate_peak_matrix/__init__.py

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

### Interval tree ###
from bx.intervals.intersection import Intersecter, Interval
def regionTree(resFrag,chrom,start,end):
    if chrom not in resFrag:
        resFrag[chrom] = Intersecter()
    resFrag[chrom].add_interval(Interval(start,end,"{0}:{1}-{2}".format(chrom, start, end)))

def regionFind(tree,chrom,start,end):
    if chrom not in tree:
        return []
    return tree[chrom].find(start,end)

peakTree = {}
peakhash = OrderedDict()
fin = readFile(sys.argv[1]) # peak file
fpeak = writeFile(sys.argv[3])
for index, line in enumerate(fin):
    tmp = line.strip().split()
    peakhash["{0}:{1}-{2}".format(tmp[0], tmp[1], tmp[2])] = index + 1
    fpeak.write("{0}:{1}-{2}\t{0}:{1}-{2}\tPeaks\n".format(tmp[0], tmp[1], tmp[2]))
    regionTree(peakTree, tmp[0], int(tmp[1]), int(tmp[2]))
fin.close()
fpeak.close()

fin = readFile(sys.argv[2])
mat = {}
totalN = 0
barcodeset = set()
for line in fin:
    tmp = line.strip().split()
    barcodeset.add(tmp[3])
    results = regionFind(peakTree, tmp[0], int(tmp[1]), int(tmp[1])+1) # start
    for kk in results:
        if tmp[3] not in mat:
            mat[tmp[3]] = {}
        if kk.value not in mat[tmp[3]]:
            mat[tmp[3]][kk.value] = 0
            totalN += 1
        mat[tmp[3]][kk.value] += 1
    results = regionFind(peakTree, tmp[0], int(tmp[2])-1, int(tmp[2])) # end
    for kk in results:
        if tmp[3] not in mat:
            mat[tmp[3]] = {}
        if kk.value not in mat[tmp[3]]:
            mat[tmp[3]][kk.value] = 0
            totalN += 1
        mat[tmp[3]][kk.value] += 1
fin.close()

fbarcode = writeFile(sys.argv[4])
barcodehash = OrderedDict()
for index, k in enumerate(sorted(barcodeset)):
    barcodehash[k] = index + 1
    fbarcode.write("{0}\n".format(k))
fbarcode.close()

fmat = writeFile(sys.argv[5])
fmat.write("%%MatrixMarket matrix coordinate integer general\n%\n")
fmat.write("{0} {1} {2}\n".format(len(peakhash), len(barcodehash), totalN))
# for k in barcodehash.keys():
#     if k not in mat:
#         continue
#     for j in peakhash.keys():
#         if j in mat[k]:
#             fmat.write("{0} {1} {2}\n".format(peakhash[j], barcodehash[k], mat[k][j]))
for k in mat.keys():
    for j in mat[k].keys():
        fmat.write("{0} {1} {2}\n".format(peakhash[j], barcodehash[k], mat[k][j]))
fmat.close()