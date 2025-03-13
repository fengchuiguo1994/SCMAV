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

fin = readFile(sys.argv[1]) # gene info file
geneorder = []
for line in fin:
    tmp = line.strip().split()
    geneorder.append(tmp[6])
fin.close()

fin = readFile(sys.argv[2]) # barcode stat
barcodeset = collections.OrderedDict()
for line in fin:
    if line.startswith("PET_total"):
        continue
    tmp = line.strip().split()
    if int(tmp[1]) >= 50:
        barcodeset[tmp[0]] = 1
fin.close()

fin = readFile(sys.argv[3]) # bedpe overlap file
myhash = {}
for line in fin:
    tmp = line.strip().split()
    if tmp[6] not in barcodeset:
        continue
    if tmp[17] not in myhash: # gene
        myhash[tmp[17]] = {}
    if tmp[6] not in myhash[tmp[17]]: # barcode
        myhash[tmp[17]][tmp[6]] = set() # hash -> gene -> barocde -> readid
    myhash[tmp[17]][tmp[6]].add(tmp[7])
fin.close()

fout = writeFile(sys.argv[4]) # output
barcodeorder = barcodeset.keys()
fout.write("{0}\n".format("\t".join(barcodeorder)))
for k in myhash.keys():
    mystr = "{0}".format(k)
    for b in barcodeorder: # loop barcode
        if b in myhash[k]:
            mystr += "\t{0}".format(len(myhash[k][b]))
        else:
            mystr += "\t0"
    fout.write("{0}\n".format(mystr))
fout.close()