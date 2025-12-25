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

if len(sys.argv) == 1:
    fin = sys.stdin
    fout = sys.stdout
else:
    if sys.argv[1] == "-":
        fin = sys.stdin
    else:
        fin = readFile(sys.argv[1])
    if len(sys.argv) == 2:
        fout = sys.stdout
    else:
        fout = writeFile(sys.argv[2])

fout.write("#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 barcode\n")

for line in fin:
    tmp = line.strip().split("\t")
    if tmp[8] == "-":
        start1 = tmp[2]
    else:
        start1 = int(tmp[1]) + 1
    if tmp[9] == "-":
        start2 = tmp[5]
    else:
        start2 = int(tmp[4]) + 1
    fout.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(".",tmp[0],start1,tmp[3],start2,tmp[8],tmp[9],tmp[10]))
fout.close()