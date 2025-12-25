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

if sys.argv[1] == "-":
    fin = sys.stdin
else:
    fin = readFile(sys.argv[1])
outprefix = sys.argv[2]
fout = writeFile("{0}.chrom.list".format(outprefix))
fheader = writeFile("{0}.header.list".format(outprefix))
filedict = {}
for line in fin:
    if line.startswith("#"):
        fheader.write(line)
        continue
    tmp = line.strip().split()
    flag = "{0}-{1}".format(tmp[1],tmp[3])
    if flag not in filedict:
        filedict[flag] = writeFile("{0}.{1}.tmp.txt".format(outprefix,flag))
    filedict[flag].write(line)
fin.close()

for k in sorted(filedict.keys()):
    filedict[k].close()
    fout.write("{0}\n".format(k.replace("-", " ")))
fout.close()
fheader.close()