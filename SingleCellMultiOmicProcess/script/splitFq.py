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

class FqRead2:
    def __init__(self,infile):
        self.fin = readFile(infile)
    def __iter__(self):
        for line in self.fin:
            rid = line.strip()
            rseq = self.fin.readline()
            rseq = rseq.strip()
            rsyb = self.fin.readline()
            rsyb = rsyb.strip()
            rqual = self.fin.readline()
            rqual = rqual.strip()
            yield rid,rseq,rsyb,rqual
        self.fin.close()

def rc(sequence):
    seq = sequence[::-1]
    trantab = str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh')
    string = seq.translate(trantab)
    return string

fin1 = FqRead2(sys.argv[1]) # barcode
fin2 = FqRead2(sys.argv[2]) # R1
fin3 = FqRead2(sys.argv[3]) # R2
nreads = int(sys.argv[4]) # number
outdir = sys.argv[5] # outdir
outprefix = sys.argv[6] # prefix
start = int(sys.argv[7]) # start
barcodelen = int(sys.argv[8]) # len
if len(sys.argv) == 10:
    if sys.argv[9] == "Reverse":
        rcflag = True
    elif sys.argv[9] == "Forward":
        rcflag = False
    else:
        sys.exit("Error")
else:
    rcflag = False

index = 0
for (rid1,rseq1,rsyb1,rqual1), (rid2,rseq2,rsyb2,rqual2), (rid3,rseq3,rsyb3,rqual3) in zip(fin1, fin2, fin3):
    if index == 0:
        fout1 = writeFile("{0}/{1}.barcode.{2}.fastq".format(outdir, outprefix, 1))
        fout2 = writeFile("{0}/{1}.R1.{2}.fastq".format(outdir, outprefix, 1))
        fout22 = writeFile("{0}/{1}.R1.{2}.tmp.fastq".format(outdir, outprefix, 1))
        fout3 = writeFile("{0}/{1}.R2.{2}.fastq".format(outdir, outprefix, 1))
    elif index % nreads == 0:
        fout1.close()
        fout2.close()
        fout3.close()
        fout1 = writeFile("{0}/{1}.barcode.{2}.fastq".format(outdir, outprefix, index//nreads + 1))
        fout2 = writeFile("{0}/{1}.R1.{2}.fastq".format(outdir, outprefix, index//nreads + 1))
        fout22 = writeFile("{0}/{1}.R1.{2}.tmp.fastq".format(outdir, outprefix, index//nreads + 1))
        fout3 = writeFile("{0}/{1}.R2.{2}.fastq".format(outdir, outprefix, index//nreads + 1))
    index += 1

    barcode = rseq1[start-1:start-1+barcodelen]
    barcodequal = rqual1[start-1:start-1+barcodelen]

    if rcflag:
        barcode = rc(barcode)
        barcodequal = barcodequal[::-1]
    fout1.write("@{0}\n{1}N\n{2}\n{3}F\n".format(rid1, barcode, rsyb1, barcodequal))
    fout2.write("{0}\n{1}\n{2}\n{3}\n".format(rid2,rseq2,rsyb2,rqual2))
    fout22.write("{0}\nNNNNNNNNN\n{1}\nFFFFFFFFF\n".format(rid2,rsyb2))
    fout3.write("{0}\n{1}\n{2}\n{3}\n".format(rid3,rseq3,rsyb3,rqual3))
fout1.close()
fout2.close()
fout22.close()
fout3.close()

fout = writeFile("{0}/{1}.readsnum".format(outdir, outprefix))
fout.write("{0}\t{1}\n".format("Total_Reads", index))
fout.close()
