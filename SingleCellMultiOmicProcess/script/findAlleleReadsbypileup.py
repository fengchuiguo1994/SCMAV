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

import pysam
def readSam(insamfile):
    if insamfile.endswith(".bam"):
        insam = pysam.AlignmentFile(insamfile,'rb')
    elif insamfile.endswith(".sam.gz"):
        insam = pysam.AlignmentFile(insamfile,'rb')
    elif insamfile.endswith(".sam"):
        insam = pysam.AlignmentFile(insamfile,'r')
    else:
        raise ValueError("the input sam/bam file is not end with sam or bam!")
    return insam
        
def writeSam(outsamfile,header):
    if outsamfile.endswith(".bam"):
        outsam = pysam.AlignmentFile(outsamfile,'wb',header=header)
    elif outsamfile.endswith(".sam"):
        outsam = pysam.AlignmentFile(outsamfile,'w',header=header)
    else:
        raise ValueError("the output sam/bam file is not end with sam or bam!")
    return outsam

snpDict = {}
fin = readFile(sys.argv[1])
for line in fin:
    tmp = line.strip().split()
    key = tmp[0] + ':' + tmp[2]
    if len(tmp) > 6:
        snpDict[key] = tmp[4:7]
    else:
        snpDict[key] = [tmp[4],tmp[5],"genome"]    
fin.close()

phase = {}
fin = readSam(sys.argv[2])
fout = writeSam(sys.argv[3],header=fin.header)
for read in fin:
    phase[read.query_name] = [[],[],[]] # P/M  block  pos

for pileupcolumn in fin.pileup(stepper="nofilter", ignore_overlaps=False, ignore_orphans=False, min_base_quality=20):
    chrom = fin.references[pileupcolumn.reference_id]
    pos = "{0}".format(pileupcolumn.pos)
    k = chrom + ':' + str(pileupcolumn.pos+1)
    if k in snpDict:
        bases = snpDict[k]
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                query_base = pileupread.alignment.query_sequence[pileupread.query_position]
                if query_base == snpDict[k][0]:
                    phase[pileupread.alignment.query_name][0].append("P")
                elif query_base == snpDict[k][1]:
                    phase[pileupread.alignment.query_name][0].append("M")
                else:
                    phase[pileupread.alignment.query_name][0].append("U")
                phase[pileupread.alignment.query_name][1].append(snpDict[k][2])
                phase[pileupread.alignment.query_name][2].append(pos)

static = {}
for read in fin.fetch():
    readID = read.query_name
    if len(phase[readID][0]) == 0:
        read.set_tag("PL","noallele")
        read.set_tag("AL","noallele")
        read.set_tag("BL","noallele")
        fout.write(read)
        if "noallele" not in static:
            static["noallele"] = 0
        static["noallele"] += 1
        continue
    read.set_tag("PL",";".join(phase[readID][2]))
    read.set_tag("AL",";".join(phase[readID][0]))
    read.set_tag("BL",";".join(phase[readID][1]))
    fout.write(read)
    if "allele" not in static:
        static["allele"] = 0
    static["allele"] += 1
fin.close()
fout.close()
for k in sorted(static.keys()):
    print("{0}: {1}".format(k,static[k]))
sys.stderr.close()
sys.stdout.close()
