import argparse
import sys

import pysam
def readSam(insamfile):
    if insamfile.endswith((".bam",".sam.gz")):
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

import gzip
def readFile(infile):
    if infile.endswith((".gz","gzip",".GZ",".GZIP")):
        fin = gzip.open(infile,'rt')
    else:
        fin = open(infile,'r')
    return fin

### Interval tree ###
from bx.intervals.intersection import Intersecter, Interval
def regionTree(resFrag,chrom,start,end,flag=None):
    if chrom not in resFrag:
        resFrag[chrom] = Intersecter()
    resFrag[chrom].add_interval(Interval(int(start),int(end),flag))

def regionFind(tree,chrom,start,end):
    if chrom not in tree:
        return []
    return tree[chrom].find(start,end)

parser = argparse.ArgumentParser(description = "  ")
parser.add_argument("-i","--bamfile",required=True,help=" the input bam file ")
parser.add_argument("-a","--allelefile",required=True,help=" the input allele file ")
parser.add_argument("-o","--output",required=True,help=" output allele reads file ")
parser.add_argument("-q", "--basequality", required=False, type=int, default=0, help=" min base quality [0] ")
args = parser.parse_args()

if __name__ == "__main__":
    alleleTree = {}
    alleledict = {}
    fin = readFile(args.allelefile)
    for line in fin:
        tmp = line.strip().split()
        if len(tmp) > 6:
            regionTree(alleleTree,tmp[0],tmp[1],tmp[2])
            alleledict["\t".join(tmp[:2])] = [tmp[4],tmp[5],tmp[6]]
        else:
            regionTree(alleleTree,tmp[0],tmp[1],tmp[2])
            alleledict["\t".join(tmp[:2])] = [tmp[4],tmp[5],"genome"]
    fin.close()

    static = {}
    fin = readSam(args.bamfile)
    fout = writeSam(args.output,fin.header)
    for read in fin:
        # print(read.qname)
        if read.is_unmapped:
            read.set_tag("PL","noallele")
            read.set_tag("AL","noallele")
            read.set_tag("BL","noallele")
            fout.write(read)
            if "unmap" not in static:
                static["unmap"] = 0
            static["unmap"] += 1
            continue
        if len(regionFind(alleleTree,read.reference_name,read.reference_start,read.reference_end)) == 0:
            read.set_tag("PL","noallele")
            read.set_tag("AL","noallele")
            read.set_tag("BL","noallele")
            fout.write(read)
            if "noallele" not in static:
                static["noallele"] = 0
            static["noallele"] += 1
            continue
        poslist = [] # pos list: PL
        PMlist = [] # allele list: AL
        blocklist = [] # block list: BL
        for readpos,refpos,refbase in read.get_aligned_pairs(matches_only=True, with_seq=True):
            readbase = read.query_sequence[readpos]
            basequality = read.query_qualities[readpos]
            if basequality < args.basequality:
                continue
            # sys.stderr.write("{0}\t{1}\n".format(readbase,refbase))
            pos = "{0}\t{1}".format(read.reference_name,refpos)
            if pos in alleledict:
                if readbase == alleledict[pos][0]:
                    PMlist.append("P")
                elif readbase == alleledict[pos][1]:
                    PMlist.append("M")
                else:
                    PMlist.append("U")
                poslist.append("{0}".format(refpos))
                blocklist.append(alleledict[pos][2])
        if len(poslist) == 0:
            read.set_tag("PL","noallele")
            read.set_tag("AL","noallele")
            read.set_tag("BL","noallele")
            fout.write(read)
            if "noallele" not in static:
                static["noallele"] = 0
            static["noallele"] += 1
            continue
        read.set_tag("PL",";".join(poslist))
        read.set_tag("AL",";".join(PMlist))
        read.set_tag("BL",";".join(blocklist))
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
