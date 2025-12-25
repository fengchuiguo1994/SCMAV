import sys
import argparse

import gzip
def readFile(infile):
    if infile.endswith((".gz","gzip")):
        fin = gzip.open(infile,'rt')
    else:
        fin = open(infile,'r')
    return fin
        
def writeFile(outfile):
    if outfile == "":
        fout = sys.stdout
    elif outfile.endswith((".gz","gzip")):
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

class SamRead:
    def __init__(self,infile):
        self.fin = readSam(infile)
    def __iter__(self):
        flag = None
        outlist = []
        for read in self.fin:
            if flag != None and flag != read.query_name:
                yield outlist
                outlist = []
            outlist.append(read)
            flag = read.query_name
        self.fin.close()
        yield outlist
    def header(self):
        return self.fin.header

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bamfile", required=True, type=str, help="The input file.")
    parser.add_argument('-o', "--outfile", required=True, type=str, help="The output file.")
    parser.add_argument('-q', "--minmapq", default="0",   type=int, help="The min mapq [default: 0]")
    # parser.add_argument('-s', "--skipNoBarcode", default="T", type=str, choices=['T', 'F'], help="T: skip the reads without barcode; F: remain the reads without barcode [default: T]")
    args = parser.parse_args()

    fin = SamRead(args.bamfile)
    fout = writeFile(args.outfile)
    mapqcutoff = args.minmapq
    fout.write("#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 num\n")

    myhash = dict()
    tmpstr = ""
    for reads in fin:
        if len(reads) == 1 or reads[0].query_name != reads[1].query_name:
            print(reads[0])
            sys.exit("just has 1 end!\n")
        if not reads[0].has_tag("CB"):
            continue
        if reads[0].is_duplicate:
            continue
        elif len(reads) > 2:
            tmp = []
            for k in reads:
                if k.is_secondary or k.is_supplementary:
                    pass
                else:
                    tmp.append(k)
            reads = tmp
        chrom1 = reads[0].reference_name
        bed1q = reads[0].mapping_quality
        if reads[0].is_reverse:
            start1 = reads[0].reference_end
            strand1 = "-"
        else:
            start1 = reads[0].reference_start + 1
            strand1 = "+"
        chrom2 = reads[1].reference_name
        bed2q = reads[1].mapping_quality
        if bed1q < bed2q:
            mapq = bed1q
        else:
            mapq = bed2q
        if mapq < mapqcutoff:
            continue
        if reads[1].is_reverse:
            start2 = reads[1].reference_end
            strand2 = "-"
        else:
            start2 = reads[1].reference_start + 1
            strand2 = "+"
        rid = reads[0].query_name
        if chrom1 > chrom2:
            tmpstr = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}".format(".",chrom2,start2,chrom1,start1,strand2,strand1,reads[0].get_tag("CB"))
            chrom2, chrom1 = chrom1, chrom2
        elif chrom1 == chrom2 and start1 > start2:
            tmpstr = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}".format(".",chrom2,start2,chrom1,start1,strand2,strand1,reads[0].get_tag("CB"))
            chrom2, chrom1 = chrom1, chrom2
        else:
            tmpstr = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}".format(".",chrom1,start1,chrom2,start2,strand1,strand2,reads[0].get_tag("CB"))
        if chrom1 not in myhash:
            myhash[chrom1] = dict()
        if chrom2 not in myhash[chrom1]:
            myhash[chrom1][chrom2] = dict()
        if tmpstr not in myhash[chrom1][chrom2]:
            myhash[chrom1][chrom2][tmpstr] = 0
        myhash[chrom1][chrom2][tmpstr] += 1

    for c1 in myhash.keys():
        for c2 in myhash[c1].keys():
            for k in myhash[c1][c2]:
                fout.write("{0}\t{1}\n".format(k, myhash[c1][c2][k]))
    fout.close()