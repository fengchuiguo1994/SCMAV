import sys
import argparse
from collections import OrderedDict

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
    args = parser.parse_args()

    fin = SamRead(args.bamfile)
    fout = writeFile(args.outfile)
    mapqcutoff = args.minmapq

    myhash = dict()
    tmpstr = ""
    read_file_counter = 0
    for reads in fin:
        if len(reads) != 2 or reads[0].query_name != reads[1].query_name:
            print(reads)
            sys.exit("error! not has two reads")
        if reads[0].is_duplicate:
            continue
        bed1c = reads[0].reference_name
        bed1s = reads[0].reference_start
        bed1e = reads[0].reference_end
        bed1q = reads[0].mapping_quality
        bed1strand = "+"
        if reads[0].is_reverse:
            bed1strand = "-"
        bed2c = reads[1].reference_name
        bed2s = reads[1].reference_start
        bed2e = reads[1].reference_end
        bed2q = reads[1].mapping_quality
        if bed1q < bed2q:
            mapq = bed1q
        else:
            mapq = bed2q
        if mapq < mapqcutoff:
            continue
        bed2strand = "+"
        if reads[1].is_reverse:
            bed2strand = "-"
        name = reads[0].query_name
        if bed1c == bed2c:
            if bed1s == bed2s:
                if bed1e > bed2e:
                    fout.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n".format(bed2c,bed2s,bed2e,bed1c,bed1s,bed1e,name,mapq,bed2strand,bed1strand))
                else:
                    fout.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n".format(bed1c,bed1s,bed1e,bed2c,bed2s,bed2e,name,mapq,bed1strand,bed2strand))
            elif bed1s > bed2s:
                fout.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n".format(bed2c,bed2s,bed2e,bed1c,bed1s,bed1e,name,mapq,bed2strand,bed1strand))
            else:
                fout.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n".format(bed1c,bed1s,bed1e,bed2c,bed2s,bed2e,name,mapq,bed1strand,bed2strand))
        elif bed1c > bed2c:
            fout.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n".format(bed2c,bed2s,bed2e,bed1c,bed1s,bed1e,name,mapq,bed2strand,bed1strand))
        else:
            fout.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n".format(bed1c,bed1s,bed1e,bed2c,bed2s,bed2e,name,mapq,bed1strand,bed2strand))
    fout.close()