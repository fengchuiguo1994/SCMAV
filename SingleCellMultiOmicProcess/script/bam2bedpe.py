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
    parser.add_argument('-d', "--dirname", required=True, type=str, help="The contact number matrix.")
    # parser.add_argument('-s', "--skipNoBarcode", default="T", type=str, choices=['T', 'F'], help="T: skip the reads without barcode; F: remain the reads without barcode [default: T]")
    args = parser.parse_args()

    fin = SamRead(args.bamfile)
    mapqcutoff = args.minmapq

    myhash = dict()
    tmpstr = ""
    read_file_counter = 0
    for reads in fin:
        if len(reads) != 2 or reads[0].query_name != reads[1].query_name:
            print(reads)
            sys.exit("error! not has two reads")
        if not reads[0].has_tag("CB"):
            continue
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
        barcode = reads[0].get_tag("CB")
        if bed1c == bed2c:
            if bed1s == bed2s:
                if bed1e > bed2e:
                    tmpstr = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}".format(bed2c,bed2s,bed2e,bed1c,bed1s,bed1e,".",mapq,bed2strand,bed1strand,barcode)
                    bed2c, bed1c = bed1c, bed2c
                else:
                    tmpstr = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}".format(bed1c,bed1s,bed1e,bed2c,bed2s,bed2e,".",mapq,bed1strand,bed2strand,barcode)
            elif bed1s > bed2s:
                tmpstr = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}".format(bed2c,bed2s,bed2e,bed1c,bed1s,bed1e,".",mapq,bed2strand,bed1strand,barcode)
                bed2c, bed1c = bed1c, bed2c
            else:
                tmpstr = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}".format(bed1c,bed1s,bed1e,bed2c,bed2s,bed2e,".",mapq,bed1strand,bed2strand,barcode)
        elif bed1c > bed2c:
            tmpstr = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}".format(bed2c,bed2s,bed2e,bed1c,bed1s,bed1e,".",mapq,bed2strand,bed1strand,barcode)
            bed2c, bed1c = bed1c, bed2c
        else:
            tmpstr = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}".format(bed1c,bed1s,bed1e,bed2c,bed2s,bed2e,".",mapq,bed1strand,bed2strand,barcode)
        if bed1c not in myhash:
            myhash[bed1c] = dict()
        if bed2c not in myhash[bed1c]:
            myhash[bed1c][bed2c] = dict()
        if tmpstr not in myhash[bed1c][bed2c]:
            myhash[bed1c][bed2c][tmpstr] = 0
            read_file_counter += 1
        myhash[bed1c][bed2c][tmpstr] += 1
    
    fout = writeFile(args.outfile)
    fout2 = writeFile("{0}.4dn".format(args.outfile))
    fout2.write("#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 barcode\n")

    fout3 = writeFile("{0}/barcodes.tsv".format(args.dirname))
    fout4 = writeFile("{0}/features.tsv".format(args.dirname))
    fout5 = writeFile("{0}/matrix.mtx".format(args.dirname))
    totalN = 0
    barcodenumdict = OrderedDict()
    numbarcodedict = OrderedDict()
    featurenumdict = OrderedDict()
    numfeaturedict = OrderedDict()
    matdict = OrderedDict()
    
    for c1 in myhash.keys():
        for c2 in myhash[c1].keys():
            for k in myhash[c1][c2]:
                fout.write("{0}\t{1}\n".format(k, myhash[c1][c2][k]))
                tmp = k.split("\t")
                if tmp[8] == "-":
                    start1 = tmp[2]
                else:
                    start1 = int(tmp[1]) + 1
                if tmp[9] == "-":
                    start2 = tmp[5]
                else:
                    start2 = int(tmp[4]) + 1
                fout2.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(".",tmp[0],start1,tmp[3],start2,tmp[8],tmp[9],tmp[10]))
                if tmp[10] not in barcodenumdict:
                    bnumber = len(barcodenumdict) + 1
                    barcodenumdict[tmp[10]] = bnumber
                    numbarcodedict[bnumber] = tmp[10]
                bnumber = barcodenumdict[tmp[10]]
                feature = "{0}-{1}".format(tmp[0], tmp[3])
                # feature = "{0}".format(tmp[0])
                # feature = "genome"
                # STAR --runMode soloCellFiltering test3 filterd3/tmp --soloCellFilter EmptyDrops_CR # three result are the same.
                if feature not in featurenumdict:
                    fnumber = len(featurenumdict) + 1
                    featurenumdict[feature] = fnumber
                    numfeaturedict[fnumber] = feature
                fnumber = featurenumdict[feature]
                if bnumber not in matdict:
                    matdict[bnumber] = dict()
                if fnumber not in matdict[bnumber]:
                    matdict[bnumber][fnumber] = 0
                    totalN += 1
                matdict[bnumber][fnumber] += 1

    fout5.write("%%MatrixMarket matrix coordinate integer general\n%\n")
    fout5.write("{0} {1} {2}\n".format(len(featurenumdict), len(barcodenumdict), totalN))
    for k in barcodenumdict.keys():
        fout3.write("{0}\n".format(k))
    for k in featurenumdict.keys():
        fout4.write("{0}\n".format(k))
    for k in matdict.keys():
        for j in matdict[k].keys():
            fout5.write("{0} {1} {2}\n".format(j, k, matdict[k][j]))
    fout.close()
    fout2.close()
    fout3.close()
    fout4.close()
    fout5.close()

    print("PET_num\t{0}\n".format(read_file_counter))
    sys.stdout.close()