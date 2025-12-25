import argparse
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

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-sb", "--STARbamfile", required=True, type=str, help="The input sam/bam file with barcode sequences")
    parser.add_argument('-bb', "--bwabamfile", required=True, type=str, help="The input sam/bam file without barcode sequences")
    parser.add_argument('-o', "--outbam", required=True, help="The output bam file. Automatically recognize sam or bam")
    parser.add_argument('-s', "--statfile", required=True, help="The stat file")
    parser.add_argument("-t", "--barcodeTranslate", default="", type=str, help="Convert barcode to the specified sequences during output [default: ]")
    args = parser.parse_args()

    barcodedict = dict()
    if args.barcodeTranslate != "":
        fin = readFile(args.barcodeTranslate)
        for line in fin:
            tmp = line.strip().split(",")
            barcodedict[tmp[1]] = tmp[0]
        fin.close()

    mydict = {}
    fin = readSam(args.STARbamfile)
    for read in fin:
        if read.has_tag("CB"):
            if read.get_tag("CB") != "-":
                if args.barcodeTranslate != "":
                    mydict[read.query_name] = barcodedict[read.get_tag("CB")]
                else:
                    mydict[read.query_name] = read.get_tag("CB")
    fin.close()

    fin = readSam(args.bwabamfile)
    fout = writeSam(args.outbam, header=fin.header)
    fstat = writeFile(args.statfile)
    index = 0
    for read in fin:
        if read.query_name in mydict:
            read.set_tag("CB", "{0}-1".format(mydict[read.query_name]))
            read.set_tag("ZB", "{0}-A".format(mydict[read.query_name]))
            fout.write(read)
            index += 1
        else:
            # read.set_tag("CB", "noBC")
            # read.set_tag("ZB", "noBC")
            pass
        # fout.write(read)
    fstat.write("{0}\t{1}\n".format("barcode_unique_mapping_reads", index//2))
    fstat.close()
    fin.close()
    fout.close()