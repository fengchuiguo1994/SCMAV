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

def detPM(mystr):
    if mystr == "noallele":
        return "noallele"
    mylist = mystr.split(";")
    if mylist.count("P") >= mylist.count("M") * 4: # 80%
        return "P"
    elif mylist.count("M") >= mylist.count("P") * 4:
        return "M"
    else:
        return "amb"

parser = argparse.ArgumentParser(description = "  ")
parser.add_argument("-i","--bamfile",required=True,help=" the input bam file ")
parser.add_argument("-o","--outputprefix",required=True,help=" output prefix ")
args = parser.parse_args()

if __name__ == "__main__":
    fin = readSam(args.bamfile)
    filehash = {}
    filehash["flt"] = writeSam("{0}.flt.bam".format(args.outputprefix),fin.header)
    filehash["M"] = writeSam("{0}.flt.M.bam".format(args.outputprefix),fin.header)
    filehash["P"] = writeSam("{0}.flt.P.bam".format(args.outputprefix),fin.header)
    filehash["no"] = writeSam("{0}.noallele.bam".format(args.outputprefix),fin.header)
    filehash["amb"] = writeSam("{0}.amb.bam".format(args.outputprefix),fin.header)
    stat = {}
    for read in fin:
        flag = detPM(read.get_tag("AL"))
        if flag == "noallele":
            filehash["no"].write(read)
            if "noallele" not in stat:
                stat["noallele"] = 0
            stat["noallele"] += 1
        elif flag == "amb":
            filehash["amb"].write(read)
            if "amb" not in stat:
                stat["amb"] = 0
            stat["amb"] += 1
        else:
            filehash["flt"].write(read)
            if "flt" not in stat:
                stat["flt"] = 0
            stat["flt"] += 1
            if flag == "M":
                filehash["M"].write(read)
            else:
                filehash["P"].write(read)
            if flag not in stat:
                stat[flag] = 0
            stat[flag] += 1
    for k in sorted(filehash.keys()):
        filehash[k].close()
    fin.close()
    for i in stat:
        print("{0}\t{1}".format(i,stat[i]))