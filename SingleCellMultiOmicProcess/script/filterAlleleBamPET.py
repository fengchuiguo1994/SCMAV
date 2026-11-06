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
    filehash["MP"] = writeSam("{0}.flt.MP.bam".format(args.outputprefix),fin.header)
    filehash["no"] = writeSam("{0}.noallele.bam".format(args.outputprefix),fin.header)
    filehash["amb"] = writeSam("{0}.amb.bam".format(args.outputprefix),fin.header)
    stat = {}
    for read1 in fin:
        read2 = next(fin)
        flag1 = detPM(read1.get_tag("AL"))
        flag2 = detPM(read2.get_tag("AL"))
        if flag1 == "noallele" and flag2 == "noallele":
            filehash["no"].write(read1)
            filehash["no"].write(read2)
            if "noallele" not in stat:
                stat["noallele"] = 0
            stat["noallele"] += 1
        elif flag1 == "amb" or flag2 == "amb":
            filehash["amb"].write(read1)
            filehash["amb"].write(read2)
            if "amb" not in stat:
                stat["amb"] = 0
            stat["amb"] += 1
        else:
            filehash["flt"].write(read1)
            filehash["flt"].write(read2)
            if "flt" not in stat:
                stat["flt"] = 0
            stat["flt"] += 1
            if flag1 == "noallele":
                if flag2 == "M":
                    filehash["M"].write(read1)
                    filehash["M"].write(read2)
                elif flag2 == "P":
                    filehash["P"].write(read1)
                    filehash["P"].write(read2)
                else:
                    print(read1)
                    sys.exit("error1")
            elif flag2 == "noallele":
                if flag1 == "M":
                    filehash["M"].write(read1)
                    filehash["M"].write(read2)
                elif flag1 == "P":
                    filehash["P"].write(read1)
                    filehash["P"].write(read2)
                else:
                    print(read1)
                    sys.exit("error2")
            elif flag1 == flag2:
                filehash[flag1].write(read1)
                filehash[flag1].write(read2)
            elif flag1 != flag2:
                filehash["MP"].write(read1)
                filehash["MP"].write(read2)
            else:
                print(read1)
                sys.exit("error3")


    for k in sorted(filehash.keys()):
        filehash[k].close()
    fin.close()
    for i in stat:
        print("{0}\t{1}".format(i,stat[i]))