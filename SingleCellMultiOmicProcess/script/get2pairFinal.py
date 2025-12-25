import sys
import copy

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
    def __init__(self,infile,mapq=30):
        self.fin = readSam(infile)
        self.mapq = mapq
    def __iter__(self):
        flag = None
        outlist1 = []
        outlist2 = []
        for read in self.fin:
            if read.is_supplementary:
                continue
            rid = read.query_name.split("_")[0]
            if flag != None and flag != rid:
                yield outlist1, outlist2
                outlist1 = []
                outlist2 = []
            if read.get_tag("SR") == "R1":
                if read.has_tag("XT") and read.get_tag("XT") == "U" and read.mapping_quality >= self.mapq:
                    outlist1.append(read)
            elif read.get_tag("SR") == "R2":
                if read.has_tag("XT") and read.get_tag("XT") == "U" and read.mapping_quality >= self.mapq:
                    outlist2.append(read)
            else:
                sys.exit("error format bam file")
            flag = rid
        self.fin.close()
        yield outlist1, outlist2
    def header(self):
        return self.fin.header

if len(sys.argv) > 3:
    fin = SamRead(sys.argv[1], int(sys.argv[3]))
else:
    fin = SamRead(sys.argv[1])
outprefix = sys.argv[2]
# fsingleTn = writeSam("{0}.single.bam".format(outprefix), header=fin.header())
fsingleTn1 = writeSam("{0}.single.withlinker.bam".format(outprefix), header=fin.header())
fsingleTn2 = writeSam("{0}.single.withoutlinker.bam".format(outprefix), header=fin.header())
# fpairTn = writeSam("{0}.pair.bam".format(outprefix), header=fin.header())
fpairTn1 = writeSam("{0}.pair.withlinker.bam".format(outprefix), header=fin.header())
fpairTn2 = writeSam("{0}.pair.withoutlinker.bam".format(outprefix), header=fin.header())

dataStat = {"R1R2":0, "R1R2_with_2End":0, "R1R2_with_0End":0, "R1R2_with_1End":0, 
            "R1":0, "R1_1tag":0, "R1_1tag_with_1End":0, "R1_1tag_with_0End":0, "R1_2tag":0, "R1_2tag_with_2End":0, "R1_2tag_with_0End":0, "R1_2tag_with_1End":0, "R2":0, "R2_1tag":0, "R2_1tag_with_1End":0, "R2_1tag_with_0End":0, "R2_2tag":0, "R2_2tag_with_2End":0, "R2_2tag_with_0End":0, "R2_2tag_with_1End":0, "None":0, "pair":0, "single":0}
for reads1, reads2 in fin:
    if len(reads1) > 0 and len(reads2) > 0:
        linker = set()
        dataStat["R1R2"] += 1
        tmp1 = [[int(k.query_name.split("_")[1]), k] for k in reads1]
        tmp1.sort(key=lambda x: x[0])
        R1read = copy.deepcopy(tmp1[0][1])
        r1tmp = R1read.query_name.split("_")
        R1read.query_name = r1tmp[0]
        R1read.set_tag("BP", "{0}_{1}".format(r1tmp[1], r1tmp[2]))
        if len(r1tmp) > 3:
            linker.update(r1tmp[3].split(";"))
        R1read.flag = R1read.flag + 64

        tmp2 = [[int(k.query_name.split("_")[1]), k] for k in reads2]
        tmp2.sort(key=lambda x: x[0])
        R2read = copy.deepcopy(tmp2[0][1])
        r2tmp = R2read.query_name.split("_")
        R2read.query_name = r2tmp[0]
        R2read.set_tag("BP", "{0}_{1}".format(r2tmp[1], r2tmp[2]))
        if len(r2tmp) > 3:
            linker.update(r2tmp[3].split(";"))
        R2read.flag = R2read.flag + 128

        if len(linker) > 0:
            R1read.set_tag("LT", "{0}".format(";".join(sorted(linker)))) # linker type
            R2read.set_tag("LT", "{0}".format(";".join(sorted(linker))))
        if r1tmp[1] == "1" and r2tmp[1] == "1":
            R1read.set_tag("FA", "R1R2_with_2End") # flag
            R2read.set_tag("FA", "R1R2_with_2End")
            R1read.set_tag("FB", "End") # flag
            R2read.set_tag("FB", "End")
            # fpairTn.write(R1read)
            # fpairTn.write(R2read)
            if R1read.has_tag("LT"):
                fpairTn1.write(R1read)
                fpairTn1.write(R2read)
            else:
                fpairTn2.write(R1read)
                fpairTn2.write(R2read)
            dataStat["R1R2_with_2End"] += 1
        elif r1tmp[1] != "1" and r2tmp[1] != "1":
            R1read.set_tag("FA", "R1R2_with_0End")
            R2read.set_tag("FA", "R1R2_with_0End")
            # fpairTn.write(R1read)
            # fpairTn.write(R2read)
            if R1read.has_tag("LT"):
                fpairTn1.write(R1read)
                fpairTn1.write(R2read)
            else:
                fpairTn2.write(R1read)
                fpairTn2.write(R2read)
            dataStat["R1R2_with_0End"] += 1
        else:
            R1read.set_tag("FA", "R1R2_with_1End")
            R2read.set_tag("FA", "R1R2_with_1End")
            if r1tmp[1] == "1":
                R1read.set_tag("FB", "End")
            else:
                R2read.set_tag("FB", "End")
            # fpairTn.write(R1read)
            # fpairTn.write(R2read)
            if R1read.has_tag("LT"):
                fpairTn1.write(R1read)
                fpairTn1.write(R2read)
            else:
                fpairTn2.write(R1read)
                fpairTn2.write(R2read)
            dataStat["R1R2_with_1End"] += 1
        dataStat["pair"] += 1
    elif len(reads1) > 0:
        dataStat["R1"] += 1
        if len(reads1) == 1:
            dataStat["R1_1tag"] += 1
            tmp = reads1[0].query_name.split("_")
            reads1[0].query_name = tmp[0]
            reads1[0].set_tag("BP", "{0}_{1}".format(tmp[1], tmp[2]))
            if len(tmp) > 3:
                reads1[0].set_tag("LT", tmp[3])
            if tmp[1] == "1":
                reads1[0].set_tag("FA", "R1_1tag_with_1End")
                reads1[0].set_tag("FB", "End")
                # fsingleTn.write(reads1[0])
                if reads1[0].has_tag("LT"):
                    fsingleTn1.write(reads1[0])
                else:
                    fsingleTn2.write(reads1[0])
                dataStat["R1_1tag_with_1End"] += 1
            else:
                reads1[0].set_tag("FA", "R1_1tag_with_0End")
                # fsingleTn.write(reads1[0])
                if reads1[0].has_tag("LT"):
                    fsingleTn1.write(reads1[0])
                else:
                    fsingleTn2.write(reads1[0])
                dataStat["R1_1tag_with_0End"] += 1
            dataStat["single"] += 1
        else:
            dataStat["R1_2tag"] += 1
            tmp = [[int(k.query_name.split("_")[1]), k] for k in reads1]
            tmp.sort(key=lambda x: x[0])
            R1read = copy.deepcopy(tmp[0][1])
            r1tmp = R1read.query_name.split("_")
            R1read.query_name = r1tmp[0]
            R1read.set_tag("BP", "{0}_{1}".format(r1tmp[1], r1tmp[2]))
            R1read.flag = R1read.flag + 64
            
            R2read = copy.deepcopy(tmp[-1][1])
            r2tmp = R2read.query_name.split("_")
            R2read.query_name = r2tmp[0]
            R2read.set_tag("BP", "{0}_{1}".format(r2tmp[1], r2tmp[2]))
            if R2read.is_reverse:
                R2read.flag = R2read.flag - 16
            else:
                R2read.flag = R2read.flag + 16
            R2read.flag = R2read.flag + 128

            if len(r1tmp) > 3:
                R1read.set_tag("LT", r1tmp[3])
                R2read.set_tag("LT", r1tmp[3])
            
            if r1tmp[1] == "1" and r2tmp[1] == "1":
                R1read.set_tag("FA", "R1_2tag_with_2End")
                R2read.set_tag("FA", "R1_2tag_with_2End")
                R1read.set_tag("FB", "End")
                R2read.set_tag("FB", "End")
                # fpairTn.write(R1read)
                # fpairTn.write(R2read)
                if R1read.has_tag("LT"):
                    fpairTn1.write(R1read)
                    fpairTn1.write(R2read)
                else:
                    fpairTn2.write(R1read)
                    fpairTn2.write(R2read)
                dataStat["R1_2tag_with_2End"] += 1
            elif r1tmp[1] != "1" and r2tmp[1] != "1":
                R1read.set_tag("FA", "R1_2tag_with_0End")
                R2read.set_tag("FA", "R1_2tag_with_0End")
                # fpairTn.write(R1read)
                # fpairTn.write(R2read)
                if R1read.has_tag("LT"):
                    fpairTn1.write(R1read)
                    fpairTn1.write(R2read)
                else:
                    fpairTn2.write(R1read)
                    fpairTn2.write(R2read)
                dataStat["R1_2tag_with_0End"] += 1
            else:
                R1read.set_tag("FA", "R1_2tag_with_1End")
                R2read.set_tag("FA", "R1_2tag_with_1End")
                if r1tmp[1] == "1":
                    R1read.set_tag("FB", "End")
                else:
                    R2read.set_tag("FB", "End")
                # fpairTn.write(R1read)
                # fpairTn.write(R2read)
                if R1read.has_tag("LT"):
                    fpairTn1.write(R1read)
                    fpairTn1.write(R2read)
                else:
                    fpairTn2.write(R1read)
                    fpairTn2.write(R2read)
                dataStat["R1_2tag_with_1End"] += 1
            dataStat["pair"] += 1
    elif len(reads2) > 0:
        dataStat["R2"] += 1
        if len(reads2) == 1:
            dataStat["R2_1tag"] += 1
            tmp = reads2[0].query_name.split("_")
            reads2[0].query_name = tmp[0]
            reads2[0].set_tag("BP", "{0}_{1}".format(tmp[1], tmp[2]))
            if len(tmp) > 3:
                reads2[0].set_tag("LT", tmp[3])
            if reads2[0].is_reverse:   # convert R2 to R1
                reads2[0].flag = reads2[0].flag - 16
            else:
                reads2[0].flag = reads2[0].flag + 16
            if tmp[1] == "1":
                reads2[0].set_tag("FA", "R2_1tag_with_1End")
                reads2[0].set_tag("FB", "End")
                # fsingleTn.write(reads2[0])
                if reads2[0].has_tag("LT"):
                    fsingleTn1.write(reads2[0])
                else:
                    fsingleTn2.write(reads2[0])
                dataStat["R2_1tag_with_1End"] += 1
            else:
                reads2[0].set_tag("FA", "R2_1tag_with_0End")
                # fsingleTn.write(reads2[0])
                if reads2[0].has_tag("LT"):
                    fsingleTn1.write(reads2[0])
                else:
                    fsingleTn2.write(reads2[0])
                dataStat["R2_1tag_with_0End"] += 1
            dataStat["single"] += 1
        else:
            dataStat["R2_2tag"] += 1
            tmp = [[int(k.query_name.split("_")[1]), k] for k in reads2]
            tmp.sort(key=lambda x: x[0])
            R1read = copy.deepcopy(tmp[-1][1])
            r1tmp = R1read.query_name.split("_")
            R1read.query_name = r1tmp[0]
            R1read.set_tag("BP", "{0}_{1}".format(r1tmp[1], r1tmp[2]))
            R1read.flag = R1read.flag + 64
            if R1read.is_reverse:
                R1read.flag = R1read.flag - 16
            else:
                R1read.flag = R1read.flag + 16
            
            R2read = copy.deepcopy(tmp[0][1])
            r2tmp = R2read.query_name.split("_")
            R2read.query_name = r2tmp[0]
            R2read.set_tag("BP", "{0}_{1}".format(r2tmp[1], r2tmp[2]))
            R2read.flag = R2read.flag + 128

            if len(r1tmp) > 3:
                R1read.set_tag("LT", r1tmp[3])
                R2read.set_tag("LT", r1tmp[3])
            
            if r1tmp[1] == "1" and r2tmp[1] == "1":
                R1read.set_tag("FA", "R2_2tag_with_2End")
                R2read.set_tag("FA", "R2_2tag_with_2End")
                R1read.set_tag("FB", "End")
                R2read.set_tag("FB", "End")
                # fpairTn.write(R1read)
                # fpairTn.write(R2read)
                if R1read.has_tag("LT"):
                    fpairTn1.write(R1read)
                    fpairTn1.write(R2read)
                else:
                    fpairTn2.write(R1read)
                    fpairTn2.write(R2read)
                dataStat["R2_2tag_with_2End"] += 1
            elif r1tmp[1] != "1" and r2tmp[1] != "1":
                R1read.set_tag("FA", "R2_2tag_with_0End")
                R2read.set_tag("FA", "R2_2tag_with_0End")
                # fpairTn.write(R1read)
                # fpairTn.write(R2read)
                if R1read.has_tag("LT"):
                    fpairTn1.write(R1read)
                    fpairTn1.write(R2read)
                else:
                    fpairTn2.write(R1read)
                    fpairTn2.write(R2read)
                dataStat["R2_2tag_with_0End"] += 1
            else:
                R1read.set_tag("FA", "R2_2tag_with_1End")
                R2read.set_tag("FA", "R2_2tag_with_1End")
                if r1tmp[1] == "1":
                    R1read.set_tag("FB", "End")
                else:
                    R2read.set_tag("FB", "End")
                # fpairTn.write(R1read)
                # fpairTn.write(R2read)
                if R1read.has_tag("LT"):
                    fpairTn1.write(R1read)
                    fpairTn1.write(R2read)
                else:
                    fpairTn2.write(R1read)
                    fpairTn2.write(R2read)
                dataStat["R2_2tag_with_1End"] += 1
            dataStat["pair"] += 1
    else: 
        dataStat["None"] += 1

# fsingleTn.close()
fsingleTn1.close()
fsingleTn2.close()
# fpairTn.close()
fpairTn1.close()
fpairTn2.close()

print("{0}\t{1}".format("R1R2", dataStat["R1R2"]))
print("{0}\t{1}".format("R1", dataStat["R1"]))
print("{0}\t{1}".format("R2", dataStat["R2"]))
print("{0}\t{1}".format("None", dataStat["None"]))

print("{0}\t{1}".format("R1_1tag", dataStat["R1_1tag"]))
print("{0}\t{1}".format("R1_2tag", dataStat["R1_2tag"]))
print("{0}\t{1}".format("R2_1tag", dataStat["R2_1tag"]))
print("{0}\t{1}".format("R2_2tag", dataStat["R2_2tag"]))

print("{0}\t{1}".format("R1R2_with_2End", dataStat["R1R2_with_2End"]))
print("{0}\t{1}".format("R1R2_with_1End", dataStat["R1R2_with_1End"]))
print("{0}\t{1}".format("R1R2_with_0End", dataStat["R1R2_with_0End"]))

print("{0}\t{1}".format("R1_1tag_with_1End", dataStat["R1_1tag_with_1End"]))
print("{0}\t{1}".format("R1_1tag_with_0End", dataStat["R1_1tag_with_0End"]))
print("{0}\t{1}".format("R1_2tag_with_2End", dataStat["R1_2tag_with_2End"]))
print("{0}\t{1}".format("R1_2tag_with_1End", dataStat["R1_2tag_with_1End"]))
print("{0}\t{1}".format("R1_2tag_with_0End", dataStat["R1_2tag_with_0End"]))

print("{0}\t{1}".format("R2_1tag_with_1End", dataStat["R2_1tag_with_1End"]))
print("{0}\t{1}".format("R2_1tag_with_0End", dataStat["R2_1tag_with_0End"]))
print("{0}\t{1}".format("R2_2tag_with_2End", dataStat["R2_2tag_with_2End"]))
print("{0}\t{1}".format("R2_2tag_with_1End", dataStat["R2_2tag_with_1End"]))
print("{0}\t{1}".format("R2_2tag_with_0End", dataStat["R2_2tag_with_0End"]))

print("{0}\t{1}".format("pair", dataStat["pair"]))
print("{0}\t{1}".format("single", dataStat["single"]))