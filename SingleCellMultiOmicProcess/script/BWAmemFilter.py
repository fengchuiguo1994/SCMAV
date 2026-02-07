import sys
import re

import gzip
def writeFile(outfile):
    if outfile.endswith((".gz","gzip")):
        fout = gzip.open(outfile,'wt')
    else:
        fout = open(outfile,'w')
    return fout

import pysam
def writeSam(outsamfile,header):
    if outsamfile.endswith(".bam"):
        outsam = pysam.AlignmentFile(outsamfile,'wb',header=header)
    elif outsamfile.endswith(".sam"):
        outsam = pysam.AlignmentFile(outsamfile,'w',header=header)
    else:
        raise ValueError("the output sam/bam file is not end with sam or bam!")
    return outsam

def detUniq(read):
    flag = False
    if read.get_tag("AS") - read.get_tag("XS") >= diff: # XA: multi-mapped, SA: (split reads)
        if read.has_tag("XA"):
            length1 = sum(int(num) for num in re.findall(r'(\d+)M', read.cigarstring))
            length2 = sum(int(num) for num in re.findall(r'(\d+)M', read.get_tag("XA").split(";")[0].split(",")[-2]))
            if length1 > length2:
                flag = True
            elif length1 == length2:
                if read.get_tag("NM") < int(read.get_tag("XA").split(";")[0].split(",")[-1]):
                    flag = True
        else:
            flag = True
    return flag

fin = pysam.AlignmentFile(sys.stdin)
fout = pysam.AlignmentFile(sys.stdout,'wb',header=fin.header)
diff = int(sys.argv[1])
mapq = int(sys.argv[2])
fstat = writeFile(sys.argv[3])

index = 0
for read1 in fin:
    read2 = next(fin)
    if read1.query_name != read2.query_name:
        print(read1.__str__())
        print(read2.__str__())
        sys.exit("error")
    if read1.mapping_quality < mapq or read2.mapping_quality < mapq:
        continue
    if detUniq(read1) and detUniq(read2):
        fout.write(read1)
        fout.write(read2)
        index += 1

fstat.write("{0}\t{1}\n".format("unique_mapping_reads", index))
fstat.close()
fin.close()
fout.close()

