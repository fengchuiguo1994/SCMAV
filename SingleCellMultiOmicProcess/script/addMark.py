import sys
import re

import pysam
def writeSam(outsamfile,header):
    if outsamfile.endswith(".bam"):
        outsam = pysam.AlignmentFile(outsamfile,'wb',header=header)
    elif outsamfile.endswith(".sam"):
        outsam = pysam.AlignmentFile(outsamfile,'w',header=header)
    else:
        raise ValueError("the output sam/bam file is not end with sam or bam!")
    return outsam

fin = pysam.AlignmentFile(sys.stdin)
fout = pysam.AlignmentFile(sys.stdout,'wb',header=fin.header)
readstrands = sys.argv[1] # R1 or R2
software = sys.argv[2] # mem or aln
if len(sys.argv) > 3:
    diff = int(sys.argv[3])
    if diff < 1:
        diff = 1
else:
    diff = 1

for read in fin:
    if software == "mem":
        read.set_tag("SC", "mem")
        if read.get_tag("AS") - read.get_tag("XS") >= diff: # XA: multi-mapped, SA: (split reads)
            if read.has_tag("XA"):
                length1 = sum(int(num) for num in re.findall(r'(\d+)M', read.cigarstring))
                length2 = sum(int(num) for num in re.findall(r'(\d+)M', read.get_tag("XA").split(";")[0].split(",")[-2]))

                if length1 > length2:
                    read.set_tag("XT", 'U', value_type="A")
                else:
                    if read.get_tag("NM") < int(read.get_tag("XA").split(";")[0].split(",")[-1]):
                        read.set_tag("XT", 'U', value_type="A")
            else:
                read.set_tag("XT", 'U', value_type="A")
    else:
        read.set_tag("SC", "aln")
        if read.is_unmapped:
            read.set_tag("XT", 'N', value_type="A")
    if readstrands == "R1":
        read.set_tag("SR", 'R1')
    elif readstrands == "R2":
        read.set_tag("SR", 'R2')
    else:
        sys.stderr.write(read.__str__())
        sys.exit("error reads")
    fout.write(read)
fin.close()
fout.close()