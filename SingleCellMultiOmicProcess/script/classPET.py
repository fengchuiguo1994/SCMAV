import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib as mpl 
mpl.rcParams['pdf.fonttype'] = 42 
mpl.rcParams['ps.fonttype'] = 42

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

fin = readSam(sys.argv[1])
distance=int(sys.argv[2]) # distance
fsingle = writeSam(sys.argv[3],header=fin.header) # self-lagation PET
fother = writeSam(sys.argv[4],header=fin.header) # other PET
fpet = writeSam(sys.argv[5],header=fin.header) # interaction PET
flog = writeFile(sys.argv[6])
singlenum = 0
othernum = 0
petnum = 0
greatdishash = {}
greatdistotal = 0
lessdishash = {}
lessdistotal = 0

for read1 in fin:
    read2 = next(fin)
    if read1.query_name != read2.query_name:
        sys.stderr.write(read1.__str__())
        sys.stderr.write(read2.__str__())
        sys.stderr.write("reads name are not match")
        sys.exit()
    if read1.reference_name == read2.reference_name:
        absisize = abs(read1.isize)
        if abs(read1.isize) >= distance and abs(read2.isize) >= distance:
            if absisize not in greatdishash:
                greatdishash[absisize] = {}
                greatdishash[absisize]["TANDEM"] = 0
                greatdishash[absisize]["FR"] = 0
                greatdishash[absisize]["RF"] = 0
                greatdishash[absisize]["total"] = 0
            greatdishash[absisize]["total"] += 1
            greatdistotal += 1
            fpet.write(read1)
            fpet.write(read2)
            petnum += 1
            if read1.is_reverse == read2.is_reverse: # TANDEM
                greatdishash[absisize]["TANDEM"] += 1
            else:
                if (read1.is_reverse == False) and (read2.is_reverse == True):
                    if read1.reference_start < read2.reference_end: # FR
                        greatdishash[absisize]["FR"] += 1
                    else:
                        greatdishash[absisize]["RF"] += 1
                else:
                    if read2.reference_start < read1.reference_end: # RF
                        greatdishash[absisize]["FR"] += 1
                    else:
                        greatdishash[absisize]["RF"] += 1
        else:
            if absisize not in lessdishash:
                lessdishash[absisize] = {}
                lessdishash[absisize]["TANDEM"] = 0
                lessdishash[absisize]["FR"] = 0
                lessdishash[absisize]["RF"] = 0
                lessdishash[absisize]["total"] = 0
            lessdishash[absisize]["total"] += 1
            lessdistotal += 1
            if read1.is_reverse == read2.is_reverse: # TANDEM
                fother.write(read1)
                fother.write(read2)
                othernum += 1
                lessdishash[absisize]["TANDEM"] += 1
            else:
                if (read1.is_reverse == False) and (read2.is_reverse == True):
                    if read1.reference_start < read2.reference_end: # FR 
                        if not read1.is_proper_pair:
                            read1.flag += 2
                        if not read2.is_proper_pair:
                            read2.flag += 2
                        fsingle.write(read1)
                        fsingle.write(read2)
                        singlenum += 1
                        lessdishash[absisize]["FR"] += 1
                    else:
                        fother.write(read1)
                        fother.write(read2)
                        othernum += 1
                        lessdishash[absisize]["RF"] += 1
                elif (read1.is_reverse == True) and (read2.is_reverse == False):
                    if read2.reference_start < read1.reference_end: # FR 
                        if not read1.is_proper_pair:
                            read1.flag += 2
                        if not read2.is_proper_pair:
                            read2.flag += 2
                        fsingle.write(read1)
                        fsingle.write(read2)
                        singlenum += 1
                        lessdishash[absisize]["FR"] += 1
                    else:
                        fother.write(read1)
                        fother.write(read2)
                        othernum += 1
                        lessdishash[absisize]["RF"] += 1
                else:
                    sys.stderr.write(read1.__str__())
                    sys.stderr.write(read2.__str__())
                    sys.exit("error1")
    else:
        fpet.write(read1)
        fpet.write(read2)
        petnum += 1

fsingle.close()
fpet.close()
fother.close()
flog.write("less_{0}_FR\t{1}\n".format(distance, singlenum))
flog.write("less_{0}_TANDEM_RF\t{1}\n".format(distance, othernum))
flog.write("great_{0}\t{1}\n".format(distance, petnum))
flog.close()

print("type\tinsertsize\tFR\tRF\tTANDEM\ttotal")
greatcutoff = greatdistotal * 0.8
greatcurrentsum = 0
greatcurrent = None
for k in sorted(greatdishash.keys()):
    print("ge{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(distance, k, greatdishash[k]["FR"], greatdishash[k]["RF"], greatdishash[k]["TANDEM"], greatdishash[k]["total"]))
    greatcurrentsum += greatdishash[k]["total"]
    if greatcurrentsum > greatcutoff and greatcurrent == None:
        greatcurrent = k

lesscutoff = lessdistotal * 0.95
lesscurrentsum = 0
lesscurrent = None
for k in sorted(lessdishash.keys()):
    print("lt{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(distance, k, lessdishash[k]["FR"], lessdishash[k]["RF"], lessdishash[k]["TANDEM"], lessdishash[k]["total"]))
    lesscurrentsum += lessdishash[k]["total"]
    if lesscurrentsum > lesscutoff and lesscurrent == None:
        lesscurrent = k

# print(greatcurrent)
# print(lesscurrent)
df = pd.DataFrame.from_dict(lessdishash, orient='index')
df = df.loc[df.index < lesscurrent]
df.sort_index(inplace=True)
# print(df)

pdf_path = "{0}.insertsize.pdf".format(sys.argv[6])
pdf = PdfPages(pdf_path)
colormap = {"FR": "red", "RF": "blue", "TANDEM": "orange"}

for col in ["FR", "RF", "TANDEM"]:
    if col == "total":
        continue
    plt.figure(figsize=(8, 4))
    plt.plot(df[col], color=colormap[col], label=col, linestyle='-', linewidth=0.5)
    plt.title(f'Insert Size (<1000) {col} Line plot', fontsize=14)
    plt.xlabel('Insert Size', fontsize=12)
    plt.ylabel('Count', fontsize=12)
    plt.legend()
    pdf.savefig(bbox_inches='tight')
    plt.close()

df = pd.DataFrame.from_dict(greatdishash, orient='index')
df = df.loc[df.index < greatcurrent]
df.sort_index(inplace=True)
# print(df)
plt.figure(figsize=(8, 4))
for col in ["FR", "RF", "TANDEM"]:
    plt.plot(df[col], color=colormap[col], label=col, linestyle='-', linewidth=1, alpha=0.1)
plt.title(f'Insert Size Line (>1000) plot', fontsize=14)
plt.xlabel('Insert Size', fontsize=12)
plt.ylabel('Count', fontsize=12)
plt.legend()
pdf.savefig(bbox_inches='tight')
plt.close()

pdf.close()
