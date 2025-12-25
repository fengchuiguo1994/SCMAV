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

mydict = {}
fin = readFile(sys.argv[1])
fout = writeFile(sys.argv[2])
bins = int(sys.argv[3]) * 1000
ex = 2*bins
libraryname = sys.argv[4]
for line in fin:
    tmp = line.strip().split()
    # barcode = tmp[6].split(":",maxsplit=1)[0]
    barcode = tmp[10]
    if barcode not in mydict:
        mydict[barcode] = {}
    if tmp[0] in ["ChrMt","ChrM","chrM","chrMt","chrm","chrmt","MT","M"]:
        continue
    if tmp[3] in ["ChrMt","ChrM","chrM","chrMt","chrm","chrmt","MT","M"]:
        continue
    if tmp[0] not in mydict[barcode]:
        mydict[barcode][tmp[0]] = {}
    if tmp[3] not in mydict[barcode][tmp[0]]:
        mydict[barcode][tmp[0]][tmp[3]] = {}
    pos1 = (int(tmp[1]) + int(tmp[2])) // ex * bins
    pos2 = (int(tmp[4]) + int(tmp[5])) // ex * bins
    if tmp[0] == tmp[3]:
        if pos1 > pos2:
            pos1, pos2 = pos2, pos1
    if pos1 not in mydict[barcode][tmp[0]][tmp[3]]:
        mydict[barcode][tmp[0]][tmp[3]][pos1] = {}
    if pos2 not in mydict[barcode][tmp[0]][tmp[3]][pos1]:
        mydict[barcode][tmp[0]][tmp[3]][pos1][pos2] = 0 
    mydict[barcode][tmp[0]][tmp[3]][pos1][pos2] += 1
fin.close()

fout.write("library\tresolution\tcell_id\tcell_name\tchrom1\tpos1\tchrom2\tpos2\tcount\n")
num = 0
for barcode in mydict:
    for chrom1 in sorted(mydict[barcode].keys()):
        for chrom2 in sorted(mydict[barcode][chrom1].keys()):
            for pos1 in sorted(mydict[barcode][chrom1][chrom2].keys()):
                for pos2 in sorted(mydict[barcode][chrom1][chrom2][pos1].keys()):
                    fout.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(libraryname,bins,barcode,barcode,chrom1,pos1,chrom2,pos2,mydict[barcode][chrom1][chrom2][pos1][pos2]))
    num += 1
fout.close()
