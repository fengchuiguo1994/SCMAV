import sys
import json
import statistics

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

def countline(infile):
    linecount = 0
    fin = readFile(infile)
    for line in fin:
        linecount += 1
    fin.close()
    return linecount

if sys.argv[-1] == "custom":
    fin = readFile(sys.argv[1]) # total
    total = int(fin.readline().strip().split()[1])
    print("Total_Reads\t{0}".format(total))
    fin.close()
    
    fin = readFile(sys.argv[2]) # barcode
    barcode = int(fin.readline().strip().split()[1])
    print("Barcode_Reads\t{0}".format(barcode))
    print("Barcode_Reads_ratio\t{0:.4f}".format(barcode/total))
    fin.close()
    
    fin = readFile(sys.argv[3]) # trim
    trim = int(fin.readline().strip().split()[1])
    print("Trim_Reads\t{0}".format(trim))
    print("Trim_Reads_ratio\t{0:.4f}".format(trim/total))
    fin.close()

    fin = readFile(sys.argv[4]) # linker
    linker = int(fin.readline().strip().split()[1])
    print("Linker_Reads\t{0}".format(linker))
    print("No_linker_reads\t{0}".format(trim - linker))
    fin.close()

    fin = readFile(sys.argv[5]) # uniq
    linkerunique = int(fin.readline().strip().split()[1])
    print("unique_mapping_reads_with_linker\t{0}".format(linkerunique))
    fin.close()
    fin = readFile(sys.argv[6]) # uniq
    nolinkerunique = int(fin.readline().strip().split()[1])
    print("unique_mapping_reads_without_linker\t{0}".format(nolinkerunique))
    fin.close()
    
    fin = readFile(sys.argv[7]) # uniq addbarcode
    linkeruniqueadd = int(fin.readline().strip().split()[1])
    print("barcode_unique_mapping_reads_with_linker\t{0}".format(linkeruniqueadd))
    fin.close()
    fin = readFile(sys.argv[8]) # uniq addbarcode
    nolinkeruniqueadd = int(fin.readline().strip().split()[1])
    print("barcode_unique_mapping_reads_without_linker\t{0}".format(nolinkeruniqueadd))
    fin.close()

    fin = readFile(sys.argv[9]) # uniq addbarcode rmPCR
    linkeruniqueaddrmPCR = int(fin.readline().strip().split()[1])
    print("barcode_unique_mapping_reads_with_linker_rmPCR\t{0}".format(linkeruniqueadd - linkeruniqueaddrmPCR))
    fin.close()
    fin = readFile(sys.argv[10]) # uniq addbarcode rmPCR
    nolinkeruniqueaddrmPCR = int(fin.readline().strip().split()[1])
    print("barcode_unique_mapping_reads_without_linker_rmPCR\t{0}".format(nolinkeruniqueadd - nolinkeruniqueaddrmPCR))
    fin.close()
    print("PCR_ratio\t{0:.4f}".format((linkeruniqueaddrmPCR+nolinkeruniqueaddrmPCR)/(linkeruniqueadd+nolinkeruniqueadd)))

    fin = readFile(sys.argv[11])
    fragment = int(fin.readline().strip().split()[1])
    print("fragment\t{0}".format(fragment))
    fin.close()
    
    print("peak\t{0}".format(countline(sys.argv[12])))

    myhash = dict()
    fin = readFile(sys.argv[13])
    fin.readline()
    fin.readline()
    aa = fin.readline().split()
    for line in fin:
        tmp = line.strip().split()
        if tmp[1] not in myhash:
            myhash[tmp[1]] = 0
        myhash[tmp[1]] += int(tmp[2])
    print("detect_cells\t{0}".format(aa[1]))
    print("median_cut_fragments_peak_per_cell\t{0}".format(statistics.median([myhash[x] for x in myhash.keys()])))
    fin.close()
