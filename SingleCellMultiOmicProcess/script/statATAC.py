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
    fin = readFile(sys.argv[1])
    total = int(fin.readline().strip().split()[1])
    print("total_reads\t{0}".format(total))
    fin.close()
    
    fin = readFile(sys.argv[2])
    barcode = int(fin.readline().strip().split()[1])
    print("barcode_reads_ratio\t{0:.4f}".format(barcode/total))
    fin.close()
    
    fin = readFile(sys.argv[3])
    unique = int(fin.readline().strip().split()[1])
    print("unique_reads_ratio\t{0:.4f}".format(unique/total))
    fin.close()
    
    fin = readFile(sys.argv[4])
    mapr = int(fin.readline().strip().split()[1])
    fin.close()
    fin = readFile(sys.argv[5])
    rmpcr = int(fin.readline().strip().split()[1])
    print("PCR_ratio\t{0:.4f}".format(rmpcr/mapr))
    fin.close()
    
    fin = readFile(sys.argv[6])
    fragment = int(fin.readline().strip().split()[1])
    print("fragment\t{0}".format(fragment))
    fin.close()
    
    print("peak\t{0}".format(countline(sys.argv[7])))

    myhash = dict()
    fin = readFile(sys.argv[8])
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

if sys.argv[-1] == "chromap":
    fin = readFile(sys.argv[1])
    for line in fin:
        if line.startswith("Number of reads:"):
            num = line.strip().split()[-1]
            total = int(num.replace(".",""))//2
        if line.startswith("Number of barcodes in whitelist:"):
            num = line.strip().split()[-1]
            barcode = int(num.replace(".",""))
        if line.startswith("Number of uniquely mapped reads:"):
            num = line.strip().split()[-1]
            unique = int(num.replace(".",""))//2
    print("total_reads\t{0}".format(total))
    print("barcode_reads_ratio\t{0:.4f}".format(barcode/total))
    print("unique_reads_ratio\t{0:.4f}".format(unique/total))
    print("PCR_ratio\tNA")
    
    fin = readFile(sys.argv[2])
    fragment = int(fin.readline().strip().split()[1])
    print("fragment\t{0}".format(fragment))
    fin.close()
    
    print("peak\t{0}".format(countline(sys.argv[3])))

    myhash = dict()
    fin = readFile(sys.argv[4])
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

if sys.argv[-1] == "10X":
    data = json.load(readFile(sys.argv[1]))
    print("total_reads\t{0}".format(data["num_fragments"]))
    print("barcode_reads_ratio\t{0:.4f}".format(data["frac_valid_barcode"]))
    print("unique_reads_ratio\t{0:.4f}".format(data["frac_reads_mapq_gt_30"]))
    print("PCR_ratio\t{0:.4f}".format(data["frac_pcr_dup_passed_filters"]))
    print("fragment\t{0}".format(data["total_passed_filter_fragments"]))
    print("peak\t{0}".format(data["total_peaks_detected"]))
    
    myhash = dict()
    fin = readFile(sys.argv[2])
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