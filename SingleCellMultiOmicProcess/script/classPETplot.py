import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib as mpl 
mpl.rcParams['pdf.fonttype'] = 42 
mpl.rcParams['ps.fonttype'] = 42
import matplotlib
matplotlib.use('agg')

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

colormap = {"FR": "red", "RF": "blue", "TANDEM": "orange"}
pdf = PdfPages("{0}.pdf".format(sys.argv[1]))

myhash = {}
fin = readFile(sys.argv[1])
for line in fin:
    tmp = line.strip().split()
    if tmp[0] == "type" or tmp[0] == "lt1000":
        continue
    nn = 10**(len(tmp[1])-2)
    size = int(tmp[1])//nn*nn
    if size not in myhash:
        myhash[size] = {}
        myhash[size]["TANDEM"] = 0
        myhash[size]["FR"] = 0
        myhash[size]["RF"] = 0
        myhash[size]["total"] = 0
        myhash[size]["step"] = nn
    myhash[size]["FR"] += int(tmp[2])
    myhash[size]["RF"] += int(tmp[3])
    myhash[size]["TANDEM"] += int(tmp[4])
    myhash[size]["total"] += int(tmp[5])
for k in sorted(myhash.keys()):
    # print("ge{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(1000, k, myhash[k]["FR"], myhash[k]["RF"], myhash[k]["TANDEM"], myhash[k]["total"], myhash[k]["step"]))
    pass
df = pd.DataFrame.from_dict(myhash, orient='index')
df.sort_index(inplace=True)
# print(df)
plt.figure(figsize=(8, 4))
for col in ["FR", "RF", "TANDEM"]:
    plt.plot(np.log2(df[col]/df["step"]+1), color=colormap[col], label=col, linestyle='-', linewidth=0.2, alpha=0.5)
plt.title(f'Insert Size Line (>1000) plot (smooth)', fontsize=14)
plt.xlabel('Insert Size', fontsize=12)
plt.ylabel('log2(Count)', fontsize=12)
plt.legend()
pdf.savefig(bbox_inches='tight')
plt.close()

dfall = pd.read_csv(sys.argv[1], sep="\t")

df = dfall[dfall["type"] == "ge1000"]
plt.figure(figsize=(8, 4))
for col in ["FR", "RF", "TANDEM"]:
    plt.plot(df["insertsize"], np.log2(df[col]+1), color=colormap[col], label=col, linestyle='-', linewidth=1, alpha=0.1)
plt.title(f'Insert Size Line (>1000) plot', fontsize=14)
plt.xlabel('Insert Size', fontsize=12)
plt.ylabel('log2(Count)', fontsize=12)
plt.legend()
pdf.savefig(bbox_inches='tight')
plt.close()

df = dfall[dfall["type"] == "lt1000"]
plt.figure(figsize=(8, 4))
for col in ["FR", "RF", "TANDEM"]:
    plt.plot(df["insertsize"], df[col], color=colormap[col], label=col, linestyle='-', linewidth=1, alpha=1)
plt.title(f'Insert Size Line (<1000) plot', fontsize=14)
plt.xlabel('Insert Size', fontsize=12)
plt.ylabel('Count', fontsize=12)
plt.legend()
pdf.savefig(bbox_inches='tight')
plt.close()

df = dfall[dfall["type"] == "lt1000"]
for col in ["FR", "RF", "TANDEM"]:
    plt.figure(figsize=(8, 4))
    plt.plot(df["insertsize"], df[col], color=colormap[col], label=col, linestyle='-', linewidth=0.5)
    plt.title(f'Insert Size (<1000) {col} Line plot', fontsize=14)
    plt.xlabel('Insert Size', fontsize=12)
    plt.ylabel('Count', fontsize=12)
    plt.legend()
    pdf.savefig(bbox_inches='tight')
    plt.close()
pdf.close()