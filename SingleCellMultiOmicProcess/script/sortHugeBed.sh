infile=$1
outfile=$2
tmpdir=$3
threads=$4

SCRIPT_DIR=$(dirname "$0")
SCRIPT_DIR=$(realpath $SCRIPT_DIR)

outdir=$(dirname $(realpath $outfile))
outname=$(basename $(realpath $outfile))

python $SCRIPT_DIR/splitBedByChrom.py $infile $tmpdir/$outname
cat $tmpdir/$outname.chrom.list | while read i; do
    sort --parallel=$threads -k2,2n -k3,3n $tmpdir/$outname.${i}.tmp.txt 
    rm $tmpdir/$outname.${i}.tmp.txt
done > $outfile
