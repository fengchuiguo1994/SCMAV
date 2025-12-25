infile=$1
outfile=$2
tmpdir=$3
threads=$4

SCRIPT_DIR=$(dirname "$0")
SCRIPT_DIR=$(realpath $SCRIPT_DIR)

outdir=$(dirname $(realpath $outfile))
outname=$(basename $(realpath $outfile))

python $SCRIPT_DIR/split4dnByChrom.py $infile $tmpdir/$outname
cat $tmpdir/$outname.header.list > $outfile
cat $tmpdir/$outname.chrom.list | while read i j; do
    sort --parallel=$threads -k3,3n -k5,5n $tmpdir/$outname.${i}-${j}.tmp.txt 
    rm $tmpdir/$outname.${i}-${j}.tmp.txt
done >> $outfile
