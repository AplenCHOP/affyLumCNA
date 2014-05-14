#!/bin/bash
#$ -cwd

for x in *.txt;
do

# derive sample id from filename
# alternatively, extract sample id from another source
base=$(basename $x ".txt")

# rename first line to include sample id
sed -i "s/Name\tChr\tPosition\tB Allele Freq\tLog R Ratio/Name\tChr\tPosition\t$base.B Allele Freq\t$base.Log R Ratio/g" $x

# remove lines with entries containing chr 0,  XY, Y, and MT
sed -i -e "/\tXY\t/d" -e "/\tY\t/d" -e "/\tMT\t/d" -e "/\t0\t0\t/d" $x

# save to new file with extension .srtd
(head -n 1 $x && tail -n +2 $x | sort -k 2,2 -k 3,3n) > $x.srtd

done
