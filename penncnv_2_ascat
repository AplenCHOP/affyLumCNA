# PennCNV to ASCAT 
# After you have obtained GC adjusted signal intensity files, we need to split the file into two, one denoting the B-Allele Frequencies, and one with the Log R Ratios. This is the format which ASCAT uses to generate segments and allele-specific copy number calls using a similar procedure as for Illumina arrays.

# Make header file
head -1 gw6.lrr_baf.txt > gw6.headers.txt

# Split_baf.sh
#!/bin/bash
#$-cwd

PAT=$1;

INFILE="/../gw6.lrr_baf.txt"
HEADER="/../gw6.headers.txt"

PAT="`cut -f${PAT} $HEADER`";
PAT="`echo $PID | tr -d ' '`";

BAF=`basename $PAT .CEL.BAlleleFreq`;
cut -f1,2,3,$i $INFILE > $BAF.baf;
sed -i "s/Name\tChr\tPosition\t$PAT/\tchrs\tpos\t$BAF/g" $BAF.baf;

LRR=`basename $PAT .CEL.LogRRatio`;
cut -f1,2,3,$i $INFILE > $LRR.lrr;
sed -i "s/Name\tChr\tPosition\t$PAT/\tchrs\tpos\t$LRR/g" $LRR.lrr;
