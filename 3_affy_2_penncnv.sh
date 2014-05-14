# adapted from the PennCNV website
# A guide to extract allele-specific signal values from raw Affymetrix CEL files using Affymetrix Power Tools software.

# Software & Library requirements
# PennCNV software
# Reference File: hapmap.quant-norm.normalization-target.txt
# Affymetrix Power Tools (APT) software
# Affymetrix GenomeWideSNP 6 library files

#
# 1. Make a file containing a list of all CEL files + location
#

ls -d $PWD/*.CEL > CEL_FILES_LIST.txt
sed -i '1s/^/cel_files/' CEL_FILES_LIST.txt

#
# 2. CEL to Allele-Specific Signals
#
# The above below will extract signal intensity values for PM probes in all the CEL files specified in the listfile. The values are then quantile normalized,  median polished, generating A & B signal intensity values for each SNP. The filehapmap.quant-norm.normalization-target.txt is is the reference quantile distribution from the HapMap3 project.

APTPROBSUM="/../apt-1.15.1-x86_64-intel-linux/bin/apt-probeset-summarize"
CDF="/../GenomeWideSNP_6.cdf"
SKETCH="/../hapmap.quant-norm.normalization-target.txt"
OUTDIR="/../"
CELFILES="/../CEL_FILES_LIST.txt"
${APTPROBSUM} \
--cdf-file ${CDF} \
--analysis quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true \
--target-sketch=${SKETCH} \
--out-dir ${OUTDIR} \
--cel-files ${CELFILES}

#
# 3. Create BirdSeed files
#
# This step generates genotyping calls from the raw Affymetrix GenomeWideSNP_6 CEL files using the Birdseed algorithm.

APTPROBGEN="/../apt-1.15.1-x86_64-intel-linux/bin/apt-probeset-genotype"
CDF="/../GenomeWideSNP_6.cdf"
BSEED="/../GenomeWideSNP_6.birdseed.models"
SPECIAL="/../GenomeWideSNP_6.specialSNPs"
OUTDIR="/../"
CELFILES="/../CEL_FILES_LIST.txt"
${APTPROBGEN} \
-c ${CDF} -a birdseed \
--read-models-birdseed ${BSEED} \
--special-snps ${SPECIAL} \
--out-dir ${OUTDIR} \
--cel-files ${CELFILES}

#
# 4. create sex file from birdseed report
#
# The file_sex file is a two-column file that annotates the sex information for each CEL file. The birdseed.report.txt file that was generated from from the previous step can bes used to construct the sexfile based on computed gender information.

OUTDIR=”/../”
fgrep male ${OUTDIR}birdseed.report.txt | cut -f 1,2 > ${OUTDIR}gender_computed.txt

#
# 4. Generate canonical genotype clustering file
#
# The following step will generate canonical genotype clusters. The affygw6.hg18.pfb file contains the annotated marker positions in hg18 (NCBI 36) human genome assembly.

GENOCLUSTER="/../affy_geno_cluster.pl"
BSEEDCALL="/../birdseed.calls.txt"
BSEEDCONF="/../birdseed.confidences.txt"
QNTNRMSUM="/../quant-norm.pm-only.med-polish.expr.summary.txt"
PFB="/../affygw6.hg18.pfb"
SEXFILE="/../gender_computed.txt"
OUTDIR="/../"
${GENOCLUSTER} ${BSEEDCALL} ${BSEEDCONF} ${QNTNRMSUM} \
-locfile ${PFB} \
-sexfile ${SEXFILE} \
-out ${OUTDIR}gw6.genocluster

#
# 5. Calculate LRR and BAF
#
# Finally, we can extract allele-specific signal intensity measures to calculate the Log R Ratio (LRR) values and the B Allele Frequency (BAF) values for each marker in each individual. The cluster file gw6.genocluster is a file which specifies the chromosome position of each SNP or CN probe.

NORMAFFY="/../normalize_affy_geno_cluster.pl"
GENOCLUST="/../gw6.genocluster"
QNTNRMSUM="/../quant-norm.pm-only.med-polish.expr.summary.txt"
PFB="/../affygw6.hg18.pfb"
BAF_LRR_ALL="/../gw6.lrr_baf.txt"
${NORMAFFY} ${GENOCLUST} ${QNTNRMSUM} \
-locfile ${PFB}
-out ${BAF_LRR_ALL}
