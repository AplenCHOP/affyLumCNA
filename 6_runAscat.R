4.	Running ASCAT
Submission Script
#!/bin/bash
#$ -cwd
# 1234 is filename prefix for the specific patient in this case

qsub -b y 'R CMD BATCH "--args patient_1234" /performAscatLum.R'

R script performAscatLum.R
args <- commandArgs(trailingOnly = TRUE)

# specify file location and names
dir <- "/myDir/cnvData/lumData/GCadjusted/splitData/"

ptlrr <- paste(dir, eval(parse(text = args[[1]])), "_PT.txt.adjusted.lrr", sep = "")
ptbaf <- paste(dir, eval(parse(text = args[[1]])), "_PT.txt.adjusted.baf", sep = "")
prlrr <- paste(dir, eval(parse(text = args[[1]])), "_PR.txt.adjusted.lrr", sep = "")
prbaf <- paste(dir, eval(parse(text = args[[1]])), "_PR.txt.adjusted.baf", sep = "")

# load in library
source("/myDir/Rscripts/ASCAT2.2/ascat.R")

# perform ASCAT_CNASTRUCT analysis
ascat.bc = ascat.loadData(ptlrr, ptbaf, prlrr, prbaf)
ascat.plotRawData(ascat.bc)
ascat.aspcf = ascat.aspcf(ascat.bc)
ascat.plotSegmentedData(ascat.aspcf)
ascat.output = ascat.runAscat(ascat.aspcf)
save.image(file = paste(eval(parse(text = args[[1]])), ".RData", sep = ""))

ascat.unmatched = ascat.loadData(ptlrr, ptbaf)
source("/myDir/Rscripts/ASCAT2.2/predictGG.R")
platform = "Illumina2.5M"
ascat.gg = ascat.predictGermlineGenotypes(ascat.unmatched, platform)
ascat.unmatched = ascat.aspcf(ascat.unmatched, ascat.gg = ascat.gg)
#ascat.plotSegmentedData(ascat.unmatched)
ascat.output.unmatched = ascat.runAscat(ascat.unmatched)
save.image(file = paste(eval(parse(text = args[[1]])), ".RData", sep = ""))
