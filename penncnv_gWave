# Genomic Wave Correction
# One of the quality assurance steps in copy number analyses involves detecting samples with genomic waves in their log ratio data. Genomic waves are phenomena whereby, despite normalization, the log ratio data appear to have a long-range wave pattern when plotted in genomic space. Waviness is hypothesized to be correlated with the GC content of the probes themselves in addition to the GC content of the region around the probes. This wave pattern wreaks havoc on copy number detection algorithms. There are some cases (tumor data) where GC content correction may introduce more bias than clean the data from it’s waviness. Therefore, we simply remove samples with extreme wave factors.

# PennCNV has a procedure for adjusting signal intensity values, without generating CNV calls. The genomic_wave.pl program in PennCNV package can be used to adjust signal intensity values, by employing the wave correction algorithm from Diskin, et al 2008. The input file must contain a field in the header line that says “Log R Ratio”.

#!/bin/bash
#$-cwd

perl -w genomic_wave.pl --listfile data/my_study_samples.txt --adjust --gcmodelfile lib/illumina_omni2.5_custom.hg19.gcModel

# genomic_wave.pl: is the file supplied by the PennCNV program,
# listfile: is a textfile that contains the filenames of all the samples that need to be evaluated.
# adjust: can be used to generate a new file with updated Log R Ratio measures. The file will have a new extension, labeled .adjusted.
# gcmodelfile: denotes the reference file containing quantitative GC wave values per SNP. These gcmodel files are platform-specific, and many can be downloaded from the PennCNV file. However, because we are using the Illumina 2.5M Omni platform, the gcmodel file is not provided, and must be created. You can email for the Illumina 2.5M Omni gcmodel file.
