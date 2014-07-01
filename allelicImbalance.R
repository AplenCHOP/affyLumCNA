#Now we call allelic imbalancase based on the output from ASCAT. 
# Any time the two alleles are not represented in equal amounts (balanced), that is defined as AI

source("/ASCAT2.1/pipeline/AI/originalBirkbak.R")

# 1. DO NOT affy and illumina segments
# They come from a different genomic build (hg 18 vs h19)
#affy.ascat.seg<-read.table(, sep="\t", T)
#lum.ascat.seg<-read.table(, sep="\t", T)
ascat.seg <- read.table("/ASCAT2.1/manrev/affyLum.AI.txt", sep = "\t", stringsAsFactors = F, T)
colnames(ascat.seg)[1] = "SampleID"
colnames(ascat.seg)[2] = "Chr"
colnames(ascat.seg)[3] = "Start"
colnames(ascat.seg)[4] = "End"
colnames(ascat.seg)[5] = "nProbes"
colnames(ascat.seg)[6] = "nA"
colnames(ascat.seg)[7] = "nB"

# remove data from X 
ascat.seg <- ascat.seg[ascat.seg$Chr != "X", ] 

# The Annotation Table
hg18 <- read.table("/ASCAT2.1/lib/chrominfo.hg18.txt",sep="\t",T)
hg19 <- read.table("/ASCAT2.1/lib/chrominfo.hg19.txt",sep="\t",T)
hg19 <- hg19[hg19$chrom != "X", ]
hg19 <- hg19[hg19$chrom != "Y", ]


# 2. Define Allelic Imbalance
AI <- c(0, 1)[match(ascat.seg[, 6] == ascat.seg[, 7], c("TRUE", "FALSE"))]
ascat.seg <- cbind(ascat.seg[, 1:5], AI, ascat.seg[, 6:7])
write.table(ascat.seg, "affylum.AI.overview.txt", sep = "\t", row.names = F, col.names = T, quote = F)
#head(ascat.seg)

# Telomere determination
minProbes <- 0
tel.ascat.seg <- no.tel(ascat.seg, chrominfo = hg19, min.probes = minProbes)
summary(tel.ascat.seg)
write.table(tel.ascat.seg, "affylum.AI.tel.txt", sep = "\t", row.names = F, col.names = T, quote = F)


tel.cna.ascat.seg <- no.tel.cna(ascat.seg, chrominfo = hg19, min.probes = minProbes)
summary(tel.cna.ascat.seg)
write.table(tel.cna.ascat.seg, "affylum.AI.tel.cna.txt", sep = "\t", row.names = F, col.names = T, quote = F)

no.cna.ascat.seg <- gaps(ascat.seg, min.probes = 50000)
summary(no.cna.ascat.seg)
#write.table(no.cna.ascat.seg, "affylum.AI.no.cna.txt", sep = "\t", row.names = F, col.names = T, quote = F)


# Now load the patient response data using the Miller-Payne grade
load("aml.mrd.RData")
aml.mrd
# output should look something like this
#P1 P3 P4 P5 P6 P7 P8 P9 P10 P11 P12 P13 P14 P15 P16 P17 P18 P20 P21 P22
#3 5 1 5 1 4 4 5 4 3 1 1 2 0 1 5 5 2 0 2

# Now load the patient germline FLT3 mutation status (1 means mutation)
load("aml.flt3.RData")
aaml.flt3
# output should look something like this
#P1 P3 P4 P5 P6 P7 P8 P9 P10 P11 P12 P13 P14 P15 P16 P17 P18 P20 P21 P22
#0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0

#Make a binary response vector for MRD, based on at least 99% remission rates
aml.mrdbin <- c(0, 1)[match(aml.mrd >= 0.1, c("FALSE", "TRUE"))]
names(aml.mrdbin) <- names(aml.mrd)
aml.mrdbin

# Fig 2A: plot the differences between mrd+ and mrd-
boxplot(no.tel.dfci1[dfci1.cont >= cont.max, 1] ~ dfci1.mpbin[dfci1.cont >= cont.max], ylab = "NtAI", xlab = "", notch = F, names = c("No response", "Response"), main = "Figure 2a", ylim = c(0, 35))
beeswarm(no.tel.dfci1[dfci1.cont >= cont.max, 1] ~ dfci1.mpbin[dfci1.cont >= cont.max], col = c(1, 2), add = T, method = "center", pwpch = c(16,17)[match(dfci1.brca[dfci1.cont >= cont.max], c(0, 1))])
legend("topleft", legend = paste("P =", signif(wilcox.test(no.tel.dfci1[dfci1.cont >= cont.max, 1] ~ dfci1.mpbin[dfci1.cont >= cont.max])$p.value, 3)), bty = "n", pch = -1)

# Calculate the ROC curve for figure 2b, and identify the optimum number of NtAI
# for separation between responders on non-responder.
# Then determine the area under the curve and the p-value for association between
# sensitivity to cisplatin and high numbers of NtAI.
R1 <- rocdemo.sca(dfci1.mpbin[dfci1.cont >= cont.max], no.tel.dfci1[dfci1.cont >= cont.max, 1])
cut.indx <- which(sqrt((1 - R1@sens)^2 + (1 - R1@spec)^2) == min(sqrt((1 - R1@sens)^2 + (1 - R1@spec)^2)))
optimum.cut <- R1@cuts[cut.indx]
optimum.cut
AUC(R1)
wilcox.test(no.tel.dfci1[dfci1.cont >= cont.max, 1] ~ dfci1.mpbin[dfci1.cont >= cont.max])

# Fig 2B
# Plot figure 2b
plot(R1, ylab = "Sensitivity", xlab = "1-specificity", lty = 1, col = 1, main = "Figure 2b", lwd = 2)
lines(c(0, 1), c(0, 1), lty = 2)
legend("bottomright", legend = paste("AUC = ", signif(AUC(R1), 2), sep = ""), pch = c(-1), col = c(-1))

# Calculate AUC confidence intervals for figure 2B
ci.auc.cis1 <- ci.auc(no.tel.dfci1[dfci1.cont >= cont.max, 1], dfci1.mpbin[dfci1.cont >= cont.max]) 
ci.auc.cis1


# 3 FIGURES
#6.2 Supplementary figure 5
# Plot supplementary figure 5

affy.telAI.seg <- telAI.fn(affy.ascat.seg, chrominfo = hg18, min.probes = min.probes.affy)
lum.telAI.seg <- telAI.fn(lum.ascat.seg, chrominfo = hg19, min.probes = min.probes.lum)
affy.names <- affy.mrd[affy.cont >= cont.max]
names(affy.names) <- names(affy.cont[affy.cont >= cont.max])
lum.names <- lum.mrd[lum.cont >= cont.max]
names(lum.names) <- names(lum.cont[lum.cont >= cont.max])
comb.tel.seg <- rbind(affy.telAI.seg, lum.telAI.seg)
comb.tel.seg <- comb.tel.seg[comb.tel.seg[, 1] %in% names(c(affy.names, lum.names)), ]
telAI.breaks <- matrix(0, 0, 4)
for (i in 1:nrow(comb.tel.seg)) {
	a <- c(2, 4)[match(comb.tel.seg[i, 1] %in% names(affy.names), c("TRUE", "FALSE"))]
	b <- c(0, 2)[match(comb.tel.seg[i, 1] %in% names(c(affy.names, lum.names)[c(affy.names, lum.names) %in% c(4:5)]), c("FALSE", "TRUE"))]
	telAI.breaks <- rbind(telAI.breaks, c(as.numeric(comb.tel.seg[i, 2]), NA, a, b))
	if (as.numeric(comb.tel.seg[i, 3]) < hg19[as.numeric(comb.tel.seg[i, 2]), 3] * 1000) {
		telAI.breaks[i, 2] <- as.numeric(comb.tel.seg[i, 4])
	}
	if (as.numeric(comb.tel.seg[i, 3]) > hg19[as.numeric(comb.tel.seg[i, 2]), 3] * 1000) {
		telAI.breaks[i, 2] <- as.numeric(comb.tel.seg[i, 3])
	}
}

CHR <- telAI.breaks[, 1]
CHR[CHR == 23] <- "X"
MapInfo <- telAI.breaks[, 2]
rownames(telAI.breaks) <- 1:nrow(telAI.breaks)


# Now Plot
par(las = 1)
chrompos <- prepareGenomePlot(data.frame(CHR, MapInfo), paintCytobands = TRUE,  organism = "hsa", sexChromosomes = T)
a <- beeswarm(chrompos[, 2] ~ chrompos[, 1], horizontal = T, vertical = F, do.plot = F, method = "square", pch = 22, pwbg = telAI.breaks[, 4], add = T, cex = 0.4)
a <- apply(a, 2, as.numeric)
rownames(a) <- 1:nrow(a)
bb <- sort(as.numeric(names(table(as.numeric(a[, 1]) - as.numeric(a[, 6])))))
bb <- bb[bb > 0][1]
a[, 1] <- a[, 1] + bb
for (i in unique(a[, 6])) {
        b <- a[a[, 6] == i, ]
        if (nrow(b[b[, 1] <= b[, 6], , drop = F]) > 0) {
                tmp <- b[as.numeric(b[, 1]) <= as.numeric(b[, 6]), , drop = F]
                tmp <- tmp[order(abs(as.numeric(tmp[, 1]) - as.numeric(tmp[, 6])), decreasing = T), , drop = F]
                tmp <- tmp[!duplicated(tmp[, 2]), , drop = F]
                for (k in 1:nrow(tmp)) {
                        ab <- as.numeric(tmp[k, 6]) - as.numeric(tmp[k, 1])
                        b[b[, 2] %in% tmp[k, 2], 1] <- b[b[, 2] %in% tmp[k, 2], 1] + ab + bb
                }
                a[rownames(b), 1] <- b[, 1]
        }
}

points(a[, 2], a[, 1], cex = 0.4, pch = 22, bg = a[, 5])
label.panel("A", xoff = 2)
legend("top", legend = c("MRD negative", "MRD positive"), col = 1, pch = 22, pt.bg = c(0, 2))

