# Required Packages
library(affy)
library(squash)
library(genefilter)
library(beeswarm)
library(gplots)
library(gdata)
library(ROC)
library(DNAcopy)
library(quantsmooth)
library(cluster)
library(jetset)
library(boot)

#
#This function puts a label on the top left corner of a gure
#

label.panel <- function(txt, xoff = 1, yoff = 2, cex = 1.2) {
	x <- grconvertX(0, from = "nfc") + (xoff * strwidth("m"))
	y <- grconvertY(1, from = "nfc") - (yoff * strheight("A"))
	text(x, y, labels = txt, font = 2, xpd = TRUE, cex = cex)
}

#
# This function re-organizes the output from running the ASCAT algorithm into a
# more compact format similar to the output from other segmentation algorithms such
# as Circular Binary Segmentation
#

organize.ascat.segments <- function(ascat.output, markers) {
	samp.names <- colnames(ascat.output$nA)
	failed.names <- ascat.output$failedarrays
	if (length(failed.names) >= 1) {
		failed.locations <- which(is.na(ascat.output$segments))
		new.samp.names <- 1:length(ascat.output$segments)
		new.samp.names[-c(failed.locations)] <- samp.names
		new.samp.names[failed.locations] <- failed.names
		samp.names <- new.samp.names
	}
	out.seg <- matrix(0, 0, ncol = 7)
	colnames(out.seg) <- c("SampleID", "Chr", "Start", "End", "nProbes", "nA", "nB")
	if (length(ascat.output$segments) != length(samp.names)) {
		stop
	}
	for (i in 1:length(ascat.output$segments)) {
		sample.segs <- ascat.output$segments[[i]]
		if (class(sample.segs) != "matrix") {
			next
		}
		if (class(sample.segs) == "matrix") {
			for (j in 1:nrow(sample.segs)) {
				tmp <- markers[sample.segs[j, 1]:sample.segs[j, 2], ]
				if (!all(tmp[, 1] == tmp[1, 1])) {
					stop("Error in segmentation")
				}
				chr <- as.character(tmp[1, 1])
				seg.start <- as.character(tmp[1, 2])
				seg.end <- as.character(tmp[nrow(tmp), 2])
				nProbes <- nrow(tmp)
				out.seg <- rbind(out.seg, c(samp.names[i], chr, seg.start, seg.end, nProbes, sample.segs[j, 3:4]))
			}
			gc()
		}
	}
	out.seg[out.seg[, 2] == "X", 2] <- 23
	out.seg <- as.data.frame(out.seg)
	out.seg[, 2:7] <- apply(out.seg[, 2:7], 2, function(x) {
		as.numeric(as.character(x))
	})
	out.seg[, 1] <- as.character(out.seg[, 1])
	return(out.seg)
}

#
# This function denes and summarizes the regions of telomeric allelic imbalance
#

no.tel <- function(seg.out, chrominfo, min.size = 0, min.probes = 1,max.size = 1e+09, cnv.check = "no", cnv.seg = NULL, cnv.gain = NULL) {
	if (class(seg.out) == "DNAcopy") {
		seg.out <- seg.out$output
	}
	tmp.segs <- seg.out[!seg.out[, 2] %in% c("MT", "Y", "24"), ]
	tmp.segs <- tmp.segs[!tmp.segs[, 5] < min.probes, ]
	tmp.segs[, 2] <- as.character(tmp.segs[, 2])
	if (nrow(tmp.segs[tmp.segs[, 2] == "X", ]) > 0) {
		tmp.segs[tmp.segs[, 2] == "X", 2] <- rep(23, nrow(tmp.segs[tmp.segs[, 2] == "X", ]))
	}
	if (cnv.check != "no") {
		if (class(cnv.seg) == "DNAcopy") {
			cnv.seg <- cnv.seg$output
		}
		tmp.cnv <- cnv.seg[!cnv.seg[, 2] %in% c("MT", "Y", "24"), ]
		tmp.cnv <- tmp.cnv[!tmp.cnv[, 5] < min.probes, ]
		tmp.cnv[, 2] <- as.character(tmp.cnv[, 2])
		if (nrow(tmp.cnv[tmp.cnv[, 2] == "X", ]) > 0) {
			tmp.cnv[tmp.cnv[, 2] == "X", 2] <- rep(23, nrow(tmp.cnv[tmp.cnv[, 2] == "X", ]))
		}
	}
	tmp.segs[tmp.segs[, 6] == 1, 6] <- 2
	for (j in 1:length(unique(tmp.segs[, 1]))) {
		tmp.sample <- tmp.segs[tmp.segs[, 1] == unique(tmp.segs[, 1])[j], ]
		for (i in 1:23) {
			tmp1 <- tmp.sample[tmp.sample[, 2] == i, , drop = F]
			if (tmp1[1, 6] == 2 & nrow(tmp1) != 1 & tmp1[1, 4] < (chrominfo[i, 3] * 1000)) {
				tmp.sample[tmp.sample[, 2] == i, 6][1] <- 1
				if (cnv.check != "no") {
					if (cnv.seg[cnv.seg[, 1] == unique(tmp.segs[, 1])[j] & cnv.seg[, 2] == i, ][1, 6] > cnv.gain) {
						tmp.sample[tmp.sample[, 2] == i, 6][1] <- 0
					}
				}
			}
			if (tmp1[nrow(tmp1), 6] == 2 & nrow(tmp1) != 1 & tmp1[nrow(tmp1), 3] > (chrominfo[i, 3] * 1000)) {
				tmp.sample[tmp.sample[, 2] == i, 6][nrow(tmp.sample[tmp.sample[, 2] == i, ])] <- 1
				if (cnv.check != "no") {
					if (cnv.seg[cnv.seg[, 1] == unique(tmp.segs[, 1])[j] & cnv.seg[, 2] == i, ][nrow(cnv.seg[cnv.seg[, 1] == unique(tmp.segs[, 1])[j] & cnv.seg[, 2] == i, ]), 6] > cnv.gain) {
						tmp.sample[tmp.sample[, 2] == i, 6][nrow(tmp.sample[tmp.sample[, 2] == i, ])] <- 0
					}
				}
			}
			if (nrow(tmp.sample[tmp.sample[, 2] == i, ]) == 1 & tmp.sample[tmp.sample[, 2] == i, 6][1] != 0) {
				tmp.sample[tmp.sample[, 2] == i, 6][1] <- 3
			}
		}
		tmp.segs[tmp.segs[, 1] == unique(tmp.segs[, 1])[j], ] <- tmp.sample
	}
	no.events <- matrix(0, nrow = length(unique(tmp.segs[, 1])), ncol = 28)
	rownames(no.events) <- unique(tmp.segs[, 1])
	colnames(no.events) <- c("Telomeric AI", "Mean size", "Interstitial AI", "Mean Size", "Wholo chr AI", 1:23)
	a <- 0
	for (i in unique(tmp.segs[, 1])) {
		a <- a + 1 
		tmp <- tmp.segs[tmp.segs[, 1] == i, ]
		tmp <- tmp[(tmp[, 4] - tmp[, 3]) > min.size, ]
		tmp <- tmp[(tmp[, 4] - tmp[, 3]) < max.size, ]
		no.events[a, 1] <- nrow(tmp[tmp[, 6] == 1, ])
		no.events[a, 2] <- mean(tmp[tmp[, 6] == 1, 4] - tmp[tmp[, 6] == 1, 3])
		no.events[a, 3] <- nrow(tmp[tmp[, 6] == 2, ])
		no.events[a, 4] <- mean(tmp[tmp[, 6] == 2, 4] - tmp[tmp[, 6] == 2, 3])
		no.events[a, 5] <- nrow(tmp[tmp[, 6] == 3, ])
		no.events[a, tmp[tmp[, 6] == 3, 2]] <- 1
	}
	return(no.events)
}

#
#This function denes and summarizes the regions of telomeric copy number changes
#

no.tel.cna <- function(seg.out, chrominfo, min.size = 0, min.probes = 1, max.size = 1e+09, gain = log2(2.5/2), loss = log2(1.5/2)) {
	if (class(seg.out) == "DNAcopy") {
		seg.out <- seg.out$output
	}
	tmp.segs <- seg.out
	tmp.segs <- tmp.segs[!tmp.segs[, 5] < min.probes, ]
	tmp.segs[, 2] <- as.character(tmp.segs[, 2])
	if (nrow(tmp.segs[tmp.segs[, 2] == "X", ]) > 0) {
		tmp.segs[tmp.segs[, 2] == "X", 2] <- rep(23, nrow(tmp.segs[tmp.segs[, 2] == "X", ]))
	}
	seg.out <- seg.out[!seg.out[, 2] %in% c("MT", "Y", "24"), ]
  	tmp.segs[tmp.segs[, 6] >= gain, 6] <- 10
	tmp.segs[tmp.segs[, 6] <= loss, 6] <- -10
	tmp.segs[tmp.segs[, 6] > loss & tmp.segs[, 6] < gain, 6] <- 0
	for (j in 1:length(unique(tmp.segs[, 1]))) { 
		tmp.sample <- tmp.segs[tmp.segs[, 1] == unique(tmp.segs[, 1])[j], ]
		for (i in 1:23) {
			tmp1 <- tmp.sample[tmp.sample[, 2] == i, , drop = F]
			if (tmp1[1, 6] == 10 & nrow(tmp1) != 1 & tmp1[1, 4] < (chrominfo[i, 3] * 1000)) {
				tmp.sample[tmp.sample[, 2] == i, 6][1] <- 1
			}
			if (tmp1[nrow(tmp1), 6] == 10 & nrow(tmp1) != 1 & tmp1[nrow(tmp1), 3] > (chrominfo[i, 3] * 1000)) {
				tmp.sample[tmp.sample[, 2] == i, 6][nrow(tmp.sample[tmp.sample[, 2] == i, ])] <- 1
			}
			if (tmp1[1, 6] == -10 & nrow(tmp1) != 1 & tmp1[1, 4] < (chrominfo[i, 3] * 1000)) {
				tmp.sample[tmp.sample[, 2] == i, 6][1] <- 2
			}
			if (tmp1[nrow(tmp1), 6] == -10 & nrow(tmp1) != 1 & tmp1[nrow(tmp1), 3] > (chrominfo[i, 3] * 1000)) {
				tmp.sample[tmp.sample[, 2] == i, 6][nrow(tmp.sample[tmp.sample[, 2] == i, ])] <- 2
			}
			if (nrow(tmp.sample[tmp.sample[, 2] == i, ]) == 10 & tmp.sample[tmp.sample[, 2] == i, 6][1] != 0) {
				tmp.sample[tmp.sample[, 2] == i, 6][1] <- 3
			}
			if (nrow(tmp.sample[tmp.sample[, 2] == i, ]) == -10 & tmp.sample[tmp.sample[, 2] == i, 6][1] != 0) {
				tmp.sample[tmp.sample[, 2] == i, 6][1] <- 4
			}
		}
		tmp.segs[tmp.segs[, 1] == unique(tmp.segs[, 1])[j], ] <- tmp.sample
	}
	no.events <- matrix(0, nrow = length(unique(tmp.segs[, 1])), ncol = 9)
	rownames(no.events) <- unique(tmp.segs[, 1])
	colnames(no.events) <- c("Total telomeric CNA", "Mean size", "Telomeric gain", "Mean size", "Telomeric loss", "Mean size", "Total whole chromosome CNA", "Whole chromosome gain", "Whole chromosome loss")
	a <- 0
	for (i in unique(tmp.segs[, 1])) {
		a <- a + 1 
		tmp <- tmp.segs[tmp.segs[, 1] == i, ]
		tmp <- tmp[(tmp[, 4] - tmp[, 3]) > min.size, ]
		tmp <- tmp[(tmp[, 4] - tmp[, 3]) < max.size, ]
		no.events[a, 1] <- nrow(tmp[tmp[, 6] == 1 | tmp[, 6] == 2, ])
		no.events[a, 2] <- mean(tmp[tmp[, 6] == 1 | tmp[, 6] == 2, 4] - tmp[tmp[, 6] == 1 | tmp[, 6] == 2, 3])
		no.events[a, 3] <- nrow(tmp[tmp[, 6] == 1, ])
		no.events[a, 4] <- mean(tmp[tmp[, 6] == 1, 4] - tmp[tmp[, 6] == 1, 3])
		no.events[a, 5] <- nrow(tmp[tmp[, 6] == 2, ])
		no.events[a, 6] <- mean(tmp[tmp[, 6] == 2, 4] - tmp[tmp[, 6] == 2, 3])
		no.events[a, 7] <- nrow(tmp[tmp[, 6] == 3 | tmp[, 6] == 4, ])
		no.events[a, 8] <- nrow(tmp[tmp[, 6] == 3, ])
		no.events[a, 9] <- nrow(tmp[tmp[, 6] == 4, ])
	}
	return(no.events)
}

#
#This function counts the total number of copy number aberrations
#

gaps <- function(cbs.out, min.size = 0, gain = log2(2.5/2), loss = log2(1.5/2), min.probes = 1) {
	if (length(gain) == 1) {
		gain <- rep(gain, length(unique(cbs.out[, 1])))
	}
	if (length(loss) == 1) {
		loss <- rep(loss, length(unique(cbs.out[, 1])))
	}
	stopifnot(length(gain) == length(unique(cbs.out[, 1])))
	stopifnot(length(loss) == length(unique(cbs.out[, 1])))
	cbs.out <- cbs.out[!(cbs.out[, 2] == 24), ]
	cbs.out <- cbs.out[!(cbs.out[, 5] <= min.probes), ]
	gaps <- matrix(rep(NA, length(unique(cbs.out[, 1]))), nrow = length(unique(cbs.out[, 1])), ncol = 7)
	rownames(gaps) <- unique(cbs.out[, 1])
	colnames(gaps) <- c("Short Aberrations", "Long Aberrations", "Short gains", "Short deletions", "Total Gains", "Total Deletions", "Total Aberrations")
	for (i in 1:nrow(gaps)) {
		a <- cbs.out[cbs.out[, 1] == rownames(gaps)[i], ]
		a <- a[a[, 6] >= gain[i] | a[, 6] <= loss[i], ]
		gaps[i, 1] <- nrow(a[a[, 4] - a[, 3] <= min.size, ])
		gaps[i, 2] <- nrow(a[a[, 4] - a[, 3] > min.size, ])
		gaps[i, 5] <- nrow(a[a[, 6] > gain[i], ])
		gaps[i, 6] <- nrow(a[a[, 6] < loss[i], ])
		b <- a[a[, 4] - a[, 3] <= min.size, ]
		gaps[i, 3] <- nrow(b[b[, 6] > gain[i], ])
		gaps[i, 4] <- nrow(b[b[, 6] < loss[i], ])
	}
	gaps[, 7] <- gaps[, 5] + gaps[, 6]
	return(gaps)
}

#
# This function identies all telomeric AI and returns the location
#

telAI.fn <- function(seg.out, chrominfo, output = "tel", min.size = 0, min.probes = 1, max.size = 1e+09) {
	if (class(seg.out) == "DNAcopy") {
		seg.out <- seg.out$output
	}
	tmp.segs <- seg.out[!seg.out[, 2] %in% c("MT", "Y", "24"), ]
	tmp.segs <- tmp.segs[!tmp.segs[, 5] < min.probes, ]
	tmp.segs[, 2] <- as.character(tmp.segs[, 2])
	if (nrow(tmp.segs[tmp.segs[, 2] == "X", ]) > 0) {
		tmp.segs[tmp.segs[, 2] == "X", 2] <- rep(23, nrow(tmp.segs[tmp.segs[, 2] == "X", ]))
	}
	tmp.segs[tmp.segs[, 6] == 1, 6] <- 2
	tmp.segs <- tmp.segs[tmp.segs[, 4] - tmp.segs[, 3] >= min.size, ]
	for (j in 1:length(unique(tmp.segs[, 1]))) {
		tmp.sample <- tmp.segs[tmp.segs[, 1] == unique(tmp.segs[, 1])[j], ]
		for (i in 1:23) {
			tmp1 <- tmp.sample[tmp.sample[, 2] == i, , drop = F]
			if (tmp1[1, 6] == 2 & nrow(tmp1) != 1 & tmp1[1, 4] < (chrominfo[i, 3] * 1000)) {
				tmp.sample[tmp.sample[, 2] == i, 6][1] <- 1
			}
			if (tmp1[nrow(tmp1), 6] == 2 & nrow(tmp1) != 1 & tmp1[nrow(tmp1), 3] > (chrominfo[i, 3] * 1000)) {
				tmp.sample[tmp.sample[, 2] == i, 6][nrow(tmp.sample[tmp.sample[, 2] == i, ])] <- 1
			}
			if (nrow(tmp.sample[tmp.sample[, 2] == i, ]) == 1 & tmp.sample[tmp.sample[, 2] == i, 6][1] != 0) {
				tmp.sample[tmp.sample[, 2] == i, 6][1] <- 3
			}
		}
		tmp.segs[tmp.segs[, 1] == unique(tmp.segs[, 1])[j], ] <- tmp.sample
	}
	if (output == "tel") {
		tmp.segs <- tmp.segs[tmp.segs[, 6] == 1, ]
	}
	if (output == "int") {
		tmp.segs <- tmp.segs[tmp.segs[, 6] == 2, ]
	}
	return(tmp.segs)
}

#
# This function takes as input breakpoint location of telomeric AI, and then test
# each breakpoint for association with a given DNA structure, supplied as a two-column
# matrix of locations
#

telBreak.CNV.fn <- function(test.breaks, test.loc, window.size = 25000) {
	break.loc <- matrix(0, 0, 2)
	for (i in 1:nrow(test.breaks)) {
		break.loc <- rbind(break.loc, c(as.numeric(test.breaks[i, 2]), NA))
		if (as.numeric(test.breaks[i, 3]) < chrominfo[as.numeric(test.breaks[i, 2]), 3] * 1000) {
			break.loc[i, 2] <- as.numeric(test.breaks[i, 4])
		}
		if (as.numeric(test.breaks[i, 3]) > chrominfo[as.numeric(test.breaks[i, 2]), 3] * 1000) {
			break.loc[i, 2] <- as.numeric(test.breaks[i, 3])
		}
	}
	test.breaksDNA <- vector()
	for (i in 1:nrow(break.loc)) {
		tmp <- test.loc[test.loc[, 1] == break.loc[i, 1], ]
		a <- c(break.loc[i, 2] - window.size, break.loc[i, 2] + window.size)
		tmp <- tmp[tmp[, 3] >= a[1], , drop = F]
		tmp <- tmp[tmp[, 2] <= a[2], , drop = F]
		test.breaksDNA[i] <- nrow(tmp)
	}
	tmp <- c(length(test.breaksDNA[test.breaksDNA > 0]), length(test.breaksDNA))
	names(tmp) <- c("Associated with", "Total")
	return(tmp)
}

#
# This functions calculate sensitivity, specicity, positive predictive value, negative
#predictive value, accuracy, and a p-value, based on predicted versus true input
#

acc.calc <- function(test, truth) {
	if (class(test) == "logical") {
		test <- c(0, 1)[match(test, c("FALSE", "TRUE"))]
	}
	if (class(truth) == "logical") {
		truth <- c(0, 1)[match(truth, c("FALSE", "TRUE"))]
	}
	TN <- length(test[test == 0 & truth == 0])
	TP <- length(test[test == 1 & truth == 1])
	FN <- length(test[test == 0 & truth == 1])
	FP <- length(test[test == 1 & truth == 0])
	c(ACC = (TP + TN)/(TP + FP + FN + TN), PPV = TP/(TP + FP), NPV = TN/(TN + FN), SENS = TP/(TP + FN), SPEC = TN/(TN + FP), P = fisher.test(table(test, truth))$p.value)
}

#
# This is a convenience function that simply takes two vectors as input, and add one as names for the other
#

name.vect <- function(x, y) {
	names(x) <- y
	return(x)
}

#
# These two functions are used to calculate condence intervals for AUC curves by
# bootstrapping 1000 times
#

AUC.stat <- function(data, indices) {
	d <- data[indices, ]
	return(AUC(rocdemo.sca(d[, 1], d[, 2])))
}

ci.auc <- function(data, pcr, subset = TRUE) {
	data <- data[subset]
	pcr <- pcr[subset]
	results <- vector()
	d <- cbind(pcr, data)
	tmp <- boot(data = d, statistic = AUC.stat, R = 1000)
	tmp <- boot.ci(tmp, type = "bca")
	results[1] <- tmp$bca[4]
	results[2] <- tmp$bca[5]
	names(results) <- c("Low", "High")
	return(results)
}

# THE ENDS #
