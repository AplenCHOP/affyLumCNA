# adjusted version
# test commands: 
# clean: load("~/dataSets/AAML_0531/cnvData/lumData/GCadjusted/RData/773376.RData")
# noisy: load("~/dataSets/AAML_0531/cnvData/lumData/GCadjusted/RData/765474.RData")
# load("ASCAT2.1/lib/omni2.5.markers.ascat.RData")
# markers <- omni2.5.markers.ascat.sorted
# source("ASCAT2.1/pipeline/AI/modifiedBirkbakLum.R")

organize.ascat.segs <- function(ascat.output, markers) {
	# init
	chromosome = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, "X", "Y")
	samp.names <- colnames(ascat.output$nA)
	failed.names <- ascat.output$failedarrays
	nSnps = 250000
	maxSeg = 150
			
	if (length(failed.names) == 0) {
		out.seg <- ascat.output$segments
        	out.seg$diff <- out.seg$endpos - out.seg$startpos
		out.seg <- out.seg[out.seg$diff != 0, ]
		# NOT NECESSARY WITH ASCATv2.2
		# Annotate row numbers to chromosome and positionout.s
		# out.seg <- matrix(0, 0, ncol = 7)
		# colnames(out.seg) <- c("sample", "chr", "startpos", "endpos", "diff", "nMajor", "nMinor")
		# sample.segs <- ascat.output$segments[[1]]
		# for (j in 1:nrow(sample.segs)) {
		#	tmp <- markers[sample.segs[j, 1]:sample.segs[j, 2], ]
		#	chr <- as.character(tmp[1, 1])
		#	seg.start <- as.character(tmp[1, 2])
		#	seg.end <- as.character(tmp[nrow(tmp), 2])
		#	diff <- nrow(tmp)
		#	out.seg <- rbind(out.seg, c(samp.names, chr, seg.start, seg.end, diff, sample.segs[j, 3:4]))
		# }
		# out.seg[out.seg[ , 2] == "X", 2] <- 23
		# out.seg <- as.data.frame(out.seg)
		# out.seg[ , 2:7] <- apply(out.seg[ , 2:7], 2, function(x) {
        	#	as.numeric(as.character(x))
		# })
		# out.seg[ , 1] <- as.character(out.seg[ , 1])
	 	# out.tmp <- out.seg
                # out.seg$cn <- log2(out.seg$nMajor + out.seg$nMinor + 0.001)
		# write.table(out.seg[, c(1, 2, 3, 4, 7, 8)], paste(samp.names, ".glad.lum.txt", sep = ""), sep = "\t", row.names = F, col.names = F, quote = F)
		
		# NOT NECESSARY WITH ASCATv2.2
		# Cleaning: Set regions with total CN > 3.5 to 2 [affy specific]
		# if (ascat.output$ploidy > 3.5) {
		#	for (i in 1:dim(out.seg)[1]) {
		#		if (out.seg$nMajor[i] > 0) {
                #                	out.seg$nMajor[i] <- out.seg$nMinor[i] - 1
		#		}
		#		if (out.seg$nMajor[i] > 0) {
                #        		out.seg$nMinor[i] <- out.seg$nMinor[i] - 1
		#		}
                #        }
                # }
		
 		# Cleaning: Set regions with CN > 2 to 2
		for (i in 1:dim(out.seg)[1]) {
			if (out.seg$nMajor[i] > 2) {
				out.seg$nMajor[i] <- 2
			}
                        if (out.seg$nMinor[i] > 2) {
                                out.seg$nMinor[i] <- 2
                        }
		}
		# out.seg$cn <- log2(out.seg$nMajor + out.seg$nMinor + 0.001)

		# Format: collapse regions if consecutive regions overlap in CN
                out.seg$volg <- 0
                out.seg$item <- 0
                for (i in 1:length(chromosome)) {
                        tmp <- out.seg[out.seg$chr == chromosome[[i]], ]
                        if (dim(tmp)[1] == 1) {
                                out.seg[out.seg$chr == chromosome[[i]], ]$volg <- 0
                                out.seg[out.seg$chr == chromosome[[i]], ]$item <- 0
                        	out.seg[out.seg$chr == chromosome[[i]], ]$diff <- out.seg[out.seg$chr == chromosome[[i]], ]$diff
			}
                        if (dim(tmp)[1] > 1){
 				out.seg[out.seg$chr == chromosome[[i]], ]$volg[1] <- 0
                                out.seg[out.seg$chr == chromosome[[i]], ]$item[1] <- 0
                                out.seg[out.seg$chr == chromosome[[i]], ]$diff[1] <- out.seg[out.seg$chr == chromosome[[i]], ]$diff[1]
                                for (j in 2:(dim(tmp)[1])) {
                                        if (tmp$nMajor[j] != tmp$nMajor[j - 1]) {
                                                out.seg[out.seg$chr == chromosome[[i]], ]$volg[j] <- out.seg[out.seg$chr == chromosome[[i]], ]$volg[j - 1] + 1
                                                out.seg[out.seg$chr == chromosome[[i]], ]$item[j] <- 0
						out.seg[out.seg$chr == chromosome[[i]], ]$diff[j] <- out.seg[out.seg$chr == chromosome[[i]], ]$diff[j]
                                        }
                                        if (tmp$nMajor[j] == tmp$nMajor[j - 1]) {
                                                if (tmp$nMinor[j] != tmp$nMinor[j - 1]) {
                                                        out.seg[out.seg$chr == chromosome[[i]], ]$volg[j] <- out.seg[out.seg$chr == chromosome[[i]], ]$volg[j - 1] + 1
                                                        out.seg[out.seg$chr == chromosome[[i]], ]$item[j] <- 0
							out.seg[out.seg$chr == chromosome[[i]], ]$diff[j] <- out.seg[out.seg$chr == chromosome[[i]], ]$diff[j]
                                                }
                                                if (tmp$nMinor[j] == tmp$nMinor[j - 1]) {
                                                        out.seg[out.seg$chr == chromosome[[i]], ]$volg[j] <- out.seg[out.seg$chr == chromosome[[i]], ]$volg[j - 1]
                                                        out.seg[out.seg$chr == chromosome[[i]], ]$item[j] <- out.seg[out.seg$chr == chromosome[[i]], ]$item[j - 1] + 1
                                                	out.seg[out.seg$chr == chromosome[[i]], ]$diff[j] <- out.seg[out.seg$chr == chromosome[[i]], ]$diff[j] + out.seg[out.seg$chr == chromosome[[i]], ]$diff[j - 1] 
						}
                                        }
                                }
                        }
                }
		out.c1 <- matrix(0, 0, ncol = 7)
                colnames(out.c1) <- c("sample", "chr", "startpos", "endpos", "diff", "nMajor", "nMinor")
                for (i in 1:(length(chromosome) - 1)) {
                        for (j in 0:(max(out.seg[out.seg$chr == chromosome[[i]] , ]$volg))) {
                                tmp <- subset(out.seg[out.seg$chr == chromosome[[i]], ], out.seg[out.seg$chr == chromosome[[i]], ]$volg == j)
                                fix <- tmp[1, ]
                                fix$endpos <- tmp$endpos[nrow(tmp)]
                                fix$diff <- sum(tmp$diff)
                                out.c1 <- rbind(out.c1, fix[ , 1:7])
                        }
                }

		# OK IF TO MANY DELETIONS THAN GIVE THEM A CLEAN PROFILE
		# if (dim(out.c1[out.c1$nMinor == 0, ])[1] >= 100 ) { 
		#	clean<-read.table("/gpfs/fs121/h/vujkovic/R/scripts/ASCAT2.1/pipeline/AI/glad.segs.template.txt",T,sep="\t", stringsAsFactors=F)
		#	write.table(clean[, c(1, 2, 3, 4, 7, 8)], paste(samp.names, ".clean.lum.txt", sep = ""), sep = "\t", row.names = F, col.names = F, quote = F)
		#}

		#if (dim(out.c1[out.c1$nMinor == 0, ]) [1] < 100 ) {
			# do the normals
			# Cleaning: Remove regions with insufficient nSnps [false positives]
                	out.c1$Flag <- as.integer(out.c1$diff < nSnps)
                	out.c2 <- matrix(0, 0, ncol = 7)
                	colnames(out.c2) <- c("sample", "chr", "startpos", "endpos", "diff", "nMajor", "nMinor")
                	for (i in 1:(length(chromosome) - 1)) {
                        	tmp <- out.c1[out.c1$chr == chromosome[[i]], ]
                        	if (sum(tmp$Flag) == 0) {
                        	        # out.c2 <- rbind(out.c2, tmp[ , 1:7])
                	        }
        	                if (sum(tmp$Flag) >= 1) {
	                                if (dim(tmp)[1] <= 2) {
                                	        if (tmp$Flag[1] == 1) {
                        	                                tmp[tmp$chr == chromosome[[i]], ]$nMajor[1] <- tmp[tmp$chr == chromosome[[i]], ]$nMajor[2]
                	                                        tmp[tmp$chr == chromosome[[i]], ]$nMinor[1] <- tmp[tmp$chr == chromosome[[i]], ]$nMinor[2]
        	                                }
	                                        if (tmp$Flag[2] == 1) {
                                                        	tmp[tmp$chr == chromosome[[i]], ]$nMajor[2] <- tmp[tmp$chr == chromosome[[i]], ]$nMajor[1]
                                                	        tmp[tmp$chr == chromosome[[i]], ]$nMinor[2] <- tmp[tmp$chr == chromosome[[i]], ]$nMinor[1]
                                        	}
                                	}
                            	    if (dim(tmp)[1] > 2) {
                                	        if (tmp$Flag[1] == 1) {
                        	                        tmp2 <- tmp[tmp$Flag == 0, ]
                	                                if (dim(tmp2)[1] == 1) {
        	                                                tmp[tmp$chr == chromosome[[i]], ]$nMajor[1] <- tmp2[tmp2$chr == chromosome[[i]], ]$nMajor[1]
	                                                        tmp[tmp$chr == chromosome[[i]], ]$nMinor[1] <- tmp2[tmp2$chr == chromosome[[i]], ]$nMinor[1]
                                        	 	}
                                        	        if (dim(tmp2)[1] > 1) {
                                	                        tmp[tmp$chr == chromosome[[i]], ]$nMajor[1] <- tmp2[tmp2$chr == chromosome[[i]], ]$nMajor[2]
                        	                                tmp[tmp$chr == chromosome[[i]], ]$nMinor[1] <- tmp2[tmp2$chr == chromosome[[i]], ]$nMinor[2]
                	                                }
        	                                }
	
                                        	for (j in 2:(dim(tmp)[1] - 1)) {
                                                	if (tmp$Flag[j] == 1) {
                                                        	if (tmp$nMajor[j - 1] == tmp$nMajor[j + 1]) {
                                                                	if (tmp$nMinor[j - 1] == tmp$nMinor[j + 1]) {
                                                                        	tmp[tmp$chr == chromosome[[i]], ]$nMajor[j] <- tmp[tmp$chr == chromosome[[i]], ]$nMajor[j + 1]
                                                                        	tmp[tmp$chr == chromosome[[i]], ]$nMinor[j] <- tmp[tmp$chr == chromosome[[i]], ]$nMinor[j + 1]
                                                                	}
                                                                	if (tmp$nMinor[j - 1] != tmp$nMinor[j + 1]) {
                                                                        	tmp[tmp$chr == chromosome[[i]], ]$nMajor[j] <- 1
                                                                        	tmp[tmp$chr == chromosome[[i]], ]$nMinor[j] <- 1
                                                                	}
                                                        	}
                                                        	if (tmp$nMajor[j - 1] != tmp$nMajor[j + 1]) {
                                                                	tmp[tmp$chr == chromosome[[i]], ]$nMajor[j] <- 1
                                                                	tmp[tmp$chr == chromosome[[i]], ]$nMinor[j] <- 1
                                                        	}
                                                	}
                                        	}
                                        	if (tmp$Flag[dim(tmp)[1]] == 1) {
                                                	tmp2 <- tmp[tmp$Flag == 0, ]
                                                	if (dim(tmp2)[1] == 1) {
                                                        	tmp[tmp$chr == chromosome[[i]], ]$nMajor[dim(tmp)[1]] <- tmp2[tmp2$chr == chromosome[[i]], ]$nMajor[1]
                                                        	tmp[tmp$chr == chromosome[[i]], ]$nMinor[dim(tmp)[1]] <- tmp2[tmp2$chr == chromosome[[i]], ]$nMinor[1]
                                                	}
							if (dim(tmp2)[1] > 1) {						
								tmp[tmp$chr == chromosome[[i]], ]$nMajor[dim(tmp)[1]] <- tmp2[tmp2$chr == chromosome[[i]], ]$nMajor[(dim(tmp2)[1] - 1)]
                                                		tmp[tmp$chr == chromosome[[i]], ]$nMinor[dim(tmp)[1]] <- tmp2[tmp2$chr == chromosome[[i]], ]$nMinor[(dim(tmp2)[1] - 1)]
							}
						}
                                	}
                        	}
                        	out.c2 <- rbind(out.c2, tmp[ , 1:7])
                	}

			# Format: collapse regions if consecutive regions overlap in CN
			out.c2$volg <- 0
			out.c2$item <- 0
			for (i in 1:length(chromosome)) {
        			tmp <- out.c2[out.c2$chr == chromosome[[i]], ]
        			if (dim(tmp)[1] == 1) {
                			out.c2[out.c2$chr == chromosome[[i]], ]$volg <- 0
                			out.c2[out.c2$chr == chromosome[[i]], ]$item <- 0
        			}
        			if (dim(tmp)[1] > 1) {
                			for (j in 2:(dim(tmp)[1])) {
                        			if (tmp$nMajor[j] != tmp$nMajor[j - 1]) {
                                			out.c2[out.c2$chr == chromosome[[i]], ]$volg[j] <- out.c2[out.c2$chr == chromosome[[i]], ]$volg[j - 1] + 1
                                			out.c2[out.c2$chr == chromosome[[i]], ]$item[j] <- 0
                        			}
                        			if (tmp$nMajor[j] == tmp$nMajor[j - 1]) {
                                			if (tmp$nMinor[j] != tmp$nMinor[j - 1]) {
                                        			out.c2[out.c2$chr == chromosome[[i]], ]$volg[j] <- out.c2[out.c2$chr == chromosome[[i]], ]$volg[j - 1] + 1
                                        			out.c2[out.c2$chr == chromosome[[i]], ]$item[j] <- 0
                                			}	
                                			if (tmp$nMinor[j] == tmp$nMinor[j - 1]) {
                                        			out.c2[out.c2$chr == chromosome[[i]], ]$volg[j] <- out.c2[out.c2$chr == chromosome[[i]], ]$volg[j - 1]
                                        			out.c2[out.c2$chr == chromosome[[i]], ]$item[j] <- out.c2[out.c2$chr == chromosome[[i]], ]$item[j - 1] + 1
                                       			}
                       				}
                			}
        			}
			}
			out.c3 <- matrix(0, 0, ncol = 7)
			colnames(out.c3) <- c("sample", "chr", "startpos", "endpos", "diff", "nMajor", "nMinor")
			for (i in 1:(length(chromosome) - 1)) {
        			for (j in 0:(max(out.c2[out.c2$chr == chromosome[[i]], ]$volg))) {
                			tmp <- subset(out.c2[out.c2$chr == chromosome[[i]], ], out.c2[out.c2$chr == chromosome[[i]], ]$volg == j)
                			fix <- tmp[1, ]
                			fix$endpos <- tmp$endpos[nrow(tmp)]
                			fix$diff <- sum(tmp$diff)
                			out.c3 <- rbind(out.c3, fix[ , 1:7])
        			}
			}
			# out.c3$cn <- log2(out.cn3$nMajor + out.c3$nMinor + 0.001)
			write.table(out.c3[, c(1, 2, 3, 4, 7, 8)], paste(samp.names, "_", nSnps/1000, "MB.glad.lum.txt", sep = ""), sep = "\t", row.names = F, col.names = F, quote = F)
	}	#}
}


