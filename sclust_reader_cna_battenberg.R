# =============================================================================
# Title: Sclust submission to the SMC-Het Challenge
# Name: sclust_reader_cna_battenberg.R
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 15/06/16
# =============================================================================

# -----------------------------------------------------------------------------
# Sclust - 2. CNA reader - Preprocessing Battenberg format
# -----------------------------------------------------------------------------
initCNA <- function(sample.cna) { 
   cna <- data.frame(matrix(0, nrow(sample.cna), 7))
   names(cna) <- c("chromosome", "start", "end", "copy_number", "major_cn", "minor_cn", "clonal_frequency")
   
   cna$chromosome <- sample.cna$chr
   cna$start <- sample.cna$startpos
   cna$end   <- sample.cna$endpos

   return(cna);
}

initCNA1 <- function(sample.cna) { 
   cna <- initCNA(sample.cna)
   
   cna$major_cn <- sample.cna$nMaj1_A
   cna$minor_cn <- sample.cna$nMin1_A
   cna$copy_number <- cna$major_cn + cna$minor_cn
   cna$clonal_frequency <- sample.cna$frac1_A
   
   return(cna);
}

initCNA2 <- function(sample.cna2) { 
   cna2 <- initCNA(sample.cna2)
   
   cna2$major_cn <- sample.cna2$nMaj2_A
   cna2$minor_cn <- sample.cna2$nMin2_A
   cna2$copy_number <- cna2$major_cn + cna2$minor_cn
   cna2$clonal_frequency <- sample.cna2$frac2_A
   
   return(cna2);
}

##
## Main
args <- commandArgs(T)
sample <- args[1]

if (length(readLines(sample)) != 0) {
   sample.cna <- read.table(sample, header=T, sep="\t", fill=T, as.is=T, comment.char="#")

   cna <- initCNA1(sample.cna)
   
   ## For clonal_frequency != 1
   sample.cna2 <- subset(sample.cna, frac1_A != 1)
   cna <- rbind(cna, initCNA2(sample.cna2))

   write.table(cna, "sclust.txt", col.names=names(cna), row.names=F, quote=F, sep="\t")
}
