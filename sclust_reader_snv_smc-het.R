# =============================================================================
# Title: Sclust for the SMC-Het Challenge
# Name: sclust_reader_snv_smc-het.R
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 15/06/16
# =============================================================================

# -----------------------------------------------------------------------------
# Sclust - 1. SNV reader - Preprocessing SMC-Het format
# -----------------------------------------------------------------------------
initSNV <- function(sample.snv) { 
   snv <- data.frame(matrix(".", nrow(sample.snv), 8))
   names(snv) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
   
   snv[,1] <- paste("chr", sample.snv$CHROM, sep="")
   snv$POS <- sample.snv$POS
   snv$ID  <- sample.snv$ID
   snv$REF <- sample.snv$REF
   snv$ALT <- sample.snv$ALT
   snv$QUAL   <- 255
   snv$FILTER <- "PASS"
   
   return(snv);
}

snvInfoDP <- function(format) {
   format1 <- unlist(strsplit(format, ":"))
   format2 <- unlist(strsplit(format1[2], ","))
   ref <- as.numeric(format2[1])
   alt <- as.numeric(format2[2])
   
   return(alt + ref)
}

snvInfoAF <- function(format) {
   format1 <- unlist(strsplit(format, ":"))
   format2 <- unlist(strsplit(format1[2], ","))
   ref <- as.numeric(format2[1])
   alt <- as.numeric(format2[2])
   
   return( round( alt / (alt + ref), 7) );
}

##
## Main
args <- commandArgs(T)
sample <- args[1]

if (length(readLines(sample)) != 0) {
   sample.snv <- read.table(sample, header=F, sep="\t", fill=T, as.is=T, comment.char="#")
   names(sample.snv) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "NORMAL", "TUMOR")
   
   snv <- initSNV(sample.snv)
   ## Read depth
   snv$INFO <- paste("DP=", mapply(v = 1:nrow(snv), function(v) snvInfoDP(sample.snv[v, 11])), sep="")
   snv$INFO <- paste(snv$INFO, "DP_N=.", sep=";")
   
   ## AF
   snv$INFO <- paste(snv$INFO, paste("AF=", mapply(v = 1:nrow(snv), function(v) snvInfoAF(sample.snv[v, 11])), sep=""), sep=";")
   snv$INFO <- paste(snv$INFO, "AF_N=.", sep=";")
   
   ## FR and TG   
   snv$INFO <- paste(snv$INFO, "FR=.", sep=";")
   snv$INFO <- paste(snv$INFO, "TG=.", sep=";")

   write.table(snv[,1:8], "sclust.vcf", col.names=names(snv[,1:8]), row.names=F, quote=F, sep="\t")
}
