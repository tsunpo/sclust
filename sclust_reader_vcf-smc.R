# =============================================================================
# Title: sclust_reader_vcf-smc.R
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 15/06/16
# =============================================================================

# -----------------------------------------------------------------------------
# Sclust 1 - VCF reader - SMC-Het format
# -----------------------------------------------------------------------------
initVCF <- function(sample.snv) { 
   vcf <- data.frame(matrix(".", nrow(sample.snv), 8))
   names(vcf) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
   
   vcf[,1] <- paste("chr", sample.snv$CHROM, sep="")
   vcf$POS <- sample.snv$POS
   vcf$ID  <- sample.snv$ID
   vcf$REF <- sample.snv$REF
   vcf$ALT <- sample.snv$ALT
   vcf$QUAL   <- 255
   vcf$FILTER <- "PASS"
   
   return(vcf);
}

vcfInfoDP <- function(format) {
   format1 <- unlist(strsplit(format, ":"))
   format2 <- unlist(strsplit(format1[2], ","))
   ref <- as.numeric(format2[1])
   alt <- as.numeric(format2[2])
   
   return(alt + ref)
}

vcfInfoAF <- function(format) {
   format1 <- unlist(strsplit(format, ":"))
   format2 <- unlist(strsplit(format1[2], ","))
   ref <- as.numeric(format2[1])
   alt <- as.numeric(format2[2])
   
   return( round( alt / (alt + ref), 7) );
}

##
args <- commandArgs(T)
sample <- args[1]

if (length(readLines(sample)) != 0) {
   sample.snv <- read.table(sample,  header=F, sep="\t", fill=T, as.is=T, comment.char="#")
   names(sample.snv) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "NORMAL", "TUMOR")
   
   vcf <- initVCF(sample.snv)
   ## Read depth
   vcf$INFO <- paste("DP=", mapply(v = 1:nrow(vcf), function(v) vcfInfoDP(sample.snv[v, 11])), sep="")
   vcf$INFO <- paste(vcf$INFO, "DP_N=.", sep=";")
   
   ## AF
   vcf$INFO <- paste(vcf$INFO, paste("AF=", mapply(v = 1:nrow(vcf), function(v) vcfInfoAF(sample.snv[v, 11])), sep=""), sep=";")
   vcf$INFO <- paste(vcf$INFO, "AF_N=.", sep=";")
   
   ## FR and TG   
   vcf$INFO <- paste(vcf$INFO, "FR=.", sep=";")
   vcf$INFO <- paste(vcf$INFO, "TG=.", sep=";")

   write.table(vcf[,1:8], "sclust.vcf", col.names=names(vcf[,1:8]), row.names=F, quote=F, sep="\t")
}
