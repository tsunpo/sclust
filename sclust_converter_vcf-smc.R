# =============================================================================
# Title: sclust_converter_vcf-smc.R
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 14/06/16
# =============================================================================
source("the_daily_package.R")

# -----------------------------------------------------------------------------
# Sclust 1 - Generate converted VCFs (SMC-Het format)
# -----------------------------------------------------------------------------
initVCF <- function(tmp) { 
   vcf <- toTable(".", 8, nrow(tmp), c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))
   
   vcf[,1] <- paste("chr", tmp$CHROM, sep="")
   vcf$POS <- tmp$POS
   vcf$ID  <- tmp$ID
   vcf$REF <- tmp$REF
   vcf$ALT <- tmp$ALT
   vcf$QUAL   <- "255"
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
   tmp        <- readTable(sample, header=F, rownames=F, sep="\t")
   names(tmp) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "NORMAL", "TUMOR")
   
   vcf <- initVCF(tmp)
   ## Read depth
   vcf$INFO <- paste("DP=", mapply(v = 1:nrow(vcf), function(v) vcfInfoDP(tmp[v, 11])), sep="")
   vcf$INFO <- paste(vcf$INFO, "DP_N=.", sep=";")
   
   ## AF
   vcf$INFO <- paste(vcf$INFO, paste("AF=", mapply(v = 1:nrow(vcf), function(v) vcfInfoAF(tmp[v, 11])), sep=""), sep=";")
   vcf$INFO <- paste(vcf$INFO, "AF_N=.", sep=";")
   
   ## FR and TG   
   vcf$INFO <- paste(vcf$INFO, "FR=.", sep=";")   ## TO-DO
   vcf$INFO <- paste(vcf$INFO, "TG=.", sep=";")

   writeTable(vcf, gzfile(paste(sample, ".gz", sep="")), colnames=T, rownames=F, sep="\t")
}
