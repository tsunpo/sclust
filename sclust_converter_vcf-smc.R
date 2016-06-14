# =============================================================================
# Title: sclust_converter_vcf-smc.R
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 14/06/16
# =============================================================================

# -----------------------------------------------------------------------------
# Sclust 1 - Generate converted VCFs (SMC-Het format)
# -----------------------------------------------------------------------------
source("the_daily_package.R")

initVCF <- function(tmp) { 
   vcf <- toTable(".", 8, nrow(tmp), c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))
   
   vcf[,1] <- paste("chr", tmp$CHROM, sep="")
   vcf$POS <- tmp$POS
   vcf$REF <- tmp$REF
   vcf$ALT <- tmp$ALT
   vcf$QUAL <- "255"
   vcf$FILTER <- "PASS"
   
   return(vcf);
}

vcfInfoDP <- function(format) {
   format1 <- unlist(strsplit(format, "t_alt_count="))
   format2 <- unlist(strsplit(format1[2], ";"))
   alt <- as.numeric(format2[1])

   format1 <- unlist(strsplit(format, "t_ref_count="))
   format2 <- unlist(strsplit(format1[2], ";"))
   ref <- as.numeric(format2[1])
   
   return(alt + ref)
}

vcfInfoAF <- function(format) {
   format1 <- unlist(strsplit(format, "t_alt_count="))
   format2 <- unlist(strsplit(format1[2], ";"))
   alt <- as.numeric(format2[1])

   format1 <- unlist(strsplit(format, "t_ref_count="))
   format2 <- unlist(strsplit(format1[2], ";"))
   ref <- as.numeric(format2[1])
   
   return( round( alt / (alt + ref), 7) );
}


##
args <- commandArgs(T)

vcfdat = read.table(args[1], sep='\t', comment.char='#')

   file <- paste(wd.in, sample, ".annotated.snv_mnv.vcf.gz", sep="")
   if (length(readLines(file)) != 0) {
      tmp        <- readTable(file, header=F, rownames=F, sep="\t")
      names(tmp) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
   
      vcf <- initVCF(tmp)
      ## Read depth
      vcf$INFO <- paste("DP=", mapply(v = 1:nrow(vcf), function(v) vcfInfoDP(tmp[v, 8])), sep="")
      vcf$INFO <- paste(vcf$INFO, "DP_N=0", sep=";")
   
      ## AF
      vcf$INFO <- paste(vcf$INFO, paste("AF=", mapply(v = 1:nrow(vcf), function(v) vcfInfoAF(tmp[v, 8])), sep=""), sep=";")
      vcf$INFO <- paste(vcf$INFO, "AF_N=0", sep=";")
   
      ## FR and TG   
      vcf$INFO <- paste(vcf$INFO, "FR=0", sep=";")   ## TO-DO
      vcf$INFO <- paste(vcf$INFO, "TG=Generic", sep=";")

      #writeTable(vcf, paste(sample.out, sample$V2, ".dkfz-snvCalling.somatic.snv_mnv.vcf", sep=""), colnames=T, rownames=F, sep="\t")
      writeTable(vcf, gzfile(paste(wd.out, sample, ".annotated.snv_mnv.vcf.gz", sep="")), colnames=T, rownames=F, sep="\t")
 