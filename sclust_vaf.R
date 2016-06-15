# =============================================================================
# Title: sclust_vaf.R
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 15/06/16
# =============================================================================

# -----------------------------------------------------------------------------
# Sclust - 3. VAF - Determining cancer cell fractions (CCFs)
# -----------------------------------------------------------------------------
initExpAF <- function(nrow) {
   expAF <- data.frame(matrix(0, nrow, 12))
   names(expAF) <- c("Mut_ID", "Chr", "Position", "Wt", "Mut", "AF_obs", "Coverage", "AF_exp", "Mut_Copies", "Mut_Copies_Raw", "Is_Subclonal_CN", "iCN")
 
   return(expAF);
}

initExpAFfromVCF <- function(sample, vcf.segment) { 
   expAF <- initExpAF(nrow(vcf.segment))
   
   expAF$Chr <- vcf.segment$CHROM[1]
   expAF$Position <- vcf.segment$POS
   expAF$Wt  <- vcf.segment$REF
   expAF$Mut <- vcf.segment$ALT
   expAF$Mut_ID <- paste(sample, paste(expAF$Chr, expAF$Position, sep=":"), "SNM", sep="_")
   
   ## VAFobs
   expAF$AF_obs <- mapply(v = 1:nrow(vcf.segment), function(v) obsVAF(vcf.segment$INFO[v]))
   expAF$Coverage <- mapply(v = 1:nrow(vcf.segment), function(v) coverage(vcf.segment$INFO[v]))

   return(expAF);
}

vcfGetSegment <- function(vcf, chromosome, start, end) {
   vcf.segment <- vcf[vcf$CHROM == paste("chr", chromosome, sep=""),]
   vcf.segment <- vcf.segment[vcf.segment $POS >= start,]
   vcf.segment <- vcf.segment[vcf.segment $POS <= end,]
   
   return(vcf.segment);
}

obsVAF <- function(format) {
   format <- unlist(strsplit(format, ";"))
   
   for (f in 1:length(format)) {
   	  value <- unlist(strsplit(format[f], "="))
      if (value[1] == "AF")
         return(value[2])
   }
}

coverage <- function(format) {
   format <- unlist(strsplit(format, ";"))
   
   for (f in 1:length(format)) {
   	  value <- unlist(strsplit(format[f], "="))
      if (value[1] == "DP")
         return(value[2])
   }
}

multiplicity <- function(p, CN, obsVAF) {
   return( round((( 2*(1-p) + CN*p ) * as.numeric(obsVAF)) / p, 7 ));
}

multiplicityHat <- function(majorCN, m) {
   mHat <- floor(m+0.5)
   
   if (mHat == 0) {
      return(1);
   } else if (mHat > majorCN) {
      return(majorCN);
   } else
      return(mHat);
}

multiplicityHatC2 <- function(theta1, m) {
   mHat <- 1
   mHatPlusTheta1 <- mHat + theta1$cellular_prevalence
   diff <- min(abs(m - mHat), abs(m - mHatPlusTheta1))
   
   for (mHat.tmp in 2:theta1$major_cn) {
      mHatPlusTheta1.tmp <- mHat.tmp + theta1$cellular_prevalence
   	  diff.tmp <- min(abs(m - mHat.tmp), abs(m - mHatPlusTheta1.tmp))
   	  
   	  if (diff.tmp < diff) {
   	  	  mHat <- mHat.tmp
   	  	  mHatPlusTheta1 <- mHatPlusTheta1.tmp
   	  	  diff <- diff.tmp
   	  }
   }
   
   return(mHatPlusTheta1)
}

expVAF <- function(p, CN, mHat) {
   return( round(mHat * p / ( 2*(1-p) + CN*p ), 7) );
}

expVAFC2 <- function(p, CN, mHatPlusTheta1, m, theta1) {
   mHat <- mHatPlusTheta1 - theta1$cellular_prevalence
   r <- mHat
   
   if (abs(m - mHatPlusTheta1) < abs(m - mHat))
      r <- mHatPlusTheta1
   	
   return( round(p * r / ( 2*(1-p) + CN*p ), 7) );
}

##
## Main
args <- commandArgs(T)

vcf <- read.table(args[1], header=F, sep="\t", fill=T, as.is=T, comment.char="#")
names(snv) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")

segments <- read.table(args[2], header=T, sep="\t", fill=T, as.is=T, comment.char="#")

purityploidy <- read.table(args[3], header=T, sep="\t", fill=T, as.is=T, comment.char="#")
purity <- purityploidy[1, 1]

## Output file *_muts_expAF.txt
expAFs <- initExpAF(0)

# -------------------------------------------------------
# 3.1 Calculate multiplicities in clonal SNVs
# -------------------------------------------------------
segments.c1 <- segments[segments$clonal_frequency == 1,]

for (y in 1:nrow(segments.c1)) {
   segment <- segments.c1[y,]
   vcf.segment <- vcfGetSegment(vcf, segment$chromosome, segment$start, segment$end)
   	  
   if (nrow(vcf.segment) != 0 && segment$copy_number != 0) {   ## ADDED 10/02/16: CN != 0
   	  expAF <- initExpAFfromVCF(sample, vcf.segment)
   	  expAF$iCN <- segment$copy_number                         ## copy_number = minor_cn + major_cn
   	  
   	  ## VAFexp (c1)
   	  expAF$Mut_Copies_Raw <- mapply(v = 1:nrow(vcf.segment), function(v) multiplicity(purity, expAF$iCN[v], expAF$AF_obs[v]))
   	  expAF$Mut_Copies     <- mapply(v = 1:nrow(vcf.segment), function(v) multiplicityHat(segment$major_cn, expAF$Mut_Copies_Raw[v]))
   	  expAF$AF_exp         <- mapply(v = 1:nrow(vcf.segment), function(v) expVAF(purity, expAF$iCN[v], expAF$Mut_Copies[v]))
   	     
   	  expAFs <- rbind(expAFs, expAF)
   }
}

# -------------------------------------------------------
# 3.2 Calculate multiplicities in subclonal SNVs
# -------------------------------------------------------
segments.c2 <- segments[segments$clonal_frequency != 1,]
      
if (nrow(segments.c2) != 0) {
   for (y in seq(1, nrow(segments.c2), 2)) {
   	  segment <- segments.c2[c(y,y+1),]
      vcf.segment <- vcfGetSegment(vcf, segment$chromosome[1], segment$start[1], segment$end[1])
   	     
   	  if (nrow(vcf.segment) != 0 && segment$copy_number[1] != 0) {   ## ADDED 10/02/16: CN != 0 (Although seems not necessary)
   	     expAF <- initExpAFfromVCF(sample, vcf.segment)
   	     expAF$Is_Subclonal_CN <- 1
   	     expAF$iCN <- (segment$copy_number[1] * segment$cellular_prevalence[1]) + (segment$copy_number[2] * segment$cellular_prevalence[2])
   	     #expAF$iCN <- 4 * 0.235247 + 3 * 0.764753
   	  
   	     ## VAFexp (c2)
   	     expAF$Mut_Copies_Raw <- mapply(v = 1:nrow(vcf.segment), function(v) multiplicity(purity, expAF$iCN[v], expAF$AF_obs[v]))
   	        
   	     theta1 <- segment[1,]
         if (theta1$copy_number < segment[2,]$copy_number)
            theta1 <- segment[2,]
   	     expAF$Mut_Copies <- mapply(v = 1:nrow(vcf.segment), function(v) multiplicityHatC2(theta1, expAF$Mut_Copies_Raw[v]))
   	     expAF$AF_exp     <- mapply(v = 1:nrow(vcf.segment), function(v) expVAFC2(purity, expAF$iCN[v], expAF$Mut_Copies[v], expAF$Mut_Copies_Raw[v], theta1))
   	     
   	     expAFs <- rbind(expAFs, expAF)
   	  }
   }
}

write.table(expAFs, "sclust_muts_expAF.txt", col.names=names(expAFs), row.names=F, quote=F, sep="\t")
