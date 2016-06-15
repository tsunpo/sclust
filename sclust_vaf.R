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

initExpAFfromSNV <- function(sample, snv.segment) { 
   expAF <- initExpAF(nrow(snv.segment))
   
   expAF$Chr <- snv.segment$CHROM[1]
   expAF$Position <- snv.segment$POS
   expAF$Wt  <- snv.segment$REF
   expAF$Mut <- snv.segment$ALT
   expAF$Mut_ID <- paste(sample, paste(expAF$Chr, expAF$Position, sep=":"), "SNM", sep="_")
   
   ## VAFobs
   expAF$AF_obs <- mapply(v = 1:nrow(snv.segment), function(v) obsVAF(snv.segment$INFO[v]))
   expAF$Coverage <- mapply(v = 1:nrow(snv.segment), function(v) coverage(snv.segment$INFO[v]))

   return(expAF);
}

snvGetSegment <- function(snv, chromosome, start, end) {
   snv.segment <- snv[snv$CHROM == paste("chr", chromosome, sep=""),]
   snv.segment <- snv.segment[snv.segment$POS >= start,]
   snv.segment <- snv.segment[snv.segment$POS <= end,]
   
   return(snv.segment);
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

cnaOrderByChromosome <- function(cna) {
   cna.sort <- data.frame(matrix(NA, 0, ncol(cna)))
   names(cna.sort) <- c("chromosome", "start", "end", "copy_number", "major_cn", "minor_cn", "clonal_frequency")
   
   chrs <- sort(unique(cna$chromosome))
   for (i in 1:length(chrs)) {
      cna.chr <- cna[cna$chromosome == chrs[i],]
      cna.chr <- cna.chr[order(cna.chr$start, decreasing=F),]
      
      cna.sort <- rbind(cna.sort, cna.chr)
   }
   
   return(cna.sort)
}

##
## Main
args <- commandArgs(T)

snv <- read.table(args[1], header=F, sep="\t", fill=T, as.is=T, comment.char="#")
names(snv) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")

cna <- read.table(args[2], header=T, sep="\t", fill=T, as.is=T, comment.char="#")

purityploidy <- read.table(args[3], header=T, sep="\t", fill=T, as.is=T, comment.char="#")
purity <- purityploidy[1, 1]

## Output file *_muts_expAF.txt
expAFs <- initExpAF(0)

# -------------------------------------------------------
# 3.1 Calculate multiplicities in clonal SNVs
# -------------------------------------------------------
cna.c1 <- cna[cna$clonal_frequency == 1,]

for (y in 1:nrow(cna.c1)) {
   segment <- cna.c1[y,]
   snv.segment <- snvGetSegment(snv, segment$chromosome, segment$start, segment$end)
   	  
   if (nrow(snv.segment) != 0 && segment$copy_number != 0) {   ## ADDED 10/02/16: CN != 0
   	  expAF <- initExpAFfromSNV(sample, snv.segment)
   	  expAF$iCN <- segment$copy_number                         ## copy_number = minor_cn + major_cn
   	  
   	  ## VAFexp (c1)
   	  expAF$Mut_Copies_Raw <- mapply(v = 1:nrow(snv.segment), function(v) multiplicity(purity, expAF$iCN[v], expAF$AF_obs[v]))
   	  expAF$Mut_Copies     <- mapply(v = 1:nrow(snv.segment), function(v) multiplicityHat(segment$major_cn, expAF$Mut_Copies_Raw[v]))
   	  expAF$AF_exp         <- mapply(v = 1:nrow(snv.segment), function(v) expVAF(purity, expAF$iCN[v], expAF$Mut_Copies[v]))
   	     
   	  expAFs <- rbind(expAFs, expAF)
   }
}

# -------------------------------------------------------
# 3.2 Calculate multiplicities in subclonal SNVs
# -------------------------------------------------------
cna.c2 <- cna[cna$clonal_frequency != 1,]
cna.c2 <- cnaOrderByChromosome(cna.c2)

if (nrow(cna.c2) != 0) {
   for (y in seq(1, nrow(cna.c2), 2)) {
   	  segment <- cna.c2[c(y,y+1),]
      snv.segment <- snvGetSegment(snv, segment$chromosome[1], segment$start[1], segment$end[1])
   	     
   	  if (nrow(snv.segment) != 0 && segment$copy_number[1] != 0) {   ## ADDED 10/02/16: CN != 0 (Although seems not necessary)
   	     expAF <- initExpAFfromSNV(sample, snv.segment)
   	     expAF$Is_Subclonal_CN <- 1
   	     expAF$iCN <- (segment$copy_number[1] * segment$cellular_prevalence[1]) + (segment$copy_number[2] * segment$cellular_prevalence[2])
   	     #expAF$iCN <- 4 * 0.235247 + 3 * 0.764753
   	  
   	     ## VAFexp (c2)
   	     expAF$Mut_Copies_Raw <- mapply(v = 1:nrow(snv.segment), function(v) multiplicity(purity, expAF$iCN[v], expAF$AF_obs[v]))
   	        
   	     theta1 <- segment[1,]
         if (theta1$copy_number < segment[2,]$copy_number)
            theta1 <- segment[2,]
   	     expAF$Mut_Copies <- mapply(v = 1:nrow(snv.segment), function(v) multiplicityHatC2(theta1, expAF$Mut_Copies_Raw[v]))
   	     expAF$AF_exp     <- mapply(v = 1:nrow(snv.segment), function(v) expVAFC2(purity, expAF$iCN[v], expAF$Mut_Copies[v], expAF$Mut_Copies_Raw[v], theta1))
   	     
   	     expAFs <- rbind(expAFs, expAF)
   	  }
   }
}

write.table(expAFs, "sclust_muts_expAF.txt", col.names=names(expAFs), row.names=F, quote=F, sep="\t")
