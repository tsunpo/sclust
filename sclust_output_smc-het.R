# =============================================================================
# Title: Sclust (for the SMC-Het Challenge)
# Name: sclust_output_smc-het.R
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 17/06/16
# =============================================================================

# -----------------------------------------------------------------------------
# Sclust - 5. Output to SMC-Het format
# -----------------------------------------------------------------------------

##
## Main
args <- commandArgs(T)

mclusters <- read.table(args[1], header=T, sep="\t", fill=T, as.is=T, comment.char="#")

assignments <- read.table(args[2], header=T, sep="\t", fill=T, as.is=T, comment.char="#")
rownames(assignments) <- assignments$Mut_ID

positions <- read.table(args[3], header=F, sep="\t", fill=T, as.is=T, comment.char="#")
positions <- as.vector(positions[,1])

## 1B
write.table(nrow(mclusters), "sclust_1B.txt", col.names=F, row.names=F, quote=F, sep="")

## 1C
c1 <- mclusters[, c(1,4,2)]
c1.tmp <- c(nrow(mclusters), length(positions) - nrow(assignments), 0)
c1 <- rbind(c1, c1.tmp)
c1[,1] <- c1[,1] + 1

for (i in 1:nrow(c1))
   if (c1[i, 3] > 1)   ## If CCF > 1, set to 1
      c1[i, 3] = 1
      
write.table(c1, "sclust_1C.txt", col.names=F, row.names=F, quote=F, sep="\t")

## 2A
a2 <- data.frame(matrix(nrow(c1), length(positions), 1))
rownames(a2) <- positions

a2.tmp <- assignments[, c("Mut_ID", "Cluster_Id")]
a2.tmp$Cluster_Id <- a2.tmp$Cluster_Id + 1
a2[rownames(a2.tmp),] <- a2.tmp$Cluster_Id

write.table(a2, "sclust_2A.txt", col.names=F, row.names=F, quote=F, sep="\t")
