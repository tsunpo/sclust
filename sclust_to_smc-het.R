# =============================================================================
# Title: Sclust for the SMC-Het Challenge
# Name: sclust_to_smc-het.R
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 17/06/16
# =============================================================================

# -----------------------------------------------------------------------------
# Sclust - 5. Report - Output to SMC-Het format
# -----------------------------------------------------------------------------

##
## Main
args <- commandArgs(T)

#mclusters <- read.table("/Users/tpyang/Work/dev/docker/concensus/sclust_mclusters.txt", header=T, sep="\t", fill=T, as.is=T, comment.char="#")
#assignments <- read.table("/Users/tpyang/Work/dev/docker/concensus/sclust_cluster_assignments.txt", header=T, sep="\t", fill=T, as.is=T, comment.char="#")
mclusters <- read.table(args[1], header=T, sep="\t", fill=T, as.is=T, comment.char="#")
assignments <- read.table(args[2], header=T, sep="\t", fill=T, as.is=T, comment.char="#")

write.table(nrow(mclusters) - 1, "sclust_1B.txt", col.names=F, row.names=F, quote=F, sep="")
write.table(mclusters[-1, c(1,4,2)], "sclust_1C.txt", col.names=F, row.names=F, quote=F, sep="\t")
write.table(subset(assignments, Cluster_Id != 0)$Cluster_Id, "sclust_2A.txt", col.names=F, row.names=F, quote=F, sep="\t")
