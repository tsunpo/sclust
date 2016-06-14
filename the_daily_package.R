# =============================================================================
# Name: The Daily Package (the_daily_package.R)
# Author: Tsun-Po Yang (tyang2@uni-koeln.de)
# Last Modified: 14/06/16
# =============================================================================

# -----------------------------------------------------------------------------
# Method: Read data from a file and return a data frame
# Tricks: 1) Row names will be assigned if an integer "column number" is provoded for 'rownames', otherwise the default is 1 if rownames=T.
#         2) Able to deal with empty file.
#         3) Able to deal with duplicate row names.
#         4) Screen preview the first 5x20 table of the data frame and prompts dimentions.
#         5) Return a vector if there is only one column in the file.
# Usage: file.df <- readTable("file.txt", header=T, rownames=T, sep="\t")
#        file.df <- readTable("file.txt", header=T, rownames=8, sep="\t")
# Last Modified: 05/01/16
# -----------------------------------------------------------------------------
readTable <- function(file, header, rownames, sep) {
	readTable0(file, header, rownames, sep, "#")
}

readTable0 <- function(file, header, rownames, sep, comment.char) {
   cat("Reading data from file [ ", file, " ] :\n", sep="")
   
   ## Deal with empty file
   if (is.na(file.info(file)$size)) {
      cat("Error in opening file :", "  no such file or directory", sep="\n")
      return(NA)   ## TODO: Should NOT return anything
   
   #} else if (file.info(file)$size == 0) {
   } else if (length(readLines(file)) == 0) {   ## Modified 06/01/16
      cat("Warning in `file.info(input)$size` :", "  no lines available in input, return NA", "[1] NA", sep="\n")
      return(NA)
   }
   
   ## Handle row names
   file.df <- NA
   if (is.logical(rownames) || is.numeric(rownames)) {
      ## ADDED 11/03/13: To handle duplicate 'row.names' which are not allowed"
      file.df <- tryCatch(read.table(file, header=header, sep=sep, fill=T, as.is=T, comment.char=comment.char, na.strings="NA"), 
                    error=function(e) read.table(file, header=header, row.names=NULL, sep=sep, fill=T, as.is=T, comment.char=comment.char, na.strings="NA"))

      ## ADDED 30/04/13: To handle a swift when column number is more than header (E.g. an extra tab at each rows)
      tryCatch(as.numeric(rownames(file.df)), 
         warning=function(w) cat("Warning in reading file :", "  there might be an extra seperator at end of each rows", "In addition: Error message:", "Returned data.frame shifted", sep="\n"))
      
   	  if (rownames == T && colnames(file.df)[1] != "row.names")   ## MODIFIED 11/03/13
   	     rownames(file.df) <- file.df[,1]
   	  
   	  else if (is.numeric(rownames))
   	     rownames(file.df) <- file.df[,rownames]
   	  
   	  else if (colnames(file.df)[1] == "row.names") {   ## ADDED 11/03/13
   	     colnames(file.df) <- colnames(file.df)[-1]
   	     file.df <- file.df[,1:dim(file.df)[2] - 1]   
   	     if (rownames)
   	        cat("Warning in `rownames=T` :", "  duplicate row names are not allowed", sep="\n")
   	  }
   	  
   } else
      file.df <- read.table(file, header=header, row.names=rownames, sep=sep, fill=T, as.is=T, comment.char=comment.char, na.strings="NA")   ## E.g. rownames=NULL/NA
      
   ## Preview data.frame or vector
   ncol <- ncol(file.df)
   nrow <- nrow(file.df)
   if (header == F && ncol == 1) {
      file.df <- as.vector(file.df[,1])
      
      if (nrow >= 20) nrow = 20
      println("  there is only one column in the file, return vector")
      println(paste("Preview first", nrow, "elements in returned vector :", sep=" "))
      print(file.df[1:nrow])
      
      println("Length :")       ## Prompts length
      print(length(file.df))
      
   } else {
   	  if (ncol >= 20) ncol = 20
      if (nrow >= 5) nrow = 5
      
      if (ncol > 20 || nrow > 5)
         println(paste("Preview first", nrow, "x", ncol, "subtable in returned data.frame :", sep=" "))  
      else
         println(paste("Preview", nrow, "x", ncol, "table in returned data.frame :", sep=" "))  
      print(file.df[1:nrow, 1:ncol])
      
      println("Dimensions :")   ## Prompts dimentions
      print(dim(file.df))
   }
   
   return(file.df)
}

## Examples
#list.v <- readTable("/lustre/scratch101/sanger/tpy/twins/imputed/610/PLINK2/8.s356.qced.noindel.duplicate.list", header=F, rownames=F, sep="\t")

#hg19 <- readTable("/nfs/team147/tpy/team147/meth/450k/smoking/Cardiogenics/expression/chr-start-end_hg19.bed", header=F, rownames=4, sep="\t")

#setwd("/nfs/team147/tpy/team147/meth/450k/smoking/Cardiogenics")
#pheno.t <- readTable("WP5_allPhenoData 8Apr11_tpy.txt", header=T, rownames=T, sep="\t")

# -----------------------------------------------------------------------------
# Method: Output table to a file
# Tricks: Adding an additional "ID" in the front of column names if rownames=T
# Usage: writeTable(table, "file.txt", colnames=T, rownames=T, sep="\t")
# Last Modified: 18/01/13
# -----------------------------------------------------------------------------
writeTable <- function(table, file, colnames, rownames, sep) {
   if (rownames == T)
   	  write.table(cbind(rownames(table), table), file=file, row.names=F, col.names=c("ID_REF", colnames(table)), sep=sep, quote=F, na="")
   else
      write.table(table, file=file, row.names=rownames, col.names=colnames, sep=sep, quote=F, na="")
}

# -----------------------------------------------------------------------------
# Method: 
# Tricks: 
# Usage:  toTable(NA, 8, 0, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))
# Last Modified: 09/10/13
# -----------------------------------------------------------------------------
toTable <- function(data, ncol, nrow, colnames) {
   table <- data.frame(matrix(data, nrow, ncol))
   names(table) <- colnames
   
   for (x in 1:ncol(table)) {
      table[,x] <- as.vector(table[,x])
   }
   
   return(table)
}
