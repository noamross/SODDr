#'Convert SODD output to a data frame
#'@export
#'@import reshape2 data.table
melt.SODD <- function(data, varnames=names(dimnames(data)),na.rm=FALSE,
                      value.name="Population") {
  DT <- data.table(reshape2:::melt.array(data=data, varnames=varnames, 
                                              na.rm=na.rm,
                                              value.name=value.name))
  newCols <- c("Species", "SizeClass", "Infected") 
  DT[, c(newCols) := as.list(strsplit(as.character(Class), "\\.")[[1]]),
     by=Class ]
  DT[, c(newCols) := lapply(.SD, factor), .SDcols=newCols]
  DT[, Class := NULL]
  setcolorder(DT,c(1,2,3,5,6,7,4))
  return(DT)
}
  

