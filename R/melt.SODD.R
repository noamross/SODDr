#'Convert SODD output to a data frame
#'@export
#'@import reshape2
melt.SODD <- function(data, varnames=names(dimnames(data)),na.rm=FALSE,
                      value.name="Population") {
  data.df <- reshape2:::melt.array(data, varnames, na.rm=na.rm, 
                                   value.name=value.name)
  cl <- which(names(data.df)=="Class")
  Classes <- do.call(rbind, strsplit(as.character(data.df$Class), "\\."))
  colnames(Classes) <- c("Species", "SizeClass", "Infected")
  data.df <- cbind(data.df[,1:(cl-1)],Classes,data.df[(cl+1):(ncol(data.df))])
  return(data.df)
}
  