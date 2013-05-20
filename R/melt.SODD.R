#'Convert SODD output to a data frame
#'@export
#'@import reshape2 data.table
melt.SODD <- function(data, varnames=names(dimnames(data)),na.rm=FALSE,
                      value.name="Population") {
  data.df <- reshape2:::melt.array(data=data, varnames=varnames, na.rm=na.rm,
                                   value.name=value.name)
  DT <- data.table(data.df)
  DT[, c("Species", "SizeClass", "Infected")  := as.list(strsplit(as.character(Class), "\\.")[[1]]), by=Class ]
  data.df <- as.data.frame(DT)
  data.df$Species <- as.factor(data.df$Species)
  data.df$SizeClass <- as.factor(data.df$SizeClass)
  data.df$Infected <- as.factor(data.df$Infected)
  data.df <- data.df[,c(1,2,4,6,7,8,5)]
  return(data.df)
}
  

library(data.table)

