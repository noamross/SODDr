#' Summarize SOD data by plot within a subsite
#' 
#' @param data Tree census data such as \code{data(SOD.data)}
#' @param subsite The subsite for which the data should be summarized
#' @param year the year for which the summary is wanted
#' @param sporulators a character vector of species which are considered sporulators
#' all other species are classified as "OTH" with no disease
#' @param size.sp Species which should be broken into size classes.  All other
#' species will only have a single sizeclass (1)
#' @param sizeclass.breaks The break points, in cm DBH, between size classes.
#' Should include lower (0) and upper (Inf) extrema.
#' @import plyr
#' @export
plot.sums <- function(data, subsite, year="2002", sporulators=c("UMCA", "LIDE"),
                      size.sp="LIDE", sizeclass.breaks=c(0,2,10,30,Inf))  {
                      
# Get a single year of live-tree data at the subsite
census <- subset(data, Year %in% year & Mortality==0 & Subsite==subsite)

# Aggregate non-sporulating species and assume no infection in these
census$Species <- factor(census$Species, 
                              levels = c(levels(census$Species), "OTH"))
census$Species[which(!(census$Species %in% sporulators))] <- "OTH"
census$Infected[which(!(census$Species %in% sporulators))] <- "S"
census <- droplevels(census)

# Assign size classes to only the LIDE trees, giving the rest SizeClass=1

census$SizeClass <- factor(1, 
                           levels=as.character(1:(length(sizeclass.breaks)-1)))
census$SizeClass[which(census$Species %in% size.sp)] <- 
  cut(census$DBH[which(census$Species %in% size.sp)], 
      breaks=sizeclass.breaks, 
      labels=as.character(1:(length(sizeclass.breaks)-1)))

# Summarize the data by plot
plot.sum <- ddply(census, c("Plot", "Species", "SizeClass", "Infected"), 
                  summarize, Count=length(Species), .drop=FALSE)

# Eliminate unneccessary classes
plot.sum <- subset(plot.sum, Species %in% size.sp | SizeClass=="1")

#Add the Site data back in
plot.sum <- merge(unique(census[,c("Site", "Subsite", "Plot",
                                        "Easting", "Northing")]), plot.sum)
return(plot.sum)

}