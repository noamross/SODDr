# This short script takes the original CSV files provided by Richasd Cobb
# And processes them for use, creating the sd data object

# Load data and combine with subsite designations
require(utils)
SOD.plots <- utils::read.csv("redwood_plot_network_full.csv", 
                             comment.char="#", sep=";")
SOD.plots <- merge(SOD.plots, utils::read.csv("plot-site list.csv", 
                                              sep=";"))
SOD.plots$Infected <- factor(ifelse(SOD.plots$Infected,"I","S"), 
                             levels=c("S","I"))
SOD.plots <- SOD.plots[,c(1,2,17,3:16)]
