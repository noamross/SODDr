require(ggplot2)
require(lubridate)
require(scales)
require(reshape2)
require(grid)
nd <- read.csv("data/noam data.csv")
ndl <- ddply(nd, .(site), summarize, Northing=max(Northing))
ndl <- ndl[order(ndl$Northing),]
nd$site <- factor(nd$site, levels=ndl$site)

nd2 <- ddply(nd, .(ID, Plot, Tree, site), summarize, Year=2002:2007, Cond=as.factor(ifelse(time.to.death+2001 < Year, "Dead", ifelse(inf+2001 < Year, "Infected", "Susceptible"))))


nd3 <- ddply(nd2, .(site, Year), summarize, Susceptible=sum(Cond=="Susceptible")/(500*length(unique(Plot))), Infected=sum(Cond=="Infected")/(500*length(unique(Plot))))

nd4 <- melt(nd3, c("site", "Year"), variable.name="Cond", value.name="Density")


nd4$Year <- as.Date(ymd(paste(nd4$Year,"1","1",sep="-")))

siteseries <- ggplot(nd4, aes(x=Year, y=Density, fill=Cond)) + 
  geom_area(alpha=1) + 
  facet_wrap(~site) + 
  scale_x_date(labels = date_format("'%y")) +
  scale_fill_manual(values=c("#68A058", "#D59946")) +
  ylab(expression(paste("Tanoak density ", (m^{-2})))) +
  theme(text=element_text(family="Lato Light", size=14),
        panel.grid.major.x=element_blank(),
                panel.grid.minor.x=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_line(colour="#ECECEC", size=0.5, linetype=1),
        axis.ticks.y=element_blank(),
        panel.background=element_blank(),
        strip.background=element_blank(),
        legend.title=element_blank(),
        legend.position=c(0.85,0.1),
        legend.key.size=unit(1, "cm"),
        legend.text=element_text(size=22),
        axis.title=element_text(size=24),
        strip.text=element_text(size=16,vjust=-.25),
        axis.text=element_text(color="black",size=13))

ggsave("~/Dropbox/openquals/siteseries.svg", siteseries, height=8, width=13.5)



  