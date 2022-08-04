##
## Plot Figures for the France Paper
##

## Libraries
library(tidyverse)

### EDA
rivers0 <- read.csv("data/MergedRiverData.csv")

## Convert time so that t=1 corresponds to one year
rivers0$year <- (rivers0$Date - min(rivers0$Date)) / 365.25 + 2010
# rivers0$altitude[which(is.na(rivers0$altitude))] = 0.0

rivers0 %>% filter(CDSTATIONM==6131550) %>% select(CDSTATIONM, Annee, Mois, Jour, year) %>% tail()
max(rivers0$year)
rivers0[which.min(rivers0$Date),c("Annee", "Mois", "Jour")]
rivers0[which.max(rivers0$Date),c("Annee", "Mois", "Jour")]

## Subset rivers data to only places that I have observed data
rivers <- rivers0 %>% 
  select(CDSTATIONM, LBSTATIONM, NO3, year, lon, lat, # altitude, 
         Date, elevation, LBDEPARTEM,
         area, p_agricole_tot, p_forest, p_urban, p_other_land_uses,
         p_sedimentay, p_igneous, IDPR, annual_specific_runoff, 
         p_waterbody, p_not_wetland, bioGeoRegion) %>% 
  drop_na() %>%
  arrange(CDSTATIONM)

## Combine PARIS and HAUTS-DE-SEINE
rivers <- rivers %>% mutate(LBDEPARTEM = fct_collapse(.f=LBDEPARTEM, PARIS=c("PARIS", "HAUTS-DE-SEINE")))

## Remove rivers with var(log(NO3))==0
kp.rivers <- rivers %>% group_by(CDSTATIONM) %>%
  summarize(varlog=var(log(NO3))) %>% filter(varlog!=0) %>%
  select(CDSTATIONM)
rivers <- rivers %>% filter(CDSTATIONM %in% kp.rivers$CDSTATIONM)

## Map and locations
france = map_data("world", region = "France")
france.plot <- ggplot() + geom_polygon(data=france, fill="white", colour = "black", aes(x = long, y = lat, group = group)) +
  coord_fixed(1.3)


france.plot.clean = france.plot + # france.plot.dep too messy 
  # geom_point(data=X, aes(x=lon, y=lat, color=(bioGeoRegion)), size=.6) +
  # labs(color="bioGeoRegion") + 
  coord_cartesian( xlim = c(-5, 7.8), ylim=c(42, 52)) +
  # guides(colour = guide_legend(override.aes = list(size=2))) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        # legend.position="none",
        panel.background=element_blank(),
        # panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        # plot.background=element_blank()
        plot.margin = unit(c(-0.03,0,-0.01,0), "null")
  )



locs <- rivers %>% group_by(CDSTATIONM) %>% 
  summarise(lat=unique(lat)[1], lon=unique(lon)[1], alt=unique(elevation)[1])

## Estimate time, amplitude and phase for each river
river.list <- split(rivers, f=rivers$CDSTATIONM)
get.coefs <- function(a.dframe){
  my.lm <- lm(log(NO3)~year+cos(2*pi*year)+sin(2*pi*year), 
              data=a.dframe)
  river.coef <- coef(my.lm)[-1]
  beta.time <- river.coef[1]
  amp <- sqrt(sum(river.coef[-1]^2))
  phase <- atan2(river.coef[3], river.coef[2])
  return(c(beta.time, amp, phase, sigma(my.lm), mean(log(a.dframe$NO3))))
}
river.coefs <- sapply(river.list, get.coefs) %>% t() 
colnames(river.coefs) <- c("Slope", "Amp", "Phase", "SD", "Meanlog")
river.coefs <- bind_cols(as.data.frame(river.coefs), locs)
river.coefs <- river.coefs %>%
  mutate(Slope=cut(Slope, breaks=quantile(Slope), include.lowest=TRUE),
         Amp=cut(Amp, breaks=quantile(Amp), include.lowest=TRUE))

library("RColorBrewer")
(blylrd4 = rev(brewer.pal(4, "RdYlBu")))


mnplot <- france.plot.clean + geom_point(data=river.coefs, aes(x=lon, y=lat, color=Meanlog), size=.25) +
  scale_color_distiller(palette="Spectral") + labs(color="Mean")
slopeplot <- france.plot.clean + geom_point(data=river.coefs, aes(x=lon, y=lat, color=Slope), size=.25) + 
  scale_color_manual(values=blylrd4) + guides(colour = guide_legend(reverse=T)) + labs(color="Trend")
ampplot <- france.plot.clean + geom_point(data=river.coefs, aes(x=lon, y=lat, color=Amp), size=.25) + 
  scale_color_manual(values=blylrd4) + labs(color="Amplitude") + guides(colour = guide_legend(reverse=T))
phaseplot <- france.plot.clean + geom_point(data=river.coefs, aes(x=lon, y=lat, color=Phase), size=.25) +
  scale_color_gradientn(limits=c(-pi,pi), colors=circular::circular.colors(10))
# scale_color_distiller(palette="Spectral")
s2plot <- france.plot.clean + geom_point(data=river.coefs, aes(x=lon, y=lat, color=SD), size=.25) +
  scale_color_distiller(palette="Spectral")


gridExtra::grid.arrange(slopeplot, ampplot, phaseplot, s2plot)

lay = rbind(c(rep(1,13), rep(2,10)), c(rep(3,13), rep(4,10)))
gridExtra::grid.arrange(mnplot, slopeplot, phaseplot, ampplot, layout_matrix=lay)
plt = gridExtra::grid.arrange(slopeplot, phaseplot, ampplot, s2plot, layout_matrix=lay)
ggsave(plt, filename = "plots/eda_estimates.pdf", width=6.5, height=4)





### Time series
s_indx = c(4155500, 6131550, 3080660, 6080000, 6049550, 4215750)
s_indx = c(4155500, 6131550, 3080660, 6080000, 5058935) # 5089080 5106850

# ggplot(rivers %>% filter(CDSTATIONM %in% s_indx), aes(x=year, y=log(NO3), color=factor(CDSTATIONM))) + 
#   geom_line() + geom_point() + facet_wrap(~LBSTATIONM)

p = ggplot(rivers %>% filter(CDSTATIONM %in% s_indx), aes(x=year, y=log(NO3), color=factor(LBSTATIONM))) + 
  geom_line(lwd=0.5) + geom_point(size=0.5) + geom_point(aes(x=2016.5, y=-0.25, color=factor(LBSTATIONM), shape=factor(LBSTATIONM)), size=5) +
  scale_shape_manual(values=LETTERS[1:5]) +
  facet_wrap(~toupper(LBSTATIONM)) + 
  scale_y_continuous(breaks=log(c(1,2,5,10,20,50)), labels=c(1,2,5,10,20,50)) +
  theme(legend.position = "none") +
  ylab("NO3 mg/L")

# ggplot(rivers %>% filter(CDSTATIONM %in% s_indx), aes(x=year, y=log(NO3))) + 
#   geom_line() + geom_point() + facet_wrap(~LBSTATIONM)


fr = map_data("world", region = "France")
pfr = ggplot(fr) + geom_polygon(fill="gray85", colour = "white", aes(x = long, y = lat, group = group)) 

pmap = pfr + geom_point(data=rivers %>% filter(CDSTATIONM %in% s_indx), 
                        aes(x=lon, y=lat, color=factor(LBSTATIONM), shape=factor(LBSTATIONM)), size=5) + 
  scale_shape_manual(values=LETTERS[1:5]) +
  theme_void() + theme(legend.position = "none")



library("grid")
print(p)
print(pmap, vp=viewport(width=0.375, height=0.50, x=0.855, y=0.25))
# just save from export pane
