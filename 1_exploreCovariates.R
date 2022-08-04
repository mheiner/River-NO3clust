library(tidyverse)

source("scripts/0_readData.R")

## Catchment area not included in clustering model, but appears in imputation.
## Data source for catchment area: https://www.data.gouv.fr/ fr/datasets/bd-alti-r-75- m-250-m-1-000-m/


X <- rivers %>% group_by(CDSTATIONM) %>%
  summarize(Department=unique(LBDEPARTEM), lat=unique(lat)[1], lon=unique(lon)[1], 
            area=mean(area),
            # altitude=unique(altitude)[1], 
            elevation=mean(elevation+2),
            bioGeoRegion=unique(bioGeoRegion), IDPR=mean(IDPR),
            annual_specific_runoff=mean(annual_specific_runoff),
            p_sedimentay=mean(p_sedimentay),
            p_waterbody=mean(p_waterbody),
            p_wetland=mean(100 - p_waterbody - p_not_wetland),
            p_agricole_tot=mean(p_agricole_tot),
            p_forest=mean(p_forest),
            p_urban=mean(p_urban),
            p_other_land_uses=mean(p_other_land_uses))

X.land <- rivers %>% group_by(CDSTATIONM) %>%
  summarize(p_agricole_tot=mean(p_agricole_tot), p_forest=mean(p_forest),
            p_urban=mean(p_urban), p_other_land_uses=mean(p_other_land_uses)) %>%
  dplyr::select(-CDSTATIONM) %>% scale(.)
V.land <- eigen(cor(X.land))$vectors[,1:(ncol(X.land)-1)]
X <- cbind(X, X.land%*%V.land)
colnames(X)[(ncol(X)-2):ncol(X)] <- paste0("V", 1:3)

head(X)
dim(X)

summary(X)

library("corrplot")
corrplot.mixed(cor(X[,-grep("lat|lon|altitude|Depart|Intercept|CDSTATIONM|V", colnames(X))] %>% 
                     mutate(log_area=log(area), 
                            log_elev=log(elevation+2), log_runoff=log(annual_specific_runoff)) %>% 
                     select(., -c("area", "elevation", "annual_specific_runoff")) %>%
                     model.matrix(~. - 1 + bioGeoRegion*log_runoff, data=.)
),
tl.pos="lt") # can't have this interaction

corrplot.mixed(cor(X[,-grep("lat|lon|altitude|Depart|Intercept|CDSTATIONM|V", colnames(X))] %>% 
                     mutate(log_area=log(area), log_elev=log(elevation+2), log_runoff=log(annual_specific_runoff)) %>% 
                     select(., -c("area", "elevation", "annual_specific_runoff")) %>%
                     model.matrix(~. - 1 + I(log_runoff*(log_elev < quantile(log_elev, 0.2))), data=.)
),
tl.pos="lt")

corrplot.mixed(cor(X[,-grep("lat|lon|altitude|Depart|Intercept|CDSTATIONM|V", colnames(X))] %>% 
                     mutate(log_area=log(area), log_elev=log(elevation+2), log_runoff=log(annual_specific_runoff)) %>% 
                     select(., -c("area", "elevation", "annual_specific_runoff")) %>%
                     model.matrix(~. - 1 + log_runoff*p_agricole_tot, data=.)
),
tl.pos="lt")


pdf(file="plots/corrplotX.pdf", width=9, height=8.4)
corrplot.mixed(cor(X[,-grep("lat|lon|altitude|Depart|Intercept|CDSTATIONM|V", colnames(X))] %>% 
                     mutate(log_area=log(area), log_elev=log(elevation+2), log_runoff=log(annual_specific_runoff)) %>% 
                     select(., -c("area", "elevation", "annual_specific_runoff")) %>%
                     model.matrix(~. - 1, data=.)
  ),
  tl.pos="lt")
dev.off()




table(X$bioGeoRegion)

france = map_data("world", region = "France")
france.plot = ggplot() + geom_polygon(data=france, fill="white", colour = "black", aes(x = long, y = lat, group = group)) +
  coord_fixed(1.3)



library("RColorBrewer")


france.plot + 
  geom_point(data=X, aes(x=lon, y=lat, color=(elevation+2)), size=.6) +
  labs(color="Elevation") + 
  scale_color_distiller(palette="RdYlBu", trans="log")

france.plot.dep = ggplot() + geom_polygon(data=map_data("france"), fill="white", colour = "black", aes(x = long, y = lat, group = group)) +
  coord_fixed(1.3)


# scale_color_gradient2(limits=c(0,2000), midpoint=100, low="forestgreen", mid="orange", high="yellow")
# scale_color_gradient2(limits=c(-0.5,1.0), midpoint=0.25, low="yellow", mid="green", high="blue")
# scale_color_gradientn(limits=c(-0.5, 1.0), colors=rev(brewer.pal(11, name="Spectral")))

hist(X$elevation) # log was used in  model
france.plot + 
  geom_point(data=X, aes(x=lon, y=lat, color=(elevation+2)), size=.6) +
  labs(color="Elevation") + 
  scale_color_distiller(palette="RdYlBu", trans="log")

plt = france.plot + # france.plot.dep too messy 
  geom_point(data=X, aes(x=lon, y=lat, color=(bioGeoRegion)), size=.6) +
  labs(color="bioGeoRegion")  + coord_cartesian( xlim = c(-5, 7.8)) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
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
        )

ggsave("plots/XbioGeoRegion.pdf", plt, height=5, width=6.5)

summary(X)

hist(X$IDPR)
france.plot + 
  geom_point(data=X, aes(x=lon, y=lat, color=(IDPR)), size=.6) +
  labs(color="IDPR") + 
  scale_color_distiller(palette="RdYlBu")


hist(X$annual_specific_runoff) # not logged in model
hist(log(X$annual_specific_runoff))
france.plot + 
  geom_point(data=X, aes(x=lon, y=lat, color=(annual_specific_runoff)), size=.6) +
  labs(color="Annual Specific Runoff") + 
  scale_color_distiller(palette="RdYlBu", trans="log")


hist(X$p_sedimentay)
france.plot + 
  geom_point(data=X, aes(x=lon, y=lat, color=(p_sedimentay)), size=.6) +
  labs(color="% Sedimentary") + 
  scale_color_distiller(palette="RdYlBu")

hist(X$p_waterbody)
france.plot + 
  geom_point(data=X, aes(x=lon, y=lat, color=(p_waterbody + 0.01)), size=.6) +
  labs(color="% Water body") + 
  scale_color_distiller(palette="RdYlBu", trans="log")

summary(X)

hist(X$p_wetland)
france.plot + 
  geom_point(data=X, aes(x=lon, y=lat, color=(p_wetland + 0.01)), size=.6) +
  labs(color="% Wetland") + 
  scale_color_distiller(palette="RdYlBu", trans="log")
france.plot + 
  geom_point(data=X, aes(x=lon, y=lat, color=(p_wetland)), size=.6) +
  labs(color="% Wetland") + 
  scale_color_distiller(palette="RdYlBu")

hist(X$p_agricole_tot)
france.plot + 
  geom_point(data=X, aes(x=lon, y=lat, color=(p_agricole_tot)), size=.6) +
  labs(color="% Agricole") + 
  scale_color_distiller(palette="RdYlBu")

hist(X$p_forest)
france.plot + 
  geom_point(data=X, aes(x=lon, y=lat, color=(p_forest)), size=.6) +
  labs(color="% Forest") + 
  scale_color_distiller(palette="RdYlBu")


hist(X$p_urban)
france.plot + 
  geom_point(data=X, aes(x=lon, y=lat, color=(p_urban + 0.01)), size=.6) +
  labs(color="% Urban") + 
  scale_color_distiller(palette="RdYlBu", trans="log")

hist(X$p_other_land_uses)
france.plot + 
  geom_point(data=X, aes(x=lon, y=lat, color=(p_other_land_uses + 0.01)), size=.6) +
  labs(color="% Other land uses") + 
  scale_color_distiller(palette="RdYlBu", trans="log")


hist(X$V1)
france.plot + 
  geom_point(data=X, aes(x=lon, y=lat, color=(V1)), size=.6) +
  labs(color="Land PC 1") + 
  scale_color_distiller(palette="RdYlBu")

hist(X$V2)
hist(((X$V2 + 1)^0.15))
france.plot + 
  geom_point(data=X, aes(x=lon, y=lat, color=((V2+1)^0.15)), size=.6) +
  labs(color="Land PC 2 (transformed)") + 
  scale_color_distiller(palette="RdYlBu")

hist(X$V3)
hist(((X$V3 + 3)^0.15))
france.plot + 
  geom_point(data=X, aes(x=lon, y=lat, color=V3), size=.6) +
  labs(color="Land PC 3") + 
  scale_color_distiller(palette="RdYlBu")



plot(log(X$annual_specific_runoff), log(X$elevation))
ggplot(X, aes(x=log(elevation+2), y=log(annual_specific_runoff), color=bioGeoRegion)) + geom_point()
ggplot(X, aes(x=log(elevation+2), y=log(annual_specific_runoff), color=p_agricole_tot)) + geom_point()

