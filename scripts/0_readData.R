rivers0 <- read.csv("data/MergedRiverData.csv")

## Convert time so that t=1 corresponds to one year
rivers0$time <- (rivers0$Date - min(rivers0$Date)) / 365.25

## Subset rivers data to only places that I have observed data
rivers <- rivers0 %>% 
  dplyr::select(CDSTATIONM, NO3, time, lon, lat, Date, elevation, LBDEPARTEM,
         area, p_agricole_tot, p_forest, p_urban, p_other_land_uses,
         p_sedimentay, p_igneous, IDPR, annual_specific_runoff, 
         p_waterbody, p_not_wetland, bioGeoRegion) %>% 
  drop_na() %>%
  arrange(CDSTATIONM) # %>% 
# filter(time<=4)

## Combine PARIS and HAUTS-DE-SEINE
rivers <- rivers %>% 
  mutate(LBDEPARTEM = fct_collapse(.f=LBDEPARTEM, PARIS=c("PARIS", "HAUTS-DE-SEINE")))

## Remove rivers with var(log(NO3))==0
kp.rivers <- rivers %>% group_by(CDSTATIONM) %>%
  summarize(varlog=var(log(NO3))) %>% filter(varlog!=0) %>%
  dplyr::select(CDSTATIONM)
rivers <- rivers %>% filter(CDSTATIONM %in% kp.rivers$CDSTATIONM)

## Outline of France
france = map_data("world", region = "France")
france.plot = ggplot() + geom_polygon(data=france, fill="white", colour = "black", aes(x = long, y = lat, group = group)) +
  coord_fixed(1.3)

## Plot the station locations
locs <- rivers %>% group_by(CDSTATIONM) %>% 
  summarise(lat=unique(lat)[1], lon=unique(lon)[1], alt=unique(elevation)[1])
# france.plot + geom_point(data=locs, aes(x=lon, y=lat))
