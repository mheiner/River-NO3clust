##
## Plot Figures for the France Paper
##

## Libraries
library(tidyverse)

## Functions that I need
dprojnorm <- function(angle, mn=c(0,0), log=FALSE){
  u <- cbind(cos(angle), sin(angle))
  if(length(mn)==2){
    ut.mu <- rowSums(scale(u, center=FALSE, scale=1/mn))
    mn <- matrix(mn, nrow=1)
  } else if(length(mn)>2) {
    ut.mu <- rowSums(u*mn)
  }
  pdfout <- -log(2*pi) - 0.5*rowSums(mn^2) +
    log(1+ut.mu*pnorm(ut.mu)/dnorm(ut.mu))
  if(log==TRUE){
    return(pdfout)
  } else {
    return(exp(pdfout))
  }
}

cb  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## Create a plot of a projected normal
mu <- c(-1,0)
plt.df <- data.frame(x=seq(-pi, pi, length=1000))
plt.df <- plt.df %>%
  mutate(PDF=dprojnorm(x, mu)) %>%
  mutate(Reg1=cut(x, breaks=4,  labels=1:4)) %>%
  mutate(Reg2=cut(x, breaks=c(-pi, sort(runif(3, -pi, pi)), pi),  labels=1:4, include.lowest=TRUE)) %>%
  mutate(Reg3=cut(x, breaks=c(-pi, sort(runif(3, -pi, pi)), pi),  labels=1:4, include.lowest=TRUE)) %>%
  mutate(lwr=0)
p1 <- ggplot() +
  geom_ribbon(data=plt.df %>% filter(Reg1==1), mapping=aes(x=x, ymin=lwr, ymax=PDF), fill=cb[1]) +
  geom_ribbon(data=plt.df %>% filter(Reg1==2), mapping=aes(x=x, ymin=lwr, ymax=PDF), fill=cb[2]) +
  geom_ribbon(data=plt.df %>% filter(Reg1==3), mapping=aes(x=x, ymin=lwr, ymax=PDF), fill=cb[3]) +
  geom_ribbon(data=plt.df %>% filter(Reg1==4), mapping=aes(x=x, ymin=lwr, ymax=PDF), fill=cb[4]) +
  geom_line(data=plt.df, aes(x=x, y=PDF), size=1.5) +
  xlim(-pi, pi)
p2 <- ggplot() +
  geom_ribbon(data=plt.df %>% filter(Reg2==1), mapping=aes(x=x, ymin=lwr, ymax=PDF), fill=cb[1]) +
  geom_ribbon(data=plt.df %>% filter(Reg2==2), mapping=aes(x=x, ymin=lwr, ymax=PDF), fill=cb[2]) +
  geom_ribbon(data=plt.df %>% filter(Reg2==3), mapping=aes(x=x, ymin=lwr, ymax=PDF), fill=cb[3]) +
  geom_ribbon(data=plt.df %>% filter(Reg2==4), mapping=aes(x=x, ymin=lwr, ymax=PDF), fill=cb[4]) +
  geom_line(data=plt.df, aes(x=x, y=PDF), size=1.5) +
  xlim(-pi, pi)
p3 <- ggplot() +
  geom_ribbon(data=plt.df %>% filter(Reg3==1), mapping=aes(x=x, ymin=lwr, ymax=PDF), fill=cb[1]) +
  geom_ribbon(data=plt.df %>% filter(Reg3==2), mapping=aes(x=x, ymin=lwr, ymax=PDF), fill=cb[2]) +
  geom_ribbon(data=plt.df %>% filter(Reg3==3), mapping=aes(x=x, ymin=lwr, ymax=PDF), fill=cb[3]) +
  geom_ribbon(data=plt.df %>% filter(Reg3==4), mapping=aes(x=x, ymin=lwr, ymax=PDF), fill=cb[4]) +
  geom_line(data=plt.df, aes(x=x, y=PDF), size=1.5) +
  xlim(-pi, pi)
gridExtra::grid.arrange(p1, p2, p3, nrow=1)



### for two different mu values
mu1 <- c(7.0, 5.0)
plt.df <- data.frame(u=seq(-pi, pi, length=1000))
plt.df <- plt.df %>%
  mutate(PDF=dprojnorm(u, mu1)) %>%
  mutate(Reg1=cut(u, breaks=4,  labels=1:4)) %>%
  mutate(Reg2=cut(u, breaks=seq(-pi, pi, by=pi/2),  labels=1:4, include.lowest=TRUE)) %>%
  mutate(Reg3=cut(u, breaks=seq(-pi, pi, by=pi/2),  labels=1:4, include.lowest=TRUE)) %>%
  mutate(lwr=0)
p1 <- ggplot() +
  geom_ribbon(data=plt.df %>% filter(Reg1==1), mapping=aes(x=u, ymin=lwr, ymax=PDF), fill=cb[2]) +
  geom_ribbon(data=plt.df %>% filter(Reg1==2), mapping=aes(x=u, ymin=lwr, ymax=PDF), fill=cb[1]) +
  geom_ribbon(data=plt.df %>% filter(Reg1==3), mapping=aes(x=u, ymin=lwr, ymax=PDF), fill=cb[3]) +
  geom_ribbon(data=plt.df %>% filter(Reg1==4), mapping=aes(x=u, ymin=lwr, ymax=PDF), fill=cb[4]) +
  geom_line(data=plt.df, aes(x=u, y=PDF), size=1.5) +
  xlim(-pi, pi)
p1

mu2 <- c(-1.0, -1.0)
plt.df <- data.frame(u=seq(-pi, pi, length=1000))
plt.df <- plt.df %>%
  mutate(PDF=dprojnorm(u, mu2)) %>%
  mutate(Reg1=cut(u, breaks=4,  labels=1:4)) %>%
  mutate(Reg2=cut(u, breaks=seq(-pi, pi, by=pi/2),  labels=1:4, include.lowest=TRUE)) %>%
  mutate(Reg3=cut(u, breaks=seq(-pi, pi, by=pi/2),  labels=1:4, include.lowest=TRUE)) %>%
  mutate(lwr=0)
p2 <- ggplot() +
  geom_ribbon(data=plt.df %>% filter(Reg1==1), mapping=aes(x=u, ymin=lwr, ymax=PDF), fill=cb[2]) +
  geom_ribbon(data=plt.df %>% filter(Reg1==2), mapping=aes(x=u, ymin=lwr, ymax=PDF), fill=cb[1]) +
  geom_ribbon(data=plt.df %>% filter(Reg1==3), mapping=aes(x=u, ymin=lwr, ymax=PDF), fill=cb[3]) +
  geom_ribbon(data=plt.df %>% filter(Reg1==4), mapping=aes(x=u, ymin=lwr, ymax=PDF), fill=cb[4]) +
  geom_line(data=plt.df, aes(x=u, y=PDF), size=1.5) +
  xlim(-pi, pi)
p2

gridExtra::grid.arrange(p1, p2,nrow=2)




### EDA
rivers0 <- read.csv("data/MergedRiverData.csv")

## Convert time so that t=1 corresponds to one year
rivers0$year <- (rivers0$Date - min(rivers0$Date)) / 365.25 + 2010
rivers0$altitude[which(is.na(rivers0$altitude))] = 0.0

rivers0 %>% filter(CDSTATIONM==6131550) %>% select(CDSTATIONM, Annee, Mois, Jour, year) %>% tail()
max(rivers0$year)
rivers0[which.min(rivers0$Date),c("Annee", "Mois", "Jour")]
rivers0[which.max(rivers0$Date),c("Annee", "Mois", "Jour")]

## Subset rivers data to only places that I have observed data
rivers <- rivers0 %>% 
  select(CDSTATIONM, LBSTATIONM, NO3, year, lon, lat, altitude, Date, elevation, LBDEPARTEM,
         area, p_agricole_tot, p_forest, p_urban, p_other_land_uses,
         p_sedimentay, p_igneous, IDPR, annual_specific_runoff, 
         p_waterbody, p_not_wetland, bioGeoRegion) %>% 
  drop_na() %>%
  arrange(CDSTATIONM)
# %>% 
#   filter(time<=4)

## Combine PARIS and HAUTS-DE-SEINE
rivers <- rivers %>% mutate(LBDEPARTEM = fct_collapse(.f=LBDEPARTEM, PARIS=c("PARIS", "HAUTS-DE-SEINE")))

## Remove rivers with var(log(NO3))==0
kp.rivers <- rivers %>% group_by(CDSTATIONM) %>%
  summarize(varlog=var(log(NO3))) %>% filter(varlog!=0) %>%
  select(CDSTATIONM)
rivers <- rivers %>% filter(CDSTATIONM %in% kp.rivers$CDSTATIONM)

str(rivers)
length(unique(rivers$CDSTATIONM))

s_indx = c(4155500, 6131550, 3080660, 6080000, 6049550, 4215750)
s_indx = c(4155500, 6131550, 3080660, 6080000, 5058935) # 5089080 5106850

# ggplot(rivers %>% filter(CDSTATIONM %in% s_indx), aes(x=year, y=log(NO3), color=factor(CDSTATIONM))) + 
#   geom_line() + geom_point() + facet_wrap(~LBSTATIONM)

p = ggplot(rivers %>% filter(CDSTATIONM %in% s_indx), aes(x=year, y=log(NO3), color=factor(LBSTATIONM))) + 
  geom_line(lwd=0.5) + geom_point(size=0.5) + geom_point(aes(x=2016.5, y=-0.25, color=factor(LBSTATIONM), shape=factor(LBSTATIONM)), size=5) +
  scale_shape_manual(values=LETTERS[1:5]) +
  facet_wrap(~toupper(LBSTATIONM)) + 
  theme(legend.position = "none")

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





### cluster estimates and sizes (saved from posterior analysis)
pm_amp_list = list()
pm_aclust.size_list = list()
pm_time_list = list()
pm_tclust.size_list = list()
modKs = c("5-5-12", "7-7-12")

for (kk in 1:length(modKs)) {
  load(paste0("results/pm_time_amp_", modKs[kk], "_7yr.rda"))
  pm_amp_list[[kk]] = pm_amp
  pm_aclust.size_list[[kk]] = pm_aclust.size
  pm_time_list[[kk]] = pm_time
  pm_tclust.size_list[[kk]] = pm_tclust.size
}

pm_amp_list
pm_aclust.size_list
pm_time_list
pm_tclust.size_list

daf = data.frame(mod=rep(1:length(modKs), sapply(pm_amp_list, length)), 
                 clst=unlist(sapply(pm_amp_list, function(x) 1:length(x))),
                 pm_amp = unlist(pm_amp_list),
                 pm_aclustsize = unlist(pm_aclust.size_list),
                 pm_time = unlist(pm_time_list),
                 pm_tclustsize = unlist(pm_tclust.size_list))
daf

summary(daf$pm_aclustsize)
summary(daf$pm_tclustsize)

pa = ggplot(data=daf, aes(x=pm_amp, y=as.factor(mod), size=pm_aclustsize, color=pm_aclustsize)) + geom_point() + 
  theme_bw() + ylab("K") + xlab(quote(alpha)) + xlim(c(0,1.5)) +
  scale_size_continuous(breaks=c(100, 400, 700, 1000, 1300)) + 
  scale_color_continuous(breaks=c(100, 400, 700, 1000, 1300)) + 
  scale_y_discrete(labels=c("1"="5", "2"="7")) + 
  guides(color=guide_legend(title="Cluster size"), size=guide_legend(title="Cluster size"))

ggsave(pa, file="plots/amp_vals.pdf", height=2, width=5)

pb = ggplot(data=daf, aes(x=pm_time, y=as.factor(mod), size=pm_tclustsize, color=pm_tclustsize)) + geom_point() + 
  theme_bw() + ylab("K") + xlab(quote(beta)) +
  scale_size_continuous(breaks=c(100, 700, 1000, 1500)) + 
  scale_color_continuous(breaks=c(100, 700, 1000, 1500)) +
  scale_y_discrete(labels=c("1"="5", "2"="7")) + 
  guides(color=guide_legend(title="Cluster size"), size=guide_legend(title="Cluster size"))

ggsave(pb, file="plots/time_vals.pdf", height=2, width=5)

