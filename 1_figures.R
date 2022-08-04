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


