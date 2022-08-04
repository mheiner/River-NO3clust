##
## Code for clustering time, amplitude and phase
## for France time series
##

## Catchment area not included in clustering model, but appears in imputation.
## Data source for catchment area: https://www.data.gouv.fr/ fr/datasets/bd-alti-r-75- m-250-m-1-000-m/

## Libraries that I need
library(tidyverse)
library(DataExplorer)
library(caret)

## Read in the Data
rivers0 <- read.csv("data/FullRiverDataV2.csv")

## Read in the covariate data
xdata <- read.csv("data/Table_descriptors_all.csv")
names(xdata)[names(xdata)=="ID"] <- "CDSTATIONM"
plot_missing(xdata)

## Impute missing p_igneous variable using RFs
my.form <- p_igneous~area+p_agricole_tot+p_forest+
  p_urban+p_wetland_veryHighConf+p_wetland_highConf+
  p_wetland_mediumConf
rf.impute <- train(form=my.form,
                   data=xdata %>% filter(complete.cases(.)),
                   method="ranger",
                   trControl=trainControl(method="repeatedcv",
                                          number=10,
                                          repeats=1))
plot(rf.impute)
xdata <- xdata %>%
  mutate(p_igneous=ifelse(is.na(p_igneous), predict(rf.impute, xdata %>% filter(is.na(p_igneous))),
                          p_igneous)) %>%
  mutate(p_sedimentay=100-p_igneous)

## Impute missing IDRP
rf.impute <- train(form=IDPR~.,
                   data=xdata %>% select(-annual_specific_runoff, -CDSTATIONM) %>%
                     filter(!is.na(IDPR)),
                   method="ranger",
                   trControl=trainControl(method="repeatedcv",
                                          number=10,
                                          repeats=1))
plot(rf.impute)
xdata <- xdata %>%
  mutate(IDPR=ifelse(is.na(IDPR), predict(rf.impute, xdata %>% filter(is.na(IDPR))),IDPR))

## Impute missing annual_specific
rf.impute <- train(form=annual_specific_runoff~.,
                   data=xdata %>% select(-CDSTATIONM) %>%
                     filter(!is.na(annual_specific_runoff)),
                   method="ranger",
                   trControl=trainControl(method="repeatedcv",
                                          number=10,
                                          repeats=1))
plot(rf.impute)
xdata <- xdata %>%
  mutate(annual_specific_runoff=ifelse(is.na(annual_specific_runoff), 
                                       predict(rf.impute, xdata %>% filter(is.na(annual_specific_runoff))),
                                       annual_specific_runoff))
plot_missing(xdata)

## Merge the two databases
rivs <- left_join(rivers0, xdata, by="CDSTATIONM")
table(rivs$p_igneous+rivs$p_sedimentay)

## Write out Full Dataset
write.csv(file="data/MergedRiverData.csv", x=rivs, row.names=FALSE)

