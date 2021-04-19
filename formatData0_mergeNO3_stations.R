########################################################
## Merge River Data with lat/lon for spatial analysis ##
########################################################

## Libraries
library(tidyverse)

## Read in the River data
river <- read.csv("data/FullRiverDataV2_NO3.csv")

## Read in station info
station <- read.csv("data/FullRiverDataV2_stationOnly.csv")

## Merge River data with station information
fullData <- left_join(x=river, y=station, by="CDSTATIONM")

head(fullData)

## Save full data frame
write.table(file="data/FullRiverDataV2.csv", x=fullData, row.names=FALSE, sep=",") 

