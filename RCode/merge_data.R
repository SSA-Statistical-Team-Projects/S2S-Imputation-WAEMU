library(SAEplus)
library(sf)
library(nngeo)
library(data.table)
library(readstata13)
library(maptools)
library(rgdal)
library(sp)
library(osmdata)
library(dplyr)
library(data.table)


### merge the datasets received from google earth engine


## Mergin Data for GNB

### adding and merging nighttimelights

GNB_GEE.dt<-  data.table::setDT((sf::st_read(dsn = "InputData", layer = "GNB_NTL_2018SepDec")))
data.table::setnames(GNB_GEE.dt, "mean", "ntl_2018SepDec")

dt <-  data.table::setDT((sf::st_read(dsn = "InputData", layer = "GNB_NTL_2019AprJun")))
GNB_GEE.dt <-dt[,c("id", "mean")][GNB_GEE.dt, on = "id"]
data.table::setnames(GNB_GEE.dt, "mean", "ntl_2019AprJun")

### adding and merging NO2
dt <-  data.table::setDT((sf::st_read(dsn = "InputData", layer = "GNB_NO2_2018SepDec")))
GNB_GEE.dt <-dt[,c("id", "mean")][GNB_GEE.dt, on = "id"]
data.table::setnames(GNB_GEE.dt, "mean", "no2_2018SepDec")

dt <-  data.table::setDT((sf::st_read(dsn = "InputData", layer = "GNB_NO2_2019AprJun")))
GNB_GEE.dt <-dt[,c("id", "mean")][GNB_GEE.dt, on = "id"]
data.table::setnames(GNB_GEE.dt, "mean", "no2_2019AprJun")


### adding and merging in Land Cover data

dt <- data.table::setDT(sf::st_read(dsn = "InputData", layer = "GNB_LC_2018"))
GNB_GEE.dt <-dt[,c("id", "mean")][GNB_GEE.dt, on = "id"]
data.table::setnames(GNB_GEE.dt, "mean", "lc_2018")

### adding and merging in Impervious surface data
dt <- data.table::setDT(sf::st_read(dsn = "InputData", layer = "GNB_IS_2018"))
GNB_GEE.dt <-dt[,c("id", "mean")][GNB_GEE.dt, on = "id"]
data.table::setnames(GNB_GEE.dt, "mean", "is_2018")



### Merge All and save

#saveRDS(GNB_GEE.dt, file = "OutputData/GNB_GEE.rds")



## Merge in Data for MLI

### adding and merging nighttimelights

MLI_GEE.dt <-  data.table::setDT((sf::st_read(dsn = "InputData", layer = "MLI_NTL_2018OctDec")))
data.table::setnames(MLI_GEE.dt, "mean", "ntl_2018OctDec")

dt <- data.table::setDT((sf::st_read(dsn = "InputData", layer = "MLI_NTL_2019AprJul")))
MLI_GEE.dt<-dt[,c("id", "mean")][MLI_GEE.dt, on = "id"]
data.table::setnames(MLI_GEE.dt, "mean", "ntl_2019AprJuL")


### adding and merging NO2
dt <-  data.table::setDT((sf::st_read(dsn = "InputData", layer = "MLI_NO2_2018OctDec")))
MLI_GEE.dt<-dt[,c("id", "mean")][MLI_GEE.dt, on = "id"]
data.table::setnames(MLI_GEE.dt, "mean", "no2_2018OctDec")

dt <-  data.table::setDT((sf::st_read(dsn = "InputData", layer = "MLI_NO2_2019AprJul")))
MLI_GEE.dt<-dt[,c("id", "mean")][MLI_GEE.dt, on = "id"]
data.table::setnames(MLI_GEE.dt, "mean", "no2_2019AprJul")


### adding and merging in Land Cover data

dt <- data.table::setDT(sf::st_read(dsn = "InputData", layer = "MLI_LC_2018"))
MLI_GEE.dt <-dt[,c("id", "mean")][MLI_GEE.dt, on = "id"]
data.table::setnames(MLI_GEE.dt, "mean", "lc_2018")


### adding and merging IS Data into the NO2 Data

dt<- data.table::setDT(sf::st_read(dsn = "InputData", layer = "MLI_IS_2018"))
MLI_GEE.dt <-dt[,c("id", "mean")][MLI_GEE.dt, on = "id"]
data.table::setnames(MLI_GEE.dt, "mean", "is_2018")



## Divide Data into parts because >100MB
#MLI_GEE  <- readRDS("OutputData/MLI_GEE.rds")


MLI_GEE_part1 <- MLI_GEE.dt[1:437636,]

MLI_GEE_part2 <- MLI_GEE.dt[437637:875272,]

MLI_GEE_part3 <- MLI_GEE.dt[875273:1312907,]

saveRDS(MLI_GEE_part1, file = "OutputData/MLI_GEE_part1.rds")
saveRDS(MLI_GEE_part2, file = "OutputData/MLI_GEE_part2.rds")
saveRDS(MLI_GEE_part3, file = "OutputData/MLI_GEE_part3.rds")




### Merging Data for TCD

TCD_GEE.dt<-  data.table::setDT((sf::st_read(dsn = "InputData", layer = "TCD_NTL_2018SepDec")))
data.table::setnames(TCD_GEE.dt, "mean", "ntl_2018SepDec")

dt <-  data.table::setDT((sf::st_read(dsn = "InputData", layer = "TCD_NTL_2019AprJun")))
TCD_GEE.dt <-dt[,c("id", "mean")][TCD_GEE.dt, on = "id"]
data.table::setnames(TCD_GEE.dt, "mean", "ntl_2019AprJun")

### adding and merging NO2
dt<-  data.table::setDT((sf::st_read(dsn = "InputData", layer = "TCD_NO2_2018SepDec")))
TCD_GEE.dt <-dt[,c("id", "mean")][TCD_GEE.dt, on = "id"]
data.table::setnames(TCD_GEE.dt, "mean", "no2_2018SepDec")

dt <-  data.table::setDT((sf::st_read(dsn = "InputData", layer = "TCD_NO2_2019AprJun")))
TCD_GEE.dt <-dt[,c("id", "mean")][TCD_GEE.dt, on = "id"]
data.table::setnames(TCD_GEE.dt, "mean", "no2_2019AprJun")


### adding and merging in Land Cover data

dt <- data.table::setDT(sf::st_read(dsn = "InputData", layer = "TCD_LC_2018"))
TCD_GEE.dt <-dt[,c("id", "mean")][TCD_GEE.dt, on = "id"]
data.table::setnames(TCD_GEE.dt, "mean", "lc_2018")

### adding and merging IS Data  into the NO2 Data

dt <- data.table::setDT(sf::st_read(dsn = "InputData", layer = "TCD_IS_2018"))
TCD_GEE.dt <-dt[,c("id", "mean")][TCD_GEE.dt, on = "id"]
data.table::setnames(TCD_GEE.dt, "mean", "is_2018")


### Merge All and save

#saveRDS(TCD_GEE.dt, file = "OutputData/TCD_GEE.rds")


## Divide Data into parts because >100MB


TCD_GEE_part1 <- TCD_GEE.dt[1:438345,]

TCD_GEE_part2 <- TCD_GEE.dt[438346:876690,]

TCD_GEE_part3 <- TCD_GEE.dt[876691:1315035,]

saveRDS(TCD_GEE_part1, file = "OutputData/TCD_GEE_part1.rds")
saveRDS(TCD_GEE_part2, file = "OutputData/TCD_GEE_part2.rds")
saveRDS(TCD_GEE_part3, file = "OutputData/TCD_GEE_part3.rds")




### Merging Data for GIN


GIN_GEE.dt<-  data.table::setDT((sf::st_read(dsn = "InputData", layer = "GIN_NTL_2018JulSep")))


data.table::setnames(GIN_GEE.dt, "mean", "ntl_2018SepDec")

dt <-  data.table::setDT((sf::st_read(dsn = "InputData", layer = "GIN_NTL_2019AprJun")))
GIN_GEE.dt <-dt[,c("id", "mean")][GIN_GEE.dt, on = "id"]
data.table::setnames(GIN_GEE.dt, "mean", "ntl_2019AprJun")

#GIN_NTL <- merge(GIN_NTL1, GIN_NTL2,
         #        by= "id", all = TRUE)


### adding and merging NO2
dt <-  data.table::setDT((sf::st_read(dsn = "InputData", layer = "GIN_NO2_2018JulSep")))
GIN_GEE.dt <-dt[,c("id", "mean")][GIN_GEE.dt, on = "id"]
data.table::setnames(GIN_GEE.dt, "mean", "no2_2018SepDec")

dt <-  data.table::setDT((sf::st_read(dsn = "InputData", layer = "GIN_NO2_2019AprJun")))
GIN_GEE.dt <-dt[,c("id", "mean")][GIN_GEE.dt, on = "id"]
data.table::setnames(GIN_GEE.dt, "mean", "no2_2019AprJun")

#GIN_NO2 <- merge(GIN_NO2_1, GIN_NO2_2,
     #            by= "id", all = TRUE)



dt <- data.table::setDT(sf::st_read(dsn = "InputData", layer = "GIN_LC_2018JulSep"))
GIN_GEE.dt <-dt[,c("id", "mean")][GIN_GEE.dt, on = "id"]
data.table::setnames(GIN_GEE.dt, "mean", "lc_2018")


### adding and merging IS Data  into the NO2 Data

dt<- data.table::setDT(sf::st_read(dsn = "InputData", layer = "GIN_IS_2018JulSep"))
GIN_GEE.dt <-dt[,c("id", "mean")][GIN_GEE.dt, on = "id"]
data.table::setnames(GIN_GEE.dt, "mean", "is_2018")

#GIN_NO2 <- merge(GIN_NO2, GIN_IS,
               #  by= "id", all = TRUE)





### Merge All

#GIN_GEE.dt <- merge (GIN_NTL, GIN_NO2,
                 #    by = "id", all = TRUE)

saveRDS(GIN_GEE.dt, file = "OutputData/GIN_GEE.rds")







