install.packages("data.table")           # Install and load data.table
library("data.table")


### merge the datasets received from google earth engine


## Mergin Data for GNB

### adding and merging nighttimelights

GNB_NTL1 <-  data.table::setDT((sf::st_read(dsn = "InputData", layer = "GNB_NTL_2018SepDec")))
data.table::setnames(GNB_NTL1, "mean", "ntl_2018SepDec")

GNB_NTL2 <-  data.table::setDT((sf::st_read(dsn = "InputData", layer = "GNB_NTL_2019AprJun")))
data.table::setnames(GNB_NTL2, "mean", "ntl_2018AprJun")

GNB_NTL <- merge(GNB_NTL1, GNB_NTL2,
      by= "id", all = TRUE)


### adding and merging NO2
GNB_NO2_1 <-  data.table::setDT((sf::st_read(dsn = "InputData", layer = "GNB_NO2_2018SepDec")))
data.table::setnames(GNB_NO2_1, "mean", "no2_2018SepDec")

GNB_NO2_2 <-  data.table::setDT((sf::st_read(dsn = "InputData", layer = "GNB_NO2_2019AprJun")))
data.table::setnames(GNB_NO2_2, "mean", "no2_2018AprJun")

GNB_NO2 <- merge(GNB_NO2_1, GNB_NO2_2,
                 by= "id", all = TRUE)

### adding and merging IS Data into the GNB_NO2

### adding and merging in Impervious surface data
GNB_IS <- data.table::setDT(sf::st_read(dsn = "InputData", layer = "GNB_IS_2018"))

data.table::setnames(GNB_IS, "mean", "GNB_IS_2018")

GNB_NO2 <- merge(GNB_NO2, GNB_IS,
                 by= "id", all = TRUE)

### Merge All and save

GNB_GEE.dt <- merge (GNB_NTL, GNB_NO2,
                     by = "id", all = TRUE)
saveRDS(GNB_GEE.dt, file = "OutputData/GNB_GEE.rds")



## Mergin Data for MLI

### adding and merging nighttimelights

MLI_NTL1 <-  data.table::setDT((sf::st_read(dsn = "InputData", layer = "MLI_NTL_2018OctDec")))
data.table::setnames(MLI_NTL1, "mean", "ntl_2018OctDec")

MLI_NTL2 <- data.table::setDT((sf::st_read(dsn = "InputData", layer = "MLI_NTL_2019AprJul")))
data.table::setnames(MLI_NTL2, "mean", "ntl_2018AprJuL")

MLI_NTL <- merge(MLI_NTL1, MLI_NTL2,
                 by= "id", all = TRUE)


### adding and merging NO2
MLI_NO2_1 <-  data.table::setDT((sf::st_read(dsn = "InputData", layer = "MLI_NO2_2018OctDec")))
data.table::setnames(MLI_NO_1, "mean", "no2_2018OctDec")

MLI_NO2_2 <-  data.table::setDT((sf::st_read(dsn = "InputData", layer = "MLI_NO2_2019AprJul")))
data.table::setnames(MLI_NO2_2, "mean", "no2_2018AprJul")

MLI_NO2 <- merge(MLI_NO2_1, MLI_NO2_2,
                 by= "id", all = TRUE)

### adding and merging IS Data into the NO2 Data

MLI_IS <- data.table::setDT(sf::st_read(dsn = "InputData", layer = "MLI_IS_2018"))
data.table::setnames(MLI_IS, "mean", "MLI_IS_2018")

MLI_NO2 <- merge(MLI_NO2, MLI_IS,
                 by= "id", all = TRUE)

### Merge All and save

MLI_GEE.dt <- merge (MLI_NTL, MLI_NO2,
                     by = "id", all = TRUE)
saveRDS(MLI_GEE.dt, file = "OutputData/MLI_GEE.rds")




### Merging Data for TCD

TCD_NTL1 <-  data.table::setDT((sf::st_read(dsn = "InputData", layer = "TCD_NTL_2018SepDec")))
data.table::setnames(TCD_NTL1, "mean", "ntl_2018SepDec")

TCD_NTL2 <-  data.table::setDT((sf::st_read(dsn = "InputData", layer = "TCD_NTL_2019AprJun")))
data.table::setnames(TCD_NTL2, "mean", "ntl_2018AprJun")

TCD_NTL <- merge(TCD_NTL1, TCD_NTL2,
                 by= "id", all = TRUE)


### adding and merging NO2
TCD_NO2_1 <-  data.table::setDT((sf::st_read(dsn = "InputData", layer = "TCD_NO2_2018SepDec")))
data.table::setnames(TCD_NO2_1, "mean", "no2_2018SepDec")

TCD_NO2_2 <-  data.table::setDT((sf::st_read(dsn = "InputData", layer = "TCD_NO2_2019AprJun")))
data.table::setnames(TCD_NO2_2, "mean", "no2_2018AprJun")

TCD_NO2 <- merge(TCD_NO2_1, TCD_NO2_2,
                 by= "id", all = TRUE)

### adding and merging IS Data  into the NO2 Data

TCD_IS <- data.table::setDT(sf::st_read(dsn = "InputData", layer = "TCD_IS_2018"))
data.table::setnames(TCD_IS, "mean", "TCD_IS_2018")

TCD_NO2 <- merge(TCD_NO2, TCD_IS,
                 by= "id", all = TRUE)

### Merge All and save

TCD_GEE.dt <- merge (TCD_NTL, TCD_NO2,
                     by = "id", all = TRUE)

saveRDS(TCD_GEE.dt, file = "OutputData/TCD_GEE.rds")




### Merging Data for GIN


GIN_NTL1 <-  data.table::setDT((sf::st_read(dsn = "InputData", layer = "GIN_NTL_2018JulSep")))
data.table::setnames(GIN_NTL1, "mean", "ntl_2018SepDec")

GIN_NTL2 <-  data.table::setDT((sf::st_read(dsn = "InputData", layer = "GIN_NTL_2019AprJun")))
data.table::setnames(GIN_NTL2, "mean", "ntl_2018AprJun")

GIN_NTL <- merge(GIN_NTL1, GIN_NTL2,
                 by= "id", all = TRUE)


### adding and merging NO2
GIN_NO2_1 <-  data.table::setDT((sf::st_read(dsn = "InputData", layer = "GIN_NO2_2018JulSep")))
data.table::setnames(GIN_NO2_1, "mean", "no2_2018SepDec")

GIN_NO2_2 <-  data.table::setDT((sf::st_read(dsn = "InputData", layer = "GIN_NO2_2019AprJun")))
data.table::setnames(GIN_NO2_2, "mean", "no2_2018AprJun")

GIN_NO2 <- merge(GIN_NO2_1, GIN_NO2_2,
                 by= "id", all = TRUE)

### adding and merging IS Data  into the NO2 Data

GIN_IS <- data.table::setDT(sf::st_read(dsn = "InputData", layer = "GIN_IS_2018JulSep"))
data.table::setnames(GIN_IS, "mean", "MLI_IS_2018")

GIN_NO2 <- merge(GIN_NO2, GIN_IS,
                 by= "id", all = TRUE)

### Merge All

GIN_GEE.dt <- merge (GIN_NTL, GIN_NO2,
                     by = "id", all = TRUE)

saveRDS(GIN_GEE.dt, file = "OutputData/GIN_GEE.rds")







