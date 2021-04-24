## environment set up
remove(list = objects()) ## clear all objects in R workspace

options(
  stringsAsFactors = F, ## tell R to treat text as text, not factors
  width = 80, ## set the maximum width of outputs as 80 characters
  scipen = 6, ## discourage R from displaying numbers in scientific notation
  mc.cores = 6, ## set number of (mac) cores used in parallel processing
  start.time= Sys.time()
)

## a script to load in all the data

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

### pull in base country-level boundary files from remote source
###### starting with GUINEA

gin.base <- sf::st_read(dsn = "GIN_2021", layer = "sous_prefectures")

gin.grid <- gengrid(dsn = "GIN_2021", layer = "sous_prefectures",
                    raster_tif = "gin_ppp_2020_UNadj_constrained.tif")
sf::st_write(obj = gin.grid$polygon_dt, dsn = "GIN_2021", layer = "gin_poly", driver = "ESRI Shapefile")

gin.base <- sf::st_make_valid(gin.base)
gin.base <- gin.base[sf::st_geometry_type(gin.base$geometry) == "MULTIPOLYGON",]
gin.base <- rmapshaper::ms_simplify(gin.base)
gin.base <- sf::st_make_valid(gin.base)

sf::st_write(gin.base, dsn = "GIN_2021", layer = "sous_prefectures_valid", driver = "ESRI Shapefile",
             append = FALSE)


gin.ntl <- SAEplus::gee_datapull(gee_boundary = "users/ifeanyiedochie/sous_prefectures",
                                 gee_polygons = "users/ifeanyiedochie/gin_poly",
                                 gee_datestart = "2018-07-01",
                                 gee_dateend = "2018-09-30",
                                 gee_desc = "GIN_NTL_2018JulSep",
                                 ldrive_dsn = "GIN_2021/GIN_NTL_2018JulSep")
gin.ntl2 <- SAEplus::gee_datapull(gee_boundary = "users/ifeanyiedochie/sous_prefectures",
                                  gee_polygons = "users/ifeanyiedochie/gin_poly",
                                  gee_datestart = "2019-04-01",
                                  gee_dateend = "2019-06-30",
                                  gee_desc = "GIN_NTL_2019AprJun",
                                  ldrive_dsn = "GIN_2021/GIN_NTL_2019AprJun")
##name the ntl indicator real quick
ginntl_julsep.dt <- sf::st_read(dsn = "GIN_2021", layer = "GIN_NTL_2018JulSep_2021_04_08_12_56_31")
ginntl_aprjun.dt <- sf::st_read(dsn = "GIN_2021", layer = "GIN_NTL_2019AprJun_2021_04_08_12_50_58")

names(ginntl_aprjun.dt)[names(ginntl_aprjun.dt) == "mean"] <- "mean_ntlaprjun19"
names(ginntl_julsep.dt)[names(ginntl_julsep.dt) == "mean"] <- "mean_ntljulsep18"


gin.no2 <- gee_datapull(gee_boundary = "users/ifeanyiedochie/sous_prefectures",
                        gee_polygons = "users/ifeanyiedochie/GIN_NTL_2018JulSep_2021_04_08_12_56_31",
                        gee_band = "tropospheric_NO2_column_number_density",
                        gee_dataname = "COPERNICUS/S5P/NRTI/L3_NO2",
                        gee_datestart = "2018-07-01",
                        gee_dateend = "2018-09-30",
                        gee_desc = "GIN_NO2_2018JulSep",
                        ldrive_dsn = "GIN_2021/GIN_NO2_2018JulSep")
ginno2_julsep.dt <- sf::st_read(dsn = "GIN_2021", layer = "GIN_NO2_2018JulSep_2021_04_08_13_38_59")
names(ginno2_julsep.dt)[names(ginno2_julsep.dt) == "mean"] <- "mean_no2julsep18"

gin.no21 <- SAEplus::gee_datapull(gee_boundary = "users/ifeanyiedochie/sous_prefectures",
                                  gee_polygons = "users/ifeanyiedochie/GIN_NTL_2019AprJun_2021_04_08_12_50_58",
                                  gee_band = "tropospheric_NO2_column_number_density",
                                  gee_dataname = "COPERNICUS/S5P/NRTI/L3_NO2",
                                  gee_datestart = "2019-04-01",
                                  gee_dateend = "2019-06-30",
                                  gee_desc = "GIN_NO2_2019AprJun",
                                  ldrive_dsn = "GIN_2021/GIN_NO2_2019AprJun")
##pull in the landcover data
SAEplus::gee_datapull(gee_boundary = "users/ifeanyiedochie/sous_prefectures",
                      gee_polygons = "users/ifeanyiedochie/GIN_NTL_2018JulSep_2021_04_08_12_56_31",
                      gee_band = c("tree-coverfraction","urban-coverfraction","grass-coverfraction",
                                   "shrub-coverfraction","crops-coverfraction","bare-coverfraction",
                                   "water-permanent-coverfraction","water-seasonal-coverfraction",
                                   "moss-coverfraction"),
                      gee_dataname = "COPERNICUS/Landcover/100m/Proba-V-C3/Global",
                      gee_datestart = "2018-07-01",
                      gee_dateend = "2018-09-30",
                      gee_desc = "GIN_LC_2018JulSep",
                      ldrive_dsn = "GIN_2021/GIN_LC_2018JulSep")

SAEplus::gee_datapull(gee_boundary = "users/ifeanyiedochie/sous_prefectures",
                      gee_polygons = "users/ifeanyiedochie/GIN_NTL_2019AprJun_2021_04_08_12_50_58",
                      gee_band = c("tree-coverfraction","urban-coverfraction","grass-coverfraction",
                                   "shrub-coverfraction","crops-coverfraction","bare-coverfraction",
                                   "water-permanent-coverfraction","water-seasonal-coverfraction",
                                   "moss-coverfraction"),
                      gee_dataname = "COPERNICUS/Landcover/100m/Proba-V-C3/Global",
                      gee_datestart = "2019-04-01",
                      gee_dateend = "2019-06-30",
                      gee_desc = "GIN_LC_2019AprJun",
                      ldrive_dsn = "GIN_2021/GIN_LC_2019AprJun")


## pull the data on impervious surface
SAEplus::gee_pullimage(gee_polygons = "users/ifeanyiedochie/GIN_NTL_2018JulSep_2021_04_08_12_56_31",
                       gee_band = "change_year_index",
                       gee_dataname = "Tsinghua/FROM-GLC/GAIA/v10",
                       gee_desc = "GIN_IS_2018JulSep",
                       ldrive_dsn = "GIN_2021/GIN_IS_2018JulSep")

# include other GEE data on CO, global human modification, gridmet drought data
SAEplus::gee_datapull(gee_boundary = "users/ifeanyiedochie/sous_prefectures",
                      gee_polygons = "users/ifeanyiedochie/GIN_NTL_2018JulSep_2021_04_08_12_56_31",
                      gee_band = c("CO_column_number_density", "H2O_column_number_density",
                                   "cloud_height"),
                      gee_dataname = "COPERNICUS/S5P/NRTI/L3_CO",
                      gee_datestart = "2018-07-01",
                      gee_dateend = "2018-09-30",
                      gee_desc = "GIN_CO_2018JulSep",
                      ldrive_dsn = "GIN_2021/GIN_CO_2018JulSep",
                      gee_crs = "WGS84")

SAEplus::gee_datapull(gee_boundary = "users/ifeanyiedochie/sous_prefectures",
                      gee_polygons = "users/ifeanyiedochie/GIN_NTL_2019AprJun_2021_04_08_12_50_58",
                      gee_band = c("tree-coverfraction","urban-coverfraction","grass-coverfraction",
                                   "shrub-coverfraction","crops-coverfraction","bare-coverfraction",
                                   "water-permanent-coverfraction","water-seasonal-coverfraction",
                                   "moss-coverfraction"),
                      gee_dataname = "COPERNICUS/Landcover/100m/Proba-V-C3/Global",
                      gee_datestart = "2019-04-01",
                      gee_dateend = "2019-06-30",
                      gee_desc = "GIN_LC_2019AprJun",
                      ldrive_dsn = "GIN_2021/GIN_LC_2019AprJun")

### merge the datasets received from google earth engine

### adding and merging nighttimelights
gin_gee.dt <- data.table::setDT(sf::st_read(dsn = "GIN_2021", layer = "GIN_NTL_2018JulSep_2021_04_08_12_56_31"))
data.table::setnames(gin_gee.dt, "mean", "ntl_2018julsep")

dt <- data.table::setDT(sf::st_read(dsn = "GIN_2021", layer = "GIN_NTL_2019AprJun_2021_04_08_12_50_58"))

gin_gee.dt <- dt[,c("id", "mean")][gin_gee.dt, on = "id"]
data.table::setnames(gin_gee.dt, "mean", "ntl_2019aprjun")

### adding and merging in NO2 data
dt <- data.table::setDT(sf::st_read(dsn = "GIN_2021", layer = "GIN_NO2_2018JulSep_2021_04_08_13_38_59"))
gin_gee.dt <- dt[,c("id", "mean")][gin_gee.dt, on = "id"]
data.table::setnames(gin_gee.dt, "mean", "no2_2018julsep")

dt <- data.table::setDT(sf::st_read(dsn = "GIN_2021", layer = "GIN_NO2_2019AprJun_2021_04_13_08_15_55"))
gin_gee.dt <- dt[,c("id", "mean")][gin_gee.dt, on = "id"]
data.table::setnames(gin_gee.dt, "mean", "no2_2019aprjun")

### adding and merging in NO2
dt <- data.table::setDT(sf::st_read(dsn = "GIN_2021", layer = "GIN_LC_2018JulSep_2021_04_08_16_52_54"))
lc.names <- c("water.perm", "urban.cove", "shrub.cove", "bare.cover", "tree.cover", "crops.cove", "grass.cove",
              "moss.cover", "water.seas")

new.names <- c()
for (i in seq_along(lc.names)){
  new.names[i] <- paste(lc.names[i], "2018julsep", sep = "_")
}

data.table::setnames(dt, lc.names, new.names) ##relabel the variable names
gin_gee.dt <- dt[,c("id", new.names),with=F][gin_gee.dt, on = "id"]


dt <- data.table::setDT(sf::st_read(dsn = "GIN_2021", layer = "GIN_LC_2019AprJun_2021_04_08_17_32_39"))
lc.names <- c("water.perm", "urban.cove", "shrub.cove", "bare.cover", "tree.cover", "crops.cove", "grass.cove",
              "moss.cover", "water.seas")

new.names <- c()
for (i in seq_along(lc.names)){
  new.names[i] <- paste(lc.names[i], "2019aprjun", sep = "_")
}

data.table::setnames(dt, lc.names, new.names) ##relabel the variable names
gin_gee.dt <- dt[,c("id", new.names),with=F][gin_gee.dt, on = "id"]

### adding and merging in Impervious surface data
dt <- data.table::setDT(sf::st_read(dsn = "GIN_2021", layer = "GIN_IS_2018JulSep_2021_04_12_17_00_26"))

data.table::setnames(dt, "mean", "iscyindex_2018julsep")

gin_gee.dt <- dt[,c("id", "iscyindex_2018julsep")][gin_gee.dt, on = "id"]


names(gin_gee.dt) <- gsub("\\.", "", names(gin_gee.dt))

gin_gee.dt <- sf::st_as_sf(gin_gee.dt, crs = "WGS84", agr = "constant")


###### Next, we try to pull in all the building data
#available.dt <- SAEplus::wpopbuilding_vcheck()
#gin.building <- wpopbuilding_pull(iso = "GIN", ldrive_dsn = "GIN_2021")


## combine all building stats
### list all building tif data in the GIN folder
### take the sums of count, area, total length and then averages for density, urban, mean area, cv_area, mean length,
### cv length,
dt <- SAEplus::gengrid(dsn = "./../S2S-REMDI/GIN_2021",
                       layer = "sous_prefectures",
                       raster_tif = "GIN_buildings_v2_0_count.tif",
                       grid_shp = T,
                       featname = "bld_count",
                       drop_Zero = F)
gin.bld.dt <- as.data.table(dt$polygon_dt)
dt <- SAEplus::gengrid(dsn = "./../S2S-REMDI/GIN_2021",
                       layer = "sous_prefectures",
              raster_tif = "GIN_buildings_v2_0_cv_area.tif",
              stats = "mean",
              grid_shp = T,
              featname = "bld_cvarea",
              drop_Zero = F)
add.dt <- as.data.table(dt$polygon_dt)
gin.bld.dt <- cbind(gin.bld.dt, add.dt[,"bld_cvarea"])
dt <- gengrid(dsn = "./../S2S-REMDI/GIN_2021",
              layer = "sous_prefectures",
              raster_tif = "GIN_buildings_v2_0_cv_length.tif",
              stats = "mean",
              grid_shp = T,
              featname = "bld_cvlength",
              drop_Zero = F)
add.dt <- as.data.table(dt$polygon_dt)
gin.bld.dt <- cbind(gin.bld.dt, add.dt[,"bld_cvlength"])
dt <- gengrid(dsn = "./../S2S-REMDI/GIN_2021",
              layer = "sous_prefectures",
              raster_tif = "GIN_buildings_v2_0_density.tif",
              stats = "mean",
              grid_shp = T,
              featname = "bld_density",
              drop_Zero = F)
add.dt <- as.data.table(dt$polygon_dt)
gin.bld.dt <- cbind(gin.bld.dt, add.dt[,"bld_density"])
dt <- gengrid(dsn = "./../S2S-REMDI/GIN_2021",
              layer = "sous_prefectures",
              raster_tif = "GIN_buildings_v2_0_mean_area.tif",
              stats = "mean",
              grid_shp = T,
              featname = "bld_meanarea",
              drop_Zero = F)
add.dt <- as.data.table(dt$polygon_dt)
gin.bld.dt <- cbind(gin.bld.dt, add.dt[,"bld_meanarea"])
dt <- gengrid(dsn = "./../S2S-REMDI/GIN_2021",
              layer = "sous_prefectures",
              raster_tif = "GIN_buildings_v2_0_total_length.tif",
              grid_shp = T,
              stats = "mean",
              featname = "bld_totallength",
              drop_Zero = F)
add.dt <- as.data.table(dt$polygon_dt)
gin.bld.dt <- cbind(gin.bld.dt, add.dt[,"bld_totallength"])
saveRDS(gin.bld.dt, "GIN_allbuilding.RDS")


## write the data into a shapefile format that can be used to do other work
sf::st_write(st_as_sf(gin.bld.dt, crs = "WGS84", agr = "constant"),
             layer = "GIN_allbuilding", dsn = "GIN_2021",
             driver = "ESRI Shapefile")

## pull process and join the osm data
# gin.osm <- osm_datapull(country = "Guinea",
#                         ldrive = "C:/Users/ifean/Documents/WorldBankWork/SAEPlus_Other")

gin.lines <- SAEplus::osm_processlines(shapefile_path = "./../S2S-REMDI/GIN_2021/GIN_allbuilding.shp",
                                       geoid_var = "id",
                                       osm_path = "C:/Users/ifean/Documents/WorldBankWork/SAEPlus_Other/Guinea_osmlines")
saveRDS(gin.lines, file = "./../S2S-REMDI/GIN_2021/GIN_lines_obj.RDS")


gin.mp <- SAEplus::osm_processmp(shapefile_path = "./../S2S-REMDI/GIN_2021/GIN_allbuilding.shp",
                                 geoid_var = "id",
                                 osm_path = "C:/Users/ifean/Documents/WorldBankWork/SAEPlus_Other/Guinea_osmmp",
                                 feature_var = "amenity")

saveRDS(gin.mp, file = "./../S2S-REMDI/GIN_2021/GIN_mp_obj")

gin.points <- SAEplus::osm_processpoints(shapefile_path = "./../S2S-REMDI/GIN_2021/GIN_allbuilding.shp",
                                         geoid_var = "ADM3_CODE",
                                         osm_path = "C:/Users/ifean/Documents/WorldBankWork/SAEPlus_Other/Guinea_osmpoints")

saveRDS(gin.points, file = "./../S2S-REMDI/GIN_2021/GIN_points_obj") #save the object as RData

##load osm objects
# gin.lines <- readRDS("./../S2S-REMDI/GIN_2021/GIN_lines_obj.RDS")
# gin.mp <- readRDS("./../S2S-REMDI/GIN_2021/GIN_mp_obj")
# gin.points <- readRDS("./../S2S-REMDI/GIN_2021/GIN_points_obj")

#transform objects from long to wide first
ginline.dt <-
data.table::dcast(gin.lines[[1]], id ~ highway,
                  value.var = c("roaddensity", "count", "length"),
                  fun.aggregate = mean)

ginmp.dt <-
  data.table::dcast(gin.mp[[1]],
                    id ~ amenity,
                    value.var = "count",
                    fun.aggregate = mean)

### relabel variable names
labs <- colnames(ginmp.dt)[!(colnames(ginmp.dt) %in% "id")]

paste_tolist <- function(X, tag = "pointcount"){
  paste(X, tag, sep = "_")
}

varrelabs <- unlist(lapply(labs, paste_tolist))

setnames(ginmp.dt, labs, varrelabs)
setnames(ginmp.dt, "NA_pointcount", "unclassified_pointcount")


ginosm.dt <- ginmp.dt[ginline.dt, on = c("id")]

ginosm.dt <- gin.bld.dt[ginosm.dt, on = "id"] ### all open street maps data merged

#### merge this with the household data
# merge with household data
test <- data.table::as.data.table(readstata13::read.dta13("./../S2S-REMDI/GIN_2021/GIN-Grappe_GPS_2018.dta"))
test2 <- data.table::as.data.table(readstata13::read.dta13("./../S2S-REMDI/GIN_2021/ehcvm_welfare_GIN2018.dta"))
ginhhgeo.dt <- test[test2, on = c("grappe", "vague")]

### include geospatial data into the household data
ginhhgeo.dt <- sf::st_as_sf(ginhhgeo.dt, coords = c("coordonnes_gps__Longitude", "coordonnes_gps__Latitude"),
                            crs = 4326, agr = "constant")


### implement join
ginosm.dt <- sf::st_as_sf(ginosm.dt, crs = 4326, agr = "constant")
gin_master.dt <- sf::st_join(ginhhgeo.dt, gin_gee.dt)
gin_master.dt <- sf::st_join(gin_master.dt, ginosm.dt)

gin_master.dt <- as.data.table(gin_master.dt)

saveRDS(gin_master.dt, file = "GIN_2021/GIN_master.RDS")


###### pulling data for CHAD
tcd.base <- sf::st_read(dsn = "TCD_2021", layer = "Dept_Tchad")

tcd.grid <- SAEplus::gengrid(dsn = "TCD_2021", layer = "Dept_Tchad",
                             raster_tif = "tcd_ppp_2020_UNadj_constrained.tif")

sf::st_write(tcd.grid$polygon_dt, dsn = "TCD_2021", layer = "tcd_poppoly", driver = "ESRI Shapefile")

tcd.ntl <- SAEplus::gee_datapull(gee_boundary = "users/ifeanyiedochie/Dept_Tchad",
                                 gee_polygons = "users/ifeanyiedochie/Dept_Tchad",
                                 gee_datestart = "2018-09-01",
                                 gee_dateend = "2018-12-31",
                                 gee_desc = "TCD_NTL_2018SepDec",
                                 ldrive_dsn = "TCD_2021/TCD_NTL_2018SepDec")

tcd.ntl2 <- SAEplus::gee_datapull(gee_boundary = "users/ifeanyiedochie/Dept_Tchad",
                        gee_polygons = "users/ifeanyiedochie/Dept_Tchad",
                        gee_datestart = "2019-04-01",
                        gee_dateend = "2019-06-31",
                        gee_desc = "TCD_NTL_2019AprJun",
                        ldrive_dsn = "TCD_2021/TCD_NTL_2019AprJun")

available.dt <- SAEplus::wpopbuilding_vcheck()
tcd.building <- SAEplus::wpopbuilding_pull(iso = "TCD", ldrive_dsn = "TCD_2021")

tcd.osm <- SAEplus::osm_datapull(country = "Chad",
                                 ldrive = "C:/Users/ifean/Documents/WorldBankWork/SAEPlus_Other")

tcd.lines <- SAEplus::osm_processlines(shapefile_path = "./../S2S-REMDI/TCD_2021/Dept_Tchad.shp",
                                       geoid_var = "NOMDEP",
                                       osm_path = "C:/Users/ifean/Documents/WorldBankWork/SAEPlus_Other/Chad_osmlines")
saveRDS(tcd.lines, file = "./../S2S-REMDI/TCD_2021/TCD_lines_obj")


tcd.mp <- SAEplus::osm_processmp(shapefile_path = "./../S2S-REMDI/TCD_2021/Dept_Tchad.shp",
                                 geoid_var = "NOMDEP",
                                 osm_path = "C:/Users/ifean/Documents/WorldBankWork/SAEPlus_Other/Chad_osmmp")

saveRDS(tcd.mp, file = "./../S2S-REMDI/TCD_2021/TCD_mp_obj")

tcd.points <- SAEplus::osm_processpoints(shapefile_path = "./../S2S-REMDI/TCD_2021/Dept_Tchad.shp",
                                         geoid_var = "NOMDEP",
                                         osm_path = "C:/Users/ifean/Documents/WorldBankWork/SAEPlus_Other/Chad_osmpoints")

saveRDS(tcd.points, file = "./../S2S-REMDI/TCD_2021/TCD_points_obj") #save the object as RData


## pulling data from the Guinea Bissau

gnb.base <- sf::st_read(dsn = "GNB_2021", layer = "afr_mli_l04_2002")

gnb.grid <- SAEplus::gengrid(dsn = "GNB_2021", layer = "afr_gnb_l04_2002",
                             raster_tif = "gnb_ppp_2020_UNadj_constrained.tif")


sf::st_write(gnb.grid$polygon_dt, dsn = "GNB_2021", layer = "gnb_poppoly", driver = "ESRI Shapefile")

gnb.ntl <- SAEplus::gee_datapull(gee_boundary = "users/ifeanyiedochie/afr_gnb_l04_2002",
                                 gee_polygons = "users/ifeanyiedochie/gnb_poppoly",
                                 gee_datestart = "2018-09-01",
                                 gee_dateend = "2018-12-31",
                                 gee_desc = "TCD_NTL_2018SepDec",
                                 ldrive_dsn = "TCD_2021/TCD_NTL_2018SepDec")

gnb.ntl2 <- SAEplus::gee_datapull(gee_boundary = "users/ifeanyiedochie/afr_gnb_l04_2002",
                                  gee_polygons = "users/ifeanyiedochie/gnb_poppoly",
                                  gee_datestart = "2019-04-01",
                                  gee_dateend = "2019-06-31",
                                  gee_desc = "TCD_NTL_2019AprJun",
                                  ldrive_dsn = "TCD_2021/TCD_NTL_2019AprJun")

available.dt <- SAEplus::wpopbuilding_vcheck()
gnb.building <- SAEplus::wpopbuilding_pull(iso = "GNB", ldrive_dsn = "GNB_2021", wpversion = "v1.1")

gnb.osm <- SAEplus::osm_datapull(country = "Guinea Bissau",
                                 ldrive = "C:/Users/ifean/Documents/WorldBankWork/SAEPlus_Other")

gnb.lines <- SAEplus::osm_processlines(shapefile_path = "GNB_2021/gnb_poppoly.shp",
                                       osm_path = "C:/Users/ifean/Documents/WorldBankWork/SAEPlus_Other/Guinea Bissau_osmlines")
saveRDS(gnb.lines, file = "./../S2S-REMDI/GNB_2021/GNB_lines_obj")


gnb.mp <- SAEplus::osm_processmp(shapefile_path = "GNB_2021/gnb_poppoly.shp",
                                 osm_path = "C:/Users/ifean/Documents/WorldBankWork/SAEPlus_Other/Guinea Bissau_osmmp")

saveRDS(gnb.mp, file = "./../S2S-REMDI/GNB_2021/GNB_mp_obj")

gnb.points <- SAEplus::osm_processpoints(shapefile_path = "GNB_2021/gnb_poppoly.shp",
                                         osm_path = "C:/Users/ifean/Documents/WorldBankWork/SAEPlus_Other/Guinea Bissau_osmpoints")

saveRDS(gnb.points, file = "./../S2S-REMDI/GNB_2021/GNB_points_obj") #save the object as RData


### start pulling the data for MLI

mli.base <- sf::st_read(dsn = "MLI_2021", layer = "MLI_adm4")

mli.grid <- SAEplus::gengrid(dsn = "MLI_2021", layer = "MLI_adm4",
                             raster_tif = "mli_ppp_2020_UNadj_constrained.tif")


sf::st_write(mli.grid$polygon_dt, dsn = "MLI_2021", layer = "mli_poppoly", driver = "ESRI Shapefile")

mli.ntl <- SAEplus::gee_datapull(gee_boundary = "users/ifeanyiedochie/MLI_adm4",
                                 gee_polygons = "users/ifeanyiedochie/MLI_adm4",
                                 gee_datestart = "2018-10-01",
                                 gee_dateend = "2018-12-31",
                                 gee_desc = "MLI_NTL_2018OctDec",
                                 ldrive_dsn = "MLI_2021/MLI_NTL_2018OctDec")

mli.ntl2 <- SAEplus::gee_datapull(gee_boundary = "users/ifeanyiedochie/MLI_adm4",
                                  gee_polygons = "users/ifeanyiedochie/MLI_adm4",
                                  gee_datestart = "2019-04-01",
                                  gee_dateend = "2019-07-31",
                                  gee_desc = "MLI_NTL_2019AprJul",
                                  ldrive_dsn = "MLI_2021/MLI_NTL_2019AprJul")

available.dt <- SAEplus::wpopbuilding_vcheck()
mli.building <- SAEplus::wpopbuilding_pull(iso = "MLI", ldrive_dsn = "MLI_2021", wpversion = "v1.1")

mli.osm <- SAEplus::osm_datapull(country = "Mali",
                                 ldrive = "C:/Users/ifean/Documents/WorldBankWork/SAEPlus_Other")

mli.lines <- SAEplus::osm_processlines(shapefile_path = "MLI_2021/MLI_adm4.shp",
                                       geoid_var = "NAME_4",
                                       osm_path = "C:/Users/ifean/Documents/WorldBankWork/SAEPlus_Other/Mali_osmlines")
saveRDS(mli.lines, file = "./../S2S-REMDI/MLI_2021/mli_lines_obj")


mli.mp <- SAEplus::osm_processmp(shapefile_path = "mli_2021/MLI_adm4.shp",
                                 geoid_var = "NAME_4",
                                 osm_path = "C:/Users/ifean/Documents/WorldBankWork/SAEPlus_Other/Guinea Bissau_osmmp")

saveRDS(mli.mp, file = "./../S2S-REMDI/mli_2021/mli_mp_obj")

mli.points <- SAEplus::osm_processpoints(shapefile_path = "mli_2021/MLI_adm4.shp",
                                         geoid_var = "NAME_4",
                                         osm_path = "C:/Users/ifean/Documents/WorldBankWork/SAEPlus_Other/Guinea Bissau_osmpoints")

saveRDS(mli.points, file = "./../S2S-REMDI/mli_2021/mli_points_obj") #save the object as RData





























