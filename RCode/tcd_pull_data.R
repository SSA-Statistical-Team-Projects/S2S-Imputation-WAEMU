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


tcd.base <- sf::st_read(dsn = "InputData", layer = "afr_tcd_l02")

tcd.grid <- SAEplus::gengrid(dsn = "InputData", layer = "afr_tcd_l02",
                             raster_tif = "tcd_ppp_2020_UNadj_constrained.tif",
                             grid_shp=T,
                             featname="population",
                             drop_Zero=F)

sf::st_write(tcd.grid$polygon_dt, dsn = "InputData", layer = "tcd_poppoly",
             driver = "ESRI Shapefile", append = FALSE)


## pull in google earth engine
#### pull in the night time lights
tcd.ntl <- SAEplus::gee_datapull(email = "dasalm20@gmail.com",
                                 gee_boundary = "users/dasalm20/afr_tcd_l02",
                                 gee_polygons = "users/dasalm20/tcd_poppoly",
                                 gee_datestart = "2018-09-01",
                                 gee_dateend = "2018-12-31",
                                 gee_desc = "TCD_NTL_2018SepDec",
                                 ldrive_dsn = "InputData/TCD_NTL_2018SepDec")

tcd.ntl2 <- SAEplus::gee_datapull(email = "dasalm20@gmail.com",
                                  gee_boundary = "users/dasalm20/afr_tcd_l02",
                                  gee_polygons = "users/dasalm20/tcd_poppoly",
                                  gee_datestart = "2019-04-01",
                                  gee_dateend = "2019-06-30",
                                  gee_desc = "TCD_NTL_2019AprJun",
                                  ldrive_dsn = "InputData/TCD_NTL_2019AprJun")

##name the ntl indicator real quick
tcdntl_sepdec.dt <- sf::st_read(dsn = "InoutData", layer = "TCD_NTL_2018SepDec")
tcdntl_aprjun.dt <- sf::st_read(dsn = "InputData", layer = "TCD_NTL_2019AprJun")
#
names(tcdntl_aprjun.dt)[names(tcdntl_aprjun.dt) == "mean"] <- "mean_ntlaprjun"
ames(tcdntl_sepdec.dt)[names(tcdntl_sepdec.dt) == "mean"] <- "mean_ntlsepdec"
#
#
#### pull in the NO2 data
tcd.no2 <- gee_datapull(email = "dasalm20@gmail.com",
                        gee_boundary = "users/dasalm20/afr_tcd_l02",
                        gee_polygons = "users/dasalm20/tcd_poppoly",
                        gee_band = "tropospheric_NO2_column_number_density",
                        gee_dataname = "COPERNICUS/S5P/NRTI/L3_NO2",
                        gee_datestart = "2018-09-01",
                        gee_dateend = "2018-12-31",
                        gee_desc = "TCD_NO2_2018SepDec",
                        ldrive_dsn = "InputData/TCD_NO2_2018SepDec")

tcdno2_sepdec.dt <- sf::st_read(dsn = "InputData", layer = "TCD_NO2_2018SepDec")
names(tcdno2_sepdec.dt)[names(tcdno2_sepdec.dt) == "mean"] <- "mean_no2sepdec"

tcd.no21 <- gee_datapull(email = "dasalm20@gmail.com",
                         gee_boundary = "users/dasalm20/afr_tcd_l02",
                         gee_polygons = "users/dasalm20/tcd_poppoly",
                         gee_band = "tropospheric_NO2_column_number_density",
                         gee_dataname = "COPERNICUS/S5P/NRTI/L3_NO2",
                         gee_datestart = "2019-04-01",
                         gee_dateend = "2019-06-30",
                         gee_desc = "TCD_NO2_2019AprJun",
                         gdrive_folder = "/SAEplus",
                         ldrive_dsn = "InputData/TCD_NO2_2019AprJun")

# pull the data on Landcover
SAEplus::gee_pullimage(email= "dasalm20@gmail.com",
                       gee_polygons = "users/dasalm20/TCD_NTL_2018SepDec",
                       gee_band = c("tree-coverfraction","urban-coverfraction","grass-coverfraction",
                                    "shrub-coverfraction","crops-coverfraction","bare-coverfraction",
                                    "water-permanent-coverfraction","water-seasonal-coverfraction",
                                    "moss-coverfraction"),
                       gee_dataname = "COPERNICUS/Landcover/100m/Proba-V-C3/Global/2019",
                       gdrive_folder = "/SAEplus2",
                       gee_desc = "TCD_LC_2018",
                       ldrive_dsn = "InputData/TCD_LC_2018")




## pull the data on impervious surface
SAEplus::gee_pullimage(email = "dasalm20@gmail.com",
                       gee_polygons = "users/dasalm20/tcd_poppoly",
                       gee_band = "change_year_index",
                       gee_dataname = "Tsinghua/FROM-GLC/GAIA/v10",
                       gee_desc = "TCD_IS_2018",
                       ldrive_dsn = "InputData/TCD_IS_2018")





# include other GEE data on CO, global human modification, gridmet drought data
SAEplus::gee_datapull(gee_boundary = "users/dasalm20/afr_gnb_l04",
                      gee_polygons = "users/dasalm20/gnb_poppoly",
                      gee_band = c("CO_column_number_density", "H2O_column_number_density",
                                   "cloud_height"),
                      gee_dataname = "COPERNICUS/S5P/NRTI/L3_CO",
                      gee_datestart = "2018-09-01",
                      gee_dateend = "2018-12-31",
                      gee_desc = "GNB_CO_2018SepDec",
                      ldrive_dsn = "SAEplus2/GNB_CO_2018SepDec",
                      gee_crs = "WGS84")

SAEplus::gee_datapull(gee_boundary = "users/dasalm20/afr_gnb_l04",
                      gee_polygons = "users/dasalm20/gnb_poppoly",
                      gee_band = c("tree-coverfraction","urban-coverfraction","grass-coverfraction",
                                   "shrub-coverfraction","crops-coverfraction","bare-coverfraction",
                                   "water-permanent-coverfraction","water-seasonal-coverfraction",
                                   "moss-coverfraction"),
                      gee_dataname = "COPERNICUS/Landcover/100m/Proba-V-C3/Global",
                      gee_datestart = "2019-04-01",
                      gee_dateend = "2019-06-30",
                      gee_desc = "GIN_LC_2019AprJun",
                      ldrive_dsn = "SAEplus2/GIN_LC_2019AprJun")







###### Next, we try to pull in all the building data
#available.dt <- SAEplus::wpopbuilding_vcheck()
#gnb.building <- wpopbuilding_pull(iso = "GNB", ldrive_dsn = "InputData")


## combine all building stats
### list all building tif data in the GIN folder
### take the sums of count, area, total length and then averages for density, urban, mean area, cv_area, mean length,
### cv length,
dt <- SAEplus::gengrid(dsn = "InputData",
                       layer = "afr_gnb_l04",
                       raster_tif = "GNB_buildings_v2_0_count.tif",
                       grid_shp = T,
                       featname = "bld_count",
                       drop_Zero = F)
gin.bld.dt <- as.data.table(dt$polygon_dt)
dt <- SAEplus::gengrid(dsn = "INputData",
                       layer = "afr_gnb_l04",
                       raster_tif = "GNB_buildings_v2_0_cv_area.tif",
                       stats = "mean",
                       grid_shp = T,
                       featname = "bld_cvarea",
                       drop_Zero = F)
add.dt <- as.data.table(dt$polygon_dt)
gin.bld.dt <- cbind(gin.bld.dt, add.dt[,"bld_cvarea"])
dt <- gengrid(dsn = "InputData",
              layer = "afr_gnb_l04",
              raster_tif = "GNB_buildings_v2_0_cv_length.tif",
              stats = "mean",
              grid_shp = T,
              featname = "bld_cvlength",
              drop_Zero = F)
add.dt <- as.data.table(dt$polygon_dt)
gin.bld.dt <- cbind(gin.bld.dt, add.dt[,"bld_cvlength"])
dt <- gengrid(dsn = "InputData",
              layer = "afr_gnb_l04",
              raster_tif = "GNB_buildings_v2_0_density.tif",
              stats = "mean",
              grid_shp = T,
              featname = "bld_density",
              drop_Zero = F)
add.dt <- as.data.table(dt$polygon_dt)
gin.bld.dt <- cbind(gin.bld.dt, add.dt[,"bld_density"])
dt <- gengrid(dsn = "InputData",
              layer = "afr_gnb_l04",
              raster_tif = "GNB_buildings_v2_0_mean_area.tif",
              stats = "mean",
              grid_shp = T,
              featname = "bld_meanarea",
              drop_Zero = F)
add.dt <- as.data.table(dt$polygon_dt)
gin.bld.dt <- cbind(gin.bld.dt, add.dt[,"bld_meanarea"])
dt <- gengrid(dsn = "InputData",
              layer = "afr_gnb_l04",
              raster_tif = "GNB_buildings_v2_0_total_length.tif",
              grid_shp = T,
              stats = "mean",
              featname = "bld_totallength",
              drop_Zero = F)
add.dt <- as.data.table(dt$polygon_dt)
gin.bld.dt <- cbind(gin.bld.dt, add.dt[,"bld_totallength"])
saveRDS(gin.bld.dt, "GNB_allbuilding.RDS")


## write the data into a shapefile format that can be used to do other work
sf::st_write(st_as_sf(gin.bld.dt, crs = "WGS84", agr = "constant"),
             layer = "GIN_allbuilding", dsn = "GIN_2021",
             driver = "ESRI Shapefile")

## pull process and join the osm data
tcd.osm <- osm_datapull(country = "Chad",
                        ldrive = "/Users/daylansalmeron/Documents/R_git_pro/S2S-Imputation-WAEMU/InputData")

gnb.lines <- SAEplus::osm_processlines(shapefile_path = "/Users/daylansalmeron/Documents/R_git_pro/S2S-Imputation-WAEMU/InputData/gnbp\_poppoly",
                                       geoid_var = "id",
                                       osm_path = "/Users/daylansalmeron/Documents/R_git_pro/S2S-Imputation-WAEMU/InputData/Guinea-Bissau_osmlines")
saveRDS(gnb.lines, file = "InputData/GIN_lines_obj.RDS")


gin.mp <- SAEplus::osm_processmp(shapefile_path = "./../S2S-REMDI/GIN_2021/GIN_allbuilding.shp",
                                 geoid_var = "id",
                                 osm_path = "C:/Users/ifean/Documents/WorldBankWork/SAEPlus_Other/Guinea_osmmp",
                                 feature_var = "amenity")

saveRDS(gin.mp, file = "./../S2S-REMDI/GIN_2021/GIN_mp_obj")

gin.points <- SAEplus::osm_processpoints(shapefile_path = "./../S2S-REMDI/GIN_2021/GIN_allbuilding.shp",
                                         geoid_var = "ADM3_CODE",
                                         osm_path = "C:/Users/ifean/Documents/WorldBankWork/SAEPlus_Other/Guinea_osmpoints")

saveRDS(gin.points, file = "./../S2S-REMDI/GIN_2021/GIN_points_obj") #save the object as RData
















