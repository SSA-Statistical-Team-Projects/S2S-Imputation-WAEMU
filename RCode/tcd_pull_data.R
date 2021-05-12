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

tcd.no21 <- gee_datapull(gee_boundary = "users/dasalm20/afr_tcd_l02",
                         gee_polygons = "users/dasalm20/gnb_poppoly",
                         gee_band = "tropospheric_NO2_column_number_density",
                         gee_dataname = "COPERNICUS/S5P/NRTI/L3_NO2",
                         gee_datestart = "2019-04-01",
                         gee_dateend = "2019-06-30",
                         gee_desc = "TCD_NO2_2019AprJun",
                         ldrive_dsn = "InputData/TCD_NO2_2019AprJun")

#### pull in the landcover data
SAEplus::gee_datapull(gee_boundary = "users/dasalm20/afr_tcd_l02",
                      gee_polygons = "users/dasalm20/tcd_poppoly",
                      gee_band = c("tree-coverfraction","urban-coverfraction","grass-coverfraction",
                                   "shrub-coverfraction","crops-coverfraction","bare-coverfraction",
                                   "water-permanent-coverfraction","water-seasonal-coverfraction",
                                   "moss-coverfraction"),
                      gee_dataname = "COPERNICUS/Landcover/100m/Proba-V-C3/Global",
                      gee_datestart = "2018-09-01",
                      gee_dateend = "2018-12-31",
                      gee_desc = "TCD_LC_2018SepDec",
                      ldrive_dsn = "InputData/TCD_LC_2018SepDec")

SAEplus::gee_datapull(gee_boundary = "users/dasalm20/afr_tcd_l02",
                      gee_polygons = "users/dasalm20/TCD_NTL_2019AprJun",
                      gee_band = c("tree-coverfraction","urban-coverfraction","grass-coverfraction",
                                   "shrub-coverfraction","crops-coverfraction","bare-coverfraction",
                                   "water-permanent-coverfraction","water-seasonal-coverfraction",
                                   "moss-coverfraction"),
                      gee_dataname = "COPERNICUS/Landcover/100m/Proba-V-C3/Global",
                      gee_datestart = "2019-04-01",
                      gee_dateend = "2019-06-30",
                      gee_desc = "TCD_LC_2019AprJun",
                      ldrive_dsn = "InputData/TCD_LC_2019AprJun")



available.dt <- SAEplus::wpopbuilding_vcheck()
tcd.building <- SAEplus::wpopbuilding_pull(iso = "TCD", ldrive_dsn = "TCD_2021", wpversion = "v1.1")

tcd.osm <- SAEplus::osm_datapull(country = "Chad",
                                 ldrive = "/Users/daylansalmeron/Documents/R_git_pro/SAEPlus_other")

tcd.lines <- SAEplus::osm_processlines(shapefile_path = "TCD_2021/tcd_poppoly.shp",
                                       osm_path = "/Users/daylansalmeron/Documents/R_git_pro/SAEPlus_other/Chad_osmlines")
saveRDS(tcd.lines, file = "./../S2S-REMDI/GNB_2021/TCD_lines_obj")


tcd.mp <- SAEplus::osm_processmp(shapefile_path = "TCD_2021/tcd_poppoly.shp",
                                 osm_path = "/Users/daylansalmeron/Documents/R_git_pro/SAEPlus_other/Chad_osmmp")

saveRDS(tcd.mp, file = "./../S2S-REMDI/GNB_2021/TCD_mp_obj")

tcd.points <- SAEplus::osm_processpoints(shapefile_path = "tcd_2021/tcd_poppoly.shp",
                                         osm_path = "/Users/daylansalmeron/Documents/R_git_pro/SAEPlus_other/Chad_osmpoints")

saveRDS(tcd.points, file = "./../S2S-REMDI/TCD_2021/TCD_points_obj") #save the object as RData











