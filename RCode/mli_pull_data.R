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


mli.base <- sf::st_read(dsn = "InputData", layer = "afr_mli_l04")

mli.grid <- SAEplus::gengrid(dsn = "InputData", layer = "afr_mli_l04",
                             raster_tif = "mli_ppp_2020_UNadj_constrained.tif",
                             grid_shp=T,
                             featname="population",
                             drop_Zero=F)

sf::st_write(mli.grid$polygon_dt, dsn = "InputData", layer = "mli_poppoly",
             driver = "ESRI Shapefile", append = FALSE)


## pull in google earth engine
#### pull in the night time lights
mli.ntl <- SAEplus::gee_datapull(email = "dasalm20@gmail.com",
                                 gee_boundary = "users/dasalm20/afr_mli_l04",
                                 gee_polygons = "users/dasalm20/mli_poppoly",
                                 gee_datestart = "2018-10-01",
                                 gee_dateend = "2018-12-31",
                                 gee_desc = "MLI_NTL_2018OctDec",
                                 ldrive_dsn = "MLI_2021/MLI_NTL_2018OctDec")


mli.ntl2 <- SAEplus::gee_datapull(email = "dasalm20@gmail.com",
                                  gee_boundary = "users/dasalm20/afr_mli_l04",
                                  gee_polygons = "users/dasalm20/mli_poppoly",
                                  gee_datestart = "2019-04-01",
                                  gee_dateend = "2019-07-31",
                                  gee_desc = "MLI_NTL_2019AprJul",
                                  ldrive_dsn = "MLI_2021/MLI_NTL_2019AprJul")

##name the ntl indicator real quick
mlintl_sepdec.dt <- sf::st_read(dsn = "InoutData", layer = "MLI_NTL_2018OctDec")
mlintl_aprjun.dt <- sf::st_read(dsn = "InputData", layer = "MLI_NTL_2019AprJul")
#
names(mlintl_aprjun.dt)[names(mlintl_aprjul.dt) == "mean"] <- "mean_ntlaprjul"
ames(mlintl_sepdec.dt)[names(mlintl_octdec.dt) == "mean"] <- "mean_ntloctdec"
#
#
#### pull in the NO2 data
mli.no2 <- gee_datapull(email = "dasalm20@gmail.com",
                        gee_boundary = "users/dasalm20/afr_mli_l04",
                        gee_polygons = "users/dasalm20/mli_poppoly",
                        gee_band = "tropospheric_NO2_column_number_density",
                        gee_dataname = "COPERNICUS/S5P/NRTI/L3_NO2",
                        gee_datestart = "2018-10-01",
                        gee_dateend = "2018-12-31",
                        gee_desc = "MLI_NO2_2018OctDec",
                        ldrive_dsn = "InputData/MLI_NO2_2018OctDec")

mlino2_octdec.dt <- sf::st_read(dsn = "InputData", layer = "MLI_NO2_2018OctDec")
names(mlino2_octdec.dt)[names(mlino2_octdec.dt) == "mean"] <- "mean_no2octdec"

mli.no21 <- gee_datapull(email = "dasalm20@gmail.com",
                         gee_boundary = "users/dasalm20/afr_mli_l04",
                         gee_polygons = "users/dasalm20/mli_poppoly",
                         gee_band = "tropospheric_NO2_column_number_density",
                         gee_dataname = "COPERNICUS/S5P/NRTI/L3_NO2",
                         gee_datestart = "2019-04-01",
                         gee_dateend = "2019-07-31",
                         gee_desc = "MLI_NO2_2019AprJul",
                         ldrive_dsn = "InputData/MLI_NO2_2019AprJul")

#### pull in the landcover data
SAEplus::gee_datapull(email = "dasalm20@gmail.com",
                      gee_boundary = "users/dasalm20/afr_mli_l04",
                      gee_polygons = "users/dasalm20/mli_poppoly",
                      gee_band = c("tree-coverfraction","urban-coverfraction","grass-coverfraction",
                                   "shrub-coverfraction","crops-coverfraction","bare-coverfraction",
                                   "water-permanent-coverfraction","water-seasonal-coverfraction",
                                   "moss-coverfraction"),
                      gee_dataname = "COPERNICUS/Landcover/100m/Proba-V-C3/Global",
                      gee_datestart = "2018-10-01",
                      gee_dateend = "2018-12-31",
                      gee_desc = "MLI_LC_2018OctDec",
                      ldrive_dsn = "InputData/MLI_LC_2018OctDec")

SAEplus::gee_datapull(email = "dasalm20@gmail.com",
                      gee_boundary = "users/dasalm20/afr_mli_l04",
                      gee_polygons = "users/dasalm20/MLI_NTL_2019AprJul",
                      gee_band = c("tree-coverfraction","urban-coverfraction","grass-coverfraction",
                                   "shrub-coverfraction","crops-coverfraction","bare-coverfraction",
                                   "water-permanent-coverfraction","water-seasonal-coverfraction",
                                   "moss-coverfraction"),
                      gee_dataname = "COPERNICUS/Landcover/100m/Proba-V-C3/Global",
                      gee_datestart = "2019-04-01",
                      gee_dateend = "2019-07-31",
                      gee_desc = "GNB_LC_2019AprJun",
                      ldrive_dsn = "InputData/MLI_LC_2019AprJuL")


## pull the data on impervious surface
SAEplus::gee_pullimage(email = "dasalm20@gmail.com",
                       gee_polygons = "users/dasalm20/mli_poppoly",
                       gee_band = "change_year_index",
                       gee_dataname = "Tsinghua/FROM-GLC/GAIA/v10",
                       gee_desc = "MLI_IS_2018",
                       ldrive_dsn = "InputData/MLI_IS_2018")






available.dt <- SAEplus::wpopbuilding_vcheck()
mli.building <- SAEplus::wpopbuilding_pull(iso = "MLI", ldrive_dsn = "MLI_2021", wpversion = "v1.1")

gnb.osm <- SAEplus::osm_datapull(country = "Mali",
                                 ldrive = "/Users/daylansalmeron/Documents/R_git_pro/SAEPlus_other")

gnb.lines <- SAEplus::osm_processlines(shapefile_path = "MLI_2021/gnb_poppoly.shp",
                                       osm_path = "/Users/daylansalmeron/Documents/R_git_pro/SAEPlus_other/Guinea_Bissau_osmlines")
saveRDS(gnb.lines, file = "./../S2S-REMDI/MLI_2021/MLI_lines_obj")


gnb.mp <- SAEplus::osm_processmp(shapefile_path = "MLI_2021/mli_poppoly.shp",
                                 osm_path = "/Users/daylansalmeron/Documents/R_git_pro/SAEPlus_other/Mali_osmmp")

saveRDS(gnb.mp, file = "./../S2S-REMDI/MLI_2021/MLI_mp_obj")

gnb.points <- SAEplus::osm_processpoints(shapefile_path = "MLI_2021/mli_poppoly.shp",
                                         osm_path = "/Users/daylansalmeron/Documents/R_git_pro/SAEPlus_other/MLI_osmpoints")

saveRDS(gnb.points, file = "./../S2S-REMDI/GNB_2021/MLI_points_obj") #save the object as RData












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


