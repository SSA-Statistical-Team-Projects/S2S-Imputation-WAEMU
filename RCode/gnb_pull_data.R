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


gnb.base <- sf::st_read(dsn = "InputData", layer = "afr_gnb_l04_2002")

gnb.grid <- SAEplus::gengrid(dsn = "InputData", layer = "afr_gnb_l04_2002",
                            raster_tif = "gnb_ppp_2020_UNadj_constrained.tif",
                            grid_shp = T,
                            featname = "population",
                           drop_Zero = F)


sf::st_write(gnb.grid$polygon_dt, dsn = "InputData", layer = "gnb_poppoly", driver = "ESRI Shapefile",
             append = F )

gnb.ntl <- SAEplus::gee_datapull(email = "dasalm20@gmail.com",
                                 gee_boundary = "users/dasalm20/afr_gnb_l04_2002",
                                 gee_polygons = "users/dasalm20/gnb_poppoly",
                                 gee_datestart = "2018-09-01",
                                 gee_dateend = "2018-12-31",
                                 gee_desc = "GNB_NTL_2018SepDec",
                                 ldrive_dsn = "GNB_2021/GNB_NTL_2018SepDec")

gnb.ntl2 <- SAEplus::gee_datapull(email = "dasalm20@gmail.com",
                                  gee_boundary = "users/dasalm20/afr_gnb_l04_2002",
                                  gee_polygons = "users/dasalm20/gnb_poppoly",
                                  gee_datestart = "2019-04-01",
                                  gee_dateend = "2019-06-30",
                                  gee_desc = "GNB_NTL_2019AprJun",
                                  ldrive_dsn = "GNB_2021/GNB_NTL_2019AprJun")

##name the ntl indicator real quick
 gibntl_sepdec.dt <- sf::st_read(dsn = "InoutData", layer = "GNB_NTL_2018SepDec")
gnbntl_aprjun.dt <- sf::st_read(dsn = "InputData", layer = "GNB_NTL_2019AprJun")
#
names(gnbntl_aprjun.dt)[names(gnbntl_aprjun.dt) == "mean"] <- "mean_ntlaprjun"
ames(gnbntl_sepdec.dt)[names(gnbntl_sepdec.dt) == "mean"] <- "mean_ntlsepdec"
#
#

gin.no2 <- gee_datapull( email = "dasalm20@gmail.com",
                         gee_boundary = "users/dasalm20/afr_gnb_l04_2002",
                         gee_polygons = "users/dasalm20/gnb_poppoly",
                         gee_band = "tropospheric_NO2_column_number_density",
                        gee_dataname = "COPERNICUS/S5P/NRTI/L3_NO2",
                        gee_datestart = "2018-09-01",
                        gee_dateend = "2018-12-31",
                        gee_desc = "GNB_NO2_2018SepDec",
                        ldrive_dsn = "GNB_2021/GNB_NO2_2018SepDec")

ginno2_sepdec.dt <- sf::st_read(dsn = "InputData", layer = "GNB_NO2_2018SepDec_2021_04_28_11_41_33")
names(ginno2_sepdec.dt)[names(ginno2_sepdec.dt) == "mean"] <- "mean_no2sepdec"

gin.no21 <- SAEplus::gee_datapull(email = "dasalm20@gmail.com",
                                  gee_boundary = "users/dasalm20/afr_gnb_l04_2002",
                                   gee_polygons = "users/dasalm20/gnb_poppoly",
                                  gee_band = "tropospheric_NO2_column_number_density",
                                    gee_dataname = "COPERNICUS/S5P/NRTI/L3_NO2",
                                   gee_datestart = "2019-04-01",
                                 gee_dateend = "2019-06-30",
                                  gee_desc = "GNB_NO2_2019AprJun",
                                   ldrive_dsn = "GNB_2021/GNB_NO2_2019AprJun")



# ##pull in the landcover data
SAEplus::gee_datapull(email = "dasalm20@gmail.com",
                      gee_boundary = "users/dasalm20/afr_gnb_l04_2002",
                      gee_polygons = "users/dasalm20/gnb_poppoly",
                      gee_band = c("tree-coverfraction","urban-coverfraction","grass-coverfraction",
                                   "shrub-coverfraction","crops-coverfraction","bare-coverfraction",
                                   "water-permanent-coverfraction","water-seasonal-coverfraction",
                                   "moss-coverfraction"),
                     gee_dataname = "COPERNICUS/Landcover/100m/Proba-V-C3/Global",
                      gee_datestart = "2018-09-01",
                      gee_dateend = "2018-12-31",
                      gee_desc = "GNB_LC_2018SepDec",
                      ldrive_dsn = "GNB_2021/GNB_LC_2018JulSep")

SAEplus::gee_datapull(email = "dasalm20@gmail.com",
                        gee_boundary = "users/dasalm20/afr_gnb_l04_2002",
                       gee_polygons = "users/dasalm20/gnb_poppoly",
                      gee_band = c("tree-coverfraction","urban-coverfraction","grass-coverfraction",
                                   "shrub-coverfraction","crops-coverfraction","bare-coverfraction",
                                   "water-permanent-coverfraction","water-seasonal-coverfraction",
                                   "moss-coverfraction"),
                      gee_dataname = "COPERNICUS/Landcover/100m/Proba-V-C3/Global",
                      gee_datestart = "2019-04-01",
                      gee_dateend = "2019-06-30",
                      gee_desc = "GNB_LC_2019AprJun",
                      ldrive_dsn = "GNB_2021/GNB_LC_2019AprJun")
#
#
# ## pull the data on impervious surface
 SAEplus::gee_pullimage(email = "dasalm20@gmail.com",
                       gee_polygons = "users/dasalm20/gnb_poppoly",
                       gee_band = "change_year_index",
                       gee_dataname = "Tsinghua/FROM-GLC/GAIA/v10",
                       gee_desc = "GNB_IS_2018SepDec",
                       ldrive_dsn = "GNB_2021/GNB_IS_2018SepDec")
#
## include other GEE data on CO, global human modification, gridmet drought data
SAEplus::gee_datapull( email = "dasalm20@gmail.com",
                     gee_boundary = "users/dasalm20/afr_gnb_l04_2002",
                     gee_polygons = "users/dasalm20/gnb_poppoly",
                      gee_band = c("CO_column_number_density", "H2O_column_number_density",
                                   "cloud_height"),
                      gee_dataname = "COPERNICUS/S5P/NRTI/L3_CO",
                      gee_datestart = "2018-09-01",
                      gee_dateend = "2018-12-31",
                      gee_desc = "GNB_CO_2018JulSep",
                      ldrive_dsn = "GNB_2021/GNB_CO_2018SepDec",
                      gee_crs = "WGS84")
#
SAEplus::gee_datapull(email = "dasalm20@gmail.com",
                      gee_boundary = "users/dasalm20/afr_gnb_l04_2002",
                       gee_polygons = "users/dasalm20/gnb_poppoly",
                      gee_band = c("tree-coverfraction","urban-coverfraction","grass-coverfraction",
                                   "shrub-coverfraction","crops-coverfraction","bare-coverfraction",
                                   "water-permanent-coverfraction","water-seasonal-coverfraction",
                                   "moss-coverfraction"),
                       gee_dataname = "COPERNICUS/Landcover/100m/Proba-V-C3/Global",
                      gee_datestart = "2019-04-01",
                      gee_dateend = "2019-06-30",
                      gee_desc = "GNB_LC_2019AprJun",
                     ldrive_dsn = "GNB_2021/GNB_LC_2019AprJun")








available.dt <- SAEplus::wpopbuilding_vcheck()
gnb.building <- SAEplus::wpopbuilding_pull(iso = "GNB", ldrive_dsn = "GNB_2021", wpversion = "v1.1")

gnb.osm <- SAEplus::osm_datapull(country = "Guinea Bissau",
                                 ldrive = "/Users/daylansalmeron/Documents/R_git_pro/SAEPlus_other")

gnb.lines <- SAEplus::osm_processlines(shapefile_path = "GNB_2021/gnb_poppoly.shp",
                                       osm_path = "/Users/daylansalmeron/Documents/R_git_pro/SAEPlus_other/Guinea_Bissau_osmlines")
saveRDS(gnb.lines, file = "./../S2S-REMDI/GNB_2021/GNB_lines_obj")


gnb.mp <- SAEplus::osm_processmp(shapefile_path = "GNB_2021/gnb_poppoly.shp",
                                 osm_path = "/Users/daylansalmeron/Documents/R_git_pro/SAEPlus_other/Guinea_Bissau_osmmp")

saveRDS(gnb.mp, file = "./../S2S-REMDI/GNB_2021/GNB_mp_obj")

gnb.points <- SAEplus::osm_processpoints(shapefile_path = "GNB_2021/gnb_poppoly.shp",
                                         osm_path = "/Users/daylansalmeron/Documents/R_git_pro/SAEPlus_other/Guinea_Bissau_osmpoints")

saveRDS(gnb.points, file = "./../S2S-REMDI/GNB_2021/GNB_points_obj") #save the object as RData

