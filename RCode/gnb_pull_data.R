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

