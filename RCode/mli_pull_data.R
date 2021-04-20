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


