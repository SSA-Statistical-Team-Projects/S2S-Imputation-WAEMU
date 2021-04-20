## environment set up
remove(list = objects()) ## clear all objects in R workspace

options(
  stringsAsFactors = F, ## tell R to treat text as text, not factors
  width = 80, ## set the maximum width of outputs as 80 characters
  scipen = 6, ## discourage R from displaying numbers in scientific notation
  mc.cores = 6, ## set number of (mac) cores used in parallel processing
  start.time= Sys.time()
)

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
