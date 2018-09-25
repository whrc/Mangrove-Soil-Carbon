## Mangrove Restoration Potential Map (contact: Thomas Worthington <taw52@cam.ac.uk>)
## WP: Update of the SOC layers and difference maps
## R code by: tom.hengl@envirometrix.net

## Software and functions ----
library(rgdal)
library(foreign)
library(raster)
library(utils)
library(snowfall)
library(raster)
library(R.utils)
library(plotKML)
#devtools::install_github("imbs-hl/ranger")
library(ranger)
system("gdalinfo --version")
# GDAL 2.3.0, released 2018/05/04
system("saga_cmd --version")
# SAGA Version: 2.3.1
source('mangroves_soilcarbon_functions.R')

## 100m resolution whole world:
r = raster("/mnt/DATA/MERIT/MERIT_dem_100m_v28_July_2017_i.tif")
ncols = ncol(r)
nrows = nrow(r)
xllcorner = extent(r)[1]
yllcorner = extent(r)[3]
xurcorner = extent(r)[2]
yurcorner = extent(r)[4]
cellsize = res(r)[1]
cellsize2 = 1/480

## Mangroves typlogy
ogrInfo("./Soil_Carbon_Extent/Typology_All_Polygons_Update.shp")
typo.tbl = read.dbf("./Soil_Carbon_Extent/Typology_All_Polygons_Update.dbf")
str(typo.tbl)
summary(typo.tbl$Class)
#Delta Estuary  Fringe  Lagoon 
#207694  171622  261848   98801
typo.tbl$Value = as.integer(typo.tbl$Class)
write.dbf(typo.tbl, "./Soil_Carbon_Extent/Typology_All_Polygons_Update.dbf")
hb.leg = typo.tbl[,c("Value","Class")]
## rasterize to 100 m:
## TAKES >1hr
#system(paste0('saga_cmd -c=64 grid_gridding 0 -INPUT \"./Soil_Carbon_Extent/Typology_All_Polygons_Update.shp\" -FIELD \"Value\" -GRID \"./Soil_Carbon_Extent/Typology_All_Polygons_100m.sgrd\" -GRID_TYPE 2 -TARGET_DEFINITION 0 -TARGET_USER_SIZE ', cellsize, ' -TARGET_USER_XMIN ', xllcorner+cellsize/2,' -TARGET_USER_XMAX ', xurcorner-cellsize/2, ' -TARGET_USER_YMIN ', yllcorner+cellsize/2,' -TARGET_USER_YMAX ', yurcorner-cellsize/2))
## output file 300GB
#system(paste0('gdal_translate ./Soil_Carbon_Extent/Typology_All_Polygons_100m.sdat ./Soil_Carbon_Extent/Typology_All_Polygons_100m.tif -ot \"Byte\" -co \"COMPRESS=DEFLATE\"'))
save.image()

## List of tiles with mangroves ----
#system(paste0('saga_cmd shapes_grid 2 -GRIDS=\"./Soil_Carbon_Extent/Typology_All_Polygons_100m.sgrd\" -POLYGONS=\"/mnt/DATA/LandGIS/models/tiles_ll_100km.shp\" -PARALLELIZED=1 -RESULT=\"./Soil_Carbon_Extent/ov_mangrove_tiles.shp\"'))
ov_mangroves = readOGR("./Soil_Carbon_Extent/ov_mangrove_tiles.shp", "ov_mangrove_tiles")
str(ov_mangroves@data)
summary(selS.t <- !is.na(ov_mangroves$Typology_Al.5))
## remove the 300GB file
#unlink("./Soil_Carbon_Extent/Typology_All_Polygons_100m.sdat")

## 1452 tiles with values
ov_mangroves = ov_mangroves[selS.t,]
saveRDS(ov_mangroves, "ov_mangroves.rds")
#ov_mangroves = readRDS("ov_mangroves.rds")
plot(ov_mangroves)
#writeOGR(ov_mangroves, "./Soil_Carbon_Extent/ov_mangroves.gpkg", "ov_mangroves", "GPKG")
writeOGR(ov_mangroves, "./Soil_Carbon_Extent/ov_mangroves.shp", "ov_mangroves", "ESRI Shapefile")
tS.sel = as.character(ov_mangroves$ID)
newS.dirs <- paste0("/data/mangroves/tiled/T", tS.sel)
x <- lapply(newS.dirs, dir.create, recursive=TRUE, showWarnings=FALSE)

vrt = "/mnt/DATA/GlobalForestChange2000-2014/first.vrt"
landsat.r = raster(vrt)
tr = res(landsat.r)
#tr = c(0.00025, 0.00025)
#tile_shape(1)
try( detach("package:snowfall", unload=TRUE), silent=TRUE)
try( detach("package:snow", unload=TRUE), silent=TRUE)

#x = list.files("/data/mangroves/tiled", pattern=glob2rx("TYPO*.tif"), full.names=TRUE, recursive=TRUE)
#unlink(x)

## Rasterize polygons ----
pol.lst = c("./Soil_Carbon_Extent/Typology_All_Polygons_Update.shp", "./Soil_Carbon_Extent/GMW_1996_2016_Union.shp")
poln.lst = c("TYPO", "GMW")
## test:
i = which(ov_mangroves$ID=="23946")
#i = which(ov_mangroves$ID=="23948")
#tile_shape(i, shape=pol.lst[1], l=tools::file_path_sans_ext(basename(pol.lst)[1]), varname=poln.lst[1], burn=FALSE)
#plot(raster("/data/mangroves/tiled/T23948/TYPO_30m_T23948.tif"))
#tile_shape(i, shape=pol.lst[2], l=tools::file_path_sans_ext(basename(pol.lst)[2]), varname=poln.lst[2])

library(snowfall)
sfInit(parallel=TRUE, cpus=64)
sfExport("tile_shape", "ov_mangroves", "pol.lst", "poln.lst")
sfLibrary(tools)
out <- sfClusterApplyLB(1:nrow(ov_mangroves@data), fun=function(i){ try( tile_shape(i, shape=pol.lst[1], l=tools::file_path_sans_ext(basename(pol.lst)[1]), varname=poln.lst[1], burn=FALSE) )} )
sfStop()

library(snowfall)
sfInit(parallel=TRUE, cpus=64)
sfExport("tile_shape", "ov_mangroves", "pol.lst", "poln.lst")
sfLibrary(tools)
out <- sfClusterApplyLB(1:nrow(ov_mangroves@data), fun=function(i){ try( tile_shape(i, shape=pol.lst[2], l=tools::file_path_sans_ext(basename(pol.lst)[2]), varname=poln.lst[2]) )} )
sfStop()

#x = list.files("/data/mangroves/tiled", pattern=glob2rx("*_250m_*.tif"), full.names=TRUE, recursive=TRUE)

## 250m grid:
library(snowfall)
sfInit(parallel=TRUE, cpus=64)
sfExport("tile_shape", "ov_mangroves", "pol.lst", "poln.lst")
sfLibrary(tools)
out <- sfClusterApplyLB(1:nrow(ov_mangroves@data), fun=function(i){ try( tile_shape(i, shape=pol.lst[1], l=tools::file_path_sans_ext(basename(pol.lst)[1]), tr=0.002083333, varname=poln.lst[1], burn=FALSE, res.name="250m") )} )
sfStop()

## clean up
#x = list.files("/data/mangroves/tiled", pattern=glob2rx("TREL10*.tif"), full.names=TRUE, recursive=TRUE)
#x = list.files("/data/mangroves/tiled", pattern=glob2rx("*SW2L14\")_30m*.tif"), full.names=TRUE, recursive=TRUE)
#x = list.files("/data/mangroves/tiled", pattern=glob2rx("FAPAR*.tif"), full.names=TRUE, recursive=TRUE)
#unlink(x)

## Tiling list of rasters ----
## already at 30 m res:
vrtR.lst = c("/mnt/DATA/PALSAR/PALSAR-2/F02DAR_2017_HV_20m.vrt", 
             "/mnt/DATA/PALSAR/PALSAR-2/F02DAR_2017_HH_20m.vrt", 
             "/mnt/DATA/GlobalForestChange2000-2014/first.vrt",
             "/mnt/DATA/GlobalForestChange2000-2014/last.vrt",
             "/mnt/DATA/GlobalForestChange2000-2014/treecover2000.vrt", 
             "/data/AW3D30/AW3D30_30m.vrt", 
             "/mnt/DATA/GlobalSurfaceWater/occurrence.vrt", 
             "/mnt/DATA/GlobalSurfaceWater/extent.vrt", 
             "/mnt/DATA/Landsat/treecover2010.vrt")
nm.lst = c(list("HH17", "HV17", c("REDL00","NIRL00","SW1L00","SW2L00"), c("REDL14","NIRL14","SW1L14","SW2L14"), "TRCL00", "AW3D30", "OCCGSW", "EXTGSW", "TREL10"))
typR.lst = c(rep("Int16", 2), rep("Byte", 7), rep("Byte", 12))
mvFlag.lst = c(rep(-32767, 3), rep(255, 6), -32767)
#library(parallel)
## test:
#i = which(ov_mangroves$ID=="23946")
#tile_tif(i, vrt=vrtR.lst[k], name=nm.lst[[k]], type=typR.lst[k], mvFlag=mvFlag.lst[k])

library(snowfall)
for(k in 1:length(nm.lst)){
  sfInit(parallel=TRUE, cpus=64)
  sfExport("ov_mangroves", "tile_tif", "vrtR.lst", "mvFlag.lst", "typR.lst", "nm.lst", "k")
  sfLibrary(rgdal)
  out <- sfClusterApplyLB(1:nrow(ov_mangroves@data), fun=function(i){ try( tile_tif(i, vrt=vrtR.lst[k], name=nm.lst[[k]], type=typR.lst[k], mvFlag=mvFlag.lst[k]) , silent = TRUE ) } )
  sfStop()
}

## FAPAR 250m ----
m.lst <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
library(snowfall)
for(j in 1:length(m.lst)){
  sfInit(parallel=TRUE, cpus=64)
  sfExport("downscale_tif", "ov_mangroves", "j", "m.lst")
  sfLibrary(rgdal)
  sfLibrary(RSAGA)
  out <- sfClusterApplyLB(1:nrow(ov_mangroves), function(i){ try( downscale_tif(i, vrt=paste0("/mnt/DATA/LandGIS/layers250m/veg_fapar_proba.v.", tolower(m.lst[j]), "_d_250m_s0..0cm_2014..2017_v1.0.tif"), name=paste0("FAPAR.", tolower(m.lst[j])), tiles=ov_mangroves, cpus=1, fill.gaps.method="resampling") , silent = TRUE )})
  sfStop()
  ## TH: Warning! this function breaks because of some readGDAL problem most likely
  #Error in checkForRemoteErrors(val) : 
  #  3 nodes produced errors; first error: Error in .local(.Object, ...) :
}

## SoilGrids250m ----
#del.lst = list.files("/data/mangroves/tiled", pattern=glob2rx("SOCS_0_30cm_30m_T*.tif"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)
library(snowfall)
sfInit(parallel=TRUE, cpus=64)
sfExport("downscale_tif", "ov_mangroves")
sfLibrary(rgdal)
sfLibrary(RSAGA)
out <- sfClusterApplyLB(1:nrow(ov_mangroves), function(i){ try( downscale_tif(i, vrt="/mnt/DATA/LDN/OCSTHA_M_30cm_300m_ll.tif", name="SOCS_0_30cm", tiles=ov_mangroves, cpus=1, fill.gaps.method="spline") , silent = TRUE )})
sfStop()

## Tidal range ----
#del.lst = list.files("/data/mangroves/tiled", pattern=glob2rx("TidalRange_30m_T*.tif"), full.names=TRUE, recursive=TRUE)
sfInit(parallel=TRUE, cpus=64)
sfExport("downscale_tif", "ov_mangroves")
sfLibrary(rgdal)
sfLibrary(RSAGA)
out <- sfClusterApplyLB(1:nrow(ov_mangroves), function(i){ try( downscale_tif(i, vrt="/data/mangroves/TidalRange/amplitude_Layer_f.tif", name="TidalRange", tiles=ov_mangroves, cpus=1, fill.gaps.method="resampling") , silent = TRUE )})
sfStop()

## clean up tiles without content:
x1 = list.files("/data/mangroves/tiled", pattern=glob2rx("TYPO_30m_T*.tif"), full.names=TRUE, recursive=TRUE)
x0 = list.files("/data/mangroves/tiled", pattern=glob2rx("FAPAR.jan_30m_T*.tif"), full.names=TRUE, recursive=TRUE)
str(which(!dirname(x1) %in% dirname(x0)))
#"/data/mangroves/tiled/T17231" "/data/mangroves/tiled/T19063" "/data/mangroves/tiled/T26040"
## 3 tiles without content
unlink(dirname(x1)[which(!dirname(x1) %in% dirname(x0))], recursive=TRUE)

## Sea surface temperature ----
## https://neo.sci.gsfc.nasa.gov/view.php?datasetId=MYD28M 
## downscale to 30 m resolution:
library(snowfall)
for(j in c(1:4)){
  sfInit(parallel=TRUE, cpus=64)
  sfExport("downscale_tif", "ov_mangroves", "j")
  sfLibrary(rgdal)
  sfLibrary(RSAGA)
  out <- sfClusterApplyLB(1:nrow(ov_mangroves), function(i){ try( downscale_tif(i, vrt=paste0("/data/mangroves/SurfaceTemp/sst_season_", j,"_f.tif"), name=paste0("SST_", j), tiles=ov_mangroves, cpus=1, fill.gaps.method="resampling") )})
  sfStop()
}

## TSM data ----
## downscale to 30 m resolution:
tsm_f.lst = list.files("./TSM_data", pattern=glob2rx("TSM_*_f.tif$"), full.names = TRUE)
library(snowfall)
for(j in 1:length(tsm_f.lst)){
  sfInit(parallel=TRUE, cpus=64)
  sfExport("downscale_tif", "ov_mangroves", "j", "tsm_f.lst")
  sfLibrary(rgdal)
  sfLibrary(RSAGA)
  out <- sfClusterApplyLB(1:nrow(ov_mangroves), function(i){ try( downscale_tif(i, vrt=tsm_f.lst[j], name=paste0("TSM_", j), tiles=ov_mangroves, cpus=1, fill.gaps.method="spline"), silent = TRUE)})
  sfStop()
}

## Soil profiles (soil profile Mangrove DB) ----
## https://docs.google.com/spreadsheets/d/1xVh1cxH1l9cVpKqbqh3Vzx7iqRq3gvne-mMiqm7QVn4/edit#gid=944439000
profs <- read.csv("./soildata/mangrove_soc_database_v10_sites.csv", skip = 1)
str(profs)
summary(is.na(profs$Latitude_Adjusted))
profs.f <- plyr::rename(profs, c("Site.name"="SOURCEID", "Longitude_Adjusted"="LONWGS84", "Latitude_Adjusted"="LATWGS84"))
profs.f$TIMESTRR <- as.Date(profs.f$Years_collected, format="%Y")
profs.f$SOURCEDB = "MangrovesDB"
profs.f$LONWGS84 = as.numeric(paste(profs.f$LONWGS84))
profs.f$LATWGS84 = as.numeric(paste(profs.f$LATWGS84))
str(profs.f)
## 'data.frame':	1812 obs. of  32 variables
length(levels(unique(as.factor(paste0("ID", profs.f$LONWGS84, profs.f$LATWGS84, sep="_")))))
## 1667
## Additional points ----
## from https://doi.org/10.1038/s41558-018-0162-5
profsA <- read.csv("./soildata/41558_2018_162_MOESM2_ESM.csv")
str(profsA)
profsA$SOURCEID = paste("ESM", 1:nrow(profsA), sep="_")
profsA$SOURCEDB = "ESM"
profsA.f <- plyr::rename(profsA[,c("SOURCEID","SOURCEDB","Longitude","Latitude")], c("Longitude"="LONWGS84", "Latitude"="LATWGS84"))
horsA <- data.frame(SOURCEID=profsA$SOURCEID, UHDICM=0, LHDICM=100, DEPTH=50, OCDENS=profsA$SOC..mg.cm.3.)
summary(horsA$OCDENS)
#View(horsA)

## horizon data:
hors <- read.csv("./soildata/mangrove_soc_database_v10_horizons.csv", skip=1)
hors.f <- plyr::rename(hors, c("Site.name"="SOURCEID", "U_depth"="UHDICM", "L_depth"="LHDICM", "OC_final"="ORCDRC", "BD_final"="BLD"))
hors.f$ORCDRC <- hors.f$ORCDRC*10
hors.f$BLD.f <- hors.f$BLD*1000
## convert depths to cm:
hors.f$UHDICM = hors.f$UHDICM*100
hors.f$LHDICM = hors.f$LHDICM*100
summary(hors.f$CD_calc)
hors.f$OCDENS = hors.f$CD_calc * 1000
hors.f$DEPTH <- hors.f$UHDICM + (hors.f$LHDICM - hors.f$UHDICM)/2
hors.fs = hor2xyd(hors.f)
summary(hors.fs$OCDENS)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.00   15.00   26.00   29.41   42.00  134.00       2
## 14780 values now
profs.ALL <- plyr::rbind.fill(profs.f[,c("SOURCEID","SOURCEDB","TIMESTRR","LONWGS84","LATWGS84")], profsA.f[,c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84")])
hors.ALL <- plyr::rbind.fill(hors.fs[,c("SOURCEID","UHDICM","LHDICM","DEPTH","BLD.f","ORCDRC","OCDENS")], horsA[,c("SOURCEID","UHDICM","LHDICM","DEPTH","OCDENS")])

## Merge everything together
SPROPS.MangrovesDB <- plyr::join(profs.ALL, hors.ALL)
summary(SPROPS.MangrovesDB$LATWGS84)
str(SPROPS.MangrovesDB[which(SPROPS.MangrovesDB$LATWGS84 > 60),])
SPROPS.MangrovesDB = SPROPS.MangrovesDB[!is.na(SPROPS.MangrovesDB$LONWGS84) & !is.na(SPROPS.MangrovesDB$LATWGS84) & !is.na(SPROPS.MangrovesDB$DEPTH),]
SPROPS.MangrovesDB = SPROPS.MangrovesDB[!SPROPS.MangrovesDB$DEPTH>1000,]
str(SPROPS.MangrovesDB)
## 'data.frame':	14907 obs. of  11 variables
hist(SPROPS.MangrovesDB$DEPTH)
summary(is.na(SPROPS.MangrovesDB$DEPTH)); summary(is.na(SPROPS.MangrovesDB$OCDENS))

coordinates(SPROPS.MangrovesDB) = ~LONWGS84+LATWGS84
proj4string(SPROPS.MangrovesDB) = CRS("+proj=longlat +datum=WGS84")
#View(SPROPS.MangrovesDB@data)
SPROPS.MangrovesDB$LOC_ID = as.factor(paste0("ID", SPROPS.MangrovesDB@coords[,1], SPROPS.MangrovesDB@coords[,2], sep="_"))
length(levels(unique(SPROPS.MangrovesDB$LOC_ID)))
## 1918 profiles in total
unlink("./soildata/mangroves_SOC_points.gpkg")
writeOGR(SPROPS.MangrovesDB, "./soildata/mangroves_SOC_points.gpkg", "mangroves_SOC_points", "GPKG")
summary(as.factor(SPROPS.MangrovesDB$SOURCEDB))
#ESM MangrovesDB 
#551       14767 
plot(SPROPS.MangrovesDB)

## Prepare RDS tiles for overlay / prediction ----
in.covs = basename(list.files("./tiled/T28356", pattern=glob2rx("*_30m_*.tif$"), full.names = TRUE))
in.covs = in.covs[-grep(pattern="OCDENS", in.covs)]
in.covs = in.covs[-grep(pattern="dSOCS", in.covs)]
in.covs = paste(sapply(in.covs, function(i){strsplit(i, "_T")[[1]][1]}))
## 39 layers

#del.lst = list.files("/data/mangroves/tiled", pattern=glob2rx("T*.rds$"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)
#new_data(i=which(ov_mangroves$ID=="28356"), mask="TYPO_30m", in.covs, tiles=ov_mangroves)
#m = readRDS("/data/mangroves/tiled/T31029/T31029.rds")
#x = readGDAL("/data/mangroves/tiled/T28356/TYPO_30m_T28356.tif")
#plot(m["TYPO_30m"])

## takes ca 1 hr...
library(snowfall)
snowfall::sfInit(parallel=TRUE, cpus=64)
snowfall::sfLibrary(rgdal)
snowfall::sfLibrary(raster)
snowfall::sfExport("new_data", "ov_mangroves", "in.covs")
out <- snowfall::sfClusterApplyLB(1:nrow(ov_mangroves), function(i){ new_data(i, mask="TYPO_30m", in.covs, tiles=ov_mangroves) })
snowfall::sfStop()
## Error in checkForRemoteErrors(val) : 
## 2 nodes produced errors; first error: NAs not permitted in row index

## Overlay points and 30 m res covs ----
tile.pol = rgdal::readOGR("/mnt/DATA/LandGIS/models/tiles_ll_100km.shp", "tiles_ll_100km")
source("/mnt/DATA/LandGIS/R/LandGIS_functions.R")
ovM <- extract.tiled(x=SPROPS.MangrovesDB, tile.pol=tile.pol, path="/data/mangroves/tiled", ID="ID", cpus=64)
rmatrix = ovM[,c("SOURCEID","DEPTH","OCDENS",in.covs)]
#names(rmatrix) = gsub("L00_30m", "_30m", names(rmatrix))
str(rmatrix)
## 'data.frame':	15318 obs. of  41 variables
summary(rmatrix$OCDENS)
summary(rmatrix$AW3D30_30m)
summary(rmatrix$SOCS_0_30cm_30m)
summary(rmatrix$SST_3_30m)
## 4165 missing values ?!
head(rmatrix[which(is.na(rmatrix$AW3D30_30m))[1],])
head(rmatrix[which(is.na(rmatrix$AW3D30_30m))[1452],])
## TH: missing values are due to the mangrove TYPO mask map
rmatrix$TYPO_30m = as.factor(rmatrix$TYPO_30m)
summary(rmatrix$TYPO_30m)
i.rm = model.matrix(~TYPO_30m-1, rmatrix)
str(i.rm)
rmatrix.f = cbind(rmatrix[!is.na(rmatrix$TYPO_30m),], data.frame(i.rm))

## Fit a ranger model ----
library(ranger)
fm.OCDENS <- as.formula(paste0("OCDENS ~ DEPTH + SOCS_0_30cm_30m + AW3D30_30m + HH17_30m + HV17_30m + TRCL00_30m + TREL10_30m + SW1L00_30m + SW2L00_30m + REDL00_30m + NIRL00_30m + SW1L14_30m + SW2L14_30m + REDL14_30m + NIRL14_30m + ", paste0("SST_",c(1:4),"_30m", collapse="+"), "+", paste0("TSM_",c(1:4),"_30m", collapse="+"), "+", paste0("FAPAR.", tolower(m.lst), "_30m", collapse = "+"), "+", paste0("TYPO_30m", 1:4, collapse = "+")))
fm.OCDENS
which(!all.vars(fm.OCDENS) %in% names(rmatrix.f))
## Note: banding artifacts in TSM_3 to TSM_5
rmatrix.f = rmatrix.f[complete.cases(rmatrix.f[,all.vars(fm.OCDENS)]),]
m.OCDENS_30m <- ranger(fm.OCDENS, rmatrix.f, num.trees = 150, importance='impurity', mtry=25)
m.OCDENS_30m
# Number of trees:                  150 
# Sample size:                      11153 
# Number of independent variables:  39 
# Mtry:                             11 
# Target node size:                 5 
# Variable importance mode:         impurity 
# Splitrule:                        variance 
# OOB prediction error (MSE):       56.00599 
# R squared (OOB):                  0.822202
# RMSE = +/- 7.5 kg/m3
## TH: THIS R-SQUARE IS NOT REALISTIC!
xl <- as.list(ranger::importance(m.OCDENS_30m))
#xl <- as.list(m.OCDENSq_30m$variable.importance)
print(t(data.frame(xl[order(unlist(xl), decreasing=TRUE)[1:10]])))
# SOCS_0_30cm_30m 1130729.04
# DEPTH            184335.41
# TSM_1_30m        176223.36
# TSM_3_30m        165466.16
# TSM_4_30m        165434.19
# TSM_2_30m        132054.48
# REDL00_30m       108308.54
# TRCL00_30m        82648.16
# REDL14_30m        81517.38
# SST_3_30m         64882.69

## Fit ensamble model ----
## TH: the right way to do it!
## https://cran.r-project.org/web/packages/SuperLearner/vignettes/Guide-to-SuperLearner.html#multicore-parallelization
## Use location of point so that the training is realistic:
## folds per tile - very important otherwise over-optimistic predictions!
proj4string(tile.pol) = proj4string(SPROPS.MangrovesDB)
id.t = over(SPROPS.MangrovesDB, tile.pol)
str(id.t)
id.d = plyr::join(rmatrix.f["SOURCEID"], data.frame(SOURCEID=SPROPS.MangrovesDB$SOURCEID, ID=id.t$ID), match="first")
str(id.d)
## Each 100 x 100 km block is used separately in training / validation.

detach("package:snowfall", unload=TRUE)
library(parallel)
library(SuperLearner)
cl <- parallel::makeCluster(64)
x <- parallel::clusterEvalQ(cl, library(SuperLearner))
## fitting with 10-fold cross-validation included
sl <- snowSuperLearner(Y = rmatrix.f$OCDENS, X = rmatrix.f[,all.vars(fm.OCDENS)[-1]], cluster = cl, SL.library = c("SL.xgboost","SL.glmnet", "SL.ranger"), id=id.d$ID, cvControl=list(V=10))
sl
# Risk      Coef
# SL.xgboost_All 175.6226 0.3885847
# SL.glmnet_All  182.9464 0.2702721
# SL.ranger_All  175.8531 0.3411432

## Cross-validation ----
cv_sl <- CV.SuperLearner(Y = rmatrix.f$OCDENS, X = rmatrix.f[,all.vars(fm.OCDENS)[-1]], parallel = cl, SL.library = c("SL.xgboost","SL.glmnet","SL.ranger"), V=10, id=id.d$ID, verbose=TRUE)
#SL.library = c("SL.xgboost","SL.ksvm", "SL.glmnet", "SL.ranger")
summary(cv_sl)
# Risk is based on: Mean Squared Error
# 
# All risk estimates are based on V =  10 
# 
# Algorithm    Ave     se    Min    Max
# Super Learner 151.60 3.6084 82.903 232.18
# Discrete SL 169.59 3.8744 91.006 220.39
# SL.xgboost_All 172.20 3.5805 97.799 269.82
# SL.ksvm_All 184.28 4.0730 86.710 291.63
# SL.glmnet_All 170.14 3.9622 95.836 227.53
# SL.ranger_All 165.25 3.8343 91.006 237.54
## RMSE = 12.3 kg/m3
stopCluster(cl)

## clean up:
#del.lst = list.files("/data/mangroves/tiled", pattern=glob2rx("OCDENS_*.tif$"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)
#del.lst = list.files("/data/mangroves/tiled", pattern=glob2rx("SOCS_0_30cm_30m_*.tif$"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)
# del.lst = list.files("/data/mangroves/tiled", pattern=glob2rx("dSOCS_0_*00cm_*.tif$"), full.names=TRUE, recursive=TRUE)
# unlink(del.lst)

## Predictions ----
#predict_x(i=which(ov_mangroves$ID=="31419"), tiles=ov_mangroves, gm=sl)
#predict_x(i=which(ov_mangroves$ID=="31029"), tiles=ov_mangroves, gm=sl)
#predict_x(i=which(ov_mangroves$ID=="31420"), tiles=ov_mangroves, gm=sl)
#predict_x(i=which(ov_mangroves$ID=="24778"), tiles=ov_mangroves, gm=sl)
#predict_x(i=which(ov_mangroves$ID=="26085"), tiles=ov_mangroves, gm=sl)

library(snowfall)
snowfall::sfInit(parallel=TRUE, cpus=24)
snowfall::sfLibrary(rgdal)
snowfall::sfLibrary(ranger)
snowfall::sfLibrary(SuperLearner)
snowfall::sfExport("predict_x", "ov_mangroves", "sl")
out <- snowfall::sfClusterApplyLB(1:nrow(ov_mangroves), fun=function(i){ try( predict_x(i, tiles=ov_mangroves, gm=sl) ) })
snowfall::sfStop()
#Error in checkForRemoteErrors(val) : 
#  17 nodes produced errors; first error: Error in `$<-.data.frame`(`*tmp*`, "DEPTH", value = 0) : 
#  replacement has 1 row, data has 0

## export mosaic:
t.lst <- list.files("/data/mangroves/tiled", pattern=glob2rx("dSOCS_0_100cm_*.tif$"), full.names=TRUE, recursive=TRUE)
#x = file.copy(t.lst, paste0("/data/tmp/mangroves_SOC30m/", basename(t.lst)), overwrite = TRUE)
x = file.copy(t.lst, paste0("/data/mangroves_SOC30m/", basename(t.lst)), overwrite = TRUE)
out.tmp <- "/data/mangroves_SOC30m/t_list.txt"
vrt.tmp <- "/data/mangroves_SOC30m/mangroves_dSOC_0_100cm_30m.vrt"
cat(paste0("/data/mangroves_SOC30m/", basename(t.lst)), sep="\n", file=out.tmp)
system(paste0('gdalbuildvrt -input_file_list ', out.tmp, ' ', vrt.tmp))
