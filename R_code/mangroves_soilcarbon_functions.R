
## expand tiles for 1-2pix ----
fill.NA.cells <- function(i, tiles, var="MNGUSG_30m", out.path="/data/mangroves/tiled/"){
  out.tif = paste0(out.path, "T", tiles@data[i,"ID"], "/", var, "f_T", tiles@data[i,"ID"], ".tif")
  if(!file.exists(out.tif)){
    m = readGDAL(paste0(out.path, "/T", tiles@data[i,"ID"], "/", var, "_T", tiles@data[i,"ID"], ".tif"), silent = TRUE)
    r = raster(m["band1"])
    rf = focal(r, w=matrix(1,3,3), fun=mean, pad=TRUE, na.rm=TRUE, NAonly=TRUE)
    m$out = as(rf, "SpatialGridDataFrame")@data[,1]
    writeGDAL(m["out"], out.tif, type="Byte", mvFlag=0, options=c("COMPRESS=DEFLATE"))
  }
}

## tile shapes ----
tile_shape = function(i, shape="./DownloadPack-WCMC-010-MangroveUSGS2011-Ver1-3/WCMC-010-MangroveUSGS2011-ver1-3.shp", l="WCMC-010-MangroveUSGS2011-ver1-3", out.dir="/data/mangroves/tiled/", tr=0.00025, varname="MNGUSG", res.name="30m"){
  out.tif = paste0(out.dir, 'T', ov_mangroves@data[i,"ID"], '/', varname, '_', res.name, '_T', ov_mangroves@data[i,"ID"], '.tif')
  if(!file.exists(out.tif)){
    te = paste(unlist(ov_mangroves@data[i,c("xl","yl","xu","yu")]), collapse = " ")
    system(paste0('gdal_rasterize ', shape, ' -l ',  l, ' -te ', te, ' -tr ', tr, ' ', tr, ' -ot Byte ', ' -burn 100 ', out.tif, ' -a_nodata 0 -co \"COMPRESS=DEFLATE\"'))
  }
}

## resample 30m images using tiles ----
tile_tif <- function(i, vrt="/mnt/cartman/GlobalForestChange2000-2014/first.vrt", name=c("REDL00","NIRL00","SW1L00","SW2L00"), out.dir="/data/mangroves/tiled/", tr=0.00025, type="Byte", mvFlag=255, fix.mask=TRUE){
  out.tif = paste0(out.dir, "T", ov_mangroves@data[i,"ID"], "/", name, "_30m_T", ov_mangroves@data[i,"ID"], ".tif")
  if(any(!file.exists(out.tif))){
    tmp.file = paste0(tempfile(), ".tif")
    te = paste(unlist(ov_mangroves@data[i,c("xl","yl","xu","yu")]), collapse = " ")
    system(paste0('gdalwarp ', vrt, ' ', tmp.file, ' -r \"near\" -te ', te, ' -tr ', tr, ' ', tr))
    r = readGDAL(tmp.file)
    if(fix.mask==TRUE){
      r$mask = readGDAL(paste0(out.dir, "T", ov_mangroves@data[i,"ID"], "/MNGUSG_30mf_T", ov_mangroves@data[i,"ID"], ".tif"))$band1
      for(k in 1:length(name)){
        r$fix = ifelse(is.na(r$mask), NA, r@data[,k])
        writeGDAL(r["fix"], drivername="GTiff", type=type, mvFlag=mvFlag, out.tif[k], options=c("COMPRESS=DEFLATE"))
      }
    } else {
      for(k in 1:length(name)){
        writeGDAL(r[k], drivername="GTiff", type=type, mvFlag=mvFlag, out.tif[k], options=c("COMPRESS=DEFLATE"))
      }
    }
    unlink(tmp.file)
    gc(); gc()
  }
}


## dowscale function ----
downscale_tif <- function(i, vrt, name, tiles, out.dir="/data/mangroves/tiled/", tr=0.00025, type="Int16", mvFlag=-32768, cpus=48, fill.gaps.coarse=FALSE, fill.gaps.method="resampling"){
  out.tif = paste0(out.dir, "T", tiles@data[i,"ID"], "/", name, "_30m_T", tiles@data[i,"ID"], ".tif")
  if(any(!file.exists(out.tif))){
    if(fill.gaps.coarse==TRUE){
      maskTile = paste0('/data/mangroves/tiled/T', ov_mangroves@data[i,"ID"], '/MNGUSG_250m_T', ov_mangroves@data[i,"ID"], '.sgrd')
    } else { 
      maskTile = paste0('/data/mangroves/tiled/T', ov_mangroves@data[i,"ID"], '/MNGUSG_30mf_T', ov_mangroves@data[i,"ID"], '.sgrd')
    }
    tmp.file = paste0(tempfile(), ".sdat")
    te = paste(unlist(tiles@data[i,c("xl","yl","xu","yu")]), collapse = " ")
    a_srs = proj4string(tiles)
    if(fill.gaps.coarse==TRUE){
      system(paste0('gdalwarp ', vrt, ' ', tmp.file, ' -of \"SAGA\" -ot \"Int16\"  -te ', te, ' -tr 0.002083333 0.002083333'))
    } else {
      system(paste0('gdalwarp ', vrt, ' ', tmp.file, ' -of \"SAGA\" -r \"cubicspline\" -ot \"Int16\" -dstnodata ', ifelse(mvFlag=="-32768", "-32767", mvFlag), ' -te ', te, ' -tr ', tr, ' ', tr))
    }
    tmp.file2 = paste0(tempfile(), ".sdat")
    ## http://saga-gis.org/saga_module_doc/2.2.0/grid_tools_7.html 
    if(fill.gaps.method=="spline"){ 
      system(paste0('saga_cmd -c=', cpus, ' grid_tools 7 -INPUT ', RSAGA::set.file.extension(tmp.file, ".sgrd"), ' -THRESHOLD 1 -MASK ', maskTile, ' -RESULT ', RSAGA::set.file.extension(tmp.file2, ".sgrd")))
    }
    if(fill.gaps.method=="resampling"){
      system(paste0('saga_cmd -c=', cpus, ' grid_tools 29 -RESAMPLING 2 -INPUT ', RSAGA::set.file.extension(tmp.file, ".sgrd"), ' -MASK ', maskTile, ' -RESULT ', RSAGA::set.file.extension(tmp.file2, ".sgrd")))
    }
    if(file.size(tmp.file2)>0){
      s = readGDAL(tmp.file2, silent=TRUE)
      if(fill.gaps.coarse==TRUE){ 
        s$mask = readGDAL(paste0(out.dir, "T", ov_mangroves@data[i,"ID"], "/MNGUSG_250m_T", ov_mangroves@data[i,"ID"], ".tif"))$band1
      } else {
        s$mask = readGDAL(paste0(out.dir, "T", ov_mangroves@data[i,"ID"], "/MNGUSG_30mf_T", ov_mangroves@data[i,"ID"], ".tif"))$band1
      }
      s$fix = ifelse(is.na(s$mask), NA, s$band1)
      if(is.numeric(s$fix)){
        if(fill.gaps.coarse==TRUE){
          tmp.file3 = paste0(tempfile(), ".tif")
          writeGDAL(s["fix"], drivername="GTiff", type=type, mvFlag=mvFlag, fname=tmp.file3, options=c("COMPRESS=DEFLATE"))
          system(paste0('gdalwarp ', tmp.file3, ' ', out.tif, ' -co \"COMPRESS=DEFLATE\" -r \"cubicspline\" -te ', te, ' -tr ', tr, ' ', tr))
          unlink(tmp.file3)
        } else {
          writeGDAL(s["fix"], drivername="GTiff", type=type, mvFlag=mvFlag, out.tif, options=c("COMPRESS=DEFLATE"))
        }
      }
    }
    unlink(tmp.file)
    unlink(tmp.file2)
    gc(); gc()
  }
}

## prepare new data ----
new_data <- function(i, in.covs, tiles, out.dir="/data/mangroves/tiled"){
  out.rds <- paste0(out.dir, "/T", tiles@data[i,"ID"], "/T", tiles@data[i,"ID"], ".rds")
  if(!file.exists(out.rds)){
    m = readGDAL(paste0(out.dir, "/T", tiles@data[i,"ID"], "/", in.covs[1], "_T", tiles@data[i,"ID"], ".tif"), silent = TRUE)
    m = as(m, "SpatialPixelsDataFrame")
    for(j in 2:length(in.covs)){
      in.tif = paste0(out.dir, "/T", tiles@data[i,"ID"], "/", in.covs[j], "_T", tiles@data[i,"ID"], ".tif")
      if(!file.exists(in.tif)){
        m@data[,j] = 0
      } else {
        m@data[,j] = readGDAL(in.tif, silent = TRUE)$band1[m@grid.index]
      }
    }
    names(m) = in.covs
    ## Filter out water bodies:
    m = m[m$NIRL00_30m>10,]
    ## Fill-in remaining missing values (can be very tricky):
    sel.mis = sapply(m@data, function(x){sum(is.na(x))>0})
    if(sum(sel.mis)>0){
      x = which(sel.mis)
      for(k in 1:length(x)){
        if(!is.factor(m@data[,x[k]])){
          if(length(grep(pattern="MTYP", names(m)[x[k]]))>0 | length(grep(pattern="GSW", names(m)[x[k]]))>0 | length(grep(pattern="mfw", names(m)[x[k]]))>0 | length(grep(pattern="TR", names(m)[x[k]]))>0 | sum(is.na(m@data[,x[k]]))==nrow(m) ){ 
            repn = rep(0, nrow(m)) 
          } else {
            ## To be on the safe side, replace missing values two times:
            #repn = quantile(m@data[,x[k]], probs=.5, na.rm=TRUE)
            r = raster::raster(m[x[k]])
            ## 1 using proximity filter:
            rf = raster::focal(r, w=matrix(1,15,15), fun=mean, pad=TRUE, na.rm=TRUE, NAonly=TRUE)
            repn = as(rf, "SpatialGridDataFrame")@data[m@grid.index,1]
            ## 2 using dominant value:
            repn = ifelse(is.na(repn), quantile(repn, probs=.5, na.rm=TRUE), repn)
          }
          m@data[,x[k]] = ifelse(is.na(m@data[,x[k]]), repn, m@data[,x[k]])
        }
      }
    }
    ## limit to complete data sets only?
    #m = m[complete.cases(m@data),]
    saveRDS(m, out.rds)
    gc(); gc()
  }  
}

## predict at all locations at 30 m ----
predict_x <- function(i, tiles, gm, year="00", depths=c(0,30,100,200), out.dir="/data/mangroves/tiled", with.se=FALSE){
  out.tif = paste0(out.dir, "/T", tiles@data[i,"ID"], "/OCDENS_", depths, "_year20", year, "_30m_T", tiles@data[i,"ID"], ".tif")
  if(any(!file.exists(out.tif)) | !file.exists(paste0(out.dir, "/T", tiles@data[i,"ID"], "/dSOCS_0_100cm_year20", year, "_30m_T", tiles@data[i,"ID"], ".tif"))){
    m <- readRDS(paste0(out.dir, "/T", tiles@data[i,"ID"], "/T", tiles@data[i,"ID"], ".rds"))
    names(m) = gsub(paste0("L",year,"_30m"), "_30m", names(m))
    if(year=="00"){ names(m) <- gsub("a2000mfw_30m", "mfw_30m", names(m)) }
    if(year=="14"){ names(m) <- gsub("a2012mfw_30m", "mfw_30m", names(m)) }
    for(j in 1:length(depths)){
      m$DEPTH = depths[j]
      if(class(gm)=="list"){
        xpl = list(NULL)
        for(i in 1:length(gm)){
          xpl[[i]] = predict(gm[[i]], m@data)$predictions
        }
        m@data[,paste0("OCDENS_",j)] = rowMeans(data.frame(xpl))
      } else {
        if(class(gm)=="randomForest"){
          m@data[,paste0("OCDENS_",j)] = predict(gm, m@data)
        } else {
          if(with.se==TRUE){
            ## Too computational at the moment
            x = predict(gm, m@data[1:1e4,], type="se")
            m@data[,paste0("OCDENS_",j)] = x$predictions 
            m@data[,paste0("OCDENS_",j,".sd")] = x$se
          } else {
            m@data[,paste0("OCDENS_",j)] = predict(gm, m@data)$predictions 
          }
        }
      }
      writeGDAL(m[paste0("OCDENS_",j)], out.tif[j], type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE")
      if(with.se==TRUE){ 
        writeGDAL(m[paste0("OCDENS_",j,".sd")], gsub("OCDENS_", "OCDENS.sd_", out.tif[j]), type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE")
      }
    }
    depthT = c(30,70,100)
    x = list(NULL)
    sx = list(NULL)
    for(k in 1:length(depthT)){
      x[[k]] = rowMeans(m@data[,c(paste0("OCDENS_",k),paste0("OCDENS_",k+1))], na.rm=TRUE)*depthT[k]/100 
      if(with.se==TRUE){
        ## propagated error (http://lectureonline.cl.msu.edu/~mmp/labs/error/e2.htm):
        sx[[k]] = 1/2*sqrt(m@data[,paste0("OCDENS_",k)]^2+m@data[,paste0("OCDENS_",k+1)]^2)
      }
    }
    m$SOCS = rowSums(as.data.frame(x[1:2]), na.rm=TRUE)*10
    writeGDAL(m["SOCS"], paste0(out.dir, "/T", tiles@data[i,"ID"], "/dSOCS_0_100cm_year20", year, "_30m_T", tiles@data[i,"ID"], ".tif"), type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE")
    if(with.se==TRUE){
      m$SOCS.sd = sqrt(sx[[1]]^2+sx[[2]]^2)*10
      writeGDAL(m["SOCS.sd"], paste0(out.dir, "/T", tiles@data[i,"ID"], "/dSOCS_sd_0_100cm_year20", year, "_30m_T", tiles@data[i,"ID"], ".tif"), type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE")
    }
    m$SOCS2 = rowSums(as.data.frame(x[1:3]), na.rm=TRUE)*10
    writeGDAL(m["SOCS2"], paste0(out.dir, "/T", tiles@data[i,"ID"], "/dSOCS_0_200cm_year20", year, "_30m_T", tiles@data[i,"ID"], ".tif"), type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE")
    if(with.se==TRUE){
      m$SOCS2.sd = sqrt(sx[[1]]^2+sx[[2]]^2+sx[[3]]^2)*10
      writeGDAL(m["SOCS2.sd"], paste0(out.dir, "/T", tiles@data[i,"ID"], "/dSOCS_sd_0_200cm_year20", year, "_30m_T", tiles@data[i,"ID"], ".tif"), type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE")
    }
  }
}
