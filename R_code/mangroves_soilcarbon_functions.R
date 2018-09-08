## Prediction of soil organic carbon stocks for world Mangroves at 30 m
## tom.hengl@gmail.com

## copy horizon values for 3D regression
hor2xyd = function(x, U="UHDICM", L="LHDICM", treshold.T=15){
  x$DEPTH <- x[,U] + (x[,L] - x[,U])/2
  x$THICK <- x[,L] - x[,U]
  sel = x$THICK < treshold.T
  ## begin and end of the horizon:
  x1 = x[!sel,]; x1$DEPTH = x1[,L]
  x2 = x[!sel,]; x2$DEPTH = x1[,U]
  y = do.call(rbind, list(x, x1, x2))
  return(y)
}

## expand tiles for 1-2pix ----
fill.NA.cells <- function(i, tiles, var="TYPO_30m", out.path="/data/mangroves/tiled/"){
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
tile_shape = function(i, shape, l, out.dir="/data/mangroves/tiled/", tr=0.00025, varname, res.name="30m", burn=TRUE, value.col="Value"){
  out.tif = paste0(out.dir, 'T', ov_mangroves@data[i,"ID"], '/', varname, '_', res.name, '_T', ov_mangroves@data[i,"ID"], '.tif')
  if(!file.exists(out.tif)){
    te = paste(unlist(ov_mangroves@data[i,c("xl","yl","xu","yu")]), collapse = " ")
    if(burn==TRUE){
      system(paste0('gdal_rasterize ', shape, ' -l ',  l, ' -te ', te, ' -tr ', tr, ' ', tr, ' -ot Byte ', ' -burn 100 ', out.tif, ' -a_nodata 0 -co \"COMPRESS=DEFLATE\"'), ignore.stdout=TRUE)
    } else {
      system(paste0('gdal_rasterize ', shape, ' -l ',  l, ' -te ', te, ' -tr ', tr, ' ', tr, ' -ot Byte ', ' -a ', value.col, ' ', out.tif, ' -a_nodata 255 -co \"COMPRESS=DEFLATE\"'), ignore.stdout=TRUE)
    }
  }
}

## resample 30m images using tiles ----
tile_tif <- function(i, vrt="/mnt/cartman/GlobalForestChange2000-2014/first.vrt", name=c("REDL00","NIRL00","SW1L00","SW2L00"), out.dir="/data/mangroves/tiled/", tr=0.00025, type="Byte", mvFlag=255, fix.mask=TRUE){
  out.tif = paste0(out.dir, "T", ov_mangroves@data[i,"ID"], "/", name, "_30m_T", ov_mangroves@data[i,"ID"], ".tif")
  if(any(!file.exists(out.tif))){
    mask.tif = paste0(out.dir, "T", ov_mangroves@data[i,"ID"], "/TYPO_30m_T", ov_mangroves@data[i,"ID"], ".tif")
    if(file.exists(mask.tif)){
      tmp.file = paste0(tempfile(), ".tif")
      te = paste(unlist(ov_mangroves@data[i,c("xl","yl","xu","yu")]), collapse = " ")
      system(paste0('gdalwarp ', vrt, ' ', tmp.file, ' -r \"near\" -te ', te, ' -tr ', tr, ' ', tr), ignore.stdout=TRUE)
      r = readGDAL(tmp.file)
      if(fix.mask==TRUE){
        r$mask = readGDAL(mask.tif, silent=TRUE)$band1
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
}


## dowscale function ----
downscale_tif <- function(i, vrt, name, tiles, out.dir="/data/mangroves/tiled/", tr=0.00025, type="Int16", mvFlag=-32768, cpus=64, fill.gaps.coarse=FALSE, fill.gaps.method="resampling"){
  out.tif = paste0(out.dir, "T", tiles@data[i,"ID"], "/", name, "_30m_T", tiles@data[i,"ID"], ".tif")
  if(any(!file.exists(out.tif))){
    if(fill.gaps.coarse==TRUE){
      maskTile = paste0('/data/mangroves/tiled/T', ov_mangroves@data[i,"ID"], '/TYPO_250m_T', ov_mangroves@data[i,"ID"], '.sgrd')
    } else { 
      maskTile = paste0('/data/mangroves/tiled/T', ov_mangroves@data[i,"ID"], '/TYPO_30m_T', ov_mangroves@data[i,"ID"], '.sgrd')
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
        s$mask = readGDAL(paste0(out.dir, "T", ov_mangroves@data[i,"ID"], "/TYPO_250m_T", ov_mangroves@data[i,"ID"], ".tif"))$band1
      } else {
        s$mask = readGDAL(paste0(out.dir, "T", ov_mangroves@data[i,"ID"], "/TYPO_30m_T", ov_mangroves@data[i,"ID"], ".tif"))$band1
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
new_data <- function(i, mask, in.covs, tiles, out.dir="/data/mangroves/tiled"){
  out.rds <- paste0(out.dir, "/T", tiles@data[i,"ID"], "/T", tiles@data[i,"ID"], ".rds")
  if(!file.exists(out.rds)){
    m = readGDAL(paste0(out.dir, "/T", tiles@data[i,"ID"], "/", mask, "_T", tiles@data[i,"ID"], ".tif"), silent = TRUE)
    m = as(m, "SpatialPixelsDataFrame")
    l.covs = in.covs[-which(in.covs %in% mask)]
    for(j in 1:length(l.covs)){
      in.tif = paste0(out.dir, "/T", tiles@data[i,"ID"], "/", l.covs[j], "_T", tiles@data[i,"ID"], ".tif")
      if(!file.exists(in.tif)){
        m@data[,j+1] = 0
      } else {
        m@data[,j+1] = readGDAL(in.tif, silent = TRUE)$band1[m@grid.index]
      }
    }
    names(m) = c(mask, l.covs)
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
predict_x <- function(i, tiles, gm, year="00", depths=c(0,30,100,200), out.dir="/data/mangroves/tiled", with.se=FALSE, factor.c="TYPO_30m", type="sl", num.threads = 3){
  out.tif = paste0(out.dir, "/T", tiles@data[i,"ID"], "/OCDENS_", depths, "_year20", year, "_30m_T", tiles@data[i,"ID"], ".tif")
  if(any(!file.exists(out.tif)) | !file.exists(paste0(out.dir, "/T", tiles@data[i,"ID"], "/dSOCS_0_100cm_year20", year, "_30m_T", tiles@data[i,"ID"], ".tif"))){
    require(rgdal)
    require(SuperLearner)
    m <- readRDS(paste0(out.dir, "/T", tiles@data[i,"ID"], "/T", tiles@data[i,"ID"], ".rds"))
    if(!is.null(factor.c)){
      m = m[!is.na(m@data[,factor.c]),]
      m@data[,factor.c] <- factor(m@data[,factor.c], levels=1:4)
      m.i = model.matrix(as.formula(paste0("~", factor.c, "-1")), m@data)
      new.data = cbind(m@data, data.frame(m.i))
    }
    for(j in 1:length(depths)){
      new.data$DEPTH = depths[j]
      if(any(class(gm)=="list")){
        xpl = list(NULL)
        for(i in 1:length(gm)){
          xpl[[i]] = predict(gm[[i]], new.data)$predictions
        }
        m@data[,paste0("OCDENS_",j)] = rowMeans(data.frame(xpl))
      } else {
        if(any(class(gm)=="randomForest")){
          m@data[,paste0("OCDENS_",j)] = predict(gm, new.data)
        } else {
          if(with.se==TRUE){
            ## Too computational at the moment
            #x = predict(gm, m@data, quantiles = c(0.05, 0.5, 0.95), all=FALSE)
            ## TH: does not work for large matrices
            seq.x = seq.int(by=1e4, from=1, to= nrow(m@data))
            x = lapply(2:length(seq.x), function(k){predict(gm, m@data[seq.x[k-1] : seq.x[k],], quantiles = c(0.05, 0.5, 0.95), all=FALSE)})
            x = do.call(rbind, x)
            m@data[,paste0("OCDENS_",j)] = x[,2]
            m@data[,paste0("OCDENS_",j,".L")] = ifelse(x[,1]<0, 0, x[,1])
            m@data[,paste0("OCDENS_",j,".U")] = x[,3]
          } else {
            if(type=="sl"){
              m@data[,paste0("OCDENS_",j)]  <- predict(gm, new.data[,gm$varNames], onlySL = TRUE, num.threads = num.threads)$pred[,1]
            } else {
              m@data[,paste0("OCDENS_",j)] <- predict(gm, new.data)$predictions 
            }
          }
        }
      }
      writeGDAL(m[paste0("OCDENS_",j)], out.tif[j], type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE")
      if(with.se==TRUE){ 
        writeGDAL(m[paste0("OCDENS_",j,".L")], gsub("OCDENS_", "OCDENS.L_", out.tif[j]), type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE")
        writeGDAL(m[paste0("OCDENS_",j,".U")], gsub("OCDENS_", "OCDENS.U_", out.tif[j]), type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE")
      }
    }
    depthT = c(30,70,100)
    x = list(NULL)
    sx = list(NULL)
    sy = list(NULL)
    for(k in 1:length(depthT)){
      x[[k]] = rowMeans(m@data[,c(paste0("OCDENS_",k),paste0("OCDENS_",k+1))], na.rm=TRUE)*depthT[k]/100 
      if(with.se==TRUE){
        ## propagated error (http://lectureonline.cl.msu.edu/~mmp/labs/error/e2.htm):
        #sx[[k]] = 1/2*sqrt(m@data[,paste0("OCDENS_",k)]^2+m@data[,paste0("OCDENS_",k+1)]^2)
        sx[[k]] = rowMeans(m@data[,c(paste0("OCDENS_",k,".L"),paste0("OCDENS_",k+1,".L"))], na.rm=TRUE)*depthT[k]/100 
        sy[[k]] = rowMeans(m@data[,c(paste0("OCDENS_",k,".U"),paste0("OCDENS_",k+1,".U"))], na.rm=TRUE)*depthT[k]/100
      }
    }
    m$SOCS = rowSums(as.data.frame(x[1:2]), na.rm=TRUE)*10
    writeGDAL(m["SOCS"], paste0(out.dir, "/T", tiles@data[i,"ID"], "/dSOCS_0_100cm_year20", year, "_30m_T", tiles@data[i,"ID"], ".tif"), type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE")
    if(with.se==TRUE){
      #m$SOCS.sd = sqrt(sx[[1]]^2+sx[[2]]^2)*10
      m$SOCS_L = rowSums(as.data.frame(sx[1:2]), na.rm=TRUE)*10
      writeGDAL(m["SOCS_L"], paste0(out.dir, "/T", tiles@data[i,"ID"], "/dSOCS_sL_0_100cm_year20", year, "_30m_T", tiles@data[i,"ID"], ".tif"), type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE")
      m$SOCS_U = rowSums(as.data.frame(sy[1:2]), na.rm=TRUE)*10
      writeGDAL(m["SOCS_U"], paste0(out.dir, "/T", tiles@data[i,"ID"], "/dSOCS_sU_0_100cm_year20", year, "_30m_T", tiles@data[i,"ID"], ".tif"), type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE")
    }
    m$SOCS2 = rowSums(as.data.frame(x[1:3]), na.rm=TRUE)*10
    writeGDAL(m["SOCS2"], paste0(out.dir, "/T", tiles@data[i,"ID"], "/dSOCS_0_200cm_year20", year, "_30m_T", tiles@data[i,"ID"], ".tif"), type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE")
    #if(with.se==TRUE){
      #m$SOCS2.sd = sqrt(sx[[1]]^2+sx[[2]]^2+sx[[3]]^2)*10
      #writeGDAL(m["SOCS2.sd"], paste0(out.dir, "/T", tiles@data[i,"ID"], "/dSOCS_sd_0_200cm_year20", year, "_30m_T", tiles@data[i,"ID"], ".tif"), type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE")
    #}
  }
}

## predict soil properties in parallel:
predict_parallelP <- function(j, sel, varn, formulaString, rmatrix, idcol, cpus, Nsub=1e4, remove_duplicates=FALSE, pars.ranger){
  s.train <- rmatrix[!sel==j,]
  if(remove_duplicates==TRUE){
    ## TH: optional - check how does model performs without the knowledge of the 3D dimension
    sel.dup = !duplicated(s.train[,idcol])
    s.train <- s.train[sel.dup,]
  }
  s.test <- rmatrix[sel==j,]
  n.l <- dim(s.test)[1]
  if(missing(Nsub)){ Nsub = length(all.vars(formulaString))*50 }
  if(Nsub>nrow(s.train)){ Nsub = nrow(s.train) }
  if(missing(pars.ranger)){
    gm <- quantregRanger(formulaString, s.train[complete.cases(s.train),])
  } else {
    pars.ranger$mtry = ifelse(pars.ranger$mtry >= length(all.vars(formulaString)), length(all.vars(formulaString))-1, pars.ranger$mtry)
    ##  mtry can not be larger than number of variables in data
    gm <- quantregRanger(formulaString, s.train[complete.cases(s.train),], pars.ranger)
  }
  sel.t = complete.cases(s.test)
  x.pred <- predict(gm, s.test[sel.t,], quantiles = c((1-.682)/2, 0.5, 1-(1-.682)/2))
  pred <- data.frame(predictions=x.pred[,2], se=(x.pred[,3]-x.pred[,1])/2) 
  ## export object
  obs.pred <- as.data.frame(list(s.test[sel.t,varn], pred$predictions, pred$se), col.names=c("Observed", "Predicted", "SE"))
  obs.pred[,idcol] <- s.test[sel.t,idcol]
  obs.pred$fold = j
  return(obs.pred)
}

cv_numeric <- function(formulaString, rmatrix, nfold, idcol, cpus=1, Log=FALSE, LLO=TRUE, pars.ranger){     
  varn = all.vars(formulaString)[1]
  message(paste0("Running ", nfold, "-fold cross validation with model re-fitting"))
  if(nfold > nrow(rmatrix)){ 
    stop("'nfold' argument must not exceed total number of points") 
  }
  if(sum(duplicated(rmatrix[,idcol]))>0.5*nrow(rmatrix)){
    if(LLO==TRUE){
      ## TH: Leave whole locations out
      ul <- paste(unique(rmatrix[,idcol]))
      sel.ul <- dismo::kfold(ul, k=nfold)
      sel <- lapply(1:nfold, function(o){ data.frame(row.names=which(rmatrix[,idcol] %in% ul[sel.ul==o]), x=rep(o, length(which(rmatrix[,idcol] %in% ul[sel.ul==o])))) })
      sel <- do.call(rbind, sel)
      sel <- sel[order(as.numeric(row.names(sel))),]
      message(paste0("Subsetting observations by unique location"))
    } else {
      sel <- dismo::kfold(rmatrix, k=nfold, by=rmatrix[,idcol])
      message(paste0("Subsetting observations by '", idcol, "'"))
    }
  } else {
    sel <- dismo::kfold(rmatrix, k=nfold)
    message(paste0("Simple subsetting of observations using kfolds"))
  }
  if(missing(cpus)){  
    cpus = nfold 
  }
  cpus.a <- parallel::detectCores(all.tests = FALSE, logical = FALSE) 
  cpus <- ifelse(cpus.a < cpus, cpus.a, cpus)
  if(cpus>1){
    snowfall::sfExport("predict_parallelP","idcol","formulaString","rmatrix","sel","varn","pars.ranger")
    snowfall::sfLibrary(package="plyr", character.only=TRUE)
    snowfall::sfLibrary(package="ranger", character.only=TRUE)
    out <- snowfall::sfLapply(1:nfold, function(j){predict_parallelP(j, sel=sel, varn=varn, formulaString=formulaString, rmatrix=rmatrix, idcol=idcol, pars.ranger=pars.ranger)})
    snowfall::sfStop()
  } else {
    out <- lapply(1:nfold, function(j){predict_parallelP(j, sel=sel, varn=varn, formulaString=formulaString, rmatrix=rmatrix, idcol=idcol, pars.ranger=pars.ranger)})
  }
  ## calculate mean accuracy:
  out <- plyr::rbind.fill(out)
  ME = mean(out$Observed - out$Predicted, na.rm=TRUE)
  MAE = mean(abs(out$Observed - out$Predicted), na.rm=TRUE)
  RMSE = sqrt(mean((out$Observed - out$Predicted)^2, na.rm=TRUE))
  ## Errors of errors:
  MAE.SE = mean(out$SE - abs(out$Observed - out$Predicted), na.rm=TRUE)
  ## https://en.wikipedia.org/wiki/Coefficient_of_determination
  R.squared = 1-var(out$Observed - out$Predicted, na.rm=TRUE)/var(out$Observed, na.rm=TRUE)
  if(Log==TRUE){
    logRMSE = sqrt(mean((log1p(out$Observed) - log1p(out$Predicted))^2, na.rm=TRUE))
    logR.squared = 1-var(log1p(out$Observed) - log1p(out$Predicted), na.rm=TRUE)/var(log1p(out$Observed), na.rm=TRUE)
    cv.r <- list(out, data.frame(ME=ME, MAE=MAE, RMSE=RMSE, MAE.SE=MAE.SE, R.squared=R.squared, logRMSE=logRMSE, logR.squared=logR.squared)) 
  } else {
    cv.r <- list(out, data.frame(ME=ME, MAE=MAE, RMSE=RMSE, MAE.SE=MAE.SE, R.squared=R.squared))
  }
  names(cv.r) <- c("CV_residuals", "Summary")
  return(cv.r)
}

## correlation plot:
pfun <- function(x,y, ...){
  panel.hexbinplot(x,y, ...)  
  panel.abline(0,1,lty=1,lw=2,col="black")
  panel.abline(0+RMSE,1,lty=2,lw=2,col="black")
  panel.abline(0-RMSE,1,lty=2,lw=2,col="black")
}
