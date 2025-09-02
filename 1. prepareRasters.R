# ------------------------------------------
# ------------------------------------------

rm(list = ls())
gc(reset=TRUE)

library(terra)
library(raster)
library(whitebox)
library(FNN)

# ---------

idwCMissing <- function(rasterLayer) {
  
  NACells.rasterLayer <- which( is.na(values(rasterLayer)) )
  NACoords.rasterLayer <- xyFromCell(rasterLayer,NACells.rasterLayer)
  
  r.matrix <- as.matrix(as.data.frame(rasterLayer, xy=T, na.rm=T))
  rawValues <- r.matrix[,3]
  
  closestPoints <- get.knnx(r.matrix[,1:2], NACoords.rasterLayer[,c("x","y")], k=3)
  interpValues <- numeric(nrow(NACoords.rasterLayer))
  
  for( i in 1:nrow(NACoords.rasterLayer) ) {
    
    distances <- closestPoints$nn.dist[i,]
    distances[distances == 0] <- 0.000000001
    index <- closestPoints$nn.index[i,]
    values <- rawValues[index]
    weights <- 1 / distances
    weights <- weights / sum(weights)
    interpValues[i] <- sum(weights * values)
    
  }
  
  rasterLayer[NACells.rasterLayer] <- interpValues
  return(rasterLayer)
  
}
  
idwC <- function(r,studyRegionRaster,studyRegionRaster.matrix) {
  
  r.matrix <- as.matrix(as.data.frame(r, xy=T, na.rm=T))
  r_resampled <- studyRegionRaster
  
  closestPoints <- get.knnx(r.matrix[,1:2], studyRegionRaster.matrix[,c("x","y")], k=3)
  rawValues <- r.matrix[,3]
  interpValues <- numeric(nrow(studyRegionRaster.matrix))
  
  for( i in 1:nrow(studyRegionRaster.matrix) ) {
    
    distances <- closestPoints$nn.dist[i,]
    distances[distances == 0] <- 0.000000001
    index <- closestPoints$nn.index[i,]
    values <- rawValues[index]
    weights <- 1 / distances
    weights <- weights / sum(weights)
    interpValues[i] <- sum(weights * values)
    
  }
  
  r_resampled[studyRegionRaster.matrix[,1]] <- interpValues
  return(r_resampled)
  
}

idwN <- function(r,studyRegionRaster,studyRegionRaster.matrix) {
  
  r.matrix <- as.matrix(as.data.frame(r, xy=T, na.rm=T))
  r_resampled <- studyRegionRaster
  
  closestPoints <- get.knnx(r.matrix[,1:2], studyRegionRaster.matrix[,c("x","y")], k=1)
  rawValues <- r.matrix[,3]
  interpValues <- numeric(nrow(studyRegionRaster.matrix))
  
  for( i in 1:nrow(studyRegionRaster.matrix) ) {
    
    index <- closestPoints$nn.index[i]
    interpValues[i] <- rawValues[index]
    
  }
  
  r_resampled[studyRegionRaster.matrix[,1]] <- interpValues
  return(r_resampled)
  
}

# ------------------------------------------
# Define study region 

# studyRegion <- vect("../Data/Study area/mpa_europe_starea_v2.shp")

# https://www.opendem.info/arc2meters.html
# 0.002 = 7.2 arc-sec = 170m
# 0.015 = 50.8 arc-sec = 1km

# res <- 0.015 # 0.002
# studyRegionRaster <- rast(ext(studyRegion), res = res, crs = "WGS84")
# studyRegionRaster <- rasterize(studyRegion,studyRegionRaster)
# writeRaster(studyRegionRaster,filename="../Data/Study area/studyRegionRaster.tif",overwrite=T)

# ----------------
# Define study region  [Based on wave exposure HR]

# Latitude: 0.001 degrees latitude ≈ 111 meters
# Longitude: 0.001 degrees longitude ≈ 111 meters * cos(latitude)
# Example: If you're at a latitude of 40 degrees:
# Longitude conversion: 0.001 degrees longitude ≈ 111 meters * cos(40) ≈ 85 meters

studyRegion <- vect("../Data/Study area/mpa_europe_starea_v2.shp")
res <- 0.015
studyRegionRaster <- rast(round(ext(studyRegion),1), res = res, crs = "WGS84")
studyRegionRaster <- rasterize(studyRegion,studyRegionRaster)

waveF <- list.files("../Data/Wave Model 1/Wave/", full.names = T)

for( w in waveF ) {
  cat( which(waveF == w) , "out of" , length(waveF), "\n")
  wR <- raster(w)
  wR.res <- terra::resample(rast(wR), studyRegionRaster, threads=16)
  if( w == waveF[1] ) {
    mosaic_raster <- wR.res
  } else {
    mosaic_raster <- mosaic(mosaic_raster, wR.res, fun = "mean")
  }
}

mosaic_raster.lr <- raster("../Data/Wave Model 2/eu5dgcoast1kmraswgs84wx4.tif")
mosaic_raster.lr <- terra::resample(rast(mosaic_raster.lr), mosaic_raster, threads=16)

cellsNAData <- Which( is.na(raster(mosaic_raster)), cells=TRUE)
cellswithLRData <- Which( ! is.na(raster(mosaic_raster.lr)), cells=TRUE)
cellswithLRData <- intersect(cellsNAData, cellswithLRData)
mosaic_raster[cellswithLRData] <- mosaic_raster.lr[cellswithLRData][,1]

correctRegion <- shapefile("../Data/Wave Model 2/correctRegion.shp")
correctRegion.1 <- Which( is.na(raster(mosaic_raster)), cells=TRUE)
correctRegion.2 <- cellFromPolygon(raster(mosaic_raster),correctRegion)
correctRegion <- intersect( correctRegion.1 , unlist(correctRegion.2) )
mosaic_raster[correctRegion] <- 64000

writeRaster(mosaic_raster,filename="../Data/Predictors/waveFetch.tif",overwrite=T)

studyRegionRaster <- mosaic_raster
studyRegionRaster[!is.na(studyRegionRaster)] <- 1
writeRaster(studyRegionRaster,filename="../Data/Study area/studyRegionRaster.tif",overwrite=T)

# ------------------------------------------
# ------------------------------------------

# studyRegionRaster <- rast("../Data/Study area/studyRegionRaster.tif")
# studyRegionRaster.matrix <- cells(!is.na(studyRegionRaster))
# nonNACells <- which(!is.na(values(studyRegionRaster)))
# coordinates <- xyFromCell(studyRegionRaster,nonNACells)
# studyRegionRaster.matrix <- as.matrix(data.frame(nonNACells,coordinates))

# ------------------------------------------
# River runoff

rasterLayer <- rast("../Data/River Runoff/runoff.tif")
rasterLayer <- crop(rasterLayer,studyRegionRaster)
rasterLayer <- idwCMissing(rasterLayer)
rasterLayer <- terra::resample(rasterLayer, studyRegionRaster, threads=8)
rasterLayer <- mask(rasterLayer,studyRegionRaster)
plot(rasterLayer)

rasterLayer_filled  <- mask(rasterLayer, studyRegionRaster)
rasterLayer_filled <- fill_gaps_idw(rasterLayer_filled, studyRegionRaster)
sum(is.na(values(rasterLayer_filled))) == sum(is.na(values(studyRegionRaster)))
plot(rasterLayer_filled)

writeRaster(rasterLayer_filled,filename="../Data/Predictors/runoff.tif",overwrite=T)
gc()

# ------------------------------------------
# Bathymetry data

r <- "../../../../Data/Spatial information/Rasters/GEBCO Bathymetry Global.tif"
r <- "../../../../Data/Spatial information/Rasters/Emodnet Bathymetry Atlantic.tif"

bathy <- rast(r)
bathy <- terra::resample(bathy, studyRegionRaster, threads=8)
slope <- terrain(bathy, v="slope", neighbors=8, unit="degrees")
tri <- terrain(bathy, v="TRI", neighbors=8)
gc()

bathy <- mask(bathy,studyRegionRaster)
slope <- mask(slope,studyRegionRaster)
tri <- mask(tri,studyRegionRaster)
gc()

bathy[bathy > 0] <- 0
sum(is.na(values(bathy))) == sum(is.na(values(studyRegionRaster)))

bathy_filled  <- mask(bathy, studyRegionRaster)
bathy_filled <- fill_gaps_idw(bathy_filled, studyRegionRaster)
sum(is.na(values(bathy_filled))) == sum(is.na(values(studyRegionRaster)))
plot(bathy_filled)
writeRaster(bathy_filled,filename="../Data/Predictors/bathymetry.tif",overwrite=T)

slope_filled  <- mask(slope, studyRegionRaster)
slope_filled <- fill_gaps_idw(slope_filled, studyRegionRaster)
sum(is.na(values(slope_filled))) == sum(is.na(values(studyRegionRaster)))
plot(slope_filled)
writeRaster(slope_filled,filename="../Data/Predictors/slope.tif",overwrite=T)

tri_filled  <- mask(tri, studyRegionRaster)
tri_filled <- fill_gaps_idw(tri_filled, studyRegionRaster)
sum(is.na(values(tri_filled))) == sum(is.na(values(studyRegionRaster)))
plot(tri_filled)
writeRaster(tri_filled,filename="../Data/Predictors/terrainRuggednessIndex.tif",overwrite=T)

# r_resampled <- idwC(bathy,studyRegionRaster,studyRegionRaster.matrix) 
# plot(r_resampled)
# writeRaster(r_resampled,filename="../Data/Predictors/bathymetry.tif",overwrite=T)

# r_resampled <- idwC(slope,studyRegionRaster,studyRegionRaster.matrix) 
# plot(r_resampled)
# writeRaster(r_resampled,filename="../Data/Predictors/slope.tif",overwrite=T)

# r_resampled <- idwC(tri,studyRegionRaster,studyRegionRaster.matrix) 
# plot(r_resampled)
# writeRaster(r_resampled,filename="../Data/Predictors/terrainRuggednessIndex.tif",overwrite=T)

# ------------------------------------------
# Distance to shore

# bathy <- rast("../Data/Predictors/bathymetry.tif")
# bathy[bathy > 0] <- NA
# bathy[!is.na(bathy)] <- 0
# bathy[is.na(bathy)] <- 1

rDist <- studyRegionRaster
rDist[!is.na(rDist)] <- 0
rDist[is.na(rDist)] <- 1

# install_whitebox()
wbt_init()
writeRaster(rDist,filename="../Data/rDist.tif",overwrite=T)
rDist <- wbt_euclidean_distance("../Data/rDist.tif","../Data/rDistProcess.tif")

rDist <- rast("../Data/rDistProcess.tif")
rDist <- terra::resample(rDist, studyRegionRaster, threads=8)
rDist <- mask(rDist,studyRegionRaster)
plot(rDist)

sum(is.na(values(rDist_filled))) == sum(is.na(values(studyRegionRaster)))
writeRaster(rDist_filled,filename="../Data/Predictors/distanceToShore.tif",overwrite=T)

# r_resampled <- idwC(rDist,studyRegionRaster,studyRegionRaster.matrix) 
# plot(r_resampled)
# writeRaster(r_resampled,filename="../Data/Predictors/distanceToShore.tif",overwrite=T)

file.remove("../Data/rDist.tif")
file.remove("../Data/rDistProcess.tif")

# ------------------------------------------
# Environmental data

files <- list.files("../Data/Climate", full.names = T)

for( r in files ) {
  
  rasterLayer <- rast(r)
  rasterLayer <- crop(rasterLayer,studyRegionRaster)
  rasterLayer <- idwCMissing(rasterLayer)
  rasterLayer <- terra::resample(rasterLayer, studyRegionRaster, threads=8)
  rasterLayer <- mask(rasterLayer,studyRegionRaster)

  # r_resampled <- idwC(rast(r),studyRegionRaster,studyRegionRaster.matrix) 
  # plot(r_resampled)
  
  name <- gsub("/Climate/","/Predictors/",r)
  name <- gsub(" Surface Mean","",name)
  name <- gsub(" BenthicMean","",name)
  name <- gsub(" BenthicMin","",name)
  name <- gsub(" 2010-2020","",name)
  name <- gsub(" 2000-2010","",name)
  
  if( ! sum(is.na(values(rasterLayer))) == sum(is.na(values(studyRegionRaster))) ) { stop("!") }

  writeRaster(rasterLayer,filename=name,overwrite=T)
  
}

# ------------------------------------------
# Seafloor geomorphic features

files <- list.files("../Data/Seafloor geomorphic features", full.names = T, pattern = "shp")
files <- files[!grepl("xml",files)]
files <- files[!grepl("Classification",files)]
files <- files[!grepl("Shelf_valleys",files)]

filesN <- list.files("../Data/Seafloor geomorphic features", full.names = F, pattern = "shp")
filesN <- filesN[!grepl("xml",filesN)]
filesN <- filesN[!grepl("Classification",filesN)]
filesN <- filesN[!grepl("Shelf_valleys",filesN)]

rlist <- studyRegionRaster
rlist[] <- NA
gc()

resultsDF <- data.frame()

for( f in 1:length(files)) {
  
  r_resampled <- vect(files[f])
  r_resampled <- rasterize(r_resampled,studyRegionRaster)
  gc()
  
  naIndices.1 <- cells(r_resampled)
  naIndices.2 <- cells(rlist)
  naIndices <- setdiff(naIndices.1 , naIndices.2)
  
  if( length(naIndices) > 0) { rlist[naIndices] <- f }
  gc()
  
  resultsDF <- rbind(resultsDF,data.frame(id=f,name=gsub(".shp","",filesN[f])))
  
}

r_filled  <- mask(rlist, studyRegionRaster)
r_filled <- fill_gaps_idw(r_filled, studyRegionRaster)
sum(is.na(values(r_filled))) == sum(is.na(values(studyRegionRaster)))
plot(r_filled)

writeRaster(r_filled,filename="../Data/Predictors/geomorphicFeatures.tif",overwrite=T)
write.csv(resultsDF,file="../Data/Predictors/geomorphicFeatures.csv", row.names = F, sep=";")
