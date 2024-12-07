# ------------------------------------------
# ------------------------------------------

rm(list = ls())
gc(reset=TRUE)

library(terra)
library(raster)
library(whitebox)
library(FNN)

# ---------

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

studyRegion <- vect("../Data/Study area/mpa_europe_starea_v2.shp")

# Latitude: 0.001 degrees latitude ≈ 111 meters
# Longitude: 0.001 degrees longitude ≈ 111 meters * cos(latitude)
# Example: If you're at a latitude of 40 degrees:
# Longitude conversion: 0.001 degrees longitude ≈ 111 meters * cos(40) ≈ 85 meters

res <- 0.01
studyRegionRaster <- rast(ext(studyRegion), res = res, crs = "WGS84")
studyRegionRaster <- rasterize(studyRegion,studyRegionRaster)
writeRaster(studyRegionRaster,filename="../Data/Study area/studyRegionRaster.tif",overwrite=T)

# ----------------
# Define study region  [Based on wave exposure HR]

# Latitude: 0.001 degrees latitude ≈ 111 meters
# Longitude: 0.001 degrees longitude ≈ 111 meters * cos(latitude)
# Example: If you're at a latitude of 40 degrees:
# Longitude conversion: 0.001 degrees longitude ≈ 111 meters * cos(40) ≈ 85 meters

studyRegion <- vect("../Data/Study area/mpa_europe_starea_v2.shp")
res <- 0.005
studyRegionRaster <- rast(ext(studyRegion), res = res, crs = "WGS84")
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

studyRegionRaster.matrix <- cells(!is.na(studyRegionRaster))
nonNACells <- which(!is.na(values(studyRegionRaster)))
coordinates <- xyFromCell(studyRegionRaster,nonNACells)
studyRegionRaster.matrix <- as.matrix(data.frame(nonNACells,coordinates))

# ------------------------------------------
# River runoff

lowResData <- raster("../Data/River Runoff/runoff.tif")
lowResData <- terra::resample(rast(lowResData), studyRegionRaster, threads=16)

r_resampled <- idwC(lowResData,studyRegionRaster,studyRegionRaster.matrix) 
plot(r_resampled)
writeRaster(r_resampled,filename="../Data/Predictors/runoff.tif",overwrite=T)

# ------------------------------------------
# Bathymetry data

r <- "../../../Data/Spatial information/Rasters/GEBCO Bathymetry Global.tif"

bathy <- rast(r)
bathy <- terra::resample(bathy, studyRegionRaster, threads=8)
slope <- terrain(bathy, v="slope", neighbors=8, unit="degrees")
tri <- terrain(bathy, v="TRI", neighbors=8)

bathy <- mask(bathy,studyRegionRaster)
slope <- mask(slope,studyRegionRaster)
tri <- mask(tri,studyRegionRaster)

r_resampled <- idwC(bathy,studyRegionRaster,studyRegionRaster.matrix) 
plot(r_resampled)
writeRaster(r_resampled,filename="../Data/Predictors/bathymetry.tif",overwrite=T)

r_resampled <- idwC(slope,studyRegionRaster,studyRegionRaster.matrix) 
plot(r_resampled)
writeRaster(r_resampled,filename="../Data/Predictors/slope.tif",overwrite=T)

r_resampled <- idwC(tri,studyRegionRaster,studyRegionRaster.matrix) 
plot(r_resampled)
writeRaster(r_resampled,filename="../Data/Predictors/terrainRuggednessIndex.tif",overwrite=T)

# ------------------------------------------
# Distance to shore

r <- "../../../Data/Spatial information/Rasters/GEBCO Bathymetry Global.tif"

bathy <- rast(r)
bathy <- crop(bathy,studyRegionRaster)
bathy[bathy > 0] <- NA
bathy[!is.na(bathy)] <- 0
bathy[is.na(bathy)] <- 1

# install_whitebox()
wbt_init()
writeRaster(bathy,filename="../Data/rDist.tif",overwrite=T)
rDist <- wbt_euclidean_distance("../Data/rDist.tif","../Data/rDistProcess.tif")

rDist <- rast("../Data/rDistProcess.tif")
rDist <- terra::resample(rDist, studyRegionRaster, threads=8)
rDist <- mask(rDist,studyRegionRaster)

r_resampled <- idwC(rDist,studyRegionRaster,studyRegionRaster.matrix) 
plot(r_resampled)
writeRaster(r_resampled,filename="../Data/Predictors/distanceToShore.tif",overwrite=T)

file.remove("../Data/rDist.tif")
file.remove("../Data/rDistProcess.tif")

# ------------------------------------------
# Environmental data

files <- list.files("../Data/Climate", full.names = T)

for( r in files ) {
  
  r_resampled <- idwC(rast(r),studyRegionRaster,studyRegionRaster.matrix) 
  plot(r_resampled)
  
  name <- gsub("/Climate/","/Predictors/",r)
  name <- gsub(" Surface Mean","",name)
  name <- gsub(" BenthicMean","",name)
  name <- gsub(" BenthicMin","",name)
  name <- gsub(" 2010-2020","",name)
  name <- gsub(" 2000-2010","",name)
  
  writeRaster(r_resampled,filename=name,overwrite=T)
  
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

resultsDF <- data.frame()

for( f in 1:length(files)) {
  
  r_resampled <- vect(files[f])
  r_resampled <- rasterize(r_resampled,studyRegionRaster)
 # r_resampled <- idwC(r_resampled,studyRegionRaster,studyRegionRaster.matrix) 
  
  naIndices.1 <- cells(r_resampled)
  naIndices.2 <- cells(rlist)
  naIndices <- setdiff(naIndices.1 , naIndices.2)
  
  if( length(naIndices) > 0) { rlist[naIndices] <- f }
  
  resultsDF <- rbind(resultsDF,data.frame(id=f,name=gsub(".shp","",filesN[f])))
  
}

rlist <- mask(rlist,studyRegionRaster)
rlist <- idwN(rlist,studyRegionRaster,studyRegionRaster.matrix)

all.equal(which(is.na(values(rlist))),which(is.na(values(studyRegionRaster))))

writeRaster(rlist,filename="../Data/Predictors/geomorphicFeatures.tif",overwrite=T)
write.csv(resultsDF,file="../Data/Predictors/geomorphicFeatures.csv", row.names = F, sep=";")
