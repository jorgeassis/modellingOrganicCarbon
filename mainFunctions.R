# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
#
# One Pipeline for Modelling the distribution of marine Species
#
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------

packages.to.use <- c(  
  "FNN"
  ,"dggridR"
  ,"exactextractr"
  ,"maptools"
  ,"ggnewscale"
  #,"h3"
  ,"sf"
  ,"patchwork"
  ,"usdm"
  ,"fasterize"
  ,"sm"
  ,"sf"
  ,"data.table"
  , "blockCV"
  
  , "ecospat"
  , "spThin"
  , "ecodist"
  , "credentials"
  , "rnaturalearth"
  , "ggplot2"
  , "mboost"
  , "blockCV"
  , "raster"
  , "gdata"
  , "dismo"
  , "gbm"
  , "sp"
  , "parallel"
  , "doParallel"
  , "biganalytics"
  , "nicheROVER"
  , "rgeos"
  
  , "devtools"
  , "ENMeval"
  , "remotes"
  ,"leaflet"
  , "monmlp") #, "bcp" , "rgdal"

packages.to.use <- unique(packages.to.use)

for(package in packages.to.use) {
  print(package)
  if( ! package %in% rownames(installed.packages()) & package == "patchwork" ) { devtools::install_github("thomasp85/patchwork") }
  if( ! package %in% rownames(installed.packages()) & package == "SDMTools" ) { install_version("SDMTools", "1.1-221") }
  if( ! package %in% rownames(installed.packages()) ) { install.packages( package ) }
  if( ! package %in% rownames(installed.packages()) ) { install.packages( package , type = "source") }
  if( ! package %in% rownames(installed.packages()) ) { stop("Error on package instalation") }
  library(package, character.only = TRUE)
}

## -----------------------

predictRasterTiled <- function(raster_input, model, tile_size_cells, output_dir, filename_prefix = "prediction_tile", n.trees = NULL, maxValue, nCores=10) {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  raster_extent <- ext(raster_input)
  min_x <- raster_extent[1]
  max_x <- raster_extent[2]
  min_y <- raster_extent[3]
  max_y <- raster_extent[4]
  
  # Calculate the number of tiles in x and y directions
  n_col_tiles <- ceiling((max_x - min_x) / tile_size_cells)
  n_row_tiles <- ceiling((max_y - min_y) / tile_size_cells)
  
  # List to store paths of the saved tile predictions
  tile_paths <- character(0)
  
  # Loop through rows and columns of tiles

  comb <- expand.grid(i=1:n_row_tiles,j=1:n_col_tiles)
  
  cl <- parallel::makeCluster(nCores)
  registerDoParallel(cl)
  
  doParallel <- foreach(c = 1:nrow(comb), .combine=rbind, .packages = c("gbm","dismo","terra")) %dopar% {
      
      i <- comb[c,1]
      j <- comb[c,2]

      # Calculate the extent for the current tile
      tile_min_x <- min_x + (j - 1) * tile_size_cells
      tile_max_x <- min(min_x + j * tile_size_cells, max_x)
      tile_min_y <- max_y - i * tile_size_cells
      tile_max_y <- max(max_y - (i - 1) * tile_size_cells, min_y)
      
      tile_extent <- ext(tile_min_x, tile_max_x, tile_min_y, tile_max_y)
      
      # Crop the input raster to the current tile
      tile <- terra::crop(rast(raster_input), tile_extent)
      
      # Check if the tile has any data (avoid predicting on empty tiles)
      if ( sum(!is.na(values(tile[[1]])[,1])) != 0) {
        # Predict using the provided model
        
        prediction_tile <- predict(tile, model, n.trees = ifelse(is.null(n.trees), model$n.trees, n.trees), type = "response")
        prediction_tile <- mask(prediction_tile,tile[[1]])
        prediction_tile[prediction_tile >= log(maxValue) ] <- log(maxValue)
        
        tile_filename <- file.path(output_dir, paste0(filename_prefix, "_", i, "_", j, "Log.tif"))
        writeRaster(prediction_tile, tile_filename, overwrite = TRUE)
      
        prediction_tile <- exp(prediction_tile)
        prediction_tile[prediction_tile >= maxValue ] <- maxValue
        
        tile_filename <- file.path(output_dir, paste0(filename_prefix, "_", i, "_", j, "Exp.tif"))
        writeRaster(prediction_tile, tile_filename, overwrite = TRUE)
        
    }
  }
  
  stopCluster(cl); rm(cl) ; gc(reset=TRUE)
  closeAllConnections()
  
  return(NULL)
}

## -----------------------

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

## -----------------------

themeMap <- 
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica", color = "#22211d"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_line(color = "black", size = 0),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "#FFFFFF", color = NA), # F3F3F3
    panel.background = element_rect(fill = "#FFFFFF", color = NA), # F3F3F3
    legend.background = element_rect(fill = "#FFFFFF", color = NA),
    legend.box.background = element_rect(fill='#FFFFFF'),
    panel.border = element_blank()
  )

themePlot <- theme(
  text = element_text(size=12) ,
  panel.background = element_rect(fill = "#F0F0F0", colour = "#F0F0F0",size = 0, linetype = "solid"),
  panel.grid.major = element_line(size = 0.1, linetype = 'solid',colour = "#EFEFEF"), 
  axis.ticks.y = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.y = element_text(margin = margin(t = 0, r = 12, b = 0, l = 0)) ,
  axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)) ,
  axis.text.x = element_text(size=11, margin = margin(t = 8, r = 0, b = 0, l = 0)),
  axis.text.y = element_text(size=11, margin = margin(t = 0, r = 8, b = 0, l = 0)))

## -----------------------

modelPlot <- function(model,data.to.plot,var) {
  
  options(warn=-1)
  
  logit2prob <- function(logit){
    odds <- exp(logit)
    prob <- odds / (1 + odds)
    return(prob)
  }
  
  ## -----------------------
  
  # data.to.plot <- dataset[,varsToPlot]
  
  data.to.plot <- data.to.plot[complete.cases(data.to.plot),]
  data.to.plot <- sapply(data.to.plot,as.numeric)
  dataset.variable <- data.to.plot[,var]
  
  data.to.plot <- apply(data.to.plot,2,mean,na.rm=T)
  data.to.plot <- do.call("rbind", replicate(100, data.to.plot, simplify = FALSE))
  data.to.plot[,var] <- seq(min(dataset.variable),max(dataset.variable),length.out=100)
  data.to.plot <- as.data.frame(data.to.plot)
  
  matrixEffect.m <- data.to.plot
  matrixEffect.m <- matrixEffect.m[model$feature_names]
  matrixEffect.m <- xgb.DMatrix(data = data.matrix(matrixEffect.m), label = rep(0,nrow(matrixEffect.m)))
  matrixEffect.m <- data.frame(Effect= predict(model,matrixEffect.m, type = "response") ) 
  
  matrixEffect.m <- logit2prob(matrixEffect.m)
  
  data.to.plot <- data.frame(Variable=data.to.plot[,var],value=matrixEffect.m)
  data.to.plot <- data.to.plot[sort(data.to.plot$Variable, index.return =T)$ix,]
  
  if(data.to.plot[nrow(data.to.plot),2] > data.to.plot[1,2]) {
    for( x in 2:nrow(data.to.plot) ) { if( data.to.plot[x,2] < data.to.plot[x-1,2] ) { data.to.plot[x,2] <- data.to.plot[x-1,2] } }
  }
  if(data.to.plot[nrow(data.to.plot),2] < data.to.plot[1,2]) {
    for( x in 2:nrow(data.to.plot) ) { if( data.to.plot[x,2] > data.to.plot[x-1,2] ) { data.to.plot[x,2] <- data.to.plot[x-1,2] } }
  }
  
  ## -----------------------
  
  return(data.to.plot)
  
}

## -----------------------

idw <- function(distances, values, power = 2, threshold = 1e-10) {
  if (length(distances) != length(values)) {
    stop("Lengths of distances and values must be equal.")
  }
  
  if (length(distances) == 0) {
    stop("Input vectors are empty.")
  }
  
  weights <- 1 / (distances^power + threshold)
  weighted_values <- values * weights
  interpolated_value <- sum(weighted_values) / sum(weights)
  
  return(interpolated_value)
}

## -----------------------

loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

## -----------------------

plotMap <- function(coordinates,radius,color) {
  
  set.seed(42)
  
  m <- leaflet()
  m <- addTiles(m)
  # m <- addMarkers(m, lng=coordinates[,1], lat=coordinates[,2], popup=paste0( "Species record ") , icon = greenLeafIcon)
  
  m <- addCircleMarkers(m, lng=coordinates[,1], lat=coordinates[,2], 
                        popup=paste0( "Species record ") , 
                        radius = radius, color = color , 
                        stroke = FALSE, fillOpacity = 0.5 )
  
  return(m)
  
}

## -----------------------

processLayers <- function(rasterLayers,occurrenceRecords,regionBuffer,minDepth,maxDepth,intertidal) {
  
  if( ! is.null(minDepth) ) { if( minDepth == "NULL" ) { minDepth <- NULL } }
  if( ! is.null(maxDepth) ) { if( maxDepth == "NULL" ) { maxDepth <- NULL } }
  
  rasters <- rasterLayers
  
  if( ! is.null(regionBuffer) ) {
    
    if (length(regionBuffer) == 1) { regionBuffer <- rep(regionBuffer,4)}
    
    final.extent.rasters <- c( min(occurrenceRecords[,1] , na.rm=T) - regionBuffer[1],
                               max(occurrenceRecords[,1] , na.rm=T) + regionBuffer[2],
                               min(occurrenceRecords[,2] , na.rm=T) - regionBuffer[3],
                               max(occurrenceRecords[,2] , na.rm=T) + regionBuffer[4] )
    
    final.extent.rasters <- c(  ifelse(final.extent.rasters[1] < -180 , -180 , final.extent.rasters[1]),
                                ifelse(final.extent.rasters[2] > 180 , 180 , final.extent.rasters[2]),
                                ifelse(final.extent.rasters[3] < -90 , -90 , final.extent.rasters[3]),
                                ifelse(final.extent.rasters[4] > 90 , 90 , final.extent.rasters[4]) )
    
    rasters <- crop(rasters, extent(final.extent.rasters))
    
  } else { final.extent.rasters <- c(-180,180,-90,90) }
  
  ## -----------------------
  
  if( !is.null(intertidal) ) {
    
    intertidalLayer <- raster(intertidal)
    intertidalLayer <- crop(intertidalLayer, rasters)
    rastersIntertidal <- mask(rasters,intertidalLayer)
    
  }
  
  ## -----------------------
  
  regions <- clump(subset(rasters,1))
  regionsUsed <- unique(raster::extract(regions,occurrenceRecords))
  regionsUsed <- regionsUsed[ !is.na(regionsUsed) ]
  
  regions[ ! regions %in% regionsUsed ] <- NA
  rasters <- mask(rasters,regions)
  
  ## -----------------------
  
  if( ! is.null(minDepth) | ! is.null(maxDepth) ) {
    
    if( is.null(minDepth) ) { minDepth <- 0 }
    if( is.null(maxDepth) ) { maxDepth <- 99998 }
    
    minDepth <- minDepth * (-1)
    maxDepth <- maxDepth * (-1)
    
    rclmat <- data.frame(from=c(-99999,maxDepth,minDepth),to=c(maxDepth,minDepth,0),value=c(NA,1,NA))
    rclmat <- rclmat[rclmat$from != rclmat$to,]
    
    bathymetry <- raster(bathymetryDataLayer)
    bathymetry <- crop(bathymetry, rasters)
    
    bathymetry <- reclassify(bathymetry, as.matrix(rclmat))
    rastersBathymetry <- mask(rasters,bathymetry)
    
  }
  
  if( exists("rastersIntertidal") & exists("rastersBathymetry") ) {
    
    rasters.i <- stack( sapply(1:length(names(rasters)),function(x) {  calc(stack(subset(rastersBathymetry,x),subset(rastersIntertidal,x)) , mean,na.rm=TRUE ) }) )
    names(rasters.i) <- names(rastersBathymetry)
    rm(rasters)
    rasters <- rasters.i
    
  }
  
  if( exists("rastersIntertidal") & ! exists("rastersBathymetry") ) { rasters <- mask(rasters,rastersIntertidal) }
  
  if( ! exists("rastersIntertidal") & exists("rastersBathymetry") ) { rasters <- rastersBathymetry }
  
  return(rasters)
  
}

## -----------------------

dropNoVariationLayers <- function(rasterLayers) {
  
  cave <- function(x) { length(unique(x)) / length(x) }
  randomLocations <- Which(!is.na(subset(rasterLayers,1)),cells=TRUE)
  randomLocations <- xyFromCell( subset(rasterLayers,1) , sample(randomLocations,min(length(randomLocations),1000),replace=FALSE))
  varRasterLayers <- which( apply( raster::extract(rasterLayers,randomLocations) ,2,cave) > 0.025)
  rasterLayers <- subset(rasterLayers,varRasterLayers)
  return(rasterLayers)
  
}

## -----------------------

relocateNACoords.id <- function(occurrenceRecords,rasterLayers,relocateType,relocateSpeciesDistance,relocateSpeciesDepth) {
  
  set.seed(42)
  
  if(relocateType == "nonDistance") {
    
    vals <- raster::extract(subset(rasterLayers,1),occurrenceRecords)
    occurrenceRecords <- occurrenceRecords[which(!is.na(vals)),]
    
  }
  
  if(relocateType == "distance") {
    
    if( relocateSpeciesDepth ) { 
      
      if( is.null(minDepth) ) { minDepth <- 0 }
      if( is.null(maxDepth) ) { maxDepth <- 99998 }
      
      minDepth.i <- minDepth * (-1)
      maxDepth.i <- maxDepth * (-1)
      
      rclmat <- data.frame(from=c(-99999,maxDepth.i,minDepth.i),to=c(maxDepth.i,minDepth.i,0),value=c(NA,1,NA))
      rclmat <- rclmat[rclmat$from != rclmat$to,]
      
      shape <- subset(rasterLayers,1)
      
      bathymetry <- raster(bathymetryDataLayer)
      bathymetry <- crop(bathymetry, shape )
      bathymetry <- mask(bathymetry,shape )
      
      shape <- reclassify(bathymetry, as.matrix(rclmat))
      
    }
    
    if( ! relocateSpeciesDepth ) { shape <- raster::subset(rasterLayers,1) } 
    
    to.relocate <- unique(which(is.na(raster::extract(shape,occurrenceRecords[,c("Lon","Lat")]))))
    coordinates.to.relocate <- occurrenceRecords[to.relocate,]
    
    if( nrow(coordinates.to.relocate) > 0 ) { 
      
      old.presences <- occurrenceRecords[ (1:nrow(occurrenceRecords))[! 1:nrow(occurrenceRecords) %in% to.relocate] ,]
      
      correct.points <- xyFromCell(shape,  Which(!is.na(shape), cells=TRUE) ) # result
      
      cat( paste0("Relocating ",length(to.relocate)," Points that were falling out of range"))
      cat( paste0("\n"))
      
      near.cells <- numeric(nrow(coordinates.to.relocate))
      
      for(p in 1:nrow(coordinates.to.relocate)) {
        
        near.cell.p <- spDistsN1( as.matrix(correct.points), as.matrix(coordinates.to.relocate[p,c("Lon","Lat")]),longlat=TRUE)
        
        if( near.cell.p[which.min(near.cell.p)] <= sqrt(sum(relocateSpeciesDistance^2,relocateSpeciesDistance^2)) ) {
          
          near.cell.p <- which.min(near.cell.p)
          
        } else {   near.cell.p <- NA }
        
        near.cells[p] <- near.cell.p
        
      }
      
      relocated <- which(!is.na(near.cells))
      
      if( length(relocated) > 0) {
        
        near.cells <- data.frame(resp=coordinates.to.relocate[relocated,"resp"],Lon=correct.points[near.cells[relocated],1],Lat=correct.points[near.cells[relocated],2])
        colnames(near.cells) <- colnames(old.presences)
        occurrenceRecords <- rbind(old.presences,near.cells)
        
      }
      
    }
    
    ## -----------------------
    
    if( nrow(coordinates.to.relocate) == 0) { 
      
      cat( paste0("None to Relocate"))
      cat( paste0("\n"))
      
    }
    
    ## -----------------------
    
  }
  
  toRemove <- raster::extract(shape,occurrenceRecords[,c("Lon","Lat")])
  toRemove <- which(is.na(toRemove))
  
  if( length(toRemove) > 0 ) { occurrenceRecords <- occurrenceRecords[-toRemove,] }
  
  return( occurrenceRecords )
  
}

## -----------------------

relocateNACoords <- function(occurrenceRecords,rasterLayers,relocateType,relocateSpeciesDistance,relocateSpeciesDepth) {
  
  set.seed(42)
  
  if(relocateType == "nonDistance") {
    
    vals <- raster::extract(subset(rasterLayers,1),occurrenceRecords)
    occurrenceRecords <- occurrenceRecords[which(!is.na(vals)),]
    
  }
  
  if(relocateType == "distance") {
    
    if( relocateSpeciesDepth ) { 
      
      if( is.null(minDepth) ) { minDepth <- 0 }
      if( is.null(maxDepth) ) { maxDepth <- 99998 }
      
      minDepth.i <- minDepth * (-1)
      maxDepth.i <- maxDepth * (-1)
      
      rclmat <- data.frame(from=c(-99999,maxDepth.i,minDepth.i),to=c(maxDepth.i,minDepth.i,0),value=c(NA,1,NA))
      rclmat <- rclmat[rclmat$from != rclmat$to,]
      
      shape <- subset(rasterLayers,1)
      
      bathymetry <- raster(bathymetryDataLayer)
      bathymetry <- crop(bathymetry, shape )
      bathymetry <- mask(bathymetry,shape )
      
      shape <- reclassify(bathymetry, as.matrix(rclmat))
      
    }
    
    if( ! relocateSpeciesDepth ) { shape <- raster::subset(rasterLayers,1) } 
    
    to.relocate <- unique(which(is.na(raster::extract(shape,occurrenceRecords))))
    coordinates.to.relocate <- occurrenceRecords[to.relocate,]
    
    if( nrow(coordinates.to.relocate) > 0 ) { 
      
      old.presences <- occurrenceRecords[ (1:nrow(occurrenceRecords))[! 1:nrow(occurrenceRecords) %in% to.relocate] ,]
      
      correct.points <- xyFromCell(shape,  Which(!is.na(shape), cells=TRUE) ) # result
      
      cat( paste0("Relocating ",length(to.relocate)," Points that were falling out of range"))
      cat( paste0("\n"))
      
      near.cells <- numeric(nrow(coordinates.to.relocate))
      
      for(p in 1:nrow(coordinates.to.relocate)) {
        
        near.cell.p <- spDistsN1( as.matrix(correct.points), as.matrix(coordinates.to.relocate[p,]),longlat=TRUE)
        
        if( near.cell.p[which.min(near.cell.p)] <= sqrt(sum(relocateSpeciesDistance^2,relocateSpeciesDistance^2)) ) {
          
          near.cell.p <- which.min(near.cell.p)
          
        } else {   near.cell.p <- NA }
        
        near.cells[p] <- near.cell.p
        
      }
      
      relocated <- which(!is.na(near.cells))
      
      if( length(relocated) > 0) {
        
        near.cells <- data.frame(Lon=correct.points[near.cells[relocated],1],Lat=correct.points[near.cells[relocated],2])
        colnames(near.cells) <- colnames(old.presences)
        occurrenceRecords <- rbind(old.presences,near.cells)
        
      }
      
    }
    
    ## -----------------------
    
    if( nrow(coordinates.to.relocate) == 0) { 
      
      cat( paste0("None to Relocate"))
      cat( paste0("\n"))
      
    }
    
    ## -----------------------
    
  }
  
  toRemove <- raster::extract(shape,occurrenceRecords)
  toRemove <- which(is.na(toRemove))
  
  if( length(toRemove) > 0 ) { occurrenceRecords <- occurrenceRecords[-toRemove,] }
  
  return( occurrenceRecords )
  
}

## -----------------------

spatialAutocorrelation <- function(rasterLayers,autocorrelationClassDistance,autocorrelationSubsetRecords,autocorrelationMaxDistance,autocorrelationSignif) {
  
  records <- xyFromCell( subset(rasterLayers,1), Which(!is.na(subset(rasterLayers,1)), cell=TRUE) )
  records <- records[sample(1:nrow(records),autocorrelationSubsetRecords,replace=F),]
  colnames(records) <- dataRecordsNames
  
  # --------
  
  if( length(names(rasterLayers)) > 1  ) { 
    
    corrPairs <- correlatedPairs(rasterLayers,speciesData=records,threhold=0.5,dataLayersMonotonocity=NULL)
    
    if( nrow(corrPairs) > 0) { rasters <- subset(rasterLayers,(1:length(names(rasterLayers)))[-which(names(rasterLayers) %in% corrPairs[,1])]) }
    if( nrow(corrPairs) == 0) { rasters <- rasterLayers }
    
  } else { rasters <- rasterLayers }
  
  # --------
  
  set.seed(1)
  presences.t <- records[sample(1:nrow(records), ifelse(autocorrelationSubsetRecords> nrow(records),nrow(records),autocorrelationSubsetRecords) ),dataRecordsNames]
  
  presences.environment <- data.frame(raster::extract(rasters,presences.t))
  to.remove <- unique(which(is.na(presences.environment),arr.ind = T)[,1])
  
  if(length(to.remove) > 0) {
    
    presences.t <- presences.t[-to.remove,]
    presences.environment <- raster::extract(rasterLayers,presences.t)
    
  }
  
  presences.environment <- presences.environment[,which(apply(presences.environment,2,var) != 0)]
  
  presences.environment.mean <- apply(presences.environment,2,mean, na.rm=T)
  presences.environment.sd <- apply(presences.environment,2,sd, na.rm=T)
  presences.environment.mean <- sweep(presences.environment,MARGIN=2,presences.environment.mean,'-')
  presences.environment <- sweep( presences.environment.mean , MARGIN=2,presences.environment.sd,'/')
  
  space <- spDists(as.matrix(presences.t),as.matrix(presences.t),longlat=TRUE)
  data <- ecodist::distance( presences.environment , method = "mahalanobis")
  
  n.class <- round(autocorrelationMaxDistance / autocorrelationClassDistance)
  
  resultsMatrix <- data.frame(classdistanceFrom=seq(0,autocorrelationMaxDistance-autocorrelationClassDistance,by=autocorrelationClassDistance),
                              classdistanceTo=seq(autocorrelationClassDistance,autocorrelationMaxDistance,by=autocorrelationClassDistance),
                              R=NA,
                              pVal=NA)
  
  for( i in 1:nrow(resultsMatrix)) {
    
    d1 = resultsMatrix[i,1]
    d2 = resultsMatrix[i,2]
    
    data.d <- c(as.matrix(data))
    space.d <- c(space)
    
    data.d[space.d < d1 | space.d > d2] <- NA
    space.d[space.d < d1 | space.d > d2] <- NA
    
    #data.d[space.d > d2] <- NA
    #space.d[space.d > d2] <- NA
    
    data.d <- data.d[!is.na(data.d)]
    space.d <- space.d[!is.na(space.d)]
    
    if(length(space.d) == 0) { next }
    if(length(unique(space.d)) == 1) { next }
    
    val <- cor.test(data.d, space.d, method=c("pearson"))
    val.corr <- val$estimate
    val.p.value <- val$p.value
    
    resultsMatrix[i,3] <- val.corr
    resultsMatrix[i,4] <- val.p.value
    
    # modelobject <- lm(space.d~data.d)
    # f <- summary(modelobject)$fstatistic
    # p <-0
    # tryCatch( p <- pf(f[1],f[2],f[3],lower.tail=FALSE) , error=function(e) { Error <<- TRUE })
    # resultsMatrix[i,3] <- summary(modelobject)$adj.r.squared
    # resultsMatrix[i,4] <- p
    
  }
  
  distance <- round( resultsMatrix[ which(resultsMatrix[,4] >= autocorrelationSignif)  , 2 ][1] )
  if( is.na(distance)) { distance <- autocorrelationMaxDistance }
  cat( paste0("\n"))
  cat( paste0("\n"))
  cat( paste0("First non-correlated distance: ",distance," km"))
  
  figure <- ggplot() +
    geom_hline(yintercept=0, size=0.1) +
    geom_line(data = resultsMatrix, aes(x = classdistanceTo, y = R), size = 0.1, color="Black", linetype = "dashed") +
    geom_point(shape=19, data = resultsMatrix[resultsMatrix$pVal <= autocorrelationSignif,], aes(x = classdistanceTo, y = R), size = 1.5, color="#6E6E6E") + 
    geom_point(shape=19, data = resultsMatrix[resultsMatrix$pVal > autocorrelationSignif,], aes(x = classdistanceTo, y = R), size = 1.5, color="Black") +
    geom_point(shape=19, data = resultsMatrix[resultsMatrix$pVal > autocorrelationSignif,][1,], aes(x = classdistanceTo, y = R), size = 1.75, color="Red") +
    ylab("Correlation (R)") + xlab("Geographic distance (Km)") + themePlot
  
  return( list(figure=figure, distance = distance) )
  
}

## -----------------------		

correlatedPairs <- function(rasterLayers,speciesData,threhold,dataLayersMonotonocity.i) { 
  
  list.of.cor <- data.frame()
  
  for( i in 1:(length(names(rasterLayers))-1)) {
    
    for( j in (i+1):length(names(rasterLayers))) {
      
      corr.val <- abs(cor( raster::extract(subset(rasterLayers,i),speciesData[,dataRecordsNames]) , raster::extract(subset(rasterLayers,j),speciesData[,dataRecordsNames]) ,use = "pairwise.complete.obs"))
      
      if(is.na(corr.val)) { corr.val <- 0 }
      
      if( !is.null(dataLayersMonotonocity.i)) {
        
        if( corr.val >= threhold & (sum(c(1,-1) %in% dataLayersMonotonocity.i[c(i,j)]) != 2)  ) {  list.of.cor <- rbind(list.of.cor,data.frame(Var.1=names(subset(rasterLayers,i)),
                                                                                                                                               Var.2=names(subset(rasterLayers,j)),
                                                                                                                                               Cor=corr.val,stringsAsFactors = FALSE)) }
      }
      
      if( is.null(dataLayersMonotonocity.i)) {
        
        if( corr.val >= threhold  ) {  list.of.cor <- rbind(list.of.cor,data.frame(Var.1=names(subset(rasterLayers,i)),
                                                                                   Var.2=names(subset(rasterLayers,j)),
                                                                                   Cor=corr.val,stringsAsFactors = FALSE)) }
      }
      
    }
    
  }
  
  return(list.of.cor)
  
}

## -----------------------

spatialThinning <- function(occurrenceRecords,minDistance,verbose) {
  
  dataThinning <- data.frame(Name="Sp",occurrenceRecords)

  coordinates.t <- thin(dataThinning,
                        lat.col = dataRecordsNames[2],
                        long.col = dataRecordsNames[1],
                        spec.col = "Name",
                        thin.par = minDistance,
                        reps = 1,
                        write.files = FALSE,
                        locs.thinned.list.return = TRUE,
                        verbose = FALSE)[[1]]
  
  if(verbose) {
    
    cat( paste0("\n"))
    cat( paste0("\n"))    
    cat( paste0("Input Records: ",nrow(occurrenceRecords)))
    cat( paste0("\n"))
    cat( paste0("Final Records: ",nrow(coordinates.t)))
    
  }
  
  # Remove from main dataset of occurrences
  colnames( coordinates.t ) <- dataRecordsNames
  
  # Remove log file
  if(length(list.files(".", full.names = TRUE, pattern = "spatial_thin_log.txt"))){
    file.remove( list.files(".", full.names = TRUE, pattern = "spatial_thin_log.txt") ) 
  }
  
  return(coordinates.t)
  
}

## -----------------------

generatePseudoAbsences <- function(rasterLayers,occurrenceRecords,n,spatialAutocorrelationVal) {
  
  extentOccurrenceRecords <- c(min(occurrenceRecords[,1]),max(occurrenceRecords[,1]),min(occurrenceRecords[,2]),max(occurrenceRecords[,2]))

  shape <- subset(rasterLayers,1)
  shape[!is.na(shape)] <- 1
  shape <- crop(shape,extent(c(extentOccurrenceRecords[1] - 20,
                               extentOccurrenceRecords[2] + 20,
                               extentOccurrenceRecords[3] - 20,
                               extentOccurrenceRecords[4] + 20)))
  
  shape <- xyFromCell( shape , Which(!is.na(shape), cell=TRUE) )
  colnames(shape) <- dataRecordsNames
  
  sink.points.poly <- as.data.frame(occurrenceRecords)
  coordinates( sink.points.poly ) <- dataRecordsNames
  proj4string( sink.points.poly ) <- CRS( "+proj=longlat +datum=WGS84" )
  
  sink.points.poly <- gBuffer( sink.points.poly, width=50 / 111.699, byid=TRUE )

  backgroundInformation.pts <- as.data.frame(shape)
  colnames( backgroundInformation.pts ) <- dataRecordsNames
  coordinates( backgroundInformation.pts ) <- dataRecordsNames
  proj4string( backgroundInformation.pts ) <- CRS( "+proj=longlat +datum=WGS84" )
  
  to.remove.id <- sp::over(backgroundInformation.pts,sink.points.poly)
  to.keep <- which(is.na(to.remove.id))
  shape <- shape[to.keep,]
  
  shape.i <- data.frame(1)
  ratioGeneration <- 1
  
  while( nrow(shape.i) < n ) {
    
    shape.i <- shape[sample(1:nrow(shape),min(c((n*ratioGeneration)+n,nrow(shape))), replace = FALSE),1:2]
    shape.i <- spatialThinning(shape.i,spatialAutocorrelationVal,verbose=FALSE)
    if( (n*ratioGeneration)+n >= nrow(shape) ) { break }
    ratioGeneration <- ratioGeneration + 1
    
  }

  shape.i <- shape.i[sample(1:nrow(shape.i), ifelse( nrow(shape.i) >= n , n , nrow(shape.i) ), replace = FALSE),1:2]
  
  return(shape.i)
  
}

## -----------------------

dataPartitioning <- function(speciesData,rasterLayers,type,k) {

  speciesData.i <- data.frame(speciesData)
  coordinates(speciesData.i) <- colnames(speciesData.i[,2:3])
  crs(speciesData.i) <- crs(rasterLayers)
  speciesData.i$Species <- rep("Sp",length(speciesData.i))
  
  distanceAutoCorrPlot <- NULL
  distanceAutoCorr <- NULL
  
  if( type=="spatialBlocksRowColumns" ) {

    sac <- spatialAutoRange(rasterLayer = rasterLayers,
                            sampleNumber = 5000,
                            doParallel = FALSE,
                            showPlots = FALSE)
    
    distanceAutoCorrPlot <- sac$plots$barchart
    
    distanceAutoCorr <- median(sac$rangeTable$range[-which.max(sac$rangeTable$range)])
    distanceAutoCorr <- min(distanceAutoCorr,round( ( abs(diff(c(extent(speciesData.i)[1],extent(speciesData.i)[2]))) * 111325 ) / k ) , round(( abs(diff(c(extent(speciesData.i)[3],extent(speciesData.i)[4]))) * 111325 ) / k ))
        
    crossValidation <- spatialBlock(speciesData = speciesData.i, # sf or SpatialPoints
                       species = "Sp", # the response column (binomial or multi-class)
                       rasterLayer = rasterLayers, # a raster for background (optional)
                       theRange = distanceAutoCorr, # size of the blocks in meters
                       k = k, # number of folds
                       selection = "random",
                       iteration = 10, # find evenly dispersed folds
                       seed=99,
                       biomod2Format = TRUE,
                       progress = FALSE,
                       showBlocks=FALSE,
                       verbose = FALSE)
    
  }

  if( type=="randomBlocks" ) {
    
    shape <- subset(rasterLayers,1)
    shape <- as(extent(shape), "SpatialPolygons")
    crs(shape) <- crs(rasterLayers)
    shape <- st_as_sf(shape)
    
    grd_lrg <- st_make_grid(shape, cellsize = c(round(max(maxCorrDistance.i*1000,111325) / 111325), round(max(maxCorrDistance.i*1000,111325) / 111325)))
    grd_lrg <- as_Spatial(grd_lrg)
    
    grd_lrgOver <- sp::over( grd_lrg , speciesData.i)
    grd_lrgOver <- which(!is.na(grd_lrgOver[,1]))
    grd_lrg <- grd_lrg[grd_lrgOver,]
    grd_lrg <- SpatialPolygonsDataFrame(grd_lrg, data.frame(ID= sample(1:k,length(grd_lrg),replace=TRUE) ), match.ID = F)
    
    crossValidation <- list()
    
    for(i in 1:k){
      
      grd_lrg.Test <- grd_lrg[grd_lrg$ID == i,]
      grd_lrg.Train <- grd_lrg[grd_lrg$ID != i,]
      
      if( nrow(grd_lrg.Test) == 0  | nrow(grd_lrg.Train) == 0) { next }
      
      grd_lrg.Test <- which(!is.na(sp::over( speciesData.i , grd_lrg.Test)[,1]))
      grd_lrg.Train <- which(!is.na(sp::over( speciesData.i , grd_lrg.Train)[,1]))
      
      crossValidation <- c(crossValidation,list(list(grd_lrg.Train, grd_lrg.Test)))
      
    }
    
    crossValidationPlot <<- grd_lrg
    
  }
  
  return( list(distanceAutoCorr=distanceAutoCorr,crossValidation=crossValidation,distanceAutoCorrPlot=distanceAutoCorrPlot ) )
  
}

## -----------------------

accuracyEstimate <- function(observed,predicted) {
  
  options(warn=-1)
  
  predicted.accuracy <- SDMTools::accuracy( observed , predicted , threshold = 100 )
  predicted.accuracy$deviance <- modEvA::Dsquared(obs = observed, pred=predicted, family = "binomial")
  predicted.accuracy$aicc <- AIC(glm(observed~predicted))
  
  # boyce
  boyce.p <- predicted[observed==1]
  boyce.p <- boyce.p[!is.na(boyce.p)]
  
  boyce.a <- predicted
  boyce.a <- boyce.a[!is.na(boyce.a)]
  
  if( length(boyce.p) == 1) { boyce.p <- c(boyce.p-0.01,boyce.p-0.005,boyce.p,boyce.p+0.005,boyce.p+0.01); }
  
  boyce.value <- ecospat.boyce(boyce.a , boyce.p ,PEplot=FALSE,rm.duplicate = TRUE)
  predicted.accuracy$boyce <- boyce.value$cor
  
  if(cvIndex == "boyce" ) { predicted.accuracy <- predicted.accuracy[which.min(abs(predicted.accuracy$threshold - min(  boyce.value$HS[ max( ifelse( length(which(boyce.value$F.ratio < 1)[-length(which(boyce.value$F.ratio < 1))]) > 0 , which(boyce.value$F.ratio < 1)[-length(which(boyce.value$F.ratio < 1))] , 1)) ]  ))),] }
  if(cvIndex == "area" ) { predicted.accuracy <- predicted.accuracy[max(which(predicted.accuracy$sensitivity > 0.95)),] }
  if(cvIndex == "tss" ) { predicted.accuracy <- predicted.accuracy[which.max(predicted.accuracy$specificity + predicted.accuracy$sensitivity),] }
  if(cvIndex == "auc") { predicted.accuracy <- predicted.accuracy[which.max(predicted.accuracy$AUC),] }
  
  # If non, auc by default
  if(nrow(predicted.accuracy) > 1) { predicted.accuracy <- predicted.accuracy[which.max(predicted.accuracy$AUC),] }
  
  predicted.accuracy <- data.frame( boyce = predicted.accuracy$boyce,
                                    threshold = predicted.accuracy$threshold ,
                                    auc = predicted.accuracy$AUC ,
                                    specificity = predicted.accuracy$specificity ,
                                    sensitivity = predicted.accuracy$sensitivity ,
                                    tss = predicted.accuracy$specificity + predicted.accuracy$sensitivity - 1 ,
                                    area = predicted.accuracy$sensitivity ,
                                    aicc = predicted.accuracy$aicc,
                                    deviance = predicted.accuracy$deviance )
  
  options(warn=0)
  
  return(predicted.accuracy)
  
}

## -----------------------

accuracyPredicted <- function(predicted.distribution,speciesData,type) {
  
  observed <- speciesData$PA
  predicted <- raster::extract(predicted.distribution,speciesData[,dataRecordsNames])
  predicted.accuracy <- accuracyEstimate(observed,predicted)
  predicted.accuracy$auc <- SDMTools::auc(observed ,predicted)

  return(predicted.accuracy)
  
}

## -----------------------

predictDistribution <- function(rasters,model,reclassToOne) {
  
  if( class(model)[1] == "MaxEnt" ) {
    
    predicted.distribution <- predict( model , rasters )
    predicted.distribution <- predicted.distribution / max(getValues(predicted.distribution),na.rm=TRUE)
    
  }   
  
  if( class(model)[1] == "gbm" ) {
    
    num.tress <- model$gbm.call$best.trees
    if(is.null(num.tress)) { num.tress <- length(model$trees) }
    predicted.distribution <- predict( rasters , model , n.trees=num.tress,type="response")
    
  }
  
  if( class(model)[1] == "mboost" ) {
    
    logit2prob <- function(logit){
      odds <- exp(logit)
      prob <- odds / (1 + odds)
      return(prob)
    }
    
    options(warn=-1)
    predicted.distribution <- predict( rasters , model)
    predicted.distribution <- logit2prob(predicted.distribution)
    
    options(warn=0)
    
  }
  
  if( class(model)[1] == "list" ) {
    
    shape <- subset(rasters,1)
    shape[!is.na(shape)] <- 1
    cellsPredict <- Which(!is.na(shape),cells=TRUE)
    cellsPredictValues <- as.matrix(rasters[cellsPredict])
    
    correctLayers <- names(rasters)[which(dataLayersMonotonocity == -1)]
    for(correctLayers.i in correctLayers) {
      cellsPredictValues[,correctLayers.i] <- cellsPredictValues[,correctLayers.i] * (-1)
    }
    
    predicted.distribution <- monmlp.predict( cellsPredictValues , model)
    shape[cellsPredict] <- as.numeric(predicted.distribution)
    predicted.distribution <- shape
    
  }
  
  if(reclassToOne) {
    predicted.distribution <- predicted.distribution + ( min(getValues(predicted.distribution),na.rm=T) * (-1))
    predicted.distribution <- predicted.distribution / max(getValues(predicted.distribution),na.rm=T)
  }
  
  return(predicted.distribution)
  
}

## -----------------------

reclassifyPredicted <- function(predicted.distribution,speciesData,method,reclassThreshold) {
  
  presences <- speciesData[speciesData$PA == 1,c("Lon","Lat")]  
  absences <- speciesData[speciesData$PA == 0,c("Lon","Lat")]  
  
  if(method == "directReclass") {
    
    
    predicted.distribution[predicted.distribution > reclassThreshold] <- 1
    predicted.distribution[predicted.distribution <= reclassThreshold] <- 0
    
  }
  
  if(method == "maxTSS") {
    
    observed <- c( rep(1,nrow(presences)) , rep(0,nrow(absences)) )
    predicted <- c( raster::extract(predicted.distribution,presences) , raster::extract(predicted.distribution,absences) )
    
    predicted.accuracy <- accuracy( observed , predicted , threshold =100 )
    predicted.accuracy <- predicted.accuracy[which.max(predicted.accuracy$specificity + predicted.accuracy$sensitivity),]
    
    predicted.distribution[predicted.distribution > predicted.accuracy$threshold] <- 1
    predicted.distribution[predicted.distribution <= predicted.accuracy$threshold] <- 0
    
  }
  
  if(method == "minAREA") {
    
    observed <- c( rep(1,nrow(presences)) , rep(0,nrow(absences)) )
    predicted <- c( raster::extract(predicted.distribution,presences) , raster::extract(predicted.distribution,absences) )
    
    predicted.accuracy <- accuracy( observed , predicted , threshold = 1000 )
    predicted.accuracy <- predicted.accuracy[predicted.accuracy$sensitivity > reclass.threshold,]
    predicted.accuracy <- predicted.accuracy[nrow(predicted.accuracy),]
    
    predicted.distribution[predicted.distribution > predicted.accuracy$threshold] <- 1
    predicted.distribution[predicted.distribution <= predicted.accuracy$threshold] <- 0
    
  }
  return(predicted.distribution)
}

## -----------------------

correctDateLineMap <- function(spatialObjectSF) {
  
  if(! TRUE %in% (class(spatialObjectSF) == "sf")) { stop("Object must be of class sf")}
  
  dggsOcean <- spatialObjectSF
  dggsOcean <- as(dggsOcean, "Spatial")
  listToRemove <- numeric()
  listToAdd <- list()
  
  for(i in 1:length(dggsOcean)) {
    
    if( extent(dggsOcean[i,])[1] < 0 & extent(dggsOcean[i,])[2] > 0 & max(extent(dggsOcean[i,])[1:2])-min(extent(dggsOcean[i,])[1:2]) > 180 ) { 
      
      listToRemove <- c(listToRemove,i)
      
      coords1 <- fortify(dggsOcean[i,])[,1:2]
      coords2 <- fortify(dggsOcean[i,])[,1:2]
      coords1[coords1[,1] < 0 ,1] <- 180
      coords2[coords2[,1] > 0 ,1] <- -180
      coords2 <- rbind(coords2,data.frame(long=-180,lat=max(coords2[,2])))
      coords2 <- rbind(coords2,data.frame(long=-180,lat=min(coords2[,2])))
      coords2 <- coords2[chull(coords2),] 
      
      coords1 <- spPolygons(as.matrix(coords1))
      coords2 <- spPolygons(as.matrix(coords2))
      
      listToAdd.i <- unionSpatialPolygons(union(coords1, coords2), c(1,1))
      listToAdd.i <- SpatialPolygonsDataFrame(listToAdd.i,as.data.frame(dggsOcean[i,]), match.ID = FALSE)
      
      if( length(listToAdd) >  0 ) { listToAdd <- bind(listToAdd, listToAdd.i) }
      if( length(listToAdd) == 0 ) { listToAdd <- listToAdd.i }
      
    }
    
    if( (extent(dggsOcean[i,])[1] > 0 & extent(dggsOcean[i,])[2] < 0) ) { stop(i) }
    
  }
  
  if( length(listToRemove) > 0) {
    
    dggsOcean <- dggsOcean[-listToRemove,]
    dggsOcean <- bind(dggsOcean,listToAdd)
    
  }

  dggsOcean <- st_as_sf(dggsOcean)
  
  return(dggsOcean)      
}


## ───────────────────────────────────────────────────────

fill_gaps_idw <- function(rlist, template,
                          idp  = 2,  # power
                          nmax = 12  # nearest neighbours
) {
  
  train_cells <- which(!is.na(values(rlist)))
  train_xy    <- xyFromCell(rlist, train_cells)
  train_val   <- values(rlist)[train_cells]
  train_pts <- vect(train_xy, crs = crs(rlist))
  train_pts <- as(train_pts, "Spatial")
  train_pts$z <- train_val
  train_pts <- as.data.frame(train_pts)
  
  gap_pts <- xyFromCell(rlist, which( is.na(rlist[][,1])  & ! is.na(template[][,1]) ))      # matrix with two columns
  
  nn <- RANN::nn2(data = as.matrix(train_pts[,c("x","y")]),query = as.matrix(gap_pts),k = 1)  
  preds <- train_pts[nn$nn.idx[,1],"z"]
  
  filled <- rlist
  filled[which( is.na(rlist[][,1])  & ! is.na(template[][,1]) )] <- preds
  return(filled)
}

