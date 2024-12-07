## ---------------------------------------------------------------------
## ---------------------------------------------------------------------
##      
##  biodiversityDS [ biodiversitydatascience.com ]
##
## ---------------------
##
##  Machine Learning Modelling [ Ver.301 ]
##
## ---------------------------------------------------------------------
## ---------------------------------------------------------------------

setwd("/Users/jorgeassis/Dropbox/Manuscripts/Modelling european blue carbon stocks/Scripts")

closeAllConnections()
rm(list=(ls()))
gc(reset=TRUE)

source("mainFunctions.R")
resultsDirectory <- "../Results/"

## -----------

# Theme

colorWorldMap <- "#DCDCDC" # #575757
fillWorldMap <- "#DCDCDC" # #575757
themeWorldPosition <- "above"
hexagonSize <- 100
projection <- CRS("+proj=robin +over")

## -----------

themeMap <- 
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica", color = "#22211d"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(), # Graticulate element_line(color = "black", size = 3),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "#FFFFFF", color = NA), # F3F3F3
    panel.background = element_rect(fill = "#FFFFFF", color = NA), # F3F3F3
    legend.background = element_rect(fill = "#FFFFFF", color = NA),
    legend.box.background = element_rect(fill='#FFFFFF'),
    panel.border = element_blank()
  )

themeMapLegendBottom <- theme(legend.position="bottom",
                              legend.margin=margin(0,0,0,0),
                              legend.key.height= unit(0.25, 'cm'),
                              legend.key.width= unit(0.75, 'cm'),
                              legend.box="horizontal",
                              legend.title.position = "top",
                              legend.box.background = element_blank())

## ------------------------

bb <- sf::st_union(sf::st_make_grid( st_bbox(c(xmin = -180, xmax = 180, ymax = 90, ymin = -90), crs = st_crs(4326)), n = 100))
bb <- st_transform(bb, projection)

worldMap <- ne_countries(scale = 10, returnclass = "sf")
worldMap <- st_transform(worldMap, projection)
worldMap <- st_buffer(worldMap, dist = 0.001)
worldMap <- st_crop(worldMap, bb)

## ------------------------------------------------------------------------------
## ------------------------------------------------------------------------------

comb <- expand.grid(res=c(5,25,50), log=c("Log",""))

for( i in 1:nrow(comb)) {
  
  hexagonSize <- comb[i,1] # in km
  
  if( comb[i,2] == "Log" ) { 
    rasterMap <- raster( "../Results/carbonPredictionLog.tif" )
    mapName <- "carbonPredictionLog"
    legendName <- "Organic Carbon (log)"
  }
  if( comb[i,2] == "" ) { 
    rasterMap <- raster( "../Results/carbonPrediction.tif" )
    mapName <- "carbonPrediction"
    legendName <- "Organic Carbon"
  } 

  rasterMapDF <- data.frame(xyFromCell(rasterMap, Which( !is.na(rasterMap) , cells=T)),val=rasterMap[Which( !is.na(rasterMap) , cells=T)])
  dggs <- dgconstruct(spacing=hexagonSize, metric=TRUE)
  rasterMapDF$cell <- dgGEO_to_SEQNUM(dggs,rasterMapDF$x, rasterMapDF$y)$seqnum
  grid <- dgcellstogrid(dggs, rasterMapDF$cell)
  grid <- st_wrap_dateline(grid,options = c("WRAPDATELINE=YES",  "DATELINEOFFSET=180"), quiet = TRUE)
  grid$val <- exact_extract(rasterMap, grid, 'mean')
  grid <- st_transform(grid, projection)
  
  convexHull <- st_buffer(grid[,1], dist = 0.001)
  convexHull <- st_union(convexHull)
  
  minLegend <- min(grid$val)
  maxLegend <- max(grid$val)
  scaleType <- "continuous"
  
  #----------------------
  
  nColors <- 6 # 7
  colorBreaks <- seq( minLegend , maxLegend, length.out=nColors)
  colorBreaks[1] <- minLegend
  colorBreaks[nColors] <- maxLegend
  colorBreaksLegend <- round(colorBreaks)
  
  myColors <- c("#9DCEEF","#F0DA15", "#E47723","#93003a","#530655") # Blue Yellow Red Purple
  myColors <- colorRampPalette(myColors)(nColors)
  
  #----------------------
  
  plot1 <- ggplot() + 
    geom_sf(data = grid, aes(fill=val), colour = NA, size = 0.05) + # "black"
    scale_colour_gradientn(name=legendName,
                           colours = myColors, 
                           breaks= colorBreaks[c(1,3,5,7)], aesthetics = "fill", 
                           labels=colorBreaksLegend[c(1,3,5,7)], 
                           limits=c(minLegend,maxLegend) ) +
    geom_sf(data = convexHull, fill=NA, colour = "#878787", size = 0.01) +
    geom_sf(data = bb,fill=NA, colour = "#B2B2B2" , linetype='solid', linewidth= 0.75 ) +
    
    geom_sf(data=worldMap, color=colorWorldMap,  fill=fillWorldMap , size=0.1) +
    
    themeMapLegendBottom +
    coord_sf(crs= "+proj=laea +x_0=0 +y_0=0 +lon_0=22.75 +lat_0=45", xlim = c(-4492893, 1523591), ylim = c(-2041034, 4257578))
  
  # themeMap
  
  pdf(file=paste0("../Results/",mapName," Res",hexagonSize,".pdf"),width=12,useDingbats=FALSE)
  print(plot1)
  dev.off()
  
}

## ------------------------------------------------------------------------------
## ------------------------------------------------------------------------------
## Classes

bioRegionalizationResultsFolder <- "/media/Hammerhead/Data/Distribution Models/Global distribution of cold water corals/bioRegionalization/"
file <- paste0(bioRegionalizationResultsFolder,"/solutionK8.RData")

rasterMap <- loadRData( file )
rasterMap[rasterMap == 0] <- NA
rasterMap <- aggregate(rasterMap,fact=((hexagonSize / 1110) / res(rasterMap))[1],fun="max", na.rm=T)
rasterMapDF <- data.frame(xyFromCell(rasterMap, Which( !is.na(rasterMap) , cells=T)),val=rasterMap[Which( !is.na(rasterMap) , cells=T)])
rasterMapDF <- rasterMapDF[rasterMapDF$val != 0,]

dggs <- dgconstruct(spacing=hexagonSize, metric=TRUE)
rasterMapDF$cell <- dgGEO_to_SEQNUM(dggs,rasterMapDF$x, rasterMapDF$y)$seqnum
grid <- dgcellstogrid(dggs, rasterMapDF$cell)
grid <- st_wrap_dateline(grid,options = c("WRAPDATELINE=YES",  "DATELINEOFFSET=180"), quiet = TRUE)

library(exactextractr)
grid$val <- exact_extract(rasterMap, grid, 'max')
grid <- st_transform(grid, projection)

convexHull <- st_buffer(grid[,1], dist = 0.001)
convexHull <- st_union(convexHull)

#----------------------

grid$val <- as.factor(grid$val)
nColors <- 8
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector

myColors <- c("#9DCEEF","#F0DA15", "#E47723","#93003a","#530655","#7FC97F","#CAB2D6","#FBB4AE") # Blue Yellow Red Purple

#----------------------

if( themeWorldPosition == "above") {
  
  plot1 <- ggplot() + 
    geom_sf(data=worldMap, color=colorWorldMap,  fill=fillWorldMap , size=0.1) +
    geom_sf(data = grid, aes(fill=val), colour = NA, size = 0.05) + # "black"
    scale_fill_manual(values = myColors) +
    geom_sf(data = convexHull, fill=NA, colour = "#878787", size = 0.01)
  
}

if( themeWorldPosition == "below") {
  
  plot1 <- ggplot() + 
    geom_sf(data = grid, aes(fill=val), colour = NA, size = 0.05) + 
    scale_colour_manual(values = myColors) +
    geom_sf(data = convexHull, fill=NA, colour = "#878787", size = 0.01) +
    geom_sf(data=worldMap, color=colorWorldMap,  fill=fillWorldMap , size=0.1)
  
}

plot1 <- plot1 +
  geom_sf(data = bb,fill=NA, colour = "#B2B2B2" , linetype='solid', linewidth= 0.75 ) +
  themeMap  + themeMapLegendBottom

pdf(file=gsub(".RData",paste0("Res",hexagonSize,".pdf"),file),width=12,useDingbats=FALSE)
print(plot1)
dev.off()

## ------------------------------------------------------------------------------
## ------------------------------------------------------------------------------

# Zero centered map

stop("FIX")

files <- list.files( stackResultsFolder, pattern = "RData", full.names = T , recursive = T)
files <- files[grepl("Reachable",files)] # Reachable unConstrained
files <- files[!grepl("extinction",files)] 
files <- files[!grepl("Changes",files)] 
files <- files[!grepl("Uncertainty",files)] 
files <- files[grepl("speciesRichness",files)] 

# maskIntertidal <- raster("Dependencies/Data/Rasters/coastLineRes005.tif")

maskDistribution <- calc(stack(sapply(files, function(x) { loadRData( x ) } )), sum)
maskDistribution[maskDistribution > 0] <- 1
maskDistribution[maskDistribution == 0] <- NA

files <- list.files( stackResultsFolder, pattern = "RData", full.names = T , recursive = T)
filesNames <- list.files( stackResultsFolder, pattern = "RData", full.names = F , recursive = T)
files <- files[grepl("Reachable",files)] # Reachable unConstrained
filesNames <- filesNames[grepl("Reachable",filesNames)] # Reachable unConstrained
files <- files[grepl("Changes",files)] 
filesNames <- filesNames[grepl("Changes",filesNames)]; filesNames

autoLegend <- FALSE
autoLegendValues <- data.frame()

for( m in 1:length(files) ) {
  
  rasterMap <- loadRData( files[m] )
  mapName <- gsub(".RData","",filesNames[m])
  
  if( exists("maskIntertidal")) { rasterMap <- raster::mask(rasterMap,maskIntertidal) }
  
  rasterMap[intersect(Which(is.na(maskDistribution), cells=TRUE) ,Which(!is.na(rasterMap), cells=TRUE))] <- NA
  rasterMap[intersect(Which(maskDistribution == 1, cells=TRUE) ,Which(is.na(rasterMap), cells=TRUE))] <- 0
  
  resolutionH3 <- 3
  rasterMapDF <- data.frame(xyFromCell(rasterMap, Which( !is.na(rasterMap) , cells=T)),val=rasterMap[Which( !is.na(rasterMap) , cells=T)])
  
  cl <- makeCluster( detectCores() / 2 )
  clusterExport(cl, c("resolutionH3","rasterMapDF") )
  hexAddress <- parApply(cl, data.frame(rasterMapDF), 1, function(x) { h3js::h3_geo_to_h3(x[[2]], x[[1]], res = resolutionH3) } )
  stopCluster(cl) 
  
  rasterMapDF <- data.frame(rasterMapDF,hex=hexAddress)
  
  cl <- makeCluster( detectCores() / 2 )
  clusterExport(cl, c("resolutionH3","rasterMapDF") )
  hexAddressHexagon <- parLapply(cl, unique(rasterMapDF$hex), function(x) { mean(rasterMapDF[rasterMapDF$hex == x , "val"],na.rm=T) } )
  stopCluster(cl) 
  closeAllConnections()
  
  rasterMapDF <- data.frame(hex=unique(rasterMapDF$hex),val=unlist(hexAddressHexagon))
  rasterMapDF.polygons <- h3jsr::cell_to_polygon(input =rasterMapDF$hex, simple = FALSE)
  rasterMapDF.polygons$hex <- as.character(rasterMapDF$hex)
  rasterMapDF.polygons$value <- as.numeric(as.character(rasterMapDF$val))
  rasterMapDF.polygons <- correctDateLineMap(rasterMapDF.polygons)
  
  minLegend <- min(rasterMapDF.polygons$value, na.rm=T)
  maxLegend <- max(rasterMapDF.polygons$value, na.rm=T)
  autoLegendValues <- rbind(autoLegendValues,data.frame(mapName=mapName,minLegend=minLegend,maxLegend=maxLegend))
  
  if ( ! autoLegend ) {
    
    autoLegendValues <- read.csv(paste0(stackResultsFolder,"/","LegendValues_2.csv"))[,-1]
    minLegend <- min(autoLegendValues[,2:3])
    maxLegend <- max(autoLegendValues[,2:3])
    
  }
  
  # Round data to the nearest even integer
  minLegend <- 2 * round(minLegend / 2)
  maxLegend <- 2 * round(maxLegend / 2)
  
  # ---------
  
  hexagons <- st_transform(rasterMapDF.polygons, projection)
  hexagons <- as_Spatial(hexagons)
  hexagons <- gBuffer(hexagons, byid=TRUE, width=0.001)
  hexagons <- crop(hexagons, as(bb, "Spatial"))
  hexagons <- st_as_sf(hexagons)
  
  plot1 <- mainGlobalMapEqual +
    
    { if( themeWorldPosition == "below") geom_polygon(data = worldMap, aes(x = long, y = lat, group = group), fill=ifelse(themeWorldMap == "dark","#575757","#CDCDCD"), colour = "#CDCDCD" , size=0.1 ) } +
    
    geom_sf(data = hexagons, aes(fill=value), colour ="black", size = 0.05) +
    scale_colour_gradient2(low = "#800033",mid = "white",high = "#078050",midpoint = 0, 
                           space = "Lab",na.value = "grey50", guide = "colourbar", aesthetics = "fill", 
                           limits=c(minLegend,maxLegend) ,
                           breaks=c(minLegend,minLegend / 2,0,maxLegend/2,maxLegend),labels=c(minLegend,minLegend / 2,0,maxLegend/2,maxLegend) ) + 
    
    { if( themeWorldPosition == "above") geom_polygon(data = worldMap, aes(x = long, y = lat, group = group), fill=ifelse(themeWorldMap == "dark","#575757","#CDCDCD"), colour = "#CDCDCD" , size=0.1 ) } +
    
    geom_sf(data = bb,fill=NA, colour = "white" , linetype='solid', size= 3 ) +
    theme(legend.position="bottom",
          legend.margin=margin(0,0,0,0),
          legend.key.height= unit(0.25, 'cm'),
          legend.key.width= unit(0.75, 'cm')) + theme(legend.title=element_blank(), legend.box.background = element_blank())
  
  pdf(file=gsub(".RData",".pdf",files[m]),width=12,useDingbats=FALSE)
  print(plot1)
  dev.off()
  
  # Class based colors
  hexagons$categories <- as.factor(cut(hexagons$value, breaks = c(seq(minLegend,0,length.out=4)[-4],-0.49,0.49,seq(0,maxLegend,length.out=4)[-1]), labels = c( paste0(round(seq(minLegend,0,length.out=4)[-4]), " to " ,round(seq(minLegend,0,length.out=4)[-1])) , "0" , paste0(round(seq(0,maxLegend,length.out=4)[-4])," to ",round(seq(0,maxLegend,length.out=4)[-1])) )  ) )
  
  library(RColorBrewer)
  myColors <- rev(c("#078050","#098252","#6fd49e","#FFFFFF","#ffb6ce","#d15775","#800033"))
  names(myColors) <- levels(hexagons$categories)
  
  plot1 <- mainGlobalMapEqual +
    
    { if( themeWorldPosition == "below") geom_polygon(data = worldMap, aes(x = long, y = lat, group = group), fill=ifelse(themeWorldMap == "dark","#575757","#CDCDCD"), colour = "#CDCDCD" , size=0.1 ) } +
    
    geom_sf(data = hexagons, aes(fill=categories), colour ="black", size = 0.05) +
    scale_color_manual(values = myColors, aesthetics="fill", na.value = "grey50") +
    
    { if( themeWorldPosition == "above") geom_polygon(data = worldMap, aes(x = long, y = lat, group = group), fill=ifelse(themeWorldMap == "dark","#575757","#CDCDCD"), colour = "#CDCDCD" , size=0.1 ) } +
    
    geom_sf(data = bb,fill=NA, colour = "white" , linetype='solid', size= 3 ) +
    theme(legend.position="bottom",
          legend.margin=margin(0,0,0,0),
          legend.key.height= unit(0.25, 'cm'),
          legend.key.width= unit(0.75, 'cm')) + theme(legend.title=element_blank(), legend.box.background = element_blank()) +
    guides(fill = guide_legend(nrow = 1))
  
  pdf(file=gsub(".RData","Classes.pdf",files[m]),width=12,useDingbats=FALSE)
  print(plot1)
  dev.off()
  
}

autoLegendValues
write.csv(autoLegendValues,file=paste0(stackResultsFolder,"/","LegendValues_2.csv"))

## ------------------------------------------------------------------------------
## ------------------------------------------------------------------------------

# Extrapolation map

load(paste0(stackResultsFolder,"/RData/","speciesData.RData"))

dataLayersFileType <- "tif"
dataLayers <- c("CoastalExposureW_Pred_Max.tif","AirTemperature Surface Pred LtMax.tif","AirTemperature Surface Pred LtMin.tif","Precipitation Surface Pred Mean.tif","distanceToDelta.tif","Slope.tif")
dataLayersName <- c("CoastalExposure","TempMax","TempMin","PrecipMean","distanceToDelta","Slope") # 

climateLayersDirectory <- "../Data/Climate/Present/"
rasterLayers <- list.files(climateLayersDirectory,pattern=dataLayersFileType,full.names = TRUE, recursive = TRUE)
rasterLayers <- stack(rasterLayers[as.vector(sapply(dataLayers,function(x) { which( grepl(x,rasterLayers)) } ))])
rasterLayers <- processLayers(rasterLayers,occurrenceRecords=speciesData[speciesData$PA==1,2:3],regionBuffer=NULL, minDepth="NULL" , maxDepth="NULL",intertidal="Dependencies/Data/Rasters/BO2CoastLine.tif")
names(rasterLayers) <- dataLayersName

calibrationClimaticRegion <- raster::extract(rasterLayers,speciesData[,2:3])

climateLayersDirectory <- "../Data/Climate/RCP26/"
rasterLayers <- list.files(climateLayersDirectory,pattern=dataLayersFileType,full.names = TRUE, recursive = TRUE)
rasterLayers <- stack(rasterLayers[as.vector(sapply(dataLayers,function(x) { which( grepl(x,rasterLayers)) } ))])
rasterLayers <- processLayers(rasterLayers,occurrenceRecords=speciesData[speciesData$PA==1,2:3],regionBuffer=NULL, minDepth="NULL" , maxDepth="NULL",intertidal="Dependencies/Data/Rasters/BO2CoastLine.tif")
names(rasterLayers) <- dataLayersName
names(rasterLayers)

rasterMap <- subset(rasterLayers,2)
rasterMap[rasterMap <= max(calibrationClimaticRegion[,names(rasterMap)]) ] <- 0
rasterMap[rasterMap > 0] <- 1

rasterMap <- crop(rasterMap,extent(c(min(speciesData[,2]),max(speciesData[,2]),min(speciesData[,3]),max(speciesData[,3]))))

resolutionH3 <- 3
rasterMapDF <- data.frame(xyFromCell(rasterMap, Which( !is.na(rasterMap) , cells=T)),val=rasterMap[Which( !is.na(rasterMap) , cells=T)])
rasterMapDF <- data.frame(rasterMapDF,hex=apply(rasterMapDF[,1:2],1,function(x) { h3js::h3_geo_to_h3(x[[2]], x[[1]], res = resolutionH3) } ))
rasterMapDF <- data.frame(hex=unique(rasterMapDF$hex),val=sapply(unique(rasterMapDF$hex),function(x) { max(rasterMapDF[rasterMapDF$hex == x , "val"],na.rm=T) } ))
rasterMapDF.polygons <- h3jsr::h3_to_polygon(rasterMapDF$hex, simple = FALSE)
rasterMapDF.polygons$hex <- as.character(rasterMapDF$hex)
rasterMapDF.polygons$value <- as.numeric(as.character(rasterMapDF$val))
rasterMapDF.polygons <- correctDateLineMap(rasterMapDF.polygons)

minLegend <- min(rasterMapDF.polygons$value); minLegend
maxLegend <- max(rasterMapDF.polygons$value); maxLegend

myColors <- c("#C3E6EF","#e31515") # Blue Yellow Red Purple #450751

#----------------------

plot3 <- mainGlobalMap +
  geom_sf(data = rasterMapDF.polygons, aes(fill=value), colour ="black", size = 0.1) + # round signif
  scale_colour_gradientn(colours = myColors, breaks= colorBreaks, aesthetics = "fill", labels=signif(colorBreaks, digits = 3), limits=c(minLegend,maxLegend) ) +
  theme(legend.position="bottom",
        legend.margin=margin(0,0,0,0),
        legend.key.height= unit(0.25, 'cm'),
        legend.key.width= unit(0.75, 'cm')) + theme(legend.title=element_blank(), legend.box.background = element_blank())
plot3

pdf(file=paste0(stackResultsFolder,"/extrapolationMaxTempRCP26.pdf"),width=12,useDingbats=FALSE)
plot3
dev.off()

## ------------------------------------------------------------------------------
## ------------------------------------------------------------------------------

# Zero centered Climate map

list.files("../Data/Climate/Baseline/")
# Precipitation Surface Pred LtMin.tif

rasterMapP <- raster("../Data/Climate/Baseline/OceanTemperature BenthicMin Ltmax 2010-2020.tif")
rasterMapF <- raster("../Data/Climate/ssp585/OceanTemperature BenthicMin Ltmax 2090-2100.tif")
rasterMap <- rasterMapF - rasterMapP
rasterMap

maskOcc <- loadRData( files[17] )
maskOcc[maskOcc != 1] <- NA
rasterMap <- mask(rasterMap,maskOcc)
rasterMap <- crop(rasterMap,maskOcc)
rasterMap

resolutionH3 <- 3
rasterMapDF <- data.frame(xyFromCell(rasterMap, Which( !is.na(rasterMap) , cells=T)),val=rasterMap[Which( !is.na(rasterMap) , cells=T)])
rasterMapDF <- data.frame(rasterMapDF,hex=apply(rasterMapDF[,1:2],1,function(x) { h3js::h3_geo_to_h3(x[[2]], x[[1]], res = resolutionH3) } ))
rasterMapDF <- data.frame(hex=unique(rasterMapDF$hex),val=sapply(unique(rasterMapDF$hex),function(x) { mean(rasterMapDF[rasterMapDF$hex == x , "val"],na.rm=T) } ))
rasterMapDF.polygons <- h3jsr::h3_to_polygon(rasterMapDF$hex, simple = FALSE)
rasterMapDF.polygons$hex <- as.character(rasterMapDF$hex)
rasterMapDF.polygons$value <- as.numeric(as.character(rasterMapDF$val))
rasterMapDF.polygons <- correctDateLineMap(rasterMapDF.polygons)

minLegend <- min(rasterMapDF.polygons$value); minLegend
maxLegend <- max(rasterMapDF.polygons$value); maxLegend

nColors <- 7

colfuncLower <- colorRampPalette(c("#9A0F0F","#FCB46D", "#FFFFFF")) # FFFEC7
colfuncUpper <- colorRampPalette(c("#FFFFFF","#82DC9D","#11702E"))
colorBreaks <- seq( min(rasterMapDF.polygons$value),max(rasterMapDF.polygons$value), length.out=nColors)
myColors <- c(colfuncLower(sum(colorBreaks <= 0))[- which(colfuncLower(sum(colorBreaks <= 0)) == "#FFFFFF")],"#FFFFFF",colfuncUpper(sum(colorBreaks >= 0))[-1])

if( length(myColors) > 14 ) {
  myColors <- myColors[  which(sapply( colorBreaks , function(x) f(abs(x)) | x == 0)) ]
  colorBreaks <- colorBreaks[  which(sapply( colorBreaks, function(x) f(abs(x)) | x == 0)) ]
}

plot4 <- mainGlobalMap +
  geom_sf(data = rasterMapDF.polygons, aes(fill=value), colour ="black", size = 0.05) + # round signif
  scale_colour_gradientn(colours = myColors, breaks= colorBreaks, aesthetics = "fill", labels=signif(colorBreaks, digits = 3), limits=c(minLegend,maxLegend) ) +
  theme(legend.position="bottom",
        legend.margin=margin(0,0,0,0),
        legend.key.height= unit(0.25, 'cm'),
        legend.key.width= unit(0.75, 'cm')) + theme(legend.title=element_blank(), legend.box.background = element_blank())
plot4

pdf(file=paste0(stackResultsFolder,"/climateAnomalyMaxTemperaturessp585.pdf"),width=12,useDingbats=FALSE)
plot4
dev.off()

## ------------------------
