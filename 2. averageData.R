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
setwd("~/Dropbox/Manuscripts/Modelling european blue carbon stocks/Scripts")
setwd("/Volumes/Dropbox/Dropbox/Manuscripts/Modelling european blue carbon stocks/Scripts")

closeAllConnections()
rm(list=(ls()))
gc(reset=TRUE)

source("mainFunctions.R")
library(readxl)

nCores <- 14

resultsDirectory <- "../Results/"
climateLayersDirectory <- "../Data/Predictors/"
dataLayersFileType <- "tif"

# ---------------

datasetRaw <- as.data.frame(read_excel("../Data/Database/Dataset 28-11-2024.xlsx", sheet=2))
datasetRaw <- datasetRaw[datasetRaw$Habitat != "Salt marsh",]
# datasetRaw <- as.data.frame(read.csv("../Data/Database/D4.5_avg10cm_core_final_ja.csv"))
# datasetRaw <- as.data.frame(readRDS("../Data/Database/bcdatawxsub1.rds"))

# Organic Carbon (%)

# ---------------

climateLayers <- stack(list.files(climateLayersDirectory,pattern=dataLayersFileType,full.names = TRUE, recursive = TRUE))
list.files(climateLayersDirectory,pattern=dataLayersFileType)
dataLayers <- c("Bathymetry","DiffuseAttenuation","Oxygen","Distance Shore","Landscape Geomorphic Features","Nitrate","Temperature Max","Temperature Min","Phosphate","Available Radiation","River Runoff","Salinity","Seawater Speed","Slope","Terrain Ruggedness","Wave Fetch")
names(climateLayers) <- dataLayers

# ---------------

locations <- unique(datasetRaw[,c("Longitude","Latitude")])
dataset <- raster::extract(climateLayers,locations)

# ---------------

dataNA.i <- which(is.na(dataset[,1]))
length(dataNA.i)

for( i in dataNA.i ) {
  
  cat( which(dataNA.i == i) , "out of" , length(dataNA.i), "\n")
  
  if( locations[i,1] < extent(climateLayers)[1] ) { locations[i,1] <- extent(climateLayers)[1] }
  if( locations[i,1] > extent(climateLayers)[2] ) { locations[i,1] <- extent(climateLayers)[2] }
  if( locations[i,2] < extent(climateLayers)[3] ) { locations[i,2] <- extent(climateLayers)[3] }
  if( locations[i,2] > extent(climateLayers)[4] ) { locations[i,2] <- extent(climateLayers)[4] }
  
  buffer <- 0.005
  repeat{
    climateLayers.i <- crop(climateLayers,c(locations[i,1] - buffer , locations[i,1] + buffer , locations[i,2] - buffer , locations[i,2] + buffer))
    climateLayersDF <- as.data.frame(climateLayers.i, xy=T, na.rm=T)
    if( nrow(climateLayersDF) > 2) { break }
    buffer <- buffer + 0.005
  }

  idw.nearest.r <- get.knnx( climateLayersDF[,c("x","y")] , locations[i,1:2], k=3 , algorithm="kd_tree" )
  idw.nearest.i <- idw.nearest.r$nn.index
  idw.nearest.d <- idw.nearest.r$nn.dist
  climateLayersDF <- climateLayersDF[idw.nearest.i,][,-c(1,2)]
  
  Landscape.Geomorphic.Mode <- Mode(climateLayersDF$Landscape.Geomorphic.Features)

  climateLayersDF <- as.vector(sapply( 1:ncol(climateLayersDF), function(x) { idw(idw.nearest.r$nn.dist, climateLayersDF[,x]) } ))
  names(climateLayersDF) <- names(climateLayers)
  climateLayersDF[which(names(climateLayersDF) == "Landscape.Geomorphic.Features")] <- Landscape.Geomorphic.Mode

  dataset[i,] <- climateLayersDF
  
}

head(dataset)
dataset <- as.data.frame(dataset)

dataset$Bathymetry <- abs(dataset$Bathymetry)
rownames(dataset) <- NULL
Landscape.Geomorphic.Features <- read.csv("../Data/Predictors/geomorphicFeatures.csv")
dataset$Landscape.Geomorphic.Features <- sapply(dataset$Landscape.Geomorphic.Features, function(x) { Landscape.Geomorphic.Features[Landscape.Geomorphic.Features$id == x , 2] })

# ---------------------

datasetRaw[which(datasetRaw$Methods == "CN-analyser"),"Methods"] <- "High-temperature combustion"
datasetRaw[which(datasetRaw$Methods == "CN-Analyser"),"Methods"] <- "High-temperature combustion"

dataset$OrganicCarbon <- NA
dataset$Method <- NA

for( l in 1:nrow(locations)) {
  
  loc <- which(datasetRaw$Longitude == locations[l,1] & datasetRaw$Latitude == locations[l,2])
  dataset[l,"OrganicCarbon"] <- mean(datasetRaw[loc,"Organic Carbon (%)"],na.rm=T)
  dataset[l,"Method"] <- modal(datasetRaw[loc,"Methods"],na.rm=T)

}

save(dataset, file="../Data/Database/finalDataBase.RData")

# ---------------------

vifResults <- vif(dataset[,-which(names(dataset) %in% c("Lon","Lat","Nitrate","Phosphate","Available.Radiation","Method","Landscape.Geomorphic.Features","OrganicCarbon"))])
vifResults
write.csv(vifResults,paste0(resultsDirectory,"pairsVIF.csv"))

# ---------------------

names(dataset) <- gsub("\\."," ",names(dataset))

library(ggstatsplot)
plot0 <- ggcorrmat(
  sig.level = 1,
  data = (dataset)[! names(dataset) %in% c("OrganicCarbon","Lon","Lat","Nitrate","Phosphate","Available.Radiation")],
  type = "parametric", # parametric for Pearson, nonparametric for Spearman's correlation
  colors = c("darkred", "white", "steelblue") # change default colors
)

pdf(file = paste0("../Results","/pairsPlot.pdf"), width=12, height=10 )
plot0
dev.off()

# ------------------------------------------------------------------------------------
