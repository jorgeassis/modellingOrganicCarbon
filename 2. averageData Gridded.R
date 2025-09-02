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

nCores <- 10

# Weighed averages for %OC (column:“avg_OC_percentW”), C-stock (“avg_CstockW”) and CAR (“carW”) for the top 10 cm of each core. 
response <- "avg_CstockW" # avg_OC_percentW avg_CstockW

resultsDirectory <- "../Results/"
climateLayersDirectory <- "../Data/Predictors/"
dataLayersFileType <- "tif"

# ---------------

resultsDirectory <- paste0(resultsDirectory,"/",response,"/")
if(!dir.exists(resultsDirectory)) { dir.create(resultsDirectory, recursive = TRUE) }

# ---------------

datasetRaw <- as.data.frame(read_excel("../Data/Database/Data_2025_without_sm_final.xlsx", sheet=1))
datasetRaw <- datasetRaw[datasetRaw$Habitat != "Salt marsh",]

# Organic Carbon (%)

# ---------------

climateLayers <- stack(list.files(climateLayersDirectory,pattern=dataLayersFileType,full.names = TRUE, recursive = TRUE))
list.files(climateLayersDirectory,pattern=dataLayersFileType)
dataLayers <- c("Bathymetry","DiffuseAttenuation","Oxygen","Distance Shore","Landscape Geomorphic Features","Nitrate","Temperature Max","Temperature Min","Phosphate","Available Radiation","River Runoff","Salinity","Seawater Speed","Slope","Terrain Ruggedness","Wave Fetch")
names(climateLayers) <- dataLayers

# ---------------

locations <- datasetRaw[,c("Longitude","Latitude")]
dataset <- data.frame()

for( i in 1:nrow(locations) ) {
  
  cat(paste0("Processing location ",i," of ",nrow(locations),"\n"))
  
  if( locations[i,1] < extent(climateLayers)[1] ) { locations[i,1] <- extent(climateLayers)[1] }
  if( locations[i,1] > extent(climateLayers)[2] ) { locations[i,1] <- extent(climateLayers)[2] }
  if( locations[i,2] < extent(climateLayers)[3] ) { locations[i,2] <- extent(climateLayers)[3] }
  if( locations[i,2] > extent(climateLayers)[4] ) { locations[i,2] <- extent(climateLayers)[4] }
  
  climateLayers.i <- as.data.frame(raster::extract(climateLayers,locations[i,1:2]))
  climateLayers.i <- data.frame(x=locations[i,1],y=locations[i,2],climateLayers.i)
  
  buffer <- 0.002
  if( sum(complete.cases(climateLayers.i)) == 0 ) {
    repeat{
      climateLayers.i <- crop(climateLayers,c(locations[i,1] - buffer , locations[i,1] + buffer , locations[i,2] - buffer , locations[i,2] + buffer))
      climateLayers.i <- as.data.frame(climateLayers.i,xy=TRUE, na.rm=TRUE)
      climateLayers.i <- climateLayers.i[complete.cases(climateLayers.i),]
      if( sum(complete.cases(climateLayers.i)) > 0) { break }
      buffer <- buffer + 0.002
    }
  }

  climateLayers.i <- climateLayers.i[complete.cases(climateLayers.i),]
  climateLayers.i$x <- climateLayers.i$x[1]
  climateLayers.i$y <- climateLayers.i$y[1]
  climateLayers.i$Landscape.Geomorphic.Features <- modal(climateLayers.i$Landscape.Geomorphic.Features)
  climateLayers.i <- apply(climateLayers.i,2,mean)
  
  dataset.i <- data.frame(Response=datasetRaw[i,response], t(as.data.frame(climateLayers.i)))
  rownames(dataset.i) <- NULL
  dataset <- rbind( dataset , dataset.i )
  
}

# ---------------

# Export temp file as RData

save(dataset, file=paste0("../Data/Database/tempDataset ",response,".RData"))

# Test for NA

dataset.test <- raster::extract(climateLayers,dataset[,c("x","y")])
dataset.test <- as.data.frame(dataset.test)
dataNA.i <- which(!complete.cases(dataset.test))
length(dataNA.i)

# ---------------

datasetRaster <- rasterize(dataset[,c("x","y")],climateLayers,field=datasetRaw[,response],fun=mean)
datasetRaster.loc <- xyFromCell(datasetRaster,Which( !is.na(datasetRaster) , cell = TRUE ))

datasetFinal <- data.frame(Response = raster::extract(datasetRaster,datasetRaster.loc),
                           x = datasetRaster.loc[,1], y=datasetRaster.loc[,2],
                           raster::extract(climateLayers,datasetRaster.loc)
                           )

sum( is.na(datasetFinal))

# ---------------

datasetFinal$Bathymetry <- abs(datasetFinal$Bathymetry)
rownames(dataset) <- NULL
Landscape.Geomorphic.Features <- read.csv("../Data/Predictors/geomorphicFeatures.csv")
datasetFinal$Landscape.Geomorphic.Features <- sapply(datasetFinal$Landscape.Geomorphic.Features, function(x) { Landscape.Geomorphic.Features[Landscape.Geomorphic.Features$id == x , 2] })

# ---------------------

save(datasetFinal, file=paste0("../Data/Database/finalDataBase ",response,".RData"))

# ---------------------

vifResults <- vif(datasetFinal[,-which(names(datasetFinal) %in% c("x","y","Nitrate","Phosphate","Available.Radiation","Method","Landscape.Geomorphic.Features","Response"))])
vifResults
write.csv(vifResults,paste0(resultsDirectory,"pairsVIF.csv"))

# ---------------------

names(datasetFinal) <- gsub("\\."," ",names(datasetFinal))

library(ggstatsplot)
plot0 <- ggcorrmat(
  sig.level = 1,
  data = (datasetFinal)[! names(datasetFinal) %in% c("Response","x","y","Nitrate","Phosphate","Available.Radiation")],
  type = "parametric", # parametric for Pearson, nonparametric for Spearman's correlation
  colors = c("darkred", "white", "steelblue") # change default colors
)

pdf(file = paste0(resultsDirectory,"/pairsPlot.pdf"), width=12, height=10 )
plot0
dev.off()

# ------------------------------------------------------------------------------------
