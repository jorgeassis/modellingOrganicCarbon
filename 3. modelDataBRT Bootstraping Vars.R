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

closeAllConnections()
rm(list=(ls()))
gc(reset=TRUE)

source("mainFunctions.R")
nCores <- 14

# Weighed averages for %OC (column:“avg_OC_percentW”), C-stock (“avg_CstockW”) and CAR (“carW”) for the top 10 cm of each core. 
response <- "avg_OC_percentW" # avg_CstockW avg_OC_percentW

resultsDirectory <- "../Results/"
resultsDirectory <- paste0(resultsDirectory,"/",response,"/")

# ------------------------------------------------------------------------------------

dataset <- loadRData(paste0("../Data/Database/finalDataBase ",response,".RData"))
maxObservedCarbon <- max(dataset$Response)
maxObservedCarbon

# Log the response
dataset$Response <- log(dataset$Response)
dataset <- dataset[dataset$Response > -4,]

# -----------------------------------------

dataset <- dataset[,c("Response","Bathymetry","DiffuseAttenuation","Oxygen","Distance.Shore","Landscape.Geomorphic.Features","Temperature.Max","Temperature.Min","River.Runoff","Salinity","Seawater.Speed","Slope","Wave.Fetch")]

names(dataset)[names(dataset) != "Response"]
dataLayersMonotonocity <- c(+1,+1,-1,-1 ,0,-1,+1,+1,+1,-1,-1,-1)
data.frame(names(dataset)[names(dataset) != "Response"],dataLayersMonotonocity)
names(dataLayersMonotonocity) <- names(dataset)[names(dataset) != "Response"]

climateLayersDirectory <- "../Data/Predictors/"
dataLayersFileType <- "tif"
climateLayers <- stack(list.files(climateLayersDirectory,pattern=dataLayersFileType,full.names = TRUE, recursive = FALSE))
list.files(climateLayersDirectory,pattern=dataLayersFileType)
dataLayers <- c("Bathymetry","DiffuseAttenuation","Oxygen","Distance Shore","Landscape Geomorphic Features","Nitrate","Temperature Max","Temperature Min","Phosphate","Available Radiation","River Runoff","Salinity","Seawater Speed","Slope","Terrain Ruggedness","Wave Fetch")
names(climateLayers) <- dataLayers

climateLayers <- climateLayers[[names(dataset)]]

# ------

climateLayers.i <- subset(climateLayers, which(names(climateLayers) == "Landscape.Geomorphic.Features"))
climateLayers.i <- as.factor(climateLayers.i)
climateLayers <- dropLayer(climateLayers, which(names(climateLayers) == "Landscape.Geomorphic.Features"))
climateLayers <- addLayer(climateLayers, climateLayers.i) 
names(climateLayers) %in% names(dataset)
gc(reset=TRUE)

Landscape.Geomorphic.Features <- read.csv("../Data/Predictors/geomorphicFeatures.csv")
dataset$Landscape.Geomorphic.Features <- as.factor(dataset$Landscape.Geomorphic.Features)

# -----------------------------------------

brtLearning <- seq(0.05, 1 , by=0.05)
brtTreeDepth <- 2:6
brtnTrees <- seq(50,500,50)
brtMinObsNode <- c(1,2,3,4,5,10,15,20,25) 
brtMinObsNode <- c(1,2,3,4,5,10,15,20,25) 
bagFraction <- c(0.5,0.75,1)

train.ratio <- 0.80
valid.ratio <- 0.1

bootRounds <- 1:7
cvRounds <- 1

# -----------------------------------------

bootAccuracyResultsDropVar <- data.frame()

contributionResults <- loadRData(paste0(resultsDirectory,"contributionResults.RData"))
contributionResults <- data.frame(Predictor=contributionResults[,1],
                                  Contribution=apply(contributionResults[,-1],1,mean),
                                  ContributionSD=apply(contributionResults[,-1],1,sd))

dropVarList <- contributionResults[which(contributionResults$Contribution > 5),"Predictor"]
dropVarList <- gsub(" ",".",dropVarList)

climateLayers.All <- climateLayers
dataset.All <- dataset
dataLayersMonotonocity.All <- dataLayersMonotonocity
bootAccuracyResults <- data.frame()

for( dropVar in dropVarList ){
  
    climateLayers <- dropLayer(climateLayers.All, which(names(climateLayers.All) == dropVar))
    dataset <- dataset.All[,which(names(dataset.All) != dropVar)]
    dataLayersMonotonocity <- dataLayersMonotonocity.All[which(names(dataLayersMonotonocity.All) != dropVar)]
    
    for( boot in bootRounds ){
    
      gc(reset=TRUE)
      cat("Boot round: ",boot,"\n")
      comb = expand.grid(cv.k = cvRounds, 
                         learning.complex=brtLearning,
                         tree.depth=brtTreeDepth,
                         trees=brtnTrees,
                         minObsNode=brtMinObsNode,
                         bagFraction=bagFraction)
      
      # Create shuffled indices
      total.rows <- nrow(dataset)
      set.seed( boot + 1 )
      shuffled.indices <- sample(1:total.rows, replace = FALSE)
      
      # Determine split boundaries
      train.index <- 1:floor(train.ratio * total.rows)
      valid.index <- (max(train.index) + 1):(max(train.index) + floor(valid.ratio * total.rows))
      test.index  <- (max(valid.index) + 1):total.rows
      
      # Subset data
      train.dataset <- dataset[shuffled.indices[train.index], ]
      valid.dataset <- dataset[shuffled.indices[valid.index], ]
      test.dataset  <- dataset[shuffled.indices[test.index], ]
      
      # -------
      
      cl <- parallel::makeCluster(20)
      registerDoParallel(cl)
      
      cv.accuracy <- foreach(c = 1:nrow(comb), .combine=rbind, .packages = c("gbm","dismo","SDMTools","ENMeval","modEvA")) %dopar% {
        
        cv <- comb[c,"cv.k"]
        l.rate <- comb[c,"learning.complex"]
        tree.c <- comb[c,"tree.depth"]
        trees <- comb[c,"trees"]
        minObsNode <- comb[c,"minObsNode"]
        bagFraction <- comb[c,"bagFraction"]
        
        model <- NULL
        tryCatch(       
          model <- gbm(formula = Response ~ ., distribution = "gaussian",
                       data = train.dataset, 
                       var.monotone = dataLayersMonotonocity, n.trees = trees,
                       interaction.depth = tree.c, shrinkage = l.rate,
                       n.minobsinnode = minObsNode,
                       bag.fraction = bagFraction, train.fraction = 1, cv.folds = 0,
                       keep.data = TRUE, verbose = FALSE, class.stratify.cv = NULL,
                       n.cores = NULL) , error=function(e) { model <<- NULL })
        
        if(!is.null(model)) {
          num.tress <- trees
          observed <- valid.dataset[,which(names(valid.dataset) == "Response")]
          predicted <- predict( model , valid.dataset[,-which(names(valid.dataset) == "Response")] , n.trees=num.tress,type="response")
          #model.deviance <- modEvA::Dsquared(glm(observed~predicted))
          model.deviance <- summary(lm(observed~predicted))$r.squared
          predicted.accuracy <- data.frame( cv.round=cv,tree.c=tree.c,l.rate=l.rate,trees=trees,minObsNode=minObsNode,bagFraction=bagFraction,accuracy=model.deviance)
          return(predicted.accuracy)
        }
        
      }
      
      stopCluster(cl); rm(cl) ; gc(reset=TRUE)
      closeAllConnections()
      
      # -------
      
      best.model <- aggregate(list(cvIndex=cv.accuracy[,"accuracy"]), by = list(tree.c=cv.accuracy$tree.c,l.rate=cv.accuracy$l.rate,trees=cv.accuracy$trees,minObsNode=cv.accuracy$minObsNode,bagFraction=cv.accuracy$bagFraction), mean)
      tree.c <- best.model[which.max(best.model$cvIndex),"tree.c"]
      l.rate <- best.model[which.max(best.model$cvIndex),"l.rate"]
      trees <- best.model[which.max(best.model$cvIndex),"trees"]
      minObsNode <- best.model[which.max(best.model$cvIndex),"minObsNode"]
      bagFraction <- best.model[which.max(best.model$cvIndex),"bagFraction"]
      
      model <- gbm(formula = Response ~ ., distribution = "gaussian",
                   data = train.dataset, 
                   var.monotone = dataLayersMonotonocity, n.trees = trees,
                   interaction.depth = tree.c, shrinkage = l.rate,
                   n.minobsinnode = minObsNode,
                   bag.fraction = bagFraction, train.fraction = 1, cv.folds = 0,
                   keep.data = TRUE, verbose = FALSE, class.stratify.cv = NULL,
                   n.cores = NULL)
      
      # -------
      
      meanAccuracyValid <- mean(cv.accuracy[cv.accuracy$tree.c == best.model[which.max(best.model$cvIndex),"tree.c"] & cv.accuracy$l.rate == best.model[which.max(best.model$cvIndex),"l.rate"] & cv.accuracy$trees == best.model[which.max(best.model$cvIndex),"trees"],"accuracy"])
      sdAccuracyValid <- sd(cv.accuracy[cv.accuracy$tree.c == best.model[which.max(best.model$cvIndex),"tree.c"] & cv.accuracy$l.rate == best.model[which.max(best.model$cvIndex),"l.rate"] & cv.accuracy$trees == best.model[which.max(best.model$cvIndex),"trees"],"accuracy"])
      
      observed <- test.dataset$Response
      predicted <- predict( model , test.dataset , n.trees=trees,type="response")
      accuracyTest <- summary(lm(observed~predicted))$r.squared
      corrTest <- cor(observed,predicted)
      RMSETest <- sqrt(mean((predicted - observed)^2))  
      
      observed <- train.dataset$Response
      predicted <- predict( model , train.dataset , n.trees=trees,type="response")
      accuracyTrain <- summary(lm(observed~predicted))$r.squared
      corrTrain <- cor(observed,predicted)
      RMSETrain <- sqrt(mean((predicted - observed)^2))  
    
      bootAccuracyResultsDropVar <- rbind(bootAccuracyResultsDropVar,data.frame(  dropVar=dropVar,RMSETrain=RMSETrain ))
      
      # ----------
      
      nCores.i <- 4
      carbonPrediction <- predictRasterTiled(climateLayers, model, tile_size_cells=10, output_dir="../Temp", filename_prefix = "prediction_tile", n.trees = model$n.trees, maxValue=maxObservedCarbon, nCores=nCores.i)
      gc(reset=TRUE)
      
      # Combining tiled predictions
      
      tile_paths <- list.files("../Temp",pattern="Log.tif",full.names = TRUE, recursive = TRUE)
      carbonMap <- rast(tile_paths[1])
      for( i in 2:length(tile_paths) ) { 
        carbonMap <- merge(carbonMap, rast(tile_paths[i])) 
        gc(reset=TRUE)
      }
      
      carbonMap[ carbonMap > log(100) ] <- log(100)
      writeRaster(carbonMap,filename=paste0(resultsDirectory,"/Predictions varRemove/fullModelVar",dropVar,"Round",boot,"Log.tif"),overwrite=T)
      
      file.remove(tile_paths)
      
      tile_paths <- list.files("../Temp",pattern="Exp.tif",full.names = TRUE, recursive = TRUE)
      carbonMap <- rast(tile_paths[1])
      for( i in 2:length(tile_paths) ) { 
        carbonMap <- merge(carbonMap, rast(tile_paths[i]) )
        gc(reset=TRUE)
      }

      carbonMap[ carbonMap > 100 ] <- 100
      writeRaster(carbonMap,filename=paste0(resultsDirectory,"/Predictions varRemove/fullModelVar",dropVar,"Round",boot,"Exp.tif"), overwrite=TRUE)
      
      file.remove(tile_paths)
      gc(reset=TRUE)
      
      # ----------
    
      save(bootAccuracyResultsDropVar, file=paste0(resultsDirectory,"/bootAccuracyResultsDropVar.RData"))

    }
}

bootAccuracyResultsDropVar

# -----------------------------------------

# Stack bootstrapped predictions

overwrite <- TRUE

for( index in c("Log","Exp")) {
  
  mainPrediction <- rast(paste0(resultsDirectory,"/Predictions/carbonMapMean",index,".tif"))
  mainPrediction <- app(mainPrediction, mean)
  
  for( dropVar in dropVarList ){
  
    if( ! overwrite & file.exists(paste0(resultsDirectory,"/Predictions varRemove/carbonMapMean",index,".tif")) ) { next }

    carbonMap <- list.files(paste0(resultsDirectory,"/Predictions varRemove/"),pattern=paste0(index,".tif"),full.names = TRUE, recursive = TRUE)
    carbonMap <- carbonMap[grepl(dropVar,carbonMap)]
    carbonMap <- rast(carbonMap)
    carbonMapMean <- app(carbonMap, mean)
    carbonMapDiff <- carbonMapMean - mainPrediction
    writeRaster(carbonMapDiff,filename=paste0(resultsDirectory,"/Predictions varRemove/carbonMapMeanVar",dropVar,index,".tif"),overwrite=T)
    
  }
}


# Delete files with Round

filestoRemove <- list.files(paste0(resultsDirectory,"/Predictions varRemove/"),pattern="Round",full.names = TRUE, recursive = TRUE)
filestoRemove
file.remove(filestoRemove)

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------