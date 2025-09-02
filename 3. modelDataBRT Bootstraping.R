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
response <- "avg_CstockW" # avg_CstockW avg_OC_percentW

resultsDirectory <- "../Results/"
resultsDirectory <- paste0(resultsDirectory,"/",response,"/")

# ------------------------------------------------------------------------------------

dataset <- loadRData(paste0("../Data/Database/finalDataBase ",response,".RData"))
maxObservedCarbon <- max(dataset$Response)
maxObservedCarbon

# Log the response
dataset$Response <- log(dataset$Response)
dataset <- dataset[dataset$Response > -4,]
locations <- dataset[,c("x","y")]

# -----------------------------------------

dataset <- dataset[,c("Response","Bathymetry","DiffuseAttenuation","Oxygen","Distance.Shore","Landscape.Geomorphic.Features","Temperature.Max","Temperature.Min","River.Runoff","Salinity","Seawater.Speed","Slope","Wave.Fetch")]

names(dataset)[names(dataset) != "Response"]
dataLayersMonotonocity <- c(+1,+1,-1,-1 ,0,-1,+1,+1,+1,-1,-1,-1)
data.frame(names(dataset)[names(dataset) != "Response"],dataLayersMonotonocity)

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

if( dir.exists(paste0(resultsDirectory,"/Models")) == FALSE ) { dir.create(paste0(resultsDirectory,"/Models"), recursive = TRUE) }
if( dir.exists(paste0(resultsDirectory,"/Predictions")) == FALSE ) { dir.create(paste0(resultsDirectory,"/Predictions"), recursive = TRUE) }
if( dir.exists(paste0(resultsDirectory,"/Figures")) == FALSE ) { dir.create(paste0(resultsDirectory,"/Figures"), recursive = TRUE) }

bootAccuracyResults <- data.frame()
bootObservedPredictedPlot <- data.frame()
bootPartialPlots <- list()

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
  
  cl <- parallel::makeCluster(nCores)
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
               data = rbind(train.dataset,valid.dataset), 
               var.monotone = dataLayersMonotonocity, n.trees = trees,
               interaction.depth = tree.c, shrinkage = l.rate,
               n.minobsinnode = minObsNode,
               bag.fraction = bagFraction, train.fraction = 1, cv.folds = 0,
               keep.data = TRUE, verbose = FALSE, class.stratify.cv = NULL,
               n.cores = NULL)
  
  save(model, file=paste0(resultsDirectory,"/Models/modelBoot",boot,".RData"))
  
  # -------
  
  meanAccuracyValid <- mean(cv.accuracy[cv.accuracy$tree.c == best.model[which.max(best.model$cvIndex),"tree.c"] & cv.accuracy$l.rate == best.model[which.max(best.model$cvIndex),"l.rate"] & cv.accuracy$trees == best.model[which.max(best.model$cvIndex),"trees"],"accuracy"])
  sdAccuracyValid <- sd(cv.accuracy[cv.accuracy$tree.c == best.model[which.max(best.model$cvIndex),"tree.c"] & cv.accuracy$l.rate == best.model[which.max(best.model$cvIndex),"l.rate"] & cv.accuracy$trees == best.model[which.max(best.model$cvIndex),"trees"],"accuracy"])
  
  observed <- test.dataset$Response
  predicted <- predict( model , test.dataset , n.trees=trees,type="response")
  accuracyTest <- summary(lm(observed~predicted))$r.squared
  corrTest <- cor(observed,predicted)
  RMSETest <- sqrt(mean((predicted - observed)^2))  

  observed <- rbind(train.dataset,valid.dataset)$Response
  predicted <- predict( model , rbind(train.dataset,valid.dataset) , n.trees=trees,type="response")
  accuracyTrain <- summary(lm(observed~predicted))$r.squared
  corrTrain <- cor(observed,predicted)
  RMSETrain <- sqrt(mean((predicted - observed)^2))  

  bootAccuracyResults <- rbind(bootAccuracyResults,data.frame(boot=boot,
                                                              meanAccuracyValid=meanAccuracyValid,
                                                              sdAccuracyValid=sdAccuracyValid,
                                                              accuracyTest=accuracyTest,
                                                              corrTest=corrTest,
                                                              RMSETest=RMSETest,
                                                              accuracyTrain=accuracyTrain,
                                                              corrTrain=corrTrain,
                                                              RMSETrain=RMSETrain))

  # -------
  
  if( boot == 1 ) {
    observed <- dataset$Response
    bootObservedPredictedPlot <- data.frame(observed=observed,matrix(0,ncol=max(bootRounds),nrow=nrow(dataset)))
  }
  bootObservedPredictedPlot[,boot+1] <- predict( model , dataset , n.trees=trees,type="response")
  
  # -------
  
  sum.brt <- summary(model) ; colnames(sum.brt) <- c("Predictor","Contribution")
  sum.brt <- sum.brt[sort(sum.brt$Predictor,index.return=TRUE)$ix,]
  sum.brt$Predictor <- gsub("\\."," ",sum.brt$Predictor)

  if( boot == 1 ) {
    bootContributionResults <- sum.brt
  } else {
    bootContributionResults <- merge(bootContributionResults,sum.brt,by="Predictor",all=T)
    colnames(bootContributionResults) <- c("Predictor",1:boot)
  }
  
  # -------
  
  summaryModel <- sum.brt
  summaryModel <- summaryModel[sort(summaryModel$Contribution,decreasing = TRUE,index.return=TRUE)$ix,]
  summaryModel <- summaryModel[!summaryModel$Predictor %in% c("genus","study"),]
  varsToPlot <- as.character(summaryModel[summaryModel$Contribution > 5,1])
  varsToPlot.n <- 0
  varsToPlotComposite <- character(0)

  partialPlots <- list()
  
  for( j in varsToPlot ) {
    
    varsToPlot.j <- plot(model,which(gsub("\\."," ",model$var.names) == j), return.grid = TRUE)
    varsToPlot.j <- list(varsToPlot.j)
    varsToPlot.j <- setNames(varsToPlot.j, j)
    partialPlots <- cbind(partialPlots, varsToPlot.j)
    
  }
  
  bootPartialPlots <- c(bootPartialPlots,list(partialPlots))
  
  # ----------
  
  carbonPrediction <- predictRasterTiled(climateLayers, model, tile_size_cells=10, output_dir="../Temp", filename_prefix = "prediction_tile", n.trees = model$n.trees, maxValue=maxObservedCarbon, nCores=nCores)
  gc(reset=TRUE)
  
  # Combining tiled predictions
  
  tile_paths <- list.files("../Temp",pattern="Log.tif",full.names = TRUE, recursive = TRUE)
  carbonMap <- rast(tile_paths[1])
  for( i in 2:length(tile_paths) ) { carbonMap <- merge(carbonMap, rast(tile_paths[i]) ) }
  writeRaster(carbonMap,filename=paste0(resultsDirectory,"/Predictions/fullModelRound",boot,"Log.tif"),overwrite=T)
  
  file.remove(tile_paths)
  rm(carbonMap)
  gc(reset=TRUE)
  
  tile_paths <- list.files("../Temp",pattern="Exp.tif",full.names = TRUE, recursive = TRUE)
  carbonMap <- rast(tile_paths[1])
  for( i in 2:length(tile_paths) ) { carbonMap <- merge(carbonMap, rast(tile_paths[i]) ) }
  writeRaster(carbonMap,filename=paste0(resultsDirectory,"/Predictions/fullModelRound",boot,"Exp.tif"),overwrite=T)
  
  file.remove(tile_paths)
  rm(carbonMap)
  gc(reset=TRUE)
  
  # ----------

  save(bootPartialPlots, file=paste0(resultsDirectory,"/partialPlots.RData"))
  save(bootContributionResults, file=paste0(resultsDirectory,"/contributionResults.RData"))
  save(bootObservedPredictedPlot, file=paste0(resultsDirectory,"/observedPredictedPlot.RData"))
  save(bootAccuracyResults, file=paste0(resultsDirectory,"accuracyResults.RData"))
  gc(reset=TRUE)
  
}

# -----------------------------------------

# Stack bootstrapped predictions

carbonMap <- list.files(resultsDirectory,pattern="Log.tif",full.names = TRUE, recursive = TRUE)
carbonMap <- carbonMap[grepl("fullModelRound",carbonMap)]
carbonMap <- rast(carbonMap)
carbonMapMean <- app(carbonMap, mean)
carbonMapSD <- app(carbonMap, sd)

carbonMapMean[ carbonMapMean > log(100) ] <- log(100)

writeRaster(carbonMap,filename=paste0(resultsDirectory,"/Predictions/carbonMapMeanLog.tif"),overwrite=T)
writeRaster(carbonMapSD,filename=paste0(resultsDirectory,"/Predictions/carbonMapSDLog.tif"),overwrite=T)

carbonMap <- list.files(resultsDirectory,pattern="Exp.tif",full.names = TRUE, recursive = TRUE)
carbonMap <- carbonMap[grepl("fullModelRound",carbonMap)]
carbonMap <- rast(carbonMap)
carbonMapMean <- app(carbonMap, mean)
carbonMapSD <- app(carbonMap, sd)

carbonMapMean[ carbonMapMean > 100 ] <- 100

writeRaster(carbonMap,filename=paste0(resultsDirectory,"/Predictions/carbonMapMeanExp.tif"),overwrite=T)
writeRaster(carbonMapSD,filename=paste0(resultsDirectory,"/Predictions/carbonMapSDExp.tif"),overwrite=T)

gc(reset=TRUE)

# -----------------------------------------

# Get accuracy scores

accuracyResults <- loadRData(paste0(resultsDirectory,"accuracyResults.RData"))
accuracyResultsMean <- apply(accuracyResults,2,mean)
accuracyResultsSD <- apply(accuracyResults,2,sd)

accuracyResults <- data.frame(index=names(accuracyResultsMean),accuracyResultsMean,accuracyResultsSD)[-1,]
accuracyResults
save(accuracyResults, file=paste0(resultsDirectory,"accuracyResultsMean.RData"))
write.csv(accuracyResults,paste0(resultsDirectory,"accuracyResultsMean.csv"))

# -----------------------------------------

# Get observed predicted plot

observedPredictedPlot <- loadRData(paste0(resultsDirectory,"/observedPredictedPlot.RData"))
observedPredictedPlot <- data.frame(observed=observedPredictedPlot[,1],predicted=apply(observedPredictedPlot[,-1],1,mean))

accuracyVal <- summary(lm(observedPredictedPlot$observed~observedPredictedPlot$predicted))$r.squared
corrVal <- cor(observedPredictedPlot$observed,observedPredictedPlot$predicted)
RMSEVal <- sqrt(mean((observedPredictedPlot$predicted - observedPredictedPlot$observed)^2))  

plot1 <- ggplot(observedPredictedPlot, aes(x=observed, y=predicted)) +
  geom_point(size=1.7,shape=21, color="black", fill="black",alpha = 0.2) +
  geom_smooth(method=lm , color="white", fill="red", se=TRUE,size=0.3) +
  ylab("Predicted log(organic carbon)") + xlab("Observed log(organic carbon)") + themePlot +
  annotate(geom="text", x=min(observedPredictedPlot$observed), y=max(observedPredictedPlot$predicted), label=paste0("Dev. explained: ", round(accuracyVal,digits=3)),size=4,color = "#22211d",hjust = 0) +
  annotate(geom="text", x=min(observedPredictedPlot$observed), y=max(observedPredictedPlot$predicted) - 0.25, label=paste0("Pearson's Correlation: ", round(corrVal,digits=3)),size=4, color = "#22211d",hjust = 0) +
  annotate(geom="text", x=min(observedPredictedPlot$observed), y=max(observedPredictedPlot$predicted) - 0.5, label=paste0("Root mean squared error: ", round(RMSEVal,digits=3)),size=4, color = "#22211d",hjust = 0) 
plot1

pdf(file = paste0(resultsDirectory,"/carbonFit.pdf"), width=9, height=7 )
plot1
dev.off()

# ------

carbonMapMean <- raster(paste0(resultsDirectory,"/Predictions/carbonMapMeanExp.tif"))

regionsFile <- "../Data/EuropeanRegions/Regional_seas_around_Europe.gpkg"
regions <- terra::vect(regionsFile)

# Get average residuals per regions
regionsCarbonMapMean <- terra::extract(rast(carbonMapMean),regions)
regionsCarbonMapMean <- regionsCarbonMapMean[complete.cases(regionsCarbonMapMean),]

# Box plot with regionsResiduals

df <- regionsCarbonMapMean
df$name <- regions$name[match(df$ID,1:length(regions$name))]
df$name <- as.factor(df$name)
names(df) <- c("ID","layer","name")

# Create the boxplot
plot1 <- ggplot(df, aes(x = name, y = layer)) +
  geom_boxplot(outlier.shape = 19,outlier.alpha = 0.2,outlier.size = 1.5) +
  labs( y = "Organic carbon") + coord_flip() + themePlot + theme( axis.title.y = element_blank())

pdf(file = paste0(resultsDirectory,"/Figures/carbonPerRegion.pdf"), width=15, height=9 )
print(plot1)
dev.off()

png(file = paste0(resultsDirectory,"/Figures/carbonPerRegion.png"), width=1500, height=1100 )
print(plot1)
dev.off()

regionsCarbonMapMean <- aggregate(df$layer, by = list(df$ID), mean, na.rm=T)
regionsCarbonMapMeanSD <- aggregate(df$layer, by = list(df$ID), sd, na.rm=T)
regionsCarbonMapMeanDF <- data.frame(regions$name[as.numeric(rownames(regionsCarbonMapMean))],mean=regionsCarbonMapMean$x,sd=regionsCarbonMapMeanSD$x)

save(regionsCarbonMapMeanDF, file=paste0(resultsDirectory,"regionsCarbonDF.RData"))
write.csv(regionsCarbonMapMeanDF,paste0(resultsDirectory,"/Figures/regionsCarbonDF.csv"))

# -----------------------------------------

# Map residuals

for( boot in bootRounds) {
  
  model <- loadRData(paste0(resultsDirectory,"/Models/modelBoot",boot,".RData"))
  predicted <- predict( model , dataset , n.trees=model$n.trees,type="response")
  predicted <- exp(predicted)
  predicted[predicted > 100] <- 100
  residuals <- lm(exp(dataset$Response) ~ predicted)
  residuals <- residuals(residuals)
  
  if( boot == 1) { residualDF <- data.frame(residuals) }
  if( boot != 1) { residualDF <- cbind(residualDF,data.frame(residuals)) }
  
}

residualDF <- data.frame(locations,meanResiduals=apply(residualDF,1,mean))
residualsMap <- rasterize(residualDF[,1:2],climateLayers[[1]],field=residualDF[,3])
writeRaster(residualsMap,filename=paste0(resultsDirectory,"/Predictions/residuals.tif"),overwrite=T)

regionsFile <- "../Data/EuropeanRegions/Regional_seas_around_Europe.gpkg"
regions <- terra::vect(regionsFile)

# Get average residuals per regions
regionsResiduals <- terra::extract(rast(residualsMap),regions)
regionsResiduals <- regionsResiduals[complete.cases(regionsResiduals),]

# Box plot with regionsResiduals

df <- regionsResiduals
df$name <- regions$name[match(df$ID,1:length(regions$name))]
df$name <- as.factor(df$name)

# Create the boxplot
plot1 <- ggplot(df, aes(x = name, y = layer)) +
  geom_boxplot(outlier.shape = 19,outlier.alpha = 0.2,outlier.size = 1.5) +
  labs( y = "Residuals (organic carbon)") + coord_flip() + themePlot + theme( axis.title.y = element_blank()) +
  geom_hline(yintercept = 0, color = "black", size = 0.25)

pdf(file = paste0(resultsDirectory,"/carbonFitResidualPerRegion.pdf"), width=15, height=9 )
print(plot1)
dev.off()

regionsResidualsMean <- aggregate(regionsResiduals$layer, by = list(regionsResiduals$ID), mean, na.rm=T)
regionsResidualsSD <- aggregate(regionsResiduals$layer, by = list(regionsResiduals$ID), sd, na.rm=T)
regionsResidualsDF <- data.frame(regions$name[as.numeric(rownames(regionsResidualsMean))],mean=regionsResidualsMean$x,sd=regionsResidualsSD$x)

save(regionsResidualsDF, file=paste0(resultsDirectory,"regionsResidualsDF.RData"))
write.csv(regionsResidualsDF,paste0(resultsDirectory,"regionsResidualsDF.csv"))

# -----------------------------------------

# Partial plots

bootPartialPlots <- loadRData(paste0(resultsDirectory,"partialPlots.RData"))
bootRounds <- 1:length(bootPartialPlots)

contributionResults <- loadRData(paste0(resultsDirectory,"contributionResults.RData"))
contributionResults <- data.frame(Predictor=contributionResults[,1],
                                      Contribution=apply(contributionResults[,-1],1,mean),
                                      ContributionSD=apply(contributionResults[,-1],1,sd))

varsToPlot <- contributionResults$Predictor[sort(contributionResults$Contribution, index.return=T)$ix[sort(contributionResults$Contribution, index.return=T)$x > 5]]
varsToPlot <- rev(varsToPlot)

for( var in 1:length(varsToPlot) ) {
      
    var.name <- varsToPlot[var]
    var.name <- gsub("\\."," ",var.name)
    
    availableVars <- character(0)
    for( boot in bootRounds) {
      availableVars <- unique(c(availableVars,sapply(1:length(bootPartialPlots[[boot]]), function(x) { names(bootPartialPlots[[boot]][[x]])[names(bootPartialPlots[[boot]][[x]]) != "y"] })))
    }
    availableVars <- gsub("\\."," ",availableVars)
    
    if( ! var.name %in% availableVars ) { next }
    
    if( class(df) != "function" ) { rm(df) }
    
    for( boot in bootRounds) {
      
      namesBoot <- sapply(1:length(bootPartialPlots[[boot]]), function(x) { names(bootPartialPlots[[boot]][[x]])[names(bootPartialPlots[[boot]][[x]]) != "y"] })
      namesBoot <- gsub("\\."," ",namesBoot)
      loc <- which(namesBoot == varsToPlot[var])
      
      if( length(loc) == 0 ) { next }
      
      df.b <- bootPartialPlots[[boot]][[loc]]
      names(df.b) <- c("x","y")
      
      if( varsToPlot[var] != "Landscape Geomorphic Features" & varsToPlot[var] != "Landscape.Geomorphic.Features" ) {
        
        df.b.model <- loess(y ~ x, data = df.b, span = 0.1,degree=0)
        
        if( ! exists("df") | class(df) == "function" ) {
          x.val <- df.b$x
          df <- data.frame(var=x.val)
          names(df) <- var.name
        }
        
        df$y <- predict(df.b.model, newdata = x.val)
        
      }
      
      if( varsToPlot[var] == "Landscape Geomorphic Features" | varsToPlot[var] == "Landscape.Geomorphic.Features" ) {

        if( ! exists("df") | class(df) == "function" ) {
          x.val <- df.b$x
          df <- data.frame(var=x.val)
          names(df) <- var.name
        }
        
        df[match(df.b[,1],df[,1]),"y"] <- df.b$y
        
      }
      
      names(df)[which(names(df) == "y")] <- boot
      
    }

    df <- df[,c(1,which(!is.na(apply(df[,-1],2,sum)))+1)]
    df <- data.frame( x  = df[,1], y  = apply(df[,-1],1,mean,na.rm=T), sd = apply(df[,-1],1,sd) )
    
    if( varsToPlot[var] == "Landscape Geomorphic Features" | varsToPlot[var] == "Landscape.Geomorphic.Features" ) {
    
      partialPlot <- ggplot(df, aes(x = x, y = y)) +
        geom_errorbar(aes(ymin = y - sd, ymax = y + sd), width = 0.2, size=0.35, color="gray", alpha=1) +
        geom_point( aes(x=x, y=y),color="Black", size=1.5) + 
        themePlot +
        xlab(paste0( var.name)) + 
        ylab("Effect on response")
      
    }
    
    if( varsToPlot[var] != "Landscape Geomorphic Features" & varsToPlot[var] != "Landscape.Geomorphic.Features" ) {
        
        partialPlot <- ggplot(df, aes(x = x, y = y)) +
          geom_ribbon(aes(ymin = y - sd, ymax = y + sd),fill = "gray", alpha = 0.3  ) +
          geom_line(color="Black", size=0.5) + themePlot +
          xlab(paste0( var.name)) + 
          ylab("Effect on response")
      
    }
    
    assign( paste0("partialPlotVar",var) , partialPlot )

}

pdf(file = paste0(resultsDirectory,"/partialPlotsFactor.pdf"), onefile=FALSE, width=9 )
partialPlotVar6
dev.off()

pdf(file = paste0(resultsDirectory,"/partialPlots.pdf"), onefile=FALSE, width=10, height=16 )
eval(parse(text = "(partialPlotVar1 + partialPlotVar2 + ylab('') ) / (partialPlotVar3 + partialPlotVar4 + ylab('') ) / (partialPlotVar5 + partialPlotVar7 + ylab('') ) / (partialPlotVar8 + ggplot() + theme(panel.background = element_rect(fill = 'white'),plot.background = element_rect(fill = 'white'))  )  ")) # # (partialPlotVar8 + ggplot() + theme(panel.background = element_rect(fill = 'white'),plot.background = element_rect(fill = 'white')) )
dev.off()

# -----------------------------------------

# Variable Contribution

bootContributionResults <- loadRData(paste0(resultsDirectory,"contributionResults.RData"))
bootContributionResults <- data.frame(Predictor=bootContributionResults[,1],
                                      Contribution=apply(bootContributionResults[,-1],1,mean),
                                      ContributionSD=apply(bootContributionResults[,-1],1,sd))

reLativeImportancePlot <- ggplot(bootContributionResults, aes(x= reorder(Predictor, Contribution) , y=Contribution)) +
  geom_bar( stat="identity", fill="#F0D461", alpha=1) +
  geom_errorbar( aes(ymin = Contribution - ContributionSD, ymax = Contribution + ContributionSD), width = 0.2, color="#333333") +
  coord_flip() + theme(
    axis.text=element_text(size=10),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0) , size=12),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0) , size=12),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "#EFEFEF", colour = "#EFEFEF",size = 0, linetype = "solid"),
    panel.grid.major = element_line(size = 0, linetype = 'solid',colour = "#EFEFEF"), 
    panel.grid.minor = element_line(size = 0, linetype = 'solid',colour = "#EFEFEF")
  ) + labs(x = "Predictor") + 
  labs(y = "Relative importance (%)") + geom_hline(aes(yintercept=5 ),color="Black", linetype="dashed", size=0.3) +
  annotate("text", y = 5 + 1 , x = 1 , label = "5%" , hjust = 0) + geom_hline(aes(yintercept=0 ),color="Gray", size=0.3)

reLativeImportancePlot

pdf(paste0(resultsDirectory,"variableContribution.pdf"), width = 12, height = 9 )
reLativeImportancePlot
dev.off()

write.csv(bootContributionResults,paste0(resultsDirectory,"variableContribution.csv"))

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------