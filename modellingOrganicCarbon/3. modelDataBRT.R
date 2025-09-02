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
nCores <- 12

# Weighed averages for %OC (column:“avg_OC_percentW”), C-stock (“avg_CstockW”) and CAR (“carW”) for the top 10 cm of each core. 
response <- "avg_OC_percentW"
resultsDirectory <- "../Results/"
resultsDirectory <- paste0(resultsDirectory,"/",response,"/")

# ------------------------------------------------------------------------------------

dataset <- loadRData(paste0("../Data/Database/finalDataBase ",response,".RData"))
dataset <- dataset[,-which(names(dataset) == "Terrain.Ruggedness")]
dataset <- dataset[,-which(names(dataset) == "Available.Radiation")]
dataset <- dataset[,-which(names(dataset) == "Nitrate")]
dataset <- dataset[,-which(names(dataset) == "Phosphate")]
dataset <- dataset[complete.cases(dataset),]
hist(dataset$Response, breaks=100)

head(dataset)
nrow(dataset)
# 16252

min(dataset$Response)
max(dataset$Response)
# 0.001303 - 21.04514

maxObservedCarbon <- max(dataset$Response)

mean(dataset$Response)
sd(dataset$Response)
# 1.802079 ± 2.019387

names(dataset)
dataset$Response <- log(dataset$Response + 0.0000000000000001)
dataset$Landscape.Geomorphic.Features <- as.factor(dataset$Landscape.Geomorphic.Features)

hist(dataset$Response, breaks=100)
dataset <- dataset[which(dataset$Response > -4.5),]
hist(dataset$Response, breaks=100)

plot0 <- ggplot(dataset, aes(x = Response)) +
  geom_histogram(binwidth = 0.1, fill = "gray", color = "black", alpha = 0.7) +
  labs(x = "Observed distribution of log(Organic Carbon)",
       y = "Frequency") + themePlot

pdf(file = paste0(resultsDirectory,"/carbonHistogram.pdf"), width=10, height=8 )
plot0
dev.off()

names(dataset)[names(dataset) != "Response"]
dataLayersMonotonocity <- c(+1,+1,-1,-1,0,-1,+1,+1,+1,-1,-1,-1)
#dataLayersMonotonocity <- c(+1,+1,-1,-1,0,+1,-1,+1,+1,+1,+1,+1,-1,-1,-1)
data.frame(names(dataset)[names(dataset) != "Response"],dataLayersMonotonocity)

brtLearning <- seq(0.1,1, by=0.1)
brtTreeDepth <- 2:6
brtnTrees <- seq(50,500,50)

comb = expand.grid(cv.k = 1:3, learning.complex=brtLearning,tree.depth=brtTreeDepth , trees=brtnTrees )
nrow(comb)

cl <- parallel::makeCluster(nCores)
registerDoParallel(cl)

cv.accuracy <- foreach(c = 1:nrow(comb), .combine=rbind, .packages = c("gbm","dismo","SDMTools","ENMeval","modEvA")) %dopar% {
  
  cv <- comb[c,"cv.k"]
  l.rate <- comb[c,"learning.complex"]
  tree.c <- comb[c,"tree.depth"]
  trees <- comb[c,"trees"]
  
  set.seed(cv-2)
  sample <- sample(1:nrow(dataset),nrow(dataset)*0.80, replace = F)
  train.dataset <- dataset[sample,]
  test.dataset <- dataset[(1:nrow(dataset)) [!(1:nrow(dataset)) %in% sample],]
  
  model <- NULL
  tryCatch(       
    model <- gbm(formula = Response ~ ., distribution = "gaussian",
                 data = train.dataset, 
                 var.monotone = dataLayersMonotonocity, n.trees = trees,
                 interaction.depth = tree.c, shrinkage = l.rate,
                 bag.fraction = 0.5, train.fraction = 1, cv.folds = 0,
                 keep.data = TRUE, verbose = FALSE, class.stratify.cv = NULL,
                 n.cores = NULL) , error=function(e) { model <<- NULL })
  
  if(!is.null(model)) {
    num.tress <- trees
    observed <- test.dataset[,which(names(test.dataset) == "Response")]
    predicted <- predict( model , test.dataset[,-which(names(test.dataset) == "Response")] , n.trees=num.tress,type="response")
    #model.deviance <- modEvA::Dsquared(glm(observed~predicted))
    model.deviance <- summary(lm(observed~predicted))$r.squared
    predicted.accuracy <- data.frame( cv.round=cv,tree.c=tree.c,l.rate=l.rate,trees=trees,accuracy=model.deviance)
    return(predicted.accuracy)
  }
  
}

stopCluster(cl); rm(cl) ; gc(reset=TRUE)
closeAllConnections()

# ------------

best.model <- aggregate(list(cvIndex=cv.accuracy[,"accuracy"]), by = list(tree.c=cv.accuracy$tree.c,l.rate=cv.accuracy$l.rate,trees=cv.accuracy$trees), mean)
best.model[which.max(best.model$cvIndex),]
mean(cv.accuracy[cv.accuracy$tree.c == best.model[which.max(best.model$cvIndex),"tree.c"] & cv.accuracy$l.rate == best.model[which.max(best.model$cvIndex),"l.rate"] & cv.accuracy$trees == best.model[which.max(best.model$cvIndex),"trees"],"accuracy"])
sd(cv.accuracy[cv.accuracy$tree.c == best.model[which.max(best.model$cvIndex),"tree.c"] & cv.accuracy$l.rate == best.model[which.max(best.model$cvIndex),"l.rate"] & cv.accuracy$trees == best.model[which.max(best.model$cvIndex),"trees"],"accuracy"])
# 0.3464143 ± 0.003965979

tree.c <- best.model[which.max(best.model$cvIndex),"tree.c"]
l.rate <- best.model[which.max(best.model$cvIndex),"l.rate"]
trees <- best.model[which.max(best.model$cvIndex),"trees"]

# ------------

model <- gbm(formula = Response ~ ., distribution = "gaussian",
             data = dataset, 
             var.monotone = dataLayersMonotonocity, n.trees = trees,
             interaction.depth = tree.c, shrinkage = l.rate,
             bag.fraction = 0.5, train.fraction = 1, cv.folds = 0,
             keep.data = TRUE, verbose = FALSE, class.stratify.cv = NULL,
             n.cores = NULL)

# ------------

observed <- dataset$Response
predicted <- predict( model , dataset , n.trees=trees,type="response")
summary(lm(observed~predicted))$r.squared
# 0.4252491

cor(observed,predicted)
# 0.6521112

sqrt(mean((predicted - observed)^2))
# 0.8552927

plot(observed,predicted)

# ------------

dataPlot <- data.frame(observed=observed,predicted=predicted)
exp(mean(apply(dataPlot,1,diff)))
exp(sd(apply(dataPlot,1,diff)))
#  1.000353 ± 2.352123 units

plot1 <- ggplot(dataPlot, aes(x=observed, y=predicted)) +
  geom_point(size=2,shape=21, color="black", fill="black",alpha = 0.2) +
  geom_smooth(method=lm , color="white", fill="red", se=TRUE,size=0.3) +
  xlim(c( min(c(observed,predicted)) , max(c(observed,predicted)) )) + 
  ylim(c( min(c(observed,predicted)) , max(c(observed,predicted)) )) +
  #geom_abline(intercept = 0, slope = 1, color = "black", size = 0.3) +
  ylab("Predicted log(organic carbon)") + xlab("Observed log(organic carbon)") + themePlot +
  annotate(geom="text", x=min(dataPlot$observed), y=max(dataPlot$predicted), label=paste0("Dev. explained: ", round(summary(lm(observed~predicted))$r.squared,digits=3)),size=4,color = "#22211d",hjust = 0) +
  annotate(geom="text", x=min(dataPlot$observed), y=max(dataPlot$predicted) - 0.2, label=paste0("Pearson's Correlation: ", round(cor(observed,predicted),digits=3)),size=4, color = "#22211d",hjust = 0) +
  annotate(geom="text", x=min(dataPlot$observed), y=max(dataPlot$predicted) - 0.4, label=paste0("Root mean squared error: ", round(sqrt(mean((predicted - observed)^2)),digits=3)),size=4, color = "#22211d",hjust = 0) 

plot1

pdf(file = paste0(resultsDirectory,"/carbonFit.pdf"), width=10, height=9 )
plot1
dev.off()

dataset$ResponsePredicted <- predicted
save(dataset,file=paste0(resultsDirectory,"/carbonDataFit.RData"))
write.csv(dataset,paste0(resultsDirectory,"/carbonDataFit.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------------
# Variable Contribution

sum.brt <- summary(model) ; colnames(sum.brt) <- c("Predictor","Contribution")
sum.brt <- sum.brt[sort(sum.brt$Contribution,index.return=TRUE)$ix,]
sum.brt$Predictor <- gsub("\\."," ",sum.brt$Predictor)

reLativeImportancePlot <- ggplot(sum.brt) +
  geom_bar( aes(x= reorder(Predictor, Contribution) , y=Contribution), stat="identity", fill="black", alpha=0.5) +
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

pdf(paste0(resultsDirectory,"variableContribution.pdf"), width = 12, height = 9 )
reLativeImportancePlot
dev.off()

sum.brt <- sum.brt[sort(sum.brt$Contribution, index.return=T, decreasing = T)$ix,]
rownames(sum.brt) <- NULL
sum.brt

write.csv(sum.brt,paste0(resultsDirectory,"variableContribution.csv"))

# -------------
# Partial dependence plots

summaryModel <- sum.brt
summaryModel <- summaryModel[sort(summaryModel$Contribution,decreasing = TRUE,index.return=TRUE)$ix,]
summaryModel <- summaryModel[!summaryModel$Predictor %in% c("genus","study"),]
varsToPlot <- as.character(summaryModel[summaryModel$Contribution > 5,1])
varsToPlot.n <- 0
varsToPlotComposite <- character(0)

for( j in 1:length(varsToPlot) ) {
  
  i <- which(model$var.names == gsub(" ",".",varsToPlot[j]))
  
  plot.i <- plot(model,i)
  plotData <- data.frame(x=plot.i$panel.args[[1]]$x,
                         y=plot.i$panel.args[[1]]$y)
  
  if( varsToPlot[j] == "Landscape Geomorphic Features" | varsToPlot[j] == "Landscape.Geomorphic.Features" ) {
    
    partialPlot <- ggplot() +
      geom_point( data=plotData,aes(x=x, y=y),color="Black", size=2) +
      themePlot + 
      xlab(paste0( gsub("\\."," ",varsToPlot[j]) ," (",round(summaryModel[j,2],digits = 2),"%)")) + ylab("Effect on response")
    
  }
  
  if( varsToPlot[j] != "Landscape Geomorphic Features" & varsToPlot[j] != "Landscape.Geomorphic.Features" ) {
    
    partialPlot <- ggplot() +
      geom_line( data=plotData,aes(x=x, y=y),color="Black", size=0.5) +
      themePlot + 
      xlab(paste0( gsub("\\."," ",varsToPlot[j]) ," (",round(summaryModel[j,2],digits = 2),"%)")) + ylab("Effect on response")
    
  }

  assign( paste0("partialPlotVar",j) , partialPlot )
  
  varsToPlot.n <- varsToPlot.n + 1
  varsToPlotComposite <- paste0(varsToPlotComposite, ifelse( (j %% 2) == 0 ,  paste0("(","partialPlotVar",j," + ") , paste0("partialPlotVar",j," + ylab('') )",ifelse(which(varsToPlot == varsToPlot[j]) != length(varsToPlot) , " /","")) ))
  
}

pdf(file = paste0(resultsDirectory,"/partialPlotsFactor.pdf"), onefile=FALSE, width=9 )
partialPlotVar5
dev.off()

pdf(file = paste0(resultsDirectory,"/partialPlots.pdf"), onefile=FALSE, width=10, height=16 )
eval(parse(text = "(partialPlotVar1 + partialPlotVar2 + ylab('') ) / (partialPlotVar3 + partialPlotVar4 + ylab('') ) / (partialPlotVar6 + partialPlotVar7 + ylab('') ) / (partialPlotVar8 + ggplot() + theme(panel.background = element_rect(fill = 'white'),plot.background = element_rect(fill = 'white')) )"))
dev.off()

# eval(parse(text = "(partialPlotVar1 + partialPlotVar2 + ylab('') ) / (partialPlotVar3 + partialPlotVar4 + ylab('') ) / (partialPlotVar5 + partialPlotVar6 + ylab('') ) / (partialPlotVar7 + partialPlotVar8 + ylab('') ) / (partialPlotVar9 + ggplot() + theme(panel.background = element_rect(fill = 'white'),plot.background = element_rect(fill = 'white')) )"))

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------

climateLayersDirectory <- "../Data/Predictors/"
dataLayersFileType <- "tif"

climateLayers <- stack(list.files(climateLayersDirectory,pattern=dataLayersFileType,full.names = TRUE, recursive = TRUE))
list.files(climateLayersDirectory,pattern=dataLayersFileType)

dataLayers <- c("Bathymetry","DiffuseAttenuation","Oxygen","Distance Shore","Landscape Geomorphic Features","Nitrate","Temperature Max","Temperature Min","Phosphate","Available Radiation","River Runoff","Salinity","Seawater Speed","Slope","Terrain Ruggedness","Wave Fetch")
names(climateLayers) <- dataLayers

climateLayers.i <- subset(climateLayers, which(names(climateLayers) == "Landscape.Geomorphic.Features"))
Landscape.Geomorphic.Features <- read.csv("../Data/Predictors/geomorphicFeatures.csv")

for( i in 1:nrow(Landscape.Geomorphic.Features)) {
  climateLayers.i[ climateLayers.i == i ] <- dataset$Landscape.Geomorphic.Features[which(dataset$Landscape.Geomorphic.Features == Landscape.Geomorphic.Features[Landscape.Geomorphic.Features$id == i,2])[1]]
}
names(climateLayers.i) <- "Landscape.Geomorphic.Features"

climateLayers <- dropLayer(climateLayers, which(names(climateLayers) == "Landscape.Geomorphic.Features"))
climateLayers <- addLayer(climateLayers, climateLayers.i) 

carbonPrediction <- predict(  climateLayers , model, n.trees=model$n.trees,type="response")
carbonPrediction[carbonPrediction >= log(maxObservedCarbon) ] <- log(maxObservedCarbon)
writeRaster(carbonPrediction,filename=paste0(resultsDirectory,"/carbonPredictionLog.tif"),overwrite=T)

carbonPrediction <- exp(carbonPrediction)
carbonPrediction[carbonPrediction >= maxObservedCarbon ] <- maxObservedCarbon
writeRaster(carbonPrediction,filename=paste0(resultsDirectory,"/carbonPrediction.tif"),overwrite=T)

# -------------------