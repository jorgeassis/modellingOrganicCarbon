### Modelling Blue Carbon with high-performing machine learning

## **Overview**
This repository contains the code, data, and analysis framework for modeling organic carbon in marine environments using high-performing machine learning techniques. The primary goal is to predict the distribution of organic carbon in sediments using a set of biologically relevant predictors and advanced machine learning methods.

The analysis employs **Boosted Regression Trees (BRT)**, a powerful algorithm known for its ability to handle complex relationships, avoid overfitting, and provide high predictive accuracy. This repository includes scripts, data descriptions, and key results related to the study.

## **Main Outputs**
- **Model performance reports**:  
  - Deviance explained for each fold.  
  - Best hyperparameter combination.  

- **Plots and visualizations**:  
  - **Partial dependency plots** for each predictor.  
  - **Relative contribution plots** showing predictor importance.  
  - **Spatial prediction maps** of organic carbon distribution.  

- **Processed datasets**:  
  - Aggregated site-level datasets with organic carbon estimates.  

## **Requirements**
- **Software**: The analysis was conducted in **R** (R Development Core Team, 2023).  
- **Required packages**:  
  - `gbm` — For Boosted Regression Trees.  
  - `dismo` — For modeling and environmental data.  
  - `raster` — For working with gridded spatial data.  
  - `sp` — For handling spatial objects.  
  - `caret` — For cross-validation and hyperparameter tuning.  

## **References**

- Assis, J., et al. (2023). Bio-ORACLE v3.0: Global environmental data layers for marine species distribution modeling.
- Elith, J., et al. (2008). Boosted regression trees: A new technique for modeling ecological data.
Gouvêa, D., et al. (2024). Best practices for machine learning in biodiversity studies.
- Lamichhane, S., et al. (2019). Comparative analysis of machine learning models for soil organic carbon prediction.
- Fekete, B., et al. (2002). Global runoff data for water resources research.

## **License**

This project is licensed under the MIT License.

