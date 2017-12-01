rm(list = ls())
source('./rbfBagging.R')
require('mlbench')
require('caret')
require('e1071')

DataList <- list()
DataList[['2Normals']] <- mlbench.2dnormals(n = 1000, cl = 2)
DataList[['Circle']] <- mlbench.circle(n = 1000)
DataList[['Spirals']] <- mlbench.spirals(n = 1000, cycles = 2, sd = 0.1)
DataList[['XOR']] <- mlbench.xor(n = 1000)

SEEDS <- 1:10
n_models_vet <- c(16, 32, 64, 128, 256)
n_rbf_centers <- c(32, 64)

resultDF <- data.frame()
for(seed in SEEDS){
    for(n_rbf_centers_param in n_rbf_centers){
        for(n_models in n_models_vet){
            
            for(dataSet in names(DataList)){
                Data <- DataList[[dataSet]]
                train_perc <- 0.5
                
                sub_sample_perc <- 0.5
                n_rbf_centers_method = 'single'
                r_val_param <- c(0.05, 0.25)
                r_method <- 'uniform'
                center_selection_method <- 'KMeans'
                
                set.seed(seed = seed)
                X <- Data$x
                Y <- sign(x = as.numeric(Data$classes)*2 - 3)
                
                
                n_obs <- length(x = Y)
                n_train_obs <- round(x = 0.5 * n_obs)
                train_pos <- sample(x = n_obs, size = n_train_obs)
                
                X_train <- X[train_pos, ]
                Y_train <- Y[train_pos]
                X_test <- X[-train_pos, ]
                Y_test <- Y[-train_pos]
                
                
                model <- randomRBFBags(X = X_train,
                                       Y = Y_train, 
                                       n_models = n_models,
                                       sub_sample_perc = sub_sample_perc, 
                                       n_rbf_centers_param = n_rbf_centers_param, 
                                       n_rbf_centers_method = n_rbf_centers_method, 
                                       r_val_param = r_val_param,
                                       r_method = r_method,
                                       center_selection_method = center_selection_method,
                                       seed = seed,
                                       iterative = FALSE)
                
                YhatM <- sign( randomRBFBags.predict(model = model, 
                                                     Xtest = X_test,
                                                     task = 'classification' ) )
                myModelResult <- caret::confusionMatrix(YhatM, Y_test)[['overall']][['Accuracy']]
                
                iterDF <- data.frame(dataSet = dataSet, 
                                     accuracy = myModelResult,
                                     seed = seed,
                                     n_rbf_centers = n_rbf_centers_param,
                                     n_models = n_models)
                resultDF <- rbind(resultDF, iterDF)
                print(resultDF)
                write.csv(x = resultDF, file = 'resultFile2.csv')
            }
        }
    }
}