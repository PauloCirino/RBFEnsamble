rm(list = ls())
source('./rbfBagging.R')
require('mlbench')
require('caret')
require('e1071')

Data <- mlbench::mlbench.2dnormals(n = 1000)
train_perc <- 0.5

sub_sample_perc <- 0.5
n_rbf_centers_param = c(10, 15)
n_rbf_centers_method = 'uniform'
r_val_param <- c(0.2, 0.8)
r_method <- 'uniform'
center_selection_method <- 'KMeans'
n_models <- 50

seed <- Sys.time()
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
                       center_selection_method = center_selection_method)

YhatM <- sign( randomRBFBags.predict(model = model, 
                                    Xtest = X_test) )
myModelResult <- caret::confusionMatrix(YhatM, Y_test)[['overall']][['Accuracy']]

svm_model   <- svm(x = X_train, y = Y_train)
YhatS <- sign( predict(svm_model, X_test) )
svmResult <- caret::confusionMatrix(YhatS, Y_test)[['overall']][['Accuracy']]

cat('myModelResult Accuracy = ', round(myModelResult * 100, 4))
cat('svmResult Accuracy = ', round(svmResult * 100, 4))

