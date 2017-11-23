require('pdist')
require('corpcor')

generateR <- function(r_val_param, r_method){
    r <- numeric()
    
    if(r_method == 'uniform'){
        r <- runif(n = 1, min = r_val_param[1], max = r_val_param[2])
    } else if(r_method == 'normal'){
        r <- rnorm(n = 1, mean = mean(r_val_param[1]), sd = r_val_param[2])
    } else if(r_method == 'random'){
        r <- sample(x = r_val_param, size = 1)
    } else{
        r <- r_val_param[1]
    }
    
    r
}

generateNCenters <- function(n_rbf_centers_param,
                             n_rbf_centers_method){
    n_rbf_centers <- numeric()
    
    if(n_rbf_centers_method == 'uniform'){
        n_rbf_centers <- sample(x = n_rbf_centers_param[1] : n_rbf_centers_param[2],
                                size = 1)
    } else if(n_rbf_centers_method == 'normal'){
        n_rbf_centers <- rnorm(n = 1, mean = mean(n_rbf_centers_param[1]),
                               sd = n_rbf_centers_param[2])
        n_rbf_centers <- round(n_rbf_centers)
    } else if(n_rbf_centers_method == 'random'){
        n_rbf_centers <- sample(x = n_rbf_centers_param, size = 1)
    } else{
        n_rbf_centers <- n_rbf_centers_param[1]
    }
    
    n_rbf_centers
}

select_RBF_centers <- function(X, Y, r,
                               center_selection_method, 
                               n_rbf_centers_param, n_rbf_centers_method){
    
    n_rbf_centers <- generateNCenters(n_rbf_centers_param = n_rbf_centers_param,
                                      n_rbf_centers_method = n_rbf_centers_method)
    
    C <- matrix(data = NA, nrow = n_rbf_centers, ncol = ncol(X))
    if(center_selection_method == 'KMeans'){
        C <- kmeans(x = X, centers = n_rbf_centers, iter.max = 20)[['centers']]
    } else if(center_selection_method == 'Random'){
        centers_pos <- sample(x = nrow(X), size = n_rbf_centers)
        C <- X[centers_pos, ]
    } else{
        centers_pos <- sample(x = nrow(X), size = n_rbf_centers)
        C <- X[centers_pos, ]
    }
    
    C <- as.matrix(C)
    C
}


rbf <- function(X, Y, C, r) {
    ### Centers
    #C <- X
    
    ### Important Numbers
    n_obs <- nrow(X)
    n_vars <- ncol(X)
    n_centers <- nrow(C)
    
    ### Projecoes RBF
    D <- as.matrix( pdist(X = X, Y = C) )
    R <- 2*r**2
    H <- cbind( 1, exp(-0.5*(D/R)))
    
    ### Matriz de Projecao
    I <- diag(n_obs)
    A <- solve(t(H) %*% H)
    #P <- (I - H) %*% solve(A) %*% t(H)
    
    ### Cálculo de W
    W <- A %*% t(H) %*% Y  # find RBF weights
    
    model <- list(W = W, C = C, r = r)
}

rbf.predict <- function(model, Xtest) {
    W <- model$W
    C <- model$C
    r <- model$r
    
    ### Projecoes RBF
    D <- as.matrix( pdist(X = Xtest, Y = C) )
    R <- 2*r**2
    H <- cbind( 1, exp(-0.5*(D/R)))
    
    Ypred <- H %*% W
    Ypred
}

randomRBFBags <- function(X, Y, 
                          n_models,
                          sub_sample_perc = 0.5,
                          n_rbf_centers_param = c(5, 10), 
                          n_rbf_centers_method = 'uniform',
                          r_val_param = c(0.1, 0.9),
                          r_method = 'uniform',
                          center_selection_method = 'KMeans',
                          seed = runif(n = 1)
                          
) {
    set.seed(seed)
    X <- as.matrix(X)
    Y <- as.numeric(Y)
    n_total_obs <- nrow(X)
    
    models_n_centers <- numeric()
    models_radius <- numeric()
    models_list <- list()
    for(model_num in 1:n_models){
        iter_pos <- sample(x = n_total_obs, 
                           size = round(n_total_obs * sub_sample_perc) )
        iter_X <- X[iter_pos, ]
        iter_Y <- Y[iter_pos]
        
        iter_r <- generateR(r_val_param = r_val_param, 
                            r_method = r_method)
        
        C <- select_RBF_centers(X = iter_X, Y = iter_Y, r = iter_r,
                                center_selection_method = center_selection_method, 
                                n_rbf_centers_param = n_rbf_centers_param,
                                n_rbf_centers_method = n_rbf_centers_method)
        
        iter_model <- rbf(X = iter_X, Y = iter_Y, 
                          C = C, r = iter_r)
        
        models_n_centers[model_num] <- nrow(C)
        models_radius[model_num] <- iter_r
        models_list[[model_num]] <- iter_model
    }
    
    list(n_models = n_models,
         sub_sample_perc = sub_sample_perc,
         r_method = r_method,
         n_rbf_centers_param = n_rbf_centers_param, 
         n_rbf_centers_method = n_rbf_centers_method,
         r_val_param = r_val_param,
         center_selection_method = center_selection_method,
         seed = seed,
         models_n_centers = models_n_centers,
         models_radius = models_radius,
         models_list = models_list)
}


randomRBFBags.predict <- function(model, Xtest, task = 'regression'){
    n_models <- model$n_models
    n_test_samples <- nrow(Xtest)
    
    Yhat <- matrix(nrow = n_test_samples, ncol = n_models)
    for(model_num in 1:n_models){
        iter_model <- model$models_list[[model_num]]
        iter_Yhat <- rbf.predict(model = iter_model, Xtest = Xtest)
        Yhat[, model_num] <- iter_Yhat
    }
    
    if(task == 'classification'){
        Yhat <- sign ( apply(sign(Yhat), 1, mean) )
    } else{
        Yhat <- apply(Yhat, 1, mean)
    }
    Yhat
}


