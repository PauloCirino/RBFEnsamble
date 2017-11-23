rbf <- function(X, Y, K=10, gamma=1.0) {
    N     <- dim(X)[1] # number of observations
    ncols <- dim(X)[2] # number of variables
    
    ### RBF CENTERS
    repeat {
        km <- kmeans(X, K)  # let's cluster K centers out of the dataset
        if (min(km$size)>0) # only accept if there are no empty clusters
            break
    }
    mus <- km$centers # the clusters points
    
    ### RBF MODEL
    Phi <- matrix(rep(NA,(K+1)*N), ncol=K+1)
    for (lin in 1:N) {
        Phi[lin,1] <- 1    # bias column
        for (col in 1:K) {
            Phi[lin,col+1] <- exp( -gamma * norm(as.matrix(X[lin,]-mus[col,]),"F")^2 )
        }
    }
    
    w <- pseudoinverse(t(Phi) %*% Phi) %*% t(Phi) %*% Y  # find RBF weights
    
    list(weights=w, centers=mus, gamma=gamma)  # return the rbf model
}

rbf.predict <- function(model, X, classification=FALSE) {
    
    gamma   <- model$gamma
    centers <- model$centers
    w       <- model$weights
    N       <- dim(X)[1]    # number of observations
    
    pred <- rep(w[1],N)  # we need to init to a value, so let's start with the bias
    
    for (j in 1:N) {  
        # find prediction for point xj
        for (k in 1:length(centers[,1])) {
            # the weight for center[k] is given by w[k+1] (because w[1] is the bias)
            pred[j] <- pred[j] + w[k+1] * exp( -gamma * norm(as.matrix(X[j,]-centers[k,]),"F")^2 )
        }
    }
    
    if (classification) {
        pred <- unlist(lapply(pred, sign))
    }
    pred
}