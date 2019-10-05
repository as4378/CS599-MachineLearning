
lasso <- function(Y, X){
 # check for class of Y and X
  if(class(X) != "data.frame"){
    X <- as.data.frame(X)
  }
  if(class(Y) != "data.frame"){
    Y <- as.data.frame(Y)
  }
  
  #standardize the predictors to have mean 0 and unit norm
  for(i in 1:ncol(X)){
    X[,i] <- (X[,i] - mean(X[,i]))/sd(X[,i])
  }
  
  #initial residual
  r = Y[,] - mean(Y[,])
  
  # initialize the initial coefficient matrix as all 0s
  beta <- matrix(0, nrow = ncol(X), ncol = 1)
  
  # find the correlation of residuals with active set of variables
  c <- cor(X, r)
  
  #TODO: next steps
}
