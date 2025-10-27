## Yule-Walker method related functions

# function for gamma_hat

gamma_hat <- function(y){
  n <- length(y)
  max_lag <- 12
  ybar <- mean(y)
  g_hat <- rep(0,max_lag+1)
  K <- min(max_lag,n-1)
  h <- 0
  while(h<=K){
    g_hat[h+1] <- 0
    g_hat[h+1] <- g_hat[h+1]+sum((y[1:(n-h)]-ybar)*(y[(h+1):n]-ybar))
    h <- h+1
  }
  g_hat <- g_hat/n
  return(g_hat)
}

# function for Durbin-Levinson algorithm 
DL <- function(y,gamma_0,gamma_n){
  y_cent <- y-mean(y)
  n <- length(gamma_n)
  alpha <- rep(0,n)
  v <- rep(0,n+1)
  Phi <- matrix(0,nrow=n+1,ncol=n)
  v[1] <- gamma_0
  Phi[2,1] <- gamma_n[1]/gamma_0
  alpha[1] <- Phi[2,1]
  k <- 1
  while(k<=(n-1)){
    v[k+1] <- v[k]*(1-Phi[1+k,1]^2)
    Phi[1+k+1,1] <- (gamma_n[k+1]-sum(gamma_n[1:k]*Phi[1+k,1:k]))/v[k+1]
    Phi[1+k+1,(k+1):2] <- Phi[1+k,k:1]-Phi[1+k+1,1]*Phi[1+k,1:k]
    alpha[k+1] <- Phi[1+k+1,1]
    k <- k+1
  }
  return(Phi)
}





