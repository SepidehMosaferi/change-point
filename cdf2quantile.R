## CDF to Quantile function

quantile_low <- function(a_n,b_n,c_n,d_n,xi_G,phi_G,p){
  G_L <- function(Q,p){-((2*pi)^-0.5)*(abs(Q)^0.5)*exp(-abs(Q)/8)-c_n*exp(a_n*abs(Q))*pnorm(-b_n*sqrt(abs(Q)))+
      (d_n-2+0.5*abs(Q))*pnorm(-0.5*sqrt(abs(Q)))-p}
  Q <- seq(-100,0,0.0001)
  LL <- G_L(Q,p)
  Data_L <- as.data.frame(cbind(Q,LL))
  result <- -abs(Data_L$Q[which.min(abs(Data_L$LL - 0))])
  return(result)
}

quantile_up <- function(a_p,b_p,c_p,d_p,xi_G,phi_G,p){
  G_U <- function(Q,p){1+xi_G*(phi_G*-0.5)*(1/sqrt(2*pi))*sqrt(Q)*exp(-(xi_G^2*(1/phi_G)*Q)/8)+
      c_p*exp(a_p*Q)*pnorm(-b_p*sqrt(Q))+(-d_p+2-0.5*xi_G^2+(1/phi_G)*Q)*
      pnorm(-0.5*xi_G*(1/sqrt(phi_G))*sqrt(Q))-p}
  Q <- seq(0,100,0.0001)
  UU <- G_U(Q,p)
  Data_U <- as.data.frame(cbind(Q,UU))
  result <- Data_U$Q[which.min(abs(Data_U$UU - 0))]
  return(result)
}

