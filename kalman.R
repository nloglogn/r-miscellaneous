M <- function(x,A,B,Q){
  A * x + B * rnorm(1,sd=Q**0.5)
}
O <- function(x,C,R){
  C * x + rnorm(1,sd=R**0.5)
}
kalman <- function(A,B,C,Q,R,K,P){
  K1 <- P * C / (C * P * C + R)
  P1 <- P - K * C * P
  P2 <- A * P1 * A + B * Q * B
  return(c(K1,P2))
}

x0 <- 1
P0 <- 0.1
A <- 1.0
B <- 0.2
C <- 1.0
Q <- 0.2
R <- 0.05
time <- 100


X <- c(x0)
Y <- c(O(x0,C,R))
x_prev <- x0
for(i in 1:time){
  x <- M(x_prev,A,B,Q)
  x_prev <- x
  X <- c(X,x)
  Y <- c(Y,O(x,C,R))
}

KF <- function(X,Y,A,B,C,Q,R,P0){
  K0 <- kalman(A,B,C,Q,R,1,P0)[1]
  gain <- c(K0,P0)
  X_fore <- c(Y[1])
  x_fore <- X_fore[1]
  X_est <- c()
  K <- c(K0)
  P <- c(P0)
  for(i in 1:time){
    proc <- kalman(A,B,C,Q,R,gain[1],gain[2])
    gain[1] <- proc[1]
    gain[2] <- proc[2]
    K <- c(K,gain[1])
    P <- c(P,gain[2])
    x_est <- x_fore + gain[1] * (Y[i] - C * x_fore)
    X_est <- c(X_est, x_est)
    x_fore <- A * x_est
    X_fore <- c(X_fore,x_fore)
  }
  return(list(X_hat=c(Y[1],X_est),Cov=P))
}

nlogL <- function(parms){
  A <- parms[1]
  B <- parms[2]
  C <- parms[3]
  Q <- parms[4]
  R <- parms[5]
  result <- KF(X,Y,A,B,C,Q,R,P0)
  v <- Y - C * result$X_hat
  V <- C * result$Cov * C + R
  likelihood <- -1/2 * sum(log(V) + v / V * v)
  PL <- c(parms,likelihood)
  names(PL) <- c("A","B","C","Q","R","logL")
  print(PL)
  return(-likelihood)
}


init <- c(A,B,C,Q,R) + abs(rnorm(1,0,0.0001))
out_NM <- optim(par=init,fn=nlogL,method="Nelder-Mead",control = list(parscale=init,ndeps = rep(1.0e-7,length(init))))
out_LBB <- optim(par=out_NM$par,fn=nlogL,method="L-BFGS-B",lower=0.,upper=Inf,control = list(parscale=out_NM$par,ndeps = rep(1.0e-7,length(init))))
par_opt <- out_LBB$par

result <- KF(X,Y,par_opt[1],par_opt[2],par_opt[3],par_opt[4],par_opt[5],P0)

plot(X,type="l",xlim=c(0,time),ylim=c(0,2),lty=2,ylab="", xlab="")
par(new=T)
plot(result$X_hat, type="l", lty=1, col=2, xlim=c(0,time),ylim=c(0,2),xlab="",ylab="")
par(new=T)
plot(Y, type="l",lty=2,col=4,xlim=c(0,time),ylim=c(0,2), ylab="system with noise",xlab="time")
legend("topleft",legend=c("true","observed","estimated"),col=c(1,4,2),lty=c(2,2,1))
