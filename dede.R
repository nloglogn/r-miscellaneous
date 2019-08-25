library(deSolve)

derivs <- function(t,y,parms){
  if(t < 0)
    lag <- 1
  else
    lag <- lagvalue(t - 1)
  
  dy <- r *y * (1.2 - lag)
  list(dy,dy=dy)
}

r <- 2

yinit <- 1

times <- seq(-1,100,0.01)

yout <- dede(y = yinit, times = times, func = derivs,
             parms <- c(r))

plot(yout, which=1, type="l", lwd=2,main="DDE")