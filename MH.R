# Metropolis-Hastings algorithm

parm_suggest <- function(x){
  return(x + rnorm(1,0,0.1))
}
get_likelihood <- function(mu){
  loglik <- sum(dnorm(samples,mu,1,log=T))
  return(loglik)
}
prob_threshold <- function(loglik_old,loglik_new){
  lik_ratio <- exp(loglik_new-loglik_old)
  P <- min(1, lik_ratio)
  return(P)
}
MH_procedure <- function(mu_old){
  mu_new <- parm_suggest(mu_old)
  loglik_old <- get_likelihood(mu_old)
  loglik_new <- get_likelihood(mu_new)
  P <- prob_threshold(loglik_old,loglik_new)
  u <- runif(1,0,1)
  if(u<P){
    return(mu_new)
  }else{
    return(mu_old)
  }
}
MH <- function(mu_init,t_mid,t_max){
  mu_set <- c(mu_init)
  mu_old <- mu_init
  for(t in 1:t_max){
    mu_new <- MH_procedure(mu_set[t])
    mu_set <- c(mu_set,mu_new)
  }
  return(mu_set[(t_mid+1):t_max])
}

set.seed(1234)

mu_true <- 2.345
sample_num <- 10000
samples <- rnorm(sample_num,mu_true,1)

t_mid <- 5000
t_max <- 10000

mu_init <- 1.0

out <- MH(mu_init,t_mid,t_max)
hist(out)