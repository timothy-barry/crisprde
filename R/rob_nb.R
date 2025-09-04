link <- function(mu){log(mu)}
invlink <- function(eta){exp(eta)}
derivlink <- function(mu){1/mu}
derivinvlink <- function(eta){exp(eta)}
varfunc <- function(mu,sigma){mu+sigma*mu^2}

tukeypsi <- function(r,c.tukey){
  ifelse(abs(r)>c.tukey,0,((r/c.tukey)^2-1)^2*r)
}

E.tukeypsi.1 <- function(mui,sigma,c.tukey){
  sqrtVmui <- sqrt(varfunc(mui,sigma))
  j1 <- max(c(ceiling(mui-c.tukey*sqrtVmui),0))
  j2 <- floor(mui+c.tukey*sqrtVmui)
  if (j1>j2){0}
  else {
    j12 <- j1:j2
    sum((((j12-mui)/(c.tukey*sqrtVmui))^2-1)^2*(j12-mui)*dnbinom(j12,mu=mui,size=1/sigma))/sqrtVmui
  }
}

E.tukeypsi.2 <- function(mui,sigma,c.tukey){
  sqrtVmui <- sqrt(varfunc(mui,sigma))
  j1 <- max(c(ceiling(mui-c.tukey*sqrtVmui),0))
  j2 <- floor(mui+c.tukey*sqrtVmui)
  if (j1>j2){0}
  else {
    j12 <- j1:j2
    sum((((j12-mui)/(c.tukey*sqrtVmui))^2-1)^2*(j12-mui)^2*dnbinom(j12,mu=mui,size=1/sigma))/sqrtVmui
  }
}

psi.sig.ML <- function(r,mu,sigma){
  digamma(r*sqrt(mu*(sigma*mu+1))+mu+1/sigma)-sigma*r*sqrt(mu/(sigma*mu+1))-digamma(1/sigma)-log(sigma*mu+1)
}

psi.sig.ML.mod <- function(j,mui,invsig){
  digamma(j+invsig)-digamma(invsig)-log(mui/invsig+1)-(j-mui)/(mui+invsig)
}

ai.sig.tukey <- function(mui,sigma,c.tukey){
  sqrtVmui <- sqrt(varfunc(mui,sigma))
  invsig <- 1/sigma
  j1 <- max(c(ceiling(mui-c.tukey*sqrtVmui),0))
  j2 <- floor(mui+c.tukey*sqrtVmui)
  if (j1>j2){0}
  else {
    j12 <- j1:j2
    sum((((j12-mui)/(c.tukey*sqrtVmui))^2-1)^2*psi.sig.ML.mod(j=j12,mui=mui,invsig=invsig)*dnbinom(x=j12,mu=mui,size=invsig))
  }
}

sig.rob.tukey <- function(sigma,y,mu,c.tukey,weights){
  r <- (y-mu)/sqrt(varfunc(mu,sigma))
  wi <- tukeypsi(r=r,c.tukey=c.tukey)/r
  sum(weights * wi * psi.sig.ML(r = r, mu = mu, sigma = sigma)) - sum(weights) * ai.sig.tukey(mu, sigma, c.tukey)
}

############
# UNIVARIATE
############
fit_rob_nb_univariate <- function(y, weights = rep(1, length(y)), c.tukey.beta=10, c.tukey.sigma=10, minsig=1e-3, maxsig=50, minmu=1e-10, maxmu=1e20, maxit=50, tol=1e-5, maxit.sig=30, tol.sig=1e-6, warn=FALSE) {
  n <- length(y)
  #-------------------------------------------------------------------
  # MLEs of both sigma and beta
  #-------------------------------------------------------------------
  fit <- MASS::glm.nb(formula = y ~ 1, weights = weights)
  sigma <- 1/fit$theta
  mu <- exp(fit$coefficients[[1]])
  eta <- link(mu)
  update.sigma <- T # at least 1 iteration of robust est, worst case = does not move from minsig/maxsig
  #-------------------------------------------------------------------
  # Robust estimates of both sigma and beta
  #-------------------------------------------------------------------
  sigma0 <- sigma+tol+1
  beta11 <- 0
  beta00 <- beta11+tol+1
  it <- 0
  while(abs(sigma-sigma0)>tol | max(abs(beta11-beta00))>tol & it<maxit) {
    sigma0 <- sigma
    beta00 <- beta11
    # estimate sigma given mu
    if (update.sigma) {
      tryit <- try(uniroot(f=sig.rob.tukey,interval=c(minsig,maxsig),tol=tol.sig,maxiter=maxit.sig,mu=mu,y=y,c.tukey=c.tukey.sigma,weights=weights),T)
      if (class(tryit)=='try-error') {
        if (warn){message(paste('warning: robust update of sigma failed at iteration ',it,', returning last value',sep=''))}
        update.sigma <- FALSE
      } else {
        sigma <- tryit$root
        if (sigma>maxsig) sigma <- maxsig
        if (sigma<minsig) sigma <- minsig
      }
    }
    # estimate mu given sigma
    beta1 <- 0
    beta0 <- beta1+tol+1
    it.mu <- 0
    while (max(abs(beta1-beta0)) > tol & it.mu < maxit) {
      beta0 <- beta1
      bi <- weights * E.tukeypsi.2(mui = mu, sigma = sigma, c.tukey = c.tukey.beta) * varfunc(mu = mu, sigma = sigma)^(-3/2) * mu^2
      ei <- (tukeypsi(r=(y - mu)/sqrt(varfunc(mu, sigma)), c.tukey = c.tukey.beta) - E.tukeypsi.1(mui = mu, sigma = sigma, c.tukey = c.tukey.beta))/
        E.tukeypsi.2(mui = mu, sigma = sigma, c.tukey = c.tukey.beta) * varfunc(mu, sigma) * derivlink(mu)
      zi <- eta + ei
      eta <- beta1 <- sum(zi * bi)/sum(bi)
      mu <- invlink(eta)
      if (mu > maxmu) mu <- maxmu
      if (mu < minmu) mu <- minmu
      it.mu <- it.mu + 1
    }
    beta11 <- beta1
    it <- it + 1
  }
  theta.est <- c(mu = exp(beta11), theta = 1/sigma)

  return(theta.est)
}

fit_nb_univariate <- function(y) {
  fit <- MASS::glm.nb(formula = y ~ 1)
  c(mu = exp(fit$coef[[1]]), theta = fit$theta)
}

###############
# ZERO-INFLATED
###############
fit_rob_zinf_nb_univariate <- function(y, c.tukey.beta = 10, c.tukey.sigma = 10, maxit = 100, tol = 1e-6) {
  n <- length(y)
  is_0 <- y == 0
  one_vect <- rep(1, n)

  # pilot estimate: fit hurdle model and back-calculate pi
  fit_hurdle <- pscl::hurdle(formula = y ~ 1, dist = "negbin")
  mu_pilot <- exp(fit_hurdle$coefficients$count[[1]])
  theta_pilot <- fit_hurdle$theta[[1]]
  p_nb_zero <- dnbinom(x = 0, size = theta_pilot, mu = mu_pilot)
  n_nonzero <- sum(y != 0)
  n_nb <- n_nonzero/(1 - p_nb_zero)
  n_zero_component <- n - n_nb
  pi_hat <- n_zero_component/n
}
