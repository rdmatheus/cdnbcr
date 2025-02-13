# Lifetime distributions =======================================================

## Reparameterized Weibull -----------------------------------------------------

## Probability density function
drweibull <- function(x, mu, sigma, log = FALSE){
  stats::dweibull(x, sigma, exp(log(mu) - lgamma(1/sigma + 1)), log = log)
}

## Cumulative distribution function
prweibull <- function(x, mu, sigma, lower.tail = TRUE, log.p = FALSE){
  stats::pweibull(x, sigma, exp(log(mu) - lgamma(1/sigma + 1)), lower.tail = lower.tail, log.p = log.p)
}

## Random generation
rrweibull <- function(n, mu, sigma){
  stats::rweibull(n, sigma, exp(log(mu) - lgamma(1/sigma + 1)))
}


