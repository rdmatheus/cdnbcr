# Competing causes count distributions =========================================

## Destructive negative binomial -----------------------------------------------

## "M" random generation
rDNB <- function(n, theta, phi){
  stats::rnbinom(n, size = 1/phi, mu = theta)
}

## Probability generating function of "D"
pgfDNB <- function(s, theta, phi, p, alpha){

  if (is.vector(s))
    s <- matrix(s, nrow = length(s))

  n <- dim(s)[1]
  d <- dim(s)[2]

  if(phi <= 0 | alpha < 0 | alpha > 1) return(NaN)

  theta <- matrix(theta, ncol = d)
      p <- matrix(p, ncol = d)

  theta <- do.call(rbind, replicate(n/dim(theta)[1], theta, simplify = FALSE))
      p <- do.call(rbind, replicate(n/dim(p)[1], p, simplify = FALSE))

 pgf <- matrix(NaN, n, d)
  id <- which(theta > 0 & 0 <= p & p <= 1, arr.ind = TRUE)

 pgf[id] <- (1 - p[id]) * ((1 + phi * theta[id] * (1 - s[id]) *  p[id] * (1 - alpha))^(-1/phi)) +
    p[id]  * ((1 + phi * theta[id] * (1 - s[id]) * (p[id] + alpha - p[id] * alpha))^(-1/phi))

 if(d == 1L) as.vector(pgf) else pgf

}

## Derivative of the probability generating function of "D"
pgf1DNB <- function(s, theta, phi, p, alpha){

  if (is.vector(s))
    s <- matrix(s, nrow = length(s))

  n <- dim(s)[1]
  d <- dim(s)[2]

  if(phi <= 0 | alpha < 0 | alpha > 1) return(NaN)

  theta <- matrix(theta, ncol = d)
  p <- matrix(p, ncol = d)

  theta <- do.call(rbind, replicate(n/dim(theta)[1], theta, simplify = FALSE))
  p <- do.call(rbind, replicate(n/dim(p)[1], p, simplify = FALSE))

  pgf1 <- matrix(NaN, n, d)
  id <- which(theta > 0 & 0 <= p & p <= 1, arr.ind = TRUE)

  pgf1[id] <- p[id] * (1 - p[id]) * (1 - alpha) * theta[id] * (1 + phi * theta[id] * p[id] * (1 - alpha) * (1 - s[id]) )^(-1 -1/phi) +
    p[id] * (p[id] + alpha - p[id] * alpha) * theta[id] * (1 + phi * theta[id] * (p[id] + alpha - p[id] * alpha) * (1 - s[id]))^(-1 -1/phi)

  if(d == 1L) as.vector(pgf1) else pgf1

}

## Indicates the presence of an extra parameter (phi)
extrapDNB <- TRUE
