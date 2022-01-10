# i-th possible concurrent cause ===============================================

## Correlated Bernoulli distribution -------------------------------------------

### Probability mass function
dcbern <- function(x, p, alpha){

  if (is.vector(x))
    x <- matrix(x, nrow = length(x))

  n <- dim(x)[1]
  d <- dim(x)[2]

  if(alpha < 0 | alpha > 1) return(NaN)

  p <- matrix(p, ncol = d)
  p <- do.call(rbind, replicate(n/dim(p)[1], p, simplify = FALSE))

  pmf <- matrix(0, n, d)
  pmf[which(p < 0 | p > 1, arr.ind = TRUE)] <- NaN

  id <- which(x == 0 | x == 1 & !is.nan(pmf), arr.ind = TRUE)

  pmf[id] <- (1 - p[id]) * stats::dbinom(x[id], 1, p[id] * (1 - alpha)) +
                  p[id]  * stats::dbinom(x, 1, p[id] + alpha - p[id] * alpha)

  if(d == 1L) as.vector(pmf) else pmf
}

### Random generation
rcbern <- function(n, p, alpha){
  W0 <- stats::rbinom(1, 1, p)
  W <- stats::rbinom(n, 1, p)
  V <- stats::rbinom(n, 1, alpha)

  (1 - V) * W + V * W0
}

# cdnbcr model =================================================================

## Probability density function ------------------------------------------------
dcdnbcr <- function(x, theta, phi, p, alpha, mu, sigma, log.p = FALSE){

  ## Destructive weighted Poisson model
  Dwp <- "DNB"
  pgf1D <- get(paste0("pgf1", Dwp))

  ## Survival time distribution
  Dist <- "rweibull"
  dDist <- get(paste0("d", Dist))
  pDist <- get(paste0("p", Dist))

  if (is.vector(x))
    x <- matrix(x, nrow = length(x))

  n <- dim(x)[1]
  d <- dim(x)[2]

  if(all(x < 0) | all(theta < 0) | phi < 0 | all(p < 0) | all(p > 1) |
     alpha < 0 | alpha > 1 | all(mu < 0) | sigma < 0) return(NaN)

  ## Parameter indexation
  theta <- matrix(theta, ncol = d)
  p <- matrix(p, ncol = d)
  mu <- matrix(mu, ncol = d)

  theta <- do.call(rbind, replicate(n/dim(theta)[1], theta, simplify = FALSE))
  p <- do.call(rbind, replicate(n/dim(p)[1], p, simplify = FALSE))
  mu <- do.call(rbind, replicate(n/dim(mu)[1], mu, simplify = FALSE))

  if(is.null(theta)) stop("theta and x do not have conforming sizes")
  if(is.null(p)) stop("p and x do not have conforming sizes")
  if(is.null(mu)) stop("mu and x do not have conforming sizes")

  pmf <- matrix(-Inf, n, d)
  pmf[which(x < 0 | theta < 0 | p < 0 | p > 1 | mu < 0, arr.ind = TRUE)] <- NaN

  id <- which(!is.nan(pmf), arr.ind = TRUE)

  ## Loglikelihood
  pmf[id] <- dDist(x[id], mu[id], sigma, log = TRUE) +
    log(pgf1D(1 - pDist(x[id], mu[id], sigma), theta[id], phi, p[id], alpha))

  if(d == 1L) pmf <- as.vector(pmf)
  if(log.p) pmf else exp(pmf)
}

## Cumulative distribution function --------------------------------------------
pcdnbcr <- function(q, theta, phi, p, alpha, mu, sigma, lower.tail = TRUE, log.p = FALSE){

  ## Destructive weighted Poisson model
  Dwp <- "DNB"
  pgfD <- get(paste0("pgf", Dwp))

  ## Survival time distribution
  Dist <- "rweibull"
  pDist <- get(paste0("p", Dist))

  if (is.vector(q))
    q <- matrix(q, nrow = length(q))

  n <- dim(q)[1]
  d <- dim(q)[2]

  if(all(q < 0) | all(theta < 0) | phi < 0 | all(p < 0) | all(p > 1) |
     alpha < 0 | alpha > 1 | all(mu < 0) | sigma < 0) return(NaN)

  ## Parameter indexation
  theta <- matrix(theta, ncol = d)
  p <- matrix(p, ncol = d)
  mu <- matrix(mu, ncol = d)

  theta <- do.call(rbind, replicate(n/dim(theta)[1], theta, simplify = FALSE))
  p <- do.call(rbind, replicate(n/dim(p)[1], p, simplify = FALSE))
  mu <- do.call(rbind, replicate(n/dim(mu)[1], mu, simplify = FALSE))

  if(is.null(theta)) stop("theta and x do not have conforming sizes")
  if(is.null(p)) stop("p and x do not have conforming sizes")
  if(is.null(mu)) stop("mu and x do not have conforming sizes")

  cdf <- matrix(0, n, d)
  cdf[which(theta < 0 | p < 0 | p > 1 | mu < 0, arr.ind = TRUE)] <- NaN

  id <- which(!is.nan(cdf), arr.ind = TRUE)

  ## Distribution function
  cdf[id] <- 1 - pgfD(1 - pDist(q[id], mu[id], sigma), theta[id], phi, p[id], alpha)

  if(!lower.tail) cdf <- 1 - cdf
  if(log.p) cdf <- log(cdf)
  if(d == 1L) cdf <- as.vector(cdf)
  cdf
}

## Random generation -----------------------------------------------------------
rcdnbcr <- function(n, theta, phi, p, alpha, mu, sigma){

  ## Destructive weighted Poisson model
  Dwp <- "DNB"
  rDwp <- get(paste0("r", Dwp))

  ## Survival time distribution
  Dist <- "rweibull"
  rDist <- get(paste0("r", Dist))

  ## Random generation
  m <- rDwp(n, theta, phi)
  xi <- mapply(rcbern, m, p, rep(alpha, n))
  D <- mapply(sum, xi)
  y <- mapply(rDist, D, mu, rep(sigma, n))

  suppressWarnings(mapply(min, y))
}

## Cure rate -------------------------------------------------------------------
cure_rate <- function(theta, phi, p, alpha){

  ## Destructive weighted Poisson model
  Dwp <- "DNB"
  pgfD <- get(paste0("pgf", Dwp))

  if (is.vector(theta))
    theta <- matrix(theta, nrow = length(theta))

  if (is.vector(p))
    p <- matrix(p, nrow = length(p))

  n1 <- dim(theta)[1]
  n2 <- dim(p)[1]

  ## Cure rate
  pgfD(rep(0, max(n1, n2)), theta, phi, p, alpha)

}
