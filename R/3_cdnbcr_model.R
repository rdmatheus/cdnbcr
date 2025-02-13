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
#' @name cdnbcr_model
#'
#' @title The Correlated Destructive Negative Binomial Cure Rate Model
#'
#' @description Implements functions for the Correlated Destructive Negative Binomial Cure Rate (CDNBCR) model,
#'     including probability density, cumulative distribution, cure rate, and random generation.
#'
#' @param x Vector of positive event times or quantiles.
#' @param q Vector of quantiles.
#' @param n Number of random values to return.
#' @param theta,phi Positive parameters associated with the expected initial number of competing
#'     causes and their dispersion.
#' @param p,alpha Parameters controlling the probability of a cause remaining active. The correlation
#'     among active causes is determined by \code{alpha}. Both \code{p} and
#'     \code{alpha} are restricted to the closed unit interval \code{[0, 1]}.
#'     When \code{alpha = 0}, the activations are assumed to be independent.
#' @param mu,sigma Strictly positive parameters of the Weibull survival distribution assumed for the
#'     time-to-event of active competing causes.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are
#'     \code{P(X <= x)}, otherwise, \code{P(X > x)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \code{p} are given as
#'     \code{log(p)}.
#'
#' @return \code{dcdnbcr} returns the probability density function, \code{pcdnbcr} gives the
#'     distribution function, \code{cure_rate} returns the cure rate of the model, and \code{rcdnbcr}
#'     generates random observations.
#'
#'
#' @details
#' The CDNBCR model describes the time until an event occurs while accounting for a cure fraction,
#'     assuming multiple initial competing causes with potential correlation. The model is
#'     destructive in the sense that some initial competing causes contributing to the occurrence
#'     of the event can be eliminated. Thus, the model is particularly useful in applications
#'     where interventions (e.g., treatments) may eliminate or inactivate some competing causes.
#'
#'     We assume that the number of initial competing causes is described by a negative binomial
#'     random variable with mean \code{theta} and variance \code{theta * (1 + phi * theta)}. Each of
#'     the initial competing causes may remain contributing to the occurrence of the event with
#'     probability \code{p}. The model introduces dependence on the destruction or not of these
#'     initial causes, controlled by the parameter \code{alpha}. Finally, let \eqn{Y_j} be the
#'     survival time due to the \eqn{j}th remaining competing cause, \eqn{j = 1, \ldots, D},
#'     where \eqn{D} is a latent random variable. We assume that \eqn{Y_j} follows a Weibull distribution
#'     with mean \code{mu} and variance \code{mu^2 * (Gamma(2/sigma + 1) / Gamma(1/sigma + 1)^2 - 1)}.
#'     The observed uncensored survival time is given by \eqn{T = \textrm{min}(Y_1, \ldots, Y_D)}.
#'
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
#' @author Diego I. Gallardo <\email{diego.gallardo.mateluna@gmail.com}>
#'
#' @examples
#'
#' # Parameters
#' theta <- 1.5
#' phi <- 1.2
#' p <- 0.6
#' alpha <- 0.3
#' mu <- 5
#' sigma <- 1.5
#'
#' n <- 1000
#' tobs <- rcdnbcr(n, theta, phi, p, alpha, mu, sigma)
#'
#' # Censoring
#' Cens <- 10
#' delta <- ifelse(tobs <= Cens, 1, 0)
#' tobs[delta == 0] <- Cens
#'
#' op <- par()$mfrow
#' par(mfrow = c(1, 2))
#' # Histogram
#' hist(tobs, prob = TRUE, main = " ", xlab = expression(T[obs]))
#' curve(dcdnbcr(x, theta, phi, p, alpha, mu, sigma), add = TRUE, col = 2, lwd = 2)
#' points(Cens, mean(1 - delta), pch = 16, col = 2, lwd = 2)
#'
#' # Survival function
#' plot(survival::survfit(survival::Surv(tobs, delta) ~ 1, se.fit = FALSE), xlab = expression(T[obs]),
#'      ylab = "Survival function")
#' curve(1 - pcdnbcr(x, theta, phi, p, alpha, mu, sigma), add = TRUE, col = 2, lwd = 2)
#'
#' # Cure rate
#' abline(h = cure_rate(theta, phi, p, alpha), col = 4)
#' legend("bottomleft", c("Empirical", "Theoretical", "Cure rate"),
#'         col = c(1, 2, 4), bty = "n", lty = 1)
#' par(mfrow = op)

## Probability density function ------------------------------------------------
#' @rdname cdnbcr_model
#' @export
dcdnbcr <- function(x, theta, phi, p, alpha, mu, sigma, log.p = FALSE){

  ## Destructive weighted Poisson model
  pgf1D <- pgf1DNB

  ## Survival time distribution
  dDist <- drweibull
  pDist <- prweibull

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

  theta <- do.call(rbind, replicate(n/dim(theta)[1], theta, simplify=FALSE))
  p <- do.call(rbind, replicate(n/dim(p)[1], p, simplify=FALSE))
  mu <- do.call(rbind, replicate(n/dim(mu)[1], mu, simplify=FALSE))

  pmf <- matrix(-Inf, n, d)
  pmf[which(x < 0 | theta < 0 | p < 0 | p > 1 | mu < 0, arr.ind = TRUE)] <- NaN

  id <- which(!is.nan(pmf), arr.ind = TRUE)

  ## Log-likelihood
  pmf[id] <- dDist(x[id], mu[id], sigma, log = TRUE) +
    log(pgf1D(1 - pDist(x[id], mu[id], sigma), theta[id], phi, p[id], alpha))

  if(d == 1L) pmf <- as.vector(pmf)
  if(log.p) pmf else exp(pmf)
}

## Cumulative distribution function --------------------------------------------
#' @rdname cdnbcr_model
#' @export
pcdnbcr <- function(q, theta, phi, p, alpha, mu, sigma, lower.tail = TRUE, log.p = FALSE){

  ## Destructive weighted Poisson model
  pgfD <- pgfDNB

  ## Survival time distribution
  pDist <- prweibull

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

  theta <- do.call(rbind, replicate(n/dim(theta)[1], theta, simplify=FALSE))
  p <- do.call(rbind, replicate(n/dim(p)[1], p, simplify=FALSE))
  mu <- do.call(rbind, replicate(n/dim(mu)[1], mu, simplify=FALSE))

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
#' @rdname cdnbcr_model
#' @export
rcdnbcr <- function(n, theta, phi, p, alpha, mu, sigma){

  ## Destructive weighted Poisson model
  rDwp <- rDNB

  ## Survival time distribution
  rDist <- rrweibull

  ## Random generation
  m <- rDwp(n, theta, phi)
  xi <- mapply(rcbern, m, p, rep(alpha, n))
  D <- mapply(sum, xi)
  y <- mapply(rDist, D, mu, rep(sigma, n))

  suppressWarnings(mapply(min, y))
}

## Cure rate -------------------------------------------------------------------
#' @rdname cdnbcr_model
#' @export
cure_rate <- function(theta, phi, p, alpha){

  ## Destructive weighted Poisson model
  pgfD <- pgfDNB

  if (is.vector(theta))
    theta <- matrix(theta, nrow = length(theta))

  if (is.vector(p))
    p <- matrix(p, nrow = length(p))

  n1 <- dim(theta)[1]
  n2 <- dim(p)[1]

  ## Cure rate
  pgfD(rep(0, max(n1, n2)), theta, phi, p, alpha)

}
