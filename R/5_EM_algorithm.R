## Options for estimation ------------------------------------------------------
#' Control Parameters for the EM Algorithm
#'
#' A list of control parameters for the Expectation-Maximization (EM) algorithm used in the
#'     estimation of the Correlated Destructive Negative Binomial Cure Rate Model.
#'
#' @param method Optimization method for the "M" steps. The default method is \code{"BFGS"}.
#' @param maxit Maximum number of EM algorithm iterations. The default is \code{10.000} iterations.
#'     If this limit is reached, a warning message is displayed.
#' @param start An optional vector of initial values for the algorithm. The expected order of the
#'     parameters is \code{(beta1, phi, beta2, alpha, mu, sigma)}, where \code{beta1} and \code{beta2}
#'     are vectors of coefficients associated with the regression structures of
#'     \code{theta} and \code{p}, respectively.
#' @param prec Numeric tolerance for convergence in the EM iterations.
#' @param ... Additional arguments passed to \code{\link[stats]{optim}}.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
#'
#' @return A list containing the specified control parameters.
#' @export
#'
control_EM <- function(method = "BFGS", maxit = 10000, start = NULL,
                       prec = 5e-5, ...)
{
  rval <- list(method = method, maxit = maxit, start = start, prec = prec)

  rval <- c(rval, list(...))

  if (!is.null(rval$fnscale))
    stop("fnscale must not be modified\n")

  rval$fnscale <- -1

  rval
}


# EM algorithm
#' EM Algorithm for Model Estimation
#'
#' An EM algorithm for the maximum likelihood estimation of the parameters of the Correlated
#'     Destructive Negative Binomial Cure Rate Model.
#'
#' @param time A vector of positive response variables (the follow up time).
#' @param delta Status indicator where \code{0} represents right censoring, and \code{1} indicates an
#'     event occurrence at the observed time.
#' @param Z1,Z2 Design matrices for the regression structures of \code{theta} and \code{p}, respectively.
#' @param theta.link,p.link Link functions for the regression models of \code{theta} and \code{p}.
#'     The default values are \code{"log"} for \code{theta.link} and \code{"logit"} for \code{p.link}.
#'     Other options are available via \code{\link[stats]{make.link}}. The choice of link functions
#'     should consider the possible values of \code{theta} (positive values) and \code{p} (unit interval [0,1]).
#' @param alpha If \code{TRUE} (default), \code{alpha} is included in the estimation process to
#'     model correlation between the initial competing causes. If \code{FALSE}, \code{alpha} is
#'     fixed at \code{0}, reducing the model to the uncorrelated case.
#' @param control A list of control parameters for the EM algorithm, specified using \code{\link{control_EM}}.
#' @param ... Additional arguments passed to \code{\link{control_EM}}.
#'
#' @return \code{EM} returns a list with the following components:
#' \describe{
#'   \item{estimates}{ A vector containing the estimated model parameters.}
#'   \item{par}{A named list of parameter estimates, that is, a list containing the elements
#'       \code{beta1}, \code{phi}, \code{beta2}, \code{alpha}, \code{mu}, and \code{sigma}.}
#'   \item{iterations}{ The number of EM iterations performed.}
#'   \item{convergence}{An integer indicating the convergence status, where \code{0} indicates
#'       successful convergence and \code{1} indicates that the maximum number of iterations has been reached.}
#'   \item{vcov}{The estimated covariance matrix of the estimated parameters.}
#'   \item{latent}{Latent variables obtained during the EM algorithm.}
#'  }
#'
#' @author Diego I. Gallardo \email{diego.gallardo.mateluna@gmail.com}
#' @author Rodrigo M. R. de Medeiros \email{rodrigo.matheus@ufrn.br}
EM <- function(time, delta, Z1 = NULL, Z2 = NULL,
               theta.link = "log", p.link = "logit", alpha = TRUE,
               control = control_EM(...), ...)
{

  # Initial definitions --------------------------------------------------------

  ## Sample size
  n <- length(time)

  ## Covariate matrices setting
  if (is.null(Z1)) Z1 <- matrix(1, nrow = n)
  if (is.null(Z2)) Z2 <- matrix(1, nrow = n)

  r1 <- NCOL(Z1); r2 <- NCOL(Z2)

  ## Survival time distribution
  pDist <- prweibull
  dDist <- drweibull

  ## phi and alpha identifiers
  phi_id <- TRUE
  alpha_id <- alpha

  if (!alpha_id) alpha_val <- 0L

  ## Control list --------------------------------------------------------------
  control_aux <- control
  method <- control$method
  maxit <- control$maxit
  prec <- control$prec
  start  <- control$start

  control$method <- control$maxit <- control$prec <-
    control$start <- NULL

  # Parameter index ------------------------------------------------------------
  param_id <- list(beta1 = 1:r1,
                   phi = if(!phi_id) NULL else r1 + as.numeric(phi_id),
                   beta2 = 1:r2 + r1 + as.numeric(phi_id),
                   alpha = if(!alpha_id) NULL else  r1 + r2 + as.numeric(phi_id) + as.numeric(alpha_id),
                   mu = r1 + r2 + as.numeric(phi_id) + as.numeric(alpha_id) + 1,
                   sigma = r1 + r2 + as.numeric(phi_id) + as.numeric(alpha_id) + 2)


  # Functions for the M steps --------------------------------------------------

  ## par = c(beta1, log(phi))
  Q1 <- function(par, latent){

    beta1 <- par[1:r1]
    phi <- exp(par[r1 + 1])

    M <- latent$M
    W <- latent$W
    LW <- latent$LW

    theta <- stats::make.link(theta.link)$linkinv(Z1%*%beta1)

    ll <- (LW - log(phi) - log(theta)) / phi - W / (phi * theta) - lgamma(1/phi)

    if (any(!is.finite(ll)))
      NaN
    else
      sum(ll)

  }

  EPS <- .Machine$double.eps^(1/3)

  ## par = c(beta2, qlogis(alpha))
  Q2 <- function(par, latent){

    beta2 <- par[1:r2]
    alpha <- if (alpha_id) stats::plogis(par[r2 + 1]) else alpha_val

    D <- latent$D
    Y <- latent$Y
    M <- latent$M
    DY <- latent$DY
    MY <- latent$MY

    p <- stats::make.link(p.link)$linkinv(Z2%*%beta2)

    if (alpha > 0){
      ll <- (D - DY) * (log(p) +  log1p(-alpha)) + DY * log(p * (1 - alpha) + alpha) +
        (M - D - MY + DY) * log1p(- p * (1 - alpha)) +
        (MY - DY) * log1p(- p * (1 - alpha) - alpha) +
        Y * log(p) + (1 - Y) * log1p(- p)
    }else{
      ll <- D * log(p) + (M - D) * log1p(- p)
    }


    if (any(!is.finite(ll)))
      NaN
    else
      sum(ll)

  }


  ## par = c(mu, log(sigma))
  Q3 <- function(par, latent){

    mu <- par[1]
    sigma <- exp(par[2])

    D <- latent$D

    ll <- D * pDist(time, mu, sigma, log = TRUE, lower.tail = FALSE)   +
      delta * (dDist(time, mu, sigma, log = TRUE) -
                 pDist(time, mu, sigma, log = TRUE, lower.tail = FALSE))

    if (any(!is.finite(ll)))
      NaN
    else
      sum(ll)
  }

  # Initial values ------------------------------------------------------------
  if (is.null(start)){

    if (p.link != "identity") {
      start <- c(0, rep(0, r1 - 1), 0,
                 rep(0, r2), if(alpha_id) 0 else NULL,
                 log(stats::median(time)), 0)
    } else {
      start <- c(0, rep(0, r1 - 1), 0,
                 rep(0.5, r2), if(alpha_id) 0 else NULL,
                 log(stats::median(time)), 0)
    }

  }

  psik <- start
  it <- 1
  dif <- 1

  while(dif > prec & it <= maxit){

    beta1 <- psik[param_id$beta1]
    beta2 <- psik[param_id$beta2]

    theta <- stats::make.link(theta.link)$linkinv(Z1%*%beta1)
    p <- stats::make.link(p.link)$linkinv(Z2%*%beta2)

    phi <- exp(psik[param_id$phi])
    alpha <- if (alpha_id) stats::plogis(psik[param_id$alpha]) else alpha_val
    mu <- psik[param_id$mu]
    sigma <- exp(psik[param_id$sigma])

    # E-step ------------------------------------------------------------------
    q0 <- p * (1 - alpha)
    q1 <- p * (1 - alpha) + alpha

    c0 <- 1 + phi * theta * q0 * pDist(time, mu, sigma)
    c1 <- 1 + phi * theta * q1 * pDist(time, mu, sigma)

    nu0 <- (1 - p) * (q0^delta) / (c0^(1/phi + delta))
    nu1 <-       p * (q1^delta) / (c1^(1/phi + delta))

    tau <- nu1 / (nu0 + nu1)

    Y <- tau
    W <- (1 + phi * delta) * theta * (tau / c1 + (1 - tau) / c0)
    LW <- digamma(1/phi + delta) + log(phi) + log(theta) - (tau * log(c1) + (1 - tau) * log(c0))
    D <- delta + (1 + phi * delta) * theta * pDist(time, mu, sigma, lower.tail = FALSE) *
      (tau * q1 / c1 + (1 - tau) * q0 / c0)
    DY <- delta * tau + (1 + phi * delta) * theta * pDist(time, mu, sigma, lower.tail = FALSE) *
      (tau * q1 / c1)
    M <- D + (1 + delta * phi) * theta * (tau * (1 - q1) / c1 + (1 - tau) * (1 - q0) / c0)
    MY <- DY + (1 + delta * phi) * theta * tau * (1 - q1) / c1

    latent <- list(Y = Y, M = M, D = D, DY = DY, MY = MY, W = W, LW = LW)

    # M-steps -----------------------------------------------------------------

    ### M1 ---------------------------------------------------------------------
    par1 <- suppressWarnings(stats::optim(c(beta1, log(phi)), Q1,
                                          latent = latent,
                                          method = method,
                                          control = control)$par)

    ### M2 ---------------------------------------------------------------------
    par2_inits <- if(alpha_id) c(beta2, stats::qlogis(alpha)) else beta2
    par2 <- suppressWarnings(stats::optim(par2_inits, Q2,
                                          latent = latent,
                                          method = method,
                                          control = control)$par)

    ### M3 ---------------------------------------------------------------------
    par3 <- suppressWarnings(stats::optim(c(mu, log(sigma)), Q3,
                                          latent = latent,
                                          method = method,
                                          control = control)$par)

    psi <- c(par1, par2, par3)

    dif <- max(abs(psi - psik))
    psik <- psi
    it <- it + 1

  }

  if (it >= maxit){
    convergence <- 1
    warning("Optimization failed to converge: maximum iterations reached.")
  } else {
    convergence <- 0
  }

  beta1 <- psi[param_id$beta1]
  beta2 <- psi[param_id$beta2]

  theta <- stats::make.link(theta.link)$linkinv(Z1%*%beta1)
  p <- stats::make.link(p.link)$linkinv(Z2%*%beta2)

  phi <- if (phi_id) exp(psi[param_id$phi]) else NULL
  alpha <- if (alpha_id) pmin(pmax(stats::plogis(psi[param_id$alpha]), EPS), 1 - EPS) else alpha_val
  mu <- psi[param_id$mu]
  sigma <- exp(psi[param_id$sigma])

  psi <- c(beta1, phi, beta2, if(alpha_id) alpha else NULL, mu, sigma)

  J <- pracma::hessian(ll, x0 = psi, time = time, delta = delta, Z1 = Z1, Z2 = Z2,
                       theta.link = theta.link, p.link = p.link, alpha = alpha_id)

  vcov <- try(solve(-J), silent = TRUE)
  error <- unique(grepl("Error", vcov))
  if (error) vcov <- matrix(NaN, length(psi), length(psi))

  # Latent variables
  q0 <- p * (1 - alpha)
  q1 <- p * (1 - alpha) + alpha

  c0 <- 1 + phi * theta * q0 * pDist(time, mu, sigma)
  c1 <- 1 + phi * theta * q1 * pDist(time, mu, sigma)

  nu0 <- (1 - p) * (q0^delta) / (c0^(1/phi + delta))
  nu1 <-       p * (q1^delta) / (c1^(1/phi + delta))

  tau <- nu1 / (nu0 + nu1)

  Y <- tau
  D <- delta + (1 + phi * delta) * theta *
    pDist(time, mu, sigma, lower.tail = FALSE) *
    (tau * q1 / c1 + (1 - tau) * q0 / c0)
  M <- D + (1 + delta * phi) * theta * (tau * (1 - q1) / c1 + (1 - tau) * (1 - q0) / c0)

  list(estimates =  psi,
       par = list(beta1 = beta1, phi = phi, beta2 = beta2, alpha = alpha, mu = mu, sigma = sigma),
       links = list(theta.link = theta.link, p.link = p.link),
       iterations = it,
       convergence = convergence,
       inits = start,
       vcov = vcov,
       latent = list(Y = c(Y), M = c(M), D = c(D)))

}

