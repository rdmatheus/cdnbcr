deriv.Q.psi=function(psi, psik, tobs, delta, Z1 = NULL, Z2 = NULL, Z3 = NULL,
                     control = control_EM(...), ...)
{
  n <- length(tobs)

  ## Destructive weighted Poisson model
  Dwp <- "DNB"

  ## Survival time distribution
  Dist <- "rweibull"
  pDist <- get(paste0("p", Dist), mode = "function", envir = parent.frame())

  ## phi identifier
  phi_id <- get(paste0("extrap", Dwp), mode = "logical", envir = parent.frame())

  ## Control list --------------------------------------------------------------
  method <- control$method
  maxit_EM <- control$maxit_EM
  prec <- control$prec
  start  <- control$start
  start_type <- control$start_type
  alpha_id <- control$alpha  # alpha identifier

  control$method <- control$hessian <- control$maxit_EM <- control$prec <-
    control$start <- control$start_type <- control$alpha <- NULL


  ## Covariate matrices setting
  if (is.null(Z1)) Z1 <- matrix(1, nrow = n)
  if (is.null(Z2)) Z2 <- matrix(1, nrow = n)
  if (is.null(Z3)) Z2 <- matrix(1, nrow = n)

  r1 <- NCOL(Z1); r2 <- NCOL(Z2); r3 <- NCOL(Z3)

  param_id <- list(beta1 = 1:r1,
                   phi = r1 + as.numeric(phi_id),
                   beta2 = 1:r2 + r1 + as.numeric(phi_id),
                   alpha = r1 + r2 + as.numeric(phi_id) + as.numeric(alpha_id),
                   beta3 = 1:r3 + r1 + r2 + as.numeric(phi_id) + as.numeric(alpha_id),
                   sigma = r1 + r2 + r3 + as.numeric(phi_id) + as.numeric(alpha_id) + 1)

  # Functions for E step -------------------------------------------------------
  nu <- function(u, psi){

    beta1 <- psi[param_id$beta1]
    beta2 <- psi[param_id$beta2]
    beta3 <- psi[param_id$beta3]

    theta <- exp(Z1%*%beta1)
    p <- stats::plogis(Z2%*%beta2)
    mu <- exp(Z3%*%beta3)

    alpha <- 0
    sigma <- psi[param_id$sigma]

    if (phi_id){
      phi <- psi[param_id$phi]
    }

    if (alpha_id){
      alpha <- psi[param_id$alpha]
    }

    q <- p * (1 - alpha) + u * alpha

    # Only for the DNB
    (p^u) * ((1 - p)^(1 - u)) * (q^delta) *
      ((1 + phi * theta * q * pDist(tobs, mu, sigma))^(-1/phi-delta))

  }

  a <- function(Y, psi){

    beta1 <- psi[param_id$beta1]
    beta2 <- psi[param_id$beta2]
    beta3 <- psi[param_id$beta3]

    theta <- exp(Z1%*%beta1)
    p <- stats::plogis(Z2%*%beta2)
    mu <- exp(Z3%*%beta3)

    alpha <- 0
    sigma <- psi[param_id$sigma]

    if (phi_id){
      phi <- psi[param_id$phi]
    }

    if (alpha_id){
      alpha <- psi[param_id$alpha]
    }

    q <- p * (1 - alpha) + Y * alpha
    qF <- 1 - q * pDist(tobs, mu, sigma)

    (1 + phi * delta) * theta * qF /
      (1 + phi * theta * q * pDist(tobs, mu, sigma))

  }

  b <- function(Y, psi){

    beta1 <- psi[param_id$beta1]
    beta2 <- psi[param_id$beta2]
    beta3 <- psi[param_id$beta3]

    theta <- exp(Z1%*%beta1)
    p <- stats::plogis(Z2%*%beta2)
    mu <- exp(Z3%*%beta3)

    alpha <- 0
    sigma <- psi[param_id$sigma]

    if (phi_id){
      phi <- psi[param_id$phi]
    }

    if (alpha_id){
      alpha <- psi[param_id$alpha]
    }

    q <- p * (1 - alpha) + Y * alpha
    qF <- 1 - q * pDist(tobs, mu, sigma)

    a(Y, psi) * q * pDist(tobs, mu, sigma, lower.tail = FALSE) / qF
  }

  cf <- function(Y, psi){

    beta1 <- psi[param_id$beta1]
    beta2 <- psi[param_id$beta2]
    beta3 <- psi[param_id$beta3]

    theta <- exp(Z1%*%beta1)
    p <- stats::plogis(Z2%*%beta2)
    mu <- exp(Z3%*%beta3)

    alpha <- 0
    sigma <- psi[param_id$sigma]

    if (phi_id){
      phi <- psi[param_id$phi]
    }

    if (alpha_id){
      alpha <- psi[param_id$alpha]
    }

    q <- p * (1 - alpha) + Y * alpha
    qF <- 1 - q * pDist(tobs, mu, sigma)

    1/(1 + phi * theta * q * pDist(tobs, mu, sigma))
  }

  beta1 <- psi[param_id$beta1]
  beta2 <- psi[param_id$beta2]
  beta3 <- psi[param_id$beta3]

  theta <- exp(Z1%*%beta1)
  p <- stats::plogis(Z2%*%beta2)
  mu <- exp(Z3%*%beta3)

  alpha <- 0
  sigma <- psi[param_id$sigma]

  if (phi_id){
    phi <- psi[param_id$phi]
  }

  if (alpha_id){
    alpha <- psi[param_id$alpha]
  }


  beta1k <- psik[param_id$beta1]
  beta2k <- psik[param_id$beta2]
  beta3k <- psik[param_id$beta3]

  thetak <- exp(Z1%*%beta1k)
  pk <- stats::plogis(Z2%*%beta2k)
  muk <- exp(Z3%*%beta3k)

  alphak <- 0
  sigmak <- psik[param_id$sigma]

  if (phi_id){
    phik <- psik[param_id$phi]
  }

  if (alpha_id){
    alphak <- psik[param_id$alpha]
  }


  # E-step ------------------------------------------------------------------
  tauk <- nu(1, psik)/(nu(0, psik) + nu(1, psik))
  Y <- tauk
  M <- delta + a(1, psik) * tauk + a(0, psik) * (1 - tauk)
  D <- delta + b(1, psik) * tauk + b(0, psik) * (1 - tauk)
  MY <- delta * tauk + tauk * a(1, psik)
  DY <- delta * tauk + tauk * b(1, psik)

  if (Dwp == "DNB"){
    W <- (1 + phik * delta) * thetak * (tauk * cf(1, psik) + (1 - tauk) * cf(0, psik))
    LW <- digamma(1/phik + delta) + log(phik) + log(thetak) -
      (tauk * (-log(cf(1, psik))) + (1 - tauk) *  (-log(cf(0, psik))))
  } else {
    W <- LW <- NULL
  }

  #latent <- list(Y = Y, M = M, D = D, DY = DY, MY = MY, W = W, LW = LW)
  t <- tobs

  dl.dtheta <- -1/(phi*theta)+W/(phi*theta^2)
  dl.dphi <- (log(phi)+log(theta)+digamma(1/phi)+W/theta-1-LW)/phi^2
  dl.dp <- (D-DY)/p+DY*(1-alpha)/(p*(1-alpha)+alpha)-(M-MY-D+DY)*(1-alpha)/(1-p*(1-alpha))-(MY-DY)/(1-p)+Y/p-(1-Y)/(1-p)
  dl.dalpha <- -(D-DY)/(1-alpha)+DY*(1-p)/(p*(1-alpha)+alpha)+(M-MY-D+DY)*p/(1-p*(1-alpha))-(MY-DY)/(1-alpha)
  dl.dmu <- D*(sigma*(t*gamma(1/sigma+1)/mu)^sigma/mu)-delta*sigma/mu
  dl.dsigma <- (-D*(t*gamma(1/sigma+1)/mu)^sigma+delta)*(log(t)+lgamma(1/sigma+1)-log(mu)-digamma(1/sigma+1)/sigma)+delta/sigma

  Jr <- function(r) rep(1,r)
  r1 <- ncol(Z1); r2 <- ncol(Z2); r3 <- ncol(Z3)
  deriv <- rep(0, length(psi))
  for(i in 1:length(t))
  {
    x.aux <- diag(c(Z1[i,,drop=FALSE],1,Z2[i,,drop=FALSE],1,Z3[i,,drop=FALSE],1))
    C.aux <- matrix(c(theta[i]*dl.dtheta[i]*Jr(r1), dl.dphi[i],
                   p[i]*(1-p[i])*dl.dp[i]*Jr(r2), dl.dalpha[i],
                   mu[i]*dl.dmu[i]*Jr(r3), dl.dsigma[i]), ncol=1)
    deriv <- deriv+x.aux%*%C.aux
  }
  c(deriv)
}

deriv.Q.psik <- function(psik, psi, tobs, delta, Z1 = NULL, Z2 = NULL, Z3 = NULL,
                      control = control_EM(...), ...)
{

  n <- length(tobs)

  ## Destructive weighted Poisson model
  Dwp <- "DNB"

  ## Survival time distribution
  Dist <- "rweibull"
  pDist <- get(paste0("p", Dist), mode = "function", envir = parent.frame())

  ## phi identifier
  phi_id <- get(paste0("extrap", Dwp), mode = "logical", envir = parent.frame())

  ## Control list --------------------------------------------------------------
  method <- control$method
  maxit_EM <- control$maxit_EM
  prec <- control$prec
  start  <- control$start
  start_type <- control$start_type
  alpha_id <- control$alpha  # alpha identifier

  control$method <- control$hessian <- control$maxit_EM <- control$prec <-
    control$start <- control$start_type <- control$alpha <- NULL

  ## Covariate matrices setting
  if (is.null(Z1)) Z1 <- matrix(1, nrow = n)
  if (is.null(Z2)) Z2 <- matrix(1, nrow = n)
  if (is.null(Z3)) Z2 <- matrix(1, nrow = n)

  r1 <- NCOL(Z1); r2 <- NCOL(Z2); r3 <- NCOL(Z3)

  param_id <- list(beta1 = 1:r1,
                   phi = r1 + as.numeric(phi_id),
                   beta2 = 1:r2 + r1 + as.numeric(phi_id),
                   alpha = r1 + r2 + as.numeric(phi_id) + as.numeric(alpha_id),
                   beta3 = 1:r3 + r1 + r2 + as.numeric(phi_id) + as.numeric(alpha_id),
                   sigma = r1 + r2 + r3 + as.numeric(phi_id) + as.numeric(alpha_id) + 1)

  # Functions for E step -------------------------------------------------------
  nu <- function(u, psi){

    beta1 <- psi[param_id$beta1]
    beta2 <- psi[param_id$beta2]
    beta3 <- psi[param_id$beta3]

    theta <- exp(Z1%*%beta1)
    p <- stats::plogis(Z2%*%beta2)
    mu <- exp(Z3%*%beta3)

    alpha <- 0
    sigma <- psi[param_id$sigma]

    if (phi_id){
      phi <- psi[param_id$phi]
    }

    if (alpha_id){
      alpha <- psi[param_id$alpha]
    }

    q <- p * (1 - alpha) + u * alpha

    # Only for the DNB
    (p^u) * ((1 - p)^(1 - u)) * (q^delta) *
      ((1 + phi * theta * q * pDist(tobs, mu, sigma))^(-1/phi-delta))

  }

  a <- function(Y, psi){

    beta1 <- psi[param_id$beta1]
    beta2 <- psi[param_id$beta2]
    beta3 <- psi[param_id$beta3]

    theta <- exp(Z1%*%beta1)
    p <- stats::plogis(Z2%*%beta2)
    mu <- exp(Z3%*%beta3)

    alpha <- 0
    sigma <- psi[param_id$sigma]

    if (phi_id){
      phi <- psi[param_id$phi]
    }

    if (alpha_id){
      alpha <- psi[param_id$alpha]
    }

    q <- p * (1 - alpha) + Y * alpha
    qF <- 1 - q * pDist(tobs, mu, sigma)

    (1 + phi * delta) * theta * qF /
      (1 + phi * theta * q * pDist(tobs, mu, sigma))

  }

  b <- function(Y, psi){

    beta1 <- psi[param_id$beta1]
    beta2 <- psi[param_id$beta2]
    beta3 <- psi[param_id$beta3]

    theta <- exp(Z1%*%beta1)
    p <- stats::plogis(Z2%*%beta2)
    mu <- exp(Z3%*%beta3)

    alpha <- 0
    sigma <- psi[param_id$sigma]

    if (phi_id){
      phi <- psi[param_id$phi]
    }

    if (alpha_id){
      alpha <- psi[param_id$alpha]
    }

    q <- p * (1 - alpha) + Y * alpha
    qF <- 1 - q * pDist(tobs, mu, sigma)

    a(Y, psi) * q * pDist(tobs, mu, sigma, lower.tail = FALSE) / qF
  }

  cf <- function(Y, psi){

    beta1 <- psi[param_id$beta1]
    beta2 <- psi[param_id$beta2]
    beta3 <- psi[param_id$beta3]

    theta <- exp(Z1%*%beta1)
    p <- stats::plogis(Z2%*%beta2)
    mu <- exp(Z3%*%beta3)

    alpha <- 0
    sigma <- psi[param_id$sigma]

    if (phi_id){
      phi <- psi[param_id$phi]
    }

    if (alpha_id){
      alpha <- psi[param_id$alpha]
    }

    q <- p * (1 - alpha) + Y * alpha
    qF <- 1 - q * pDist(tobs, mu, sigma)

    1/(1 + phi * theta * q * pDist(tobs, mu, sigma))
  }

  beta1 <- psi[param_id$beta1]
  beta2 <- psi[param_id$beta2]
  beta3 <- psi[param_id$beta3]

  theta <- exp(Z1%*%beta1)
  p <- stats::plogis(Z2%*%beta2)
  mu <- exp(Z3%*%beta3)

  alpha <- 0
  sigma <- psi[param_id$sigma]

  if (phi_id){
    phi <- psi[param_id$phi]
  }

  if (alpha_id){
    alpha <- psi[param_id$alpha]
  }


  beta1k <- psik[param_id$beta1]
  beta2k <- psik[param_id$beta2]
  beta3k <- psik[param_id$beta3]

  thetak <- exp(Z1%*%beta1k)
  pk <- stats::plogis(Z2%*%beta2k)
  muk <- exp(Z3%*%beta3k)

  alphak <- 0
  sigmak <- psik[param_id$sigma]

  if (phi_id){
    phik <- psik[param_id$phi]
  }

  if (alpha_id){
    alphak <- psik[param_id$alpha]
  }


  # E-step ------------------------------------------------------------------
  tauk <- nu(1, psik)/(nu(0, psik) + nu(1, psik))
  Y <- tauk
  M <- delta + a(1, psik) * tauk + a(0, psik) * (1 - tauk)
  D <- delta + b(1, psik) * tauk + b(0, psik) * (1 - tauk)
  MY <- delta * tauk + tauk * a(1, psik)
  DY <- delta * tauk + tauk * b(1, psik)

  if (Dwp == "DNB"){
    W <- (1 + phik * delta) * thetak * (tauk * cf(1, psik) + (1 - tauk) * cf(0, psik))
    LW <- digamma(1/phik + delta) + log(phik) + log(thetak) -
      (tauk * (-log(cf(1, psik))) + (1 - tauk) *  (-log(cf(0, psik))))
  } else {
    W <- LW <- NULL
  }

  #latent <- list(Y = Y, M = M, D = D, DY = DY, MY = MY, W = W, LW = LW)
  t <- tobs

  dl.dtheta <- -1/(phi*theta)+W/(phi*theta^2)
  dl.dphi <- (log(phi)+log(theta)+digamma(1/phi)+W/theta-1-LW)/phi^2
  dl.dp <- (D-DY)/p+DY*(1-alpha)/(p*(1-alpha)+alpha)-(M-MY-D+DY)*(1-alpha)/(1-p*(1-alpha))-(MY-DY)/(1-p)+Y/p-(1-Y)/(1-p)
  dl.dalpha <- -(D-DY)/(1-alpha)+DY*(1-p)/(p*(1-alpha)+alpha)+(M-MY-D+DY)*p/(1-p*(1-alpha))-(MY-DY)/(1-alpha)
  dl.dmu <- D*(sigma*(t*gamma(1/sigma+1)/mu)^sigma/mu)-delta*sigma/mu
  dl.dsigma <- (-D*(t*gamma(1/sigma+1)/mu)^sigma+delta)*(log(t)+lgamma(1/sigma+1)-log(mu)-digamma(1/sigma+1)/sigma)+delta/sigma

  Jr <- function(r) rep(1,r)
  r1 <- ncol(Z1); r2 <- ncol(Z2); r3 <- ncol(Z3)
  deriv <- rep(0, length(psi))
  for(i in 1:length(t))
  {
    x.aux <- diag(c(Z1[i,,drop=FALSE],1,Z2[i,,drop=FALSE],1,Z3[i,,drop=FALSE],1))
    C.aux <- matrix(c(theta[i]*dl.dtheta[i]*Jr(r1), dl.dphi[i],
                   p[i]*(1-p[i])*dl.dp[i]*Jr(r2), dl.dalpha[i],
                   mu[i]*dl.dmu[i]*Jr(r3), dl.dsigma[i]), ncol=1)
    deriv <- deriv+x.aux%*%C.aux
  }
  c(deriv)
}









