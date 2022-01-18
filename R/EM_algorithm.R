## Options for estimation ------------------------------------------------------
control_EM <- function(method = "BFGS", maxit_EM = 1000, hessian = TRUE,
                       alpha = TRUE, start = NULL, start_type = "opt1",
                       prec = 1e-3, ...)
{
  rval <- list(method = method, maxit_EM = maxit_EM, hessian = hessian,
               alpha = alpha, start = start, start_type = start_type,
               prec = prec)

  rval <- c(rval, list(...))

  if (!is.null(rval$fnscale))
    warning("fnscale must not be modified\n")

  rval$fnscale <- -1

  rval
}


# EM algorithm
EM <- function(tobs, delta, Z1 = NULL, Z2 = NULL, Z3 = NULL,
               se_type = c("oakes", "louis", "hessian", "null"),
               control = control_EM(...), ...)
{
  # Initial definitions --------------------------------------------------------

  ## Sample size
  n <- length(tobs)

  ## Covariate matrices setting
  if (is.null(Z1)) Z1 <- matrix(1, nrow = n)
  if (is.null(Z2)) Z2 <- matrix(1, nrow = n)
  if (is.null(Z3)) Z2 <- matrix(1, nrow = n)

  r1 <- NCOL(Z1); r2 <- NCOL(Z2); r3 <- NCOL(Z3)

  ## Destructive weighted Poisson model
  Dwp <- "DNB"

  ## Survival time distribution
  Dist <- "rweibull"
  pDist <- get(paste0("p", Dist), mode = "function", envir = parent.frame())
  dDist <- get(paste0("d", Dist), mode = "function", envir = parent.frame())


  ## phi identifier
  phi_id <- get(paste0("extrap", Dwp), mode = "logical", envir = parent.frame())

  ## Control list --------------------------------------------------------------
  control_aux <- control
  method <- control$method
  maxit_EM <- control$maxit_EM
  hessian <- control$hessian
  prec <- control$prec
  start  <- control$start
  start_type <- control$start_type
  alpha_id <- control$alpha  # alpha identifier

  control$method <- control$hessian <- control$maxit_EM <- control$prec <-
    control$start <- control$start_type <- control$alpha <- control$hessian <- NULL

  # Initial values ------------------------------------------------------------
  if (is.null(start)){
    start <- inits(tobs, delta, Z1, Z2, Z3, alpha_id, start_type)
  }

  # Parameter index ------------------------------------------------------------
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


  # Functions for the M steps --------------------------------------------------

  ## par = c(beta1, phi)
  Q1 <- function(par, latent){

    beta1 <- par[1:r1]
    if (phi_id) phi <- par[r1 + 1] else phi <- NULL

    M <- latent$M
    W <- latent$W
    LW <- latent$LW

    theta <- exp(Z1%*%beta1)

    if (any(!is.finite(theta)) | !is.finite(phi) | any(theta < 0) | phi < 0){
      NaN
    }else{
      sum((LW - log(phi) - log(theta)) / phi -
            W / (phi * theta) - lgamma(1/phi))
    }

  }

  ## par = c(beta2, alpha)
  Q2 <- function(par, latent){

    beta2 <- par[1:r2]

    D <- latent$D
    Y <- latent$Y
    M <- latent$M
    DY <- latent$DY
    MY <- latent$MY

    p <- stats::plogis(Z2%*%beta2)
    alpha <- 0

    if (alpha_id){
      alpha <- par[r2 + 1]
    }

    if (any(!is.finite(p)) | !is.finite(alpha) | any(p < 0) | any(p > 1) |
        alpha < 0 | alpha > 1){
      NaN
    }else{

      if (alpha != 0){
        sum((D - DY) * (log(p) +  log1p(-alpha)) + DY * log(p * (1 - alpha) + alpha) +
              (M - D - MY + DY) * log1p(- p * (1 - alpha)) +
              (MY - DY) * log1p(- p * (1 - alpha) - alpha) +
              Y * log(p) + (1 - Y) * log1p(- p))
      }else{
        sum(D * log(p) + (M - D) * log1p(- p))
      }

    }
  }

  ## par = c(beta3, sigma)
  Q3 <- function(par, latent){

    beta3 <- par[1:r3]
    sigma <- par[r3 + 1]

    D <- latent$D

    mu <- exp(Z3%*%beta3)

    if (any(!is.finite(mu)) | !is.finite(sigma) | any(mu < 0) | sigma < 0){
      NaN
    }else{
      ll <- D * pDist(tobs, mu, sigma, log = TRUE, lower.tail = FALSE)   +
        delta * (dDist(tobs, mu, sigma, log = TRUE) -
                   pDist(tobs, mu, sigma, log = TRUE, lower.tail = FALSE))

      if (any(!is.finite(ll)))
        NaN
      else
        sum(ll)
    }
  }


  psi0 <- start
  it <- dif <- 1
  while(dif > prec & it <= maxit_EM){

    beta1 <- psi0[param_id$beta1]
    beta2 <- psi0[param_id$beta2]
    beta3 <- psi0[param_id$beta3]

    theta <- exp(Z1%*%beta1)
    p <- stats::plogis(Z2%*%beta2)
    mu <- exp(Z3%*%beta3)

    alpha <- 0
    sigma <- psi0[param_id$sigma]

    if (phi_id){
      phi <- psi0[param_id$phi]
    } else{
      phi <- NULL
    }

    if (alpha_id){
      alpha <- psi0[param_id$alpha]
    }

    psi0_aux <- c(beta1, phi, beta2, alpha, beta3, sigma)
    # E-step ------------------------------------------------------------------
    tau <- nu(1, psi0_aux)/(nu(0, psi0_aux) + nu(1, psi0_aux))
    Y <- tau
    M <- delta + a(1, psi0_aux) * tau + a(0, psi0_aux) * (1 - tau)
    D <- delta + b(1, psi0_aux) * tau + b(0, psi0_aux) * (1 - tau)
    MY <- delta * tau + tau * a(1, psi0_aux)
    DY <- delta * tau + tau * b(1, psi0_aux)

    if (Dwp == "DNB"){
      W <- (1 + phi * delta) * theta * (tau * cf(1, psi0_aux) + (1 - tau) * cf(0, psi0_aux))
      LW <- digamma(1/phi + delta) + log(phi) + log(theta) -
        (tau * (-log(cf(1, psi0_aux))) + (1 - tau) *  (-log(cf(0, psi0_aux))))
    } else {
      W <- LW <- NULL
    }

    latent <- list(Y = Y, M = M, D = D, DY = DY, MY = MY, W = W, LW = LW)

    # M-steps -----------------------------------------------------------------

    ### M1 ---------------------------------------------------------------------
    par1 <- suppressWarnings(stats::optim(c(beta1, phi), Q1,
                                          latent = latent,
                                          method = method,
                                          control = control)$par)

    ### M2 ---------------------------------------------------------------------
    par2_inits <- if(alpha_id) c(beta2, alpha) else beta2
    par2 <- suppressWarnings(stats::optim(par2_inits, Q2,
                                          latent = latent,
                                          method = method,
                                          control = control)$par)

    ### M3 ---------------------------------------------------------------------
    par3 <- suppressWarnings(stats::optim(c(beta3, sigma), Q3,
                                          latent = latent,
                                          method = method,
                                          control = control)$par)

    psi <- c(par1, par2, par3)

    dif <- max(abs(psi0 - psi))

    psi0 <- psi
    it <- it + 1
  }

  # beta1 <- psi0[param_id$beta1]
  # beta2 <- psi0[param_id$beta2]
  # beta3 <- psi0[param_id$beta3]
  #
  # alpha <- 0
  # sigma <- psi0[param_id$sigma]
  #
  # if (phi_id){
  #   phi <- psi0[param_id$phi]
  # } else{
  #   phi <- NULL
  # }
  #
  # if (alpha_id){
  #   alpha <- psi0[param_id$alpha]
  # }
  #
  # psi0_aux <- c(beta1, phi, beta2, alpha, beta3, sigma)

  if (it >= maxit_EM){
    convergence <- 1
    warning(cat("optimization failed to converge: maximum iterations reached\n"))
  } else {
    convergence <- 0
  }

  se_type <- match.arg(se_type, c("oakes", "louis", "hessian", "null"))

  switch (se_type,
    oakes = {
      a1 <- pracma::jacobian(deriv.Q.psi, x0 = psi, psik = psi, tobs = tobs,
                             delta = delta, Z1 = Z1, Z2 = Z2, Z3 = Z3,
                             control = control_aux)
      a2 <- pracma::jacobian(deriv.Q.psik, x0 = psi, psi = psi, tobs = tobs,
                             delta = delta, Z1 = Z1, Z2 = Z2, Z3 = Z3,
                             control = control_aux)

      hessian <- a1 + a2
      vcov <- try(solve(-hessian))
      #se <- try(sqrt(diag(vcov)))
      error <- unique(grepl("Error", vcov))
      if (error) vcov <- matrix(NaN, length(psi), length(psi))
    },

    louis = {
      hessian <- Louis(psi, tobs, delta, Z1, Z2, Z3, alpha_id)
      vcov <- try(solve(-hessian))
      #se <- try(sqrt(diag(vcov)))
      error <- unique(grepl("Error", vcov))
      if (error) vcov <- matrix(NaN, length(psi), length(psi))
    },

    hessian = {
      hessian <- pracma::hessian(llike, x0 = psi, tobs = tobs, delta = delta,
                                 Z1 = Z1, Z2 = Z2, Z3 = Z3)
      vcov <- try(solve(-hessian))
      #se <- try(sqrt(diag(vcov)))
      error <- unique(grepl("Error", vcov))
      if (error) vcov <- matrix(NaN, length(psi), length(psi))
    },

    null = {
      se <- NULL
    }
  )

  #if(hessian){
  #  hessian <- -pracma::hessian(ll, x0 = psi0_aux,
  #                              tobs = tobs, delta = delta, Z1 = Z1, Z2 = Z2, Z3 = Z3,
  #                              alpha = alpha_id)
  #} else{
  #  hessian <- NULL
  #}

  list(estimates =  psi,
       iterations = it,
       convergence = convergence,
       inits = start,
       vcov = vcov)

}

