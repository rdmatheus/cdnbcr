Louis <- function(psi, tobs, delta,
                  Z1 = NULL, Z2 = NULL, Z3 = NULL,
                  alpha = TRUE){

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

  ## phi and alpha identifiers
  phi_id <- get(paste0("extrap", Dwp), mode = "logical", envir = parent.frame())
  alpha_id <- alpha

  # Parameter index ------------------------------------------------------------
  param_id <- list(beta1 = 1:r1,
                   phi = r1 + as.numeric(phi_id),
                   beta2 = 1:r2 + r1 + as.numeric(phi_id),
                   alpha = r1 + r2 + as.numeric(phi_id) + as.numeric(alpha_id),
                   beta3 = 1:r3 + r1 + r2 + as.numeric(phi_id) + as.numeric(alpha_id),
                   sigma = r1 + r2 + r3 + as.numeric(phi_id) + as.numeric(alpha_id) + 1)

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

  q0 <- p * (1 - alpha)
  q1 <- p * (1 - alpha) +  alpha
  qF0 <- 1 - q0 * pDist(tobs, mu, sigma)
  qF1 <- 1 - q1 * pDist(tobs, mu, sigma)
  S <- pDist(tobs, mu, sigma, lower.tail = FALSE)
  tau <- nu(1, psi)/(nu(0, psi) + nu(1, psi))

  # Latent variables -----------------------------------------------------------
  Y <- tau

  M <- delta + a(1, psi) * tau + a(0, psi) * (1 - tau)
  M2 <- 2 * delta * M - delta + (1 + phi * delta) * theta *
    (     tau  * qF1 * cf(1, psi) +
            (1 - tau) * qF0 * cf(0, psi) + (1 + 1/phi + delta) * phi * theta *
            (     tau  * (qF1^2) * cf(1, psi)^2 +
                    (1 - tau) * (qF0^2) * cf(0, psi)^2))

  D <- delta + b(1, psi) * tau + b(0, psi) * (1 - tau)
  D2 <- 2 * delta * D - delta + (1 + phi * delta) * theta * S *
    (tau * q1 * cf(1, psi) + (1 - tau) * q0 * cf(0, psi) +
       S * phi * theta * (1 + 1/phi + delta) *
       tau * (q1^2) * cf(1, psi)^2 + (1 - tau) * (q0^2) * cf(0, psi)^2)

  MY <- delta * tau + tau * a(1, psi)
  M2Y <- 2 * delta * MY - delta * tau + (1 + phi * delta) * theta *
    (tau * qF1 * cf(1, psi) +
       (1 + 1/phi + delta) * phi * theta * tau * (qF1^2) * cf(1, psi)^2)

  DY <- delta * tau + tau * b(1, psi)
  D2Y <- 2 * delta * DY - delta * tau + (1 + phi * delta) * theta * S *
    (tau * q1 * cf(1, psi) +
       S * phi * theta * (1 + 1/phi + delta) * tau * (q1^2) * cf(1, psi)^2)

  MD <- delta * M + delta * D - delta +
    S * (1 + phi * delta) * theta * (     tau  * q1 * cf(1, psi) +
                                            (1 - tau) * q0 + cf(0, psi)) +
    S * (1/phi + delta) * (1 + 1/phi + delta) * ((phi * theta)^2) *
    (     tau  * q1 * qF1 * cf(1, psi)^2 +
            (1 - tau) * q0 * qF0 * cf(0, psi)^2)
  MDY <- delta * MY + delta * DY - delta * tau +
    S * (1 + phi * delta) * theta * tau  * q1 * cf(1, psi)  +
    S * (1/phi + delta) * (1 + 1/phi + delta) * ((phi * theta)^2) *
    tau * q1 * qF1 * cf(1, psi)^2

  if (Dwp == "DNB"){

    ## W
    W <- (1 + phi * delta) * theta * (tau * cf(1, psi) + (1 - tau) * cf(0, psi))
    W2 <- (1/phi + delta) * ((phi * theta)^2) * (1 + 1/phi + delta) *
      (tau * cf(1, psi)^2  + (1 - tau) * cf(0, psi)^2)
    LW <- digamma(1/phi + delta) + log(phi) + log(theta) -
      (tau * (-log(cf(1, psi))) + (1 - tau) *  (-log(cf(0, psi))))
    LW2 <- trigamma(1/phi + delta) + digamma(1/phi + delta) +
      2 * digamma(1/phi + delta) * (     tau  * log(cf(1, psi) * phi * theta) +
                                           (1 - tau) * log(cf(0, psi) * phi * theta)) +
      tau * (log(cf(1, psi) * phi * theta)^2) + (1 - tau) * (log(cf(0, psi) * phi * theta)^2)
    WLW <- (1 + phi * delta) * theta *
      (digamma(1/phi + delta + 1) * (tau * cf(1, psi)  + (1 - tau) * cf(0, psi)) +
         (     tau  * cf(1, psi) * log(phi * theta * cf(1, psi)) +
                 (1 - tau) * cf(0, psi) * log(phi * theta * cf(0, psi))))

    ## Y and W
    YW <- (1 + phi * delta) * theta * tau * cf(1, psi)
    YLW <- digamma(1/phi + delta) * tau + tau * log(cf(1, psi) * phi * theta)

    ## M and W
    MW <- delta * W + (1/phi + delta) * ((phi * theta)^2) * (1 + 1/phi + delta) *
      (tau * qF1 * (cf(1, psi)^2) + (1 - tau) * qF0 * (cf(0, psi)^2))
    MLW <- delta * LW + (1 + phi * delta) * theta *
      (digamma(1/phi + delta + 1) * (     tau  * qF1 * cf(1, psi) +
                                            (1 - tau) * qF0 * cf(0, psi)) +
         tau  * qF1 * cf(1, psi) * log(cf(1, psi) * phi * theta) +
         (1 - tau) * qF0 * cf(0, psi) * log(cf(0, psi) * phi * theta))

    ## M, Y, and W
    MYW <- delta * YW + (1/phi + delta) * ((phi * theta)^2) * (1 + 1/phi + delta) *
      tau * qF1 * cf(1, psi)^2
    MYLW <- delta * YLW + (1 + phi * delta) * theta *
      (digamma(1/phi + delta + 1) * tau * qF1 * cf(1, psi) +
         tau * qF1 * cf(1, psi) * log(cf(1, psi) * phi * theta))

    ## D and W
    DW <- delta * W + S * (1/phi + delta) * ((phi * theta)^2) * (1 + 1/phi + delta) *
      (tau * q1 * (cf(1, psi)^2) + (1 - tau) * q0 * (cf(0, psi)^2))
    DLW <- delta * LW + S * (1 + phi * delta) * theta *
      (digamma(1/phi + delta + 1) * (     tau  * q1 * cf(1, psi) +
                                            (1 - tau) * q0 * cf(0, psi)) +
         (     tau  * q1 * cf(1, psi) * log(cf(1, psi) * phi * theta) +
                 (1 - tau) * q0 * cf(0, psi) * log(cf(0, psi) * phi * theta)))

    ## D, Y, and W
    DYW <- delta * YW + S * (1/phi + delta) * ((phi * theta)^2) * (1 + 1/phi + delta) *
      tau * q1 * (cf(1, psi)^2)
    DYLW <- delta * YLW + S * (1 + phi * delta) * theta *
      (digamma(1/phi + delta + 1) * tau * q1 * cf(1, psi) +
         tau * q1 * cf(1, psi) * log(cf(1, psi) * phi * theta))

  } else {
    W <- LW <- NULL
  }

  # functions g ----------------------------------------------------------------
  g1 <- (log(phi) + log(theta) + digamma(1/phi) - 1)/(phi^2)
  g2 <- (1/p + (1 - alpha)/(1 - p * (1 - alpha)))
  g3 <- (-(1-alpha)/(1 - p * (1 - alpha)))
  g4 <- (1/(p * (1 - p)))
  g5 <- (-1/p +
           (1 - alpha) / (p * (1 - alpha) + alpha) +
           (1 - alpha) / (1 - p * (1 - alpha)) +
           (1 - alpha) / (1 - p * (1 - alpha) - alpha))
  g6 <- (1 - alpha) * (1 / (1 - p * (1 - alpha)) -
                         1 / (1 - p * (1 - alpha) - alpha))
  g7 <- (-1 / (1 - alpha) - p / (1 - p * (1 - alpha)))
  g8 <- (p / (1 - p * (1 - alpha)))
  g9 <- (1 / (1 - alpha) + (1 - p) / (p * (1 - alpha) + alpha) +
           p / (1 - p * (1 - alpha)) + (1 - p) / (1 - p * (1 - alpha) - alpha))
  g10 <- (-p / (1 - p * (1 - alpha)) - (1 - p) / (1 - p * (1 - alpha) - alpha))
  g11 <- (sigma / mu) * ((tobs * gamma(1/sigma + 1) / mu)^sigma)
  g12 <- (- ((tobs * gamma(1/sigma + 1) / mu)^sigma) *
            (log(tobs) + log1p(1/sigma) - log(mu) - digamma(1/sigma + 1) / sigma))
  g13 <- delta * (log(tobs) + lgamma(1/sigma + 1) - log(mu) - digamma(1/sigma + 1) / sigma + 1 / sigma)

  # First-order derivatives ----------------------------------------------------

  ## l1 ------------------------------------------------------------------------
  dtheta <- -1/(phi * theta) + W/(phi * theta^2)
  dphi <-  g1 - LW/(phi^2) + W/(theta * phi^2)

  ## l2 ------------------------------------------------------------------------
  dp <- D * g2 + M * g3 + Y * g4 - 1/(1 - p) + DY * g5 + MY * g6
  dalpha <- D * g7 + M * g8 + DY * g9 + MY * g10

  ## l3 ------------------------------------------------------------------------
  dmu <- D * g11 - delta * sigma / mu
  dsigma <- D * g12 + g13

  # Second-order derivatives ---------------------------------------------------

  ## l1 ------------------------------------------------------------------------
  d2theta <- 1 / (phi * theta^2) - 2 * W / (phi * theta^3)
  d2phi <- (-2 * log(phi) - 2 * log(theta) + 2 * LW + 3 - 2 * W / theta -
              2 * digamma(1 / phi) - trigamma(1 / phi) / phi) / (phi^3)
  dtheta.dphi <- 1 / (theta * phi^2) - W / ((phi^2) * (theta^2))

  ## l2 ------------------------------------------------------------------------
  d2p <- - (D - DY) / (p^2) - DY * ((1 - alpha)^2) / ((p * (1 - alpha) + alpha)^2) -
    (M - D - MY + DY) * ((1 - alpha)^2) / ((1 - p * (1 - alpha))^2) -
    (MY - DY) * ((1 - alpha)^2) / ((1 - p * (1 - alpha) - alpha)^2) -
    Y / (p^2) - (1 - Y) / ((1 - p)^2)
  d2alpha <- - (D - DY) / ((1 - alpha)^2) - DY * ((1 - p)^2) / ((p * (1 - alpha) + alpha)^2) -
    (M - D - MY + DY) * (p^2) / ((1 - p * (1 - alpha))^2) -
    (MY - DY) * ((1 - p)^2) / ((1 - p * (1 - alpha) - alpha)^2)
  dp.dalpha <- - DY / ((p * (1 - alpha) + alpha)^2) +
    (M - D - MY + DY) / ((1 - p * (1 - alpha))^2)

  ## l3 ------------------------------------------------------------------------
  d2mu <- - D * sigma * (sigma + 1) * ((tobs * gamma(1/sigma+1) / mu)^sigma) / (mu ^2) +
    delta * sigma / (mu^2)
  d2sigma <- - D * ((tobs * gamma(1/sigma+1) / mu)^sigma) *
    ((log(tobs) + lgamma(1/sigma+1) - log(mu) - digamma(1/sigma + 1)/sigma)^2) +
    (- D * ((tobs * gamma(1/sigma+1) / mu)^sigma) + delta) *
    (-2 * digamma(1/sigma + 1) / (sigma^2) - trigamma(1/sigma + 1)/sigma) -
    delta / (sigma^2)
  dmu.dsigma <- D * sigma * ((tobs * gamma(1/sigma+1) / mu)^sigma) *
    (log(tobs) + lgamma(1/sigma + 1) - log(mu) - digamma(1/sigma + 1)/sigma) / mu +
    (D * ((tobs * gamma(1/sigma+1) / mu)^sigma) - delta) / mu




  # e functions ----------------------------------------------------------------
  e11 <- (1 - 2 * W / theta + (W/theta)^2) / ((phi * theta)^2)
  e12 <- -g1/(phi * theta) + W * (g1 - 1/(phi^2)) / (phi * theta^2) +
    LW / (phi^3 * theta) - WLW / (phi^3 * theta^2) + W2 / (phi^3 * theta^3)
  e13 <- -(D * g2 + M * g3 + Y * g4 - 1 / (1 - p) + DY * g5 + MY * g6) / (phi * theta) +
    (DW * g2 + MW * g3 + YW * g4 - W / (1 - p) + DYW * g5 + MYW * g6) / (phi * theta^2)
  e14 <- -(D * g7 + M * g8 + DY * g9 + MY * g10) / (phi * theta) +
    (DW * g7 + MW * g8 + DYW * g9 + MYW * g10) / (phi * theta^2)
  e15 <- -(D * g11 - delta * sigma / mu) / (phi * theta) +
    (DW * g11 - delta * W * sigma / mu) / (phi * theta^2)
  e16 <- -(D * g12 + g13) / (phi * theta) + (DW * g12 + g13) / (phi * theta^2)


  e22 <- g1^2 + LW2 / (phi^4) + W2 / (phi^4 * theta^2) - 2 * LW * g1 / (phi^2) +
    W * g1 / (phi^2 * theta) - WLW / (phi^4 * theta)
  e23 <- g1 * (D * g2 + M * g3 + Y * g4 - 1 / (1 - p) + DY * g5 + MY * g6) -
    (DLW * g2 + MLW * g3 + YLW * g4 - LW / (1 - p) + DYLW * g5 + MYLW * g6) / (phi^2) +
    (DW * g2 + MW * g3 + YW * g4 - W / (1 - p) + DYW * g5 + MYW * g6) / (phi^2 * theta)
  e24 <- g1 * (D * g7 + M * g8 + DY * g9 + MY * g10) -
    (DLW * g7 + MLW * g8 + DYLW * g9 + MYLW * g10) / (phi^2) +
    (DW * g7 + MW * g8 + DYW * g9 + MYW * g10) / (phi^2 * theta)
  e25 <- g1 * (D * g11 - delta * sigma / mu) - (DLW * g11 - delta * sigma * LW / mu)/(phi^2) +
    (DW * g11 - delta * sigma * W / mu) / (phi^2 * theta)
  e26 <- g1 * (D * g12 + g13) - (DLW * g12 + LW * g13)/(phi^2) +
    (DW * g12 + W * g13)/(phi^2 * theta)


  e33 <- D2 * (g2^2) + M2 * (g3^2) + Y * g4 * (g4 - 1 / (1 - p)) + 1 / ((1 - p)^2) +
    D2Y * g5 * (g2 + g5) + M2Y * g6 * (g3 + g6) + MD * g2 * g3 +
    DY * (g2 * g4 + g4 * g5 - g5 / (1 - p)) +
    MDY * (g2 * g6 + g3 * g5 + g5 * g6) +
    MY * (g3 * g4 + g4 * g6 - g6 / (1 - p)) - (D * g2 + M * g3) / (1 - p)
  e34 <- D2 * g2 * g7 + M2 * g3 * g8 - (D * g7 + M * g8) / (1 - p) +
    MD * (g2 * g8 + g3 * g7) + D2Y * (g2 * g9 + g5 * g7 + g5 * g9) +
    MDY * (g2 * g10 + g3 * g9 + g5 * g8 + g6 * g7 + g6 * g9 + g5 * g10) +
    M2Y * (g3 * g10 + g6 * g8 + g6 * g10) +
    DY * (g4 * g7 + g4 * g9 - g9 / (1 - p)) +
    MY * (g4 * g8 + g4 * g10 - g10 / (1 - p))
  e35 <- g11 * (D2 * g2 + MD * g3 + D2Y * g5 + MDY * g6) -
    (delta * sigma / mu) * (M * g3 + Y * g4 -  1 / (1 - p) + MY * g6) +
    DY * (g4 * g11 - delta * sigma * g5 / mu) -
    D * (1 / (1 - p) + delta * sigma * g2 / mu)
  e36 <- g12 * (D2 * g2 + MD  * g3 + D2Y * g5 + MDY * g6) -
    (delta * sigma / mu) * (M * g3 + Y * g4 - 1 / (1 - p) + MY * g6) +
    DY * (g4 * g12 + g5 * g13) - D * (1 / (1 - p) - g2 * g13)

  e44 <- D2 * (g7^2) + M2 * (g8^2) + MD * g7 * g8 + M2Y * ((g10^2) + g8 * g10) +
    D2Y * ((g9^2) + g7 * g9) + MDY * (g7 * g10 + g8 * g9 + g9 * g10)
  e45 <- g11 * (D2 * g7 + MD * g8 + D2Y * g9 + MDY * g10) -
    (delta * sigma / mu) * (D * g7 + M * g8 + DY * g9 + MY * g10)
  e46 <- g12 * (D2 * g7 + MD * g8 + D2Y * g9 + MDY * g10) +
    g13 * (D * g7 + M * g8 + DY * g9 + MY * g10)

  e55 <- D2 * (g11^2) - 2 * D * delta * sigma * g11 / mu + delta * (sigma^2) / (mu^2)
  e56 <- D2 * g11 * g12 + D * (g11 * g13 - delta * sigma * g12 / mu) -
    delta * sigma * g13 / mu

  e66 <- D2 * (g12^2) - 2 * D * g12 * g13 + g13^2

  B <- SiSi <- 0
  Si <- list()
  for(i in 1:n){
    Zi <- diag(c(Z1[i, ], 1, Z2[i, ], 1, Z3[i, ], 1))
    Ci <- c(rep(theta[i] * dtheta[i], r1), dphi[i],
            rep(p[i] * (1 - p[i]) * dp[i], r2), dalpha[i],
            rep(mu[i] * dmu[i], r3), dsigma[i])

    A1i <- rbind(cbind(theta[i] * (d2theta[i] + dtheta[i]) * diag(r1),
                       rep(dtheta.dphi[i], r1)),
                 c(rep(dtheta.dphi[i], r1), d2phi[i]))
    A2i <- rbind(cbind(p[i] * (1 - p[i]) * (d2p[i] + (1 - 2 * p[i]) * dp[i]) * diag(r2),
                       rep(p[i] * (1 - p[i]) * dp.dalpha[i], r2)),
                 c(rep(p[i] * (1 - p[i]) * dp.dalpha[i], r2), d2alpha[i]))
    A3i <- rbind(cbind(mu[i] * (mu[i] * d2mu[i] + dmu[i]) * diag(r3),
                       rep(mu[i] * dmu.dsigma[i], r3)),
                 c(rep(mu[i] * dmu.dsigma[i], r3), d2sigma[i]))

    A <- rbind(cbind(A1i, matrix(0, r1 + 1, r2 + r3 + 2)),
               cbind(matrix(0, r2 + 1, r1 + 1), A2i, matrix(0, r2 + 1, r3 + 1)),
               cbind(matrix(0, r3 + 1, r1 + 1), matrix(0, r3 + 1, r2 + 1), A3i))

    Ei <- rbind(cbind(e11[i] * diag(r1), rep(e12[i], r1), e13[i] * matrix(1, r1, r2),
                      rep(e14[i], r1), e15[i] * matrix(1, r1, r3), rep(e16[i], r1)),
                c(rep(e12[i], r1), e22[i], rep(e23[i], r2), e24[i],
                  rep(e25[i], r3), e26[i]),
                cbind(e13[i] * matrix(1, r2, r1), rep(e23[i], r2),
                      e33[i] * diag(1, r2), rep(e34[i], r2),
                      e35[i] * matrix(1, r2, r3), rep(e36[i], r2)),
                c(rep(e14[i], r1), e24[i], rep(e34[i], r2), e44[i],
                  rep(e45[i], r3), e46[i]),
                cbind(e15[i] * matrix(1, r3, r1), rep(e25[i], r3),
                      e35[i] * matrix(1, r3, r2), rep(e45[i], r3),
                      e55[i] * diag(1, r3), rep(e56[i], r3)),
                c(rep(e16[i], r1), e26[i], rep(e36[i], r2), e46[i],
                  rep(e56[i], r3), e66[i]))

    # E[B | Dobs, psi] ---------------------------------------------------------
    B <- B + t(Zi)%*%A%*%Zi

    # E[Si Si' | Dobs, psi] ----------------------------------------------------
    SiSi <- SiSi + t(Zi)%*%Ei%*%Zi

    # E[Si | Dobs, psi] --------------------------------------------------------
    Si[[i]] <- t(Zi)%*%Ci

  }


  SiSj <- 0
  for(i in 1:n){
    for(j in 1:n){
      if (i != j){
        SiSj <- SiSj + Si[[i]]%*%t(Si[[j]])
      }
    }
  }


  B - SiSi - SiSj

}

