ll <- function(psi, time, delta, Z1 = Z1, Z2 = Z2,
               theta.link = "log", p.link = "logit", alpha)
{

  n <- length(time)

  ## Survival time distribution
  pDist <- prweibull
  dDist <- drweibull

  ## phi and alpha identifiers
  phi_id <- TRUE
  alpha_id <- alpha

  if (!alpha_id) alpha_val <- 0L

  ## Covariate matrices setting
  if (is.null(Z1)) Z1 <- matrix(1, nrow = n)
  if (is.null(Z2)) Z2 <- matrix(1, nrow = n)

  r1 <- NCOL(Z1)
  r2 <- NCOL(Z2)

  # Parameter index ------------------------------------------------------------
  param_id <- list(beta1 = 1:r1,
                   phi = if(!phi_id) NULL else r1 + as.numeric(phi_id),
                   beta2 = 1:r2 + r1 + as.numeric(phi_id),
                   alpha = if(!alpha_id) NULL else  r1 + r2 + as.numeric(phi_id) + as.numeric(alpha_id),
                   mu = r1 + r2 + as.numeric(phi_id) + as.numeric(alpha_id) + 1,
                   sigma = r1 + r2 + as.numeric(phi_id) + as.numeric(alpha_id) + 2)

  beta1 <- psi[param_id$beta1]
  beta2 <- psi[param_id$beta2]

  theta <- stats::make.link(theta.link)$linkinv(Z1%*%beta1)
  p <- stats::make.link(p.link)$linkinv(Z2%*%beta2)

  phi <- if (phi_id) psi[param_id$phi] else NULL
  alpha <- if (alpha_id) psi[param_id$alpha] else alpha_val
  mu <- psi[param_id$mu]
  sigma <- psi[param_id$sigma]

  s <- pDist(time, mu, sigma, lower.tail = FALSE)
  f <- dDist(time, mu, sigma)

  Spop <- 1 - pcdnbcr(time, theta, phi, p, alpha, mu, sigma)
  log_fpop <- dcdnbcr(time, theta, phi, p, alpha, mu, sigma, log.p = TRUE)

  sum(delta * log_fpop + (1 - delta) * log(Spop))

}
