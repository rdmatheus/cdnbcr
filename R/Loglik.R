GM.1 <- function(s,theta,phi) theta*(1+phi*theta*(1-s))^(-phi^(-1)-1)
llike <- function(psi, tobs, delta, Z1 = Z1, Z2 = Z2, Z3 = Z3, alpha = TRUE)
{

  n <- length(tobs)

  ## Destructive weighted Poisson model
  Dwp <- "DNB"

  ## Survival time distribution
  Dist <- "rweibull"
  pDist <- get(paste0("p", Dist), mode = "function", envir = parent.frame())
  dDist <- get(paste0("d", Dist), mode = "function", envir = parent.frame())

  ## phi and alpha identifiers
  phi_id <- get(paste0("extrap", Dwp), mode = "logical", envir = parent.frame())
  alpha_id <- alpha

  ## Covariate matrices setting
  if (is.null(Z1)) Z1 <- matrix(1, nrow = n)
  if (is.null(Z2)) Z2 <- matrix(1, nrow = n)
  if (is.null(Z3)) Z2 <- matrix(1, nrow = n)

  r1 <- NCOL(Z1); r2 <- NCOL(Z2); r3 <- NCOL(Z3)

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


  s <- pDist(tobs, mu, sigma, lower.tail=FALSE)
  f <- dDist(tobs, mu, sigma)
  #Spop <- (1-p)*(1+phi*theta*p*(1-alpha)*(1-s))^(-phi^(-1))+
  #  p*(1+phi*theta*(p+alpha-p*alpha)*(1-s))^(-phi^(-1))
  #fpop <- f*(p*(1-p)*(1-alpha)*GM.1(1-p*(1-alpha)*(1-s),theta,phi)+
  #           p*(p+alpha-p*alpha)*GM.1(1-(p+alpha-p*alpha)*(1-s),theta,phi))
  #-sum(delta*log(fpop)+(1-delta)*log(Spop))
  Spop <- 1 - pcdnbcr(tobs, theta, phi, p, alpha, mu, sigma)
  log_fpop <- dcdnbcr(tobs, theta, phi, p, alpha, mu, sigma, log.p = TRUE)
  -sum(delta * log_fpop + (1 - delta) * log(Spop))
}
