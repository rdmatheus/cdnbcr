# Initial values
inits <- function(tobs, delta, Z1 = NULL, Z2 = NULL, Z3 = NULL,
                  alpha = TRUE, type = c("opt1", "opt2", "opt3"))
{

  # Initial definitions --------------------------------------------------------

  ## Sample size
  n <- length(tobs)

  ## Covariate matrices setting
  if (is.null(Z1)) Z1 <- matrix(1, nrow = n)
  if (is.null(Z2)) Z2 <- matrix(1, nrow = n)
  if (is.null(Z3)) Z3 <- matrix(1, nrow = n)

  r1 <- NCOL(Z1); r2 <- NCOL(Z2); r3 <- NCOL(Z3)

  Dwp <- "DNB"

  phi_id <- get(paste0("extrap", Dwp)) # phi identifier
  alpha_id <- alpha                    # alpha identifier

  ## Type of initial values
  type <- match.arg(type, c("opt1", "opt2", "opt3"))

  passoE<-function(psi, phi = NULL)
  {

    beta1 <- matrix(psi[1:r1], ncol = 1)
    beta2 <- matrix(psi[(r1 + 1):(r1 + r2)], ncol = 1)

    alpha <- psi[(r1 + r2 + 1)]
    nu <- exp(psi[r1 + r2 + 2])

    S <- exp(-exp(alpha)*tobs^nu)
    dF <- 1 - S

    theta <- exp(Z1%*%beta1)
    p <- stats::plogis(Z2%*%beta2)

    M <- c(delta + (1 + phi * delta) * theta * (1 - p * dF)/(1 + phi * theta * p * dF))
    D <- c(delta + (1 + phi * delta) * theta * p * S / (1 + phi * theta * p * dF))

    return(list(M = M, D = D))

  }

  logveroEM.lambdab <- function(lambdab, D)
  {
    alpha <- lambdab[1]
    nu <- exp(lambdab[2])

    S <- exp(-exp(alpha)*tobs^nu)
    logf <- log(nu) + (nu - 1) * log(tobs) + alpha - exp(alpha) * tobs^nu

    -sum((D - delta) * log(S) + delta * logf)
  }


  logveroEM.betab1<-function(beta1, M, phi = NULL)
  {
    beta1 <- matrix(beta1, ncol = 1)
    theta <- exp(Z1%*%beta1)

    -sum(M*(log(phi*theta)-log(1+phi*theta))-(1/phi)*log(1+phi*theta))

  }


  logveroEM.betab2<-function(beta2, M, D)
  {
    beta2<-matrix(beta2,ncol=1)
    p <- stats::plogis(Z2%*%beta2)
    ll <- -sum(D * log(p) + (M - D) * log(1 - p))
    ll
  }


  psi <- c(rep(0, r1 + r2), -log(mean(3.43)), 0)
  if (Dwp != "DLBP") phi <- diff(range(tobs))/2 else phi <- NULL
  iter <- 1; ind <- 0
  while(ind == 0 & iter <= 50)
  {
    aux <- passoE(psi, phi)
    M <- aux$M
    D <- aux$D
    beta1 <- psi[1:r1]
    beta2 <- psi[(r1 + 1):(r1 + r2)]
    lambdab <- psi[(r1 + r2 + 1):(r1 + r2 + 2)]
    maximo1 <- stats::optim(lambdab, logveroEM.lambdab, method = "BFGS", D = D,hessian = FALSE)
    lambdab <- maximo1$par
    maximo2 <- stats::optim(beta1, logveroEM.betab1, method = "BFGS", M = M, phi = phi, hessian = FALSE)
    beta1 <- maximo2$par
    maximo3 <- stats::optim(beta2, logveroEM.betab2,method = "BFGS", M = M, D = D, hessian = FALSE)
    beta2 <- maximo3$par
    psi.aux <- c(beta1, beta2, lambdab)
    ind <- ifelse(max(abs(psi - psi.aux)) < 0.001, 1, 0)
    psi <- psi.aux
    iter <- iter + 1
  }

  switch (type,
          opt1 = {
            beta10 <- psi[1:r1]
            beta20 <- psi[1:r2 + r1]
            if (phi_id){
              theta <- exp(Z1%*%beta10)
              phi0 <- 0.5 * (stats::var(D) + stats::median(theta))/(stats::median(theta)^2)
            } else {
              phi0 <- NULL
            }

          },

          opt2 = {
            aux <- survival::survfit(survival::Surv(tobs, delta) ~ 1, se.fit = FALSE)
            q0 <- min(aux$surv)
            beta10 <- c(log(2) + log(1/q0 - 1), rep(0, r1 - 1))
            beta20 <- rep(0, r2)
            if (phi_id) phi0 <- 1 else phi0 <- NULL
          },

          opt3 = {
            aux <- survival::survfit(survival::Surv(tobs, delta) ~ 1, se.fit = FALSE)
            q0 <- min(aux$surv)
            beta10 <- rep(0, r1)
            beta20 <- c(stats::qlogis(1/q0 - 1), rep(0, r2 - 1))
            if (phi_id) phi0 <- 1 else phi0 <- NULL
          }
  )

  if (alpha_id)  alpha0 <- 0.5 else alpha0 <- NULL
  beta30 <- c(log(exp(-psi[(r1 + r2 + 1)]/psi[(r1 + r2 + 2)]) *
                    gamma(1 + 1/psi[(r1 + r2 + 2)])), rep(0, r3 - 1))
  sigma0 <- exp(psi[r1 + r2 + 2])

  name_phi <- name_alpha <- NULL
  if (alpha_id) name_alpha <- "alpha"
  if (phi_id)  name_phi <- "phi"

  psi0 <- c(beta10, phi0, beta20, alpha0, beta30, sigma0)

  names(psi0) <- c(paste("beta1", 1:r1, sep = ""), name_phi,
                   paste("beta2", 1:r2, sep = ""), name_alpha,
                   paste("beta3", 1:r3, sep = ""), "sigma")

  psi0

}
