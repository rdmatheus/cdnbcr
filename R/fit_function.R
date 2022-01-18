cdnbcr <- function(formula, data, subset, na.action,
                   se_type = c("oakes", "louis", "hessian", "null"),
                   control = control_EM(...),
                   y = FALSE,
                   x = FALSE,
                   ...)
{

  Dwp <- "DNB"

  ## Call
  cl <- match.call()

  if (missing(formula))
    stop("a formula argument is required")
  if (missing(data))
    data <- environment(formula)

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE

  ## Formula
  oformula <- stats::as.formula(formula)
  formula <- Formula::as.Formula(formula)

  if (length(formula)[2L] == 1L) {
    formula <- Formula::as.Formula(formula(formula), ~ 1, ~1)
  }

  if (length(formula)[2L] == 2L) {
    formula <- Formula::as.Formula(formula(formula), ~1)
  }

  if (length(formula)[2L] > 3L) {
    formula <- Formula::Formula(formula(formula, rhs = 1:3))
    warning("formula must not have more than three RHS parts")
  }

  mf$formula <- formula

  ## Evaluate model.frame
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  ## Extract response and model matrices
  Y <- stats::model.response(mf)
  Z1 <- stats::model.matrix(formula, mf, rhs = 1)
  Z2 <- stats::model.matrix(formula, mf, rhs = 2)
  Z3 <- stats::model.matrix(formula, mf, rhs = 3)

  ## Some conditions
  if (length(Y) < 1)
    stop("empty model")
  if (!inherits(Y, "Surv"))
    stop("response must be a survival object")

  ## Response and censure indicator
  tobs <- Y[, 1]
  delta <- Y[, 2]

  ## Sample settings
  n <- length(tobs)
  r1 <- NCOL(Z1); r2 <- NCOL(Z2); r3 <- NCOL(Z3)

  # phi and alpha identifiers
  phi_id <- get(paste0("extrap", Dwp), mode = "logical", envir = parent.frame())
  alpha_id <- control$alpha

  # SE type
  se_type <- match.arg(se_type, c("oakes", "louis", "hessian", "null"))

  ## Fit
  opt <- suppressWarnings(EM(tobs, delta, Z1, Z2, Z3, se_type, control))

  ## EM convergence status
  convergence <- opt$convergence == 0

  ## Starting values
  inits <- opt$inits

  ## Iterations
  iterations <- opt$iterations

  ## Estimates
  estimates <- opt$estimates

  # Parameter index ------------------------------------------------------------
  param_id <- list(beta1 = 1:r1,
                   phi = r1 + as.numeric(phi_id),
                   beta2 = 1:r2 + r1 + as.numeric(phi_id),
                   alpha = r1 + r2 + as.numeric(phi_id) + as.numeric(alpha_id),
                   beta3 = 1:r3 + r1 + r2 + as.numeric(phi_id) + as.numeric(alpha_id),
                   sigma = r1 + r2 + r3 + as.numeric(phi_id) + as.numeric(alpha_id) + 1)

  beta1 <- estimates[param_id$beta1]
  beta2 <- estimates[param_id$beta2]
  beta3 <- estimates[param_id$beta3]

  theta <- exp(Z1%*%beta1)
  p <- stats::plogis(Z2%*%beta2)
  mu <- exp(Z3%*%beta3)

  alpha <- 0
  sigma <- estimates[param_id$sigma]

  if (phi_id){
    phi <- estimates[param_id$phi]
  }

  if (alpha_id){
    alpha <- estimates[param_id$alpha]
  }

  ## Covariance matrix
  vcov <- opt$vcov

  ## Quantile residuals
  residuals <- stats::qnorm(pcdnbcr(tobs, theta, phi, p, alpha, mu, sigma))

  ## Log-likelihood
  logLik <- -llike(c(beta1, phi, beta2, alpha, beta3, sigma), tobs, delta, Z1, Z2, Z3, alpha_id)

  theta <- exp(Z1%*%beta1)
  p <- stats::plogis(Z2%*%beta2)
  mu <- exp(Z3%*%beta3)

  if (!phi_id){
    phi <- NULL
  }

  if (!alpha_id){
    alpha <- NULL
  }

  ## set up return value
  out <- list(coefficients = list(competing = beta1,
                                  activation = beta2,
                                  mean = beta3),
              nuisance = list(phi = phi,
                              alpha = alpha,
                              sigma = sigma),
              fitted = list(competing = theta,
                            activation = p,
                            mean = mu),
              vcov = vcov,
              residuals = residuals,
              logLik = logLik,
              convergence = convergence,
              inits = inits,
              control = control,
              it = opt$iterations,
              nobs = n,
              df.null = n - 2,
              df.residual = n - r1 - r2 - r3 - sum(c(1 + alpha_id + phi_id)))



  ## Further model information
  out$call <- cl
  out$formula <- formula
  if(y) out$t <- Y
  if(x) out$x <- list(competing = Z1,
                      activation = Z2,
                      mean = Z3)

  class(out) <- "cdnbcr"
  out
}
