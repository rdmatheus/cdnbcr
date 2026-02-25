#' @name cdnbcr
#'
#' @title Correlated Destructive Negative Binomial Cure Rate Model
#'
#' @description Fits the correlated destructive negative binomial cure rate (CDNBCR) model to right-censored
#' survival data using maximum likelihood estimation via the expectation–maximization (EM)
#' algorithm.
#'
#' @param formula A model formula following the syntax of the \code{Formula} package. The left-hand
#'     side must be a \code{\link[survival]{Surv}} object specifying the observed follow-up time and
#'     censoring indicator. The right-hand side may contain two regression components separated by
#'     the \code{|} operator: the first for \code{theta} (expected number of initial competing causes)
#'     and the second for \code{p} (activation probability of an initial competing cause).
#'     If only one component is provided, it is used for \code{theta} and \code{p} is modeled with
#'     an intercept only.
#' @param data A data frame containing the variables used in the model.
#' @param subset An optional vector specifying a subset of observations to be used in the model.
#' @param na.action A function indicating how to handle missing values. See \code{\link[stats]{model.frame}}.
#' @param theta.link,p.link Link functions for the regression submodels of \code{theta} and \code{p}.
#'     Defaults are \code{"log"} for \code{theta.link} and \code{"logit"} for \code{p.link}.
#'     Alternative links can be specified using \code{\link[stats]{make.link}}.
#' @param alpha Logical; if \code{TRUE} (default), the correlation parameter \code{alpha} is estimated.
#'     If \code{FALSE}, \code{alpha} is fixed at zero, yielding the uncorrelated destructive model.
#' @param y Logical; if \code{TRUE} (default), the response vector is returned.
#' @param x Logical; if \code{TRUE} (default), the model matrices are returned. For \code{print()},
#'    \code{x} is a fitted model object of class \code{"cdnbcr"}.
#' @param control A list of control parameters for the EM algorithm, as returned by
#'     \code{\link{control_EM}}.
#' @param ... Additional arguments passed to \code{\link{control_EM}}.
#' @param digits a non-null value for digits specifies the minimum number of significant digits to
#'     be printed in values.
#'
#' @details
#' The CDNBCR model assumes that the number of initial competing causes follows a negative binomial
#' distribution with mean \eqn{\theta > 0} and dispersion parameter \eqn{\phi > 0}. Each initial
#' cause may remain active with probability \eqn{p \in (0, 1)}, and the activation indicators
#' may be correlated through the parameter \eqn{\alpha \in [0, 1]}. The event time is defined
#' as the minimum of the latent failure times associated with the remaining active causes.
#' Individuals with no active causes are considered cured and will never experience the event
#' of interest. Regression structures can be specified for both the expected number of initial
#' competing causes (\eqn{\theta}) and the activation probability (\eqn{p}).
#'
#' The CDNBCR model reduces to the (uncorrelated) destructive negative binomial cure rate model
#' (Rodrigues et al., 2011) when \eqn{\alpha = 0}. In practice, this is obtained by setting
#' \code{alpha = FALSE} in \code{cdnbcr()}.
#'
#' The response must be specified using \code{Surv(time, status)}, where \code{time} is the observed
#' follow-up time and \code{status} is the event indicator (1 = event, 0 = right-censored).
#' Regression structures can be specified separately for \eqn{\theta} and \eqn{p} using the
#' \code{|} operator. For example,
#'
#' \preformatted{
#' Surv(time, status) ~ x1 + x2 | z1 + z2
#' }
#'
#' specifies a model in which \eqn{\theta} depends on covariates \code{x1} and \code{x2}, while
#' \eqn{p} depends on covariates \code{z1} and \code{z2}.
#'
#' @return An object of class \code{"cdnbcr"} with components:
#' \describe{
#'   \item{coefficients}{A list with the estimated regression coefficients for the \code{theta}
#'        and \code{p} submodels.}
#'   \item{phi, alpha, mu, sigma}{Maximum likelihood estimates of the additional model parameters.}
#'   \item{fitted}{Fitted values of \code{theta} and \code{p} for each observation.}
#'   \item{links}{A named list with the link functions used in the \code{theta} and \code{p}
#'        regression submodels.}
#'   \item{vcov}{Estimated covariance matrix associated with the parameter estimates.}
#'   \item{logLik}{Log-likelihood of the fitted model.}
#'   \item{nobs}{Number of observations used in the fit.}
#'   \item{df.null}{Residual degrees of freedom of the null model.}
#'   \item{df.residual}{Residual degrees of freedom of the fitted model.}
#'   \item{convergence}{Convergence code of the EM algorithm (0 = successful, 1 = maximum iterations reached).}
#'   \item{inits}{Initial values used in the EM algorithm.}
#'   \item{control}{Control parameters used in the EM algorithm.}
#'   \item{iterations}{Number of EM iterations.}
#'   \item{latent}{A list with components \code{M}, \code{D}, and \code{Y} containing
#'        posterior expectations of the latent variables used in the EM algorithm. Here,
#'        \code{M} denotes the initial number of competing causes, \code{D} is the number of
#'        remaining active competing causes, and \code{Y} is a latent Bernoulli variable that governs
#'        the activation regime of the destructive mechanism (and induces dependence through
#'        \code{alpha}).}
#'   \item{call}{Matched function call.}
#'   \item{formula}{Model formula.}
#'   \item{terms}{Terms objects for the regression submodels.}
#'   \item{y}{Response vector (if \code{y = TRUE}).}
#'   \item{x}{Model matrices (if \code{x = TRUE}).}
#' }
#'
#' @seealso
#' \code{\link{summary.cdnbcr}} for detailed model summaries,
#' \code{\link{residuals.cdnbcr}} to extract Cox--Snell residuals (Cox and Snell, 1968),
#' \code{\link{plot.cdnbcr}} for diagnostic plots based on Cox--Snell residuals, and
#' \code{\link{predict.cdnbcr}} for prediction under the CDNBCR model.
#' Additional methods for \code{"cdnbcr"} objects are documented in \code{\link{cdnbcr-methods}}.
#'
#' @references
#' Cox, D. R., and Snell, E. J. (1968). A general definition of residuals.
#'     \emph{Journal of the Royal Statistical Society B}, \bold{30}, 248--265.
#'
#' De Medeiros, R. M. R., Bourguignon, M., Gómez, Y. M., and Gallardo, D. I. (2026).
#'     A correlated approach to cancer cell counting in cure rate models.
#'
#' Rodrigues, J., De Castro, M., Balakrishnan, N., and Cancho, V. G. (2011).
#'     Destructive weighted Poisson cure rate models. \emph{Lifetime Data Analysis}, \bold{17}, 333--346
#'
#' @examples
#' \donttest{## Loading survival package
#' library(survival)
#'
#' ## Dataset: e1690 (see ?e1690)
#' head(e1690)
#'
#' ## Correlated destructive fit
#' fit <- cdnbcr(formula = Surv(time, status) ~ nodeII + nodeIII + nodeIV - 1 |
#'                 sex + trt + thickness + age, data = e1690)
#' fit
#'
#' ## Uncorrelated destructive fit
#' fit0 <- cdnbcr(formula = Surv(time, status) ~ nodeII + nodeIII + nodeIV - 1 |
#'                  sex + trt + thickness + age, data = e1690, alpha = FALSE)
#' fit0}
#' @export
#' @author Diego I. Gallardo \email{diego.gallardo.mateluna@gmail.com}
#' @author Rodrigo M. R. de Medeiros \email{rodrigo.matheus@ufrn.br}
#'
cdnbcr <- function(formula, data, subset, na.action, theta.link = "log", p.link = "logit",
                   alpha = TRUE, control = control_EM(...), y = TRUE, x = TRUE, ...)
{

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
  formula <- Formula::as.Formula(formula)

  if (length(formula)[2L] == 1L) {
    formula <- Formula::as.Formula(formula(formula), ~ 1)
  }

  if (length(formula)[2L] > 2L) {
    formula <- Formula::Formula(formula(formula, rhs = 1:2))
    warning("formula must not have more than two RHS parts")
  }

  mf$formula <- formula

  ## Evaluate model.frame
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  ## Extract response and model matrices
  Y <- stats::model.response(mf)
  Z1 <- stats::model.matrix(formula, mf, rhs = 1L)
  Z2 <- stats::model.matrix(formula, mf, rhs = 2L)

  ## Some conditions
  if (length(Y) < 1)
    stop("empty model")
  if (!inherits(Y, "Surv"))
    stop("response must be a survival object")

  ## Response and censure indicator
  time <- Y[, 1]
  delta <- Y[, 2]

  ## Sample settings
  n <- length(time)
  r1 <- NCOL(Z1)
  r2 <- NCOL(Z2)

  # phi and alpha identifiers
  phi_id <- TRUE
  alpha_id <- alpha

  if (!alpha_id) alpha_val <- 0L

  ## Fit
  opt <- suppressWarnings(EM(time, delta, Z1, Z2, theta.link, p.link, alpha, control))

  ## EM convergence status
  convergence <- opt$convergence

  ## Starting values
  inits <- opt$inits
  names(inits) <- c(colnames(Z1), "phi", colnames(Z2), if(alpha_id) "alpha" else NULL, "mu", "sigma")

  ## Iterations
  iterations <- opt$iterations

  ## Estimates
  estimates <- opt$estimates

  # Parameter index ------------------------------------------------------------
  param_id <- list(beta1 = 1:r1,
                   phi = if(!phi_id) NULL else r1 + as.numeric(phi_id),
                   beta2 = 1:r2 + r1 + as.numeric(phi_id),
                   alpha = if(!alpha_id) NULL else  r1 + r2 + as.numeric(phi_id) + as.numeric(alpha_id),
                   mu = r1 + r2 + as.numeric(phi_id) + as.numeric(alpha_id) + 1,
                   sigma = r1 + r2 + as.numeric(phi_id) + as.numeric(alpha_id) + 2)

  beta1 <- estimates[param_id$beta1]
  beta2 <- estimates[param_id$beta2]

  names(beta1) <- colnames(Z1)
  names(beta2) <- colnames(Z2)

  theta <- stats::make.link(theta.link)$linkinv(Z1%*%beta1)
  p <- stats::make.link(p.link)$linkinv(Z2%*%beta2)

  phi <- estimates[param_id$phi]
  alpha <- if (alpha_id) estimates[param_id$alpha] else alpha_val
  mu <- estimates[param_id$mu]
  sigma <- estimates[param_id$sigma]

  ## Covariance matrix
  vcov <- opt$vcov
  colnames(vcov) <- rownames(vcov) <- c(colnames(Z1), "phi",
                                        colnames(Z2), if (alpha_id) "alpha" else NULL,
                                        "mu", "sigma")

  ## Log-likelihood
  logLik <- ll(estimates, time, delta, Z1, Z2, theta.link, p.link, alpha_id)

  if (!alpha_id) alpha <- NULL

  ## set up return value
  out <- list(coefficients = list(theta = beta1, p = beta2),
              phi = phi,
              alpha = alpha,
              mu = mu,
              sigma = sigma,
              fitted = list(theta = c(theta), p = c(p)),
              links = list(theta.link = theta.link, p.link = p.link),
              vcov = vcov,
              logLik = logLik,
              nobs = n,
              df.null = n - sum(c(1 + alpha_id + phi_id)),
              df.residual = n - r1 - r2 - 2 - sum(c(alpha_id + phi_id)),
              convergence = convergence,
              inits = inits,
              control = control,
              iterations = iterations,
              latent = opt$latent)



  ## Further model information
  mt <- stats::terms(formula, data = data)
  out$call <- cl
  out$formula <- formula
  out$terms <- list(theta = stats::terms(formula, data = data, rhs = 1L),
                    p = stats::terms(formula, data = data, rhs = 2L),
                    full = mt)
  if(y) out$t <- Y
  if(x) out$x <- list(theta = Z1, p = Z2)

  class(out) <- "cdnbcr"
  out
}

# Print
#' @rdname cdnbcr
#' @export
print.cdnbcr <- function(x, digits = getOption("digits"), ...)
{
  r <- x$nobs - x$df.residual

  cat("Call:\n")
  print(x$call)

  cat("\nRegression model for theta (", x$links$theta.link, " link):\n", sep = "")
  print(round(x$coefficients$theta, digits))

  cat("\nRegression model for p (", x$links$p.link, " link):\n", sep = "")
  print(round(x$coefficients$p, digits))

  cat("\nDispersion, dependence, and baseline parameters:\n")
  if (is.null(x$alpha)) {
    par <- cbind(phi = x$phi, mu = x$mu, sigma = x$sigma)
  } else {
    par <- cbind(phi = x$phi, alpha = x$alpha, mu = x$mu, sigma = x$sigma)
  }

  rownames(par) <- ""
  print(round(par, digits))

  cat("\n---\n")
  cat("Log-likelihood: ", round(x$logLik, digits), "\n", sep = "")
  cat("AIC: ", round(2 * (r - x$logLik), digits),
      "   BIC: ", round(log(x$nobs) * r - 2 * x$logLik, digits), "\n", sep = "")
  cat("Number of EM iterations: ", x$EM_iterations, "\n", sep = "")

  invisible(x)
}

