#' Correlated Destructive Negative Binomial Cure Rate Model
#'
#' Implements an Expectation-Maximization (EM) algorithm for maximum likelihood estimation in the
#'    Correlated Destructive Negative Binomial Cure Rate Model.
#'
#' @param formula A symbolic description of the model. It must include a \code{\link[survival]{Surv}}
#'     object as the response on the left-hand side of the \code{~} operator. The right-hand side
#'     of the formula can contain covariates for two sub-models: one for \code{theta} (expected number
#'     of initial competing causes) and one for \code{p} (probability of activation of an initial
#'     competing cause). Covariates are separated by the \code{|} symbol. See details below.
#' @param data A data frame containing the variables specified in the formula.
#' @param subset An optional vector specifying a subset of observations to be used in the model.
#' @param na.action A function to handle missing values, see \code{\link[stats]{model.frame}}.
#' @param theta.link,p.link Link functions for the regression submodels of \code{theta} and \code{p}.
#'     Defaults are \code{"log"} for \code{theta.link} and \code{"logit"} for \code{p.link}.
#'     Other options are available via \code{\link[stats]{make.link}}.
#' @param alpha Logical; if \code{TRUE} (default), \code{alpha} is estimated to model correlation
#'     among initial competing causes. If \code{FALSE}, \code{alpha} is fixed at 0, reducing the
#'     model to the uncorrelated case.
#' @param y Logical; if \code{TRUE}, the response vector is returned.
#' @param x Logical; if \code{TRUE}, the model matrices are returned.
#' @param control A list of control parameters for the EM algorithm, specified using \code{\link{control_EM}}.
#' @param ... Additional arguments passed to \code{\link{control_EM}}.
#'
#' @details The \code{formula} argument follows the syntax of the \code{Formula} package
#'     (Zeileis and Croissant, 2010) that allows specifying regression structures for more than one
#'      parameter in the model. For example, suppose the data are formed as follows:
#'
#'  \tabular{cccccccccccc}{
#'  \bold{time}   \tab \tab \bold{status} \tab \tab \bold{x1} \tab \tab \bold{x2} \tab \tab \bold{z1} \tab \tab \bold{z2} \cr
#'  time_1  \tab \tab status_1 \tab \tab x1_1 \tab \tab x2_1 \tab \tab z1_1 \tab \tab z2_1 \cr
#'  time_2  \tab \tab status_2 \tab \tab x1_2 \tab \tab x2_2 \tab \tab z1_2 \tab \tab z2_2 \cr
#'  ...  \tab \tab ... \tab \tab... \tab \tab ... \tab \tab ... \tab \tab ...  \cr
#'  time_n  \tab \tab status_n \tab \tab x1_n \tab \tab x2_n \tab \tab z1_n  \tab \tab z2_n \cr
#'  }
#' where, for each element of the sample, \code{time} is the follow-up time; \code{status}
#' indicates whether the observed time refers to the event of interest (\code{status = 1}) or
#' whether it was right-censored (\code{status = 0}); and \code{x1, x2, z1, z2} are
#' explanatory variables. Then, \code{formula = Surv(time, status) ~ x1 + x2 | z1 + z2}
#' specifies a CDNBCR model with two regression structures, one for \code{theta} in terms of the
#' explanatory variables \code{x1} and \code{x2} and another for \code{p} in terms of the explanatory
#' variables \code{z1} and \code{z2}.
#'
#'
#' @return The function \code{cdnbcr} returns an object of class \code{"cdnbcr"}, which consists of a list
#'     with the following components:
#' \describe{
#'   \item{coefficients}{A list containing the vectors of the regression coefficients \code{beta1}
#'        and \code{beta2}, associated with the regression structures of \code{theta} and \code{p},
#'         respectively.
#'   }
#'   \item{phi, alpha, mu, sigma}{Estimates of the additional parameters of the model.}
#'   \item{fitted}{Fitted values for \code{theta} and \code{p}}
#'   \item{residuals}{Vector of Cox-Snell residuals.}
#'   \item{vcov}{The estimated covariance matrix of the estimated parameters.}
#'   \item{logLik}{Log-likelihood of the fitted model.}
#'   \item{nobs}{Number of observations.}
#'   \item{df.null}{Residual degrees of freedom in the null model.}
#'   \item{df.residual}{Residual degrees of freedom in the fitted model.}
#'   \item{convergence}{Integer code indicating convergence status (0 = successful, 1 = max iterations reached).}
#'   \item{inits}{Initial values used in the EM algorithm.}
#'   \item{control}{List of control parameters used in the EM algorithm.}
#'   \item{EM_iterations}{Number of EM iterations performed.}
#'   \item{call}{Function call.}
#'   \item{formula}{Model formula used.}
#'   \item{terms}{List of terms objects for the two regression frameworks.}
#'   \item{y}{Response vector, if \code{y = TRUE}.}
#'   \item{x}{Model matrices, if \code{x = TRUE}.}
#'  }
#'
#' @examples
#'  \dontrun{
#'   ## Loading survival package
#'   library(survival)
#'
#'   ## Data (see ?e1690)
#'   head(e1690)
#'
#'   ## Correlated destructive fit
#'   fit <- cdnbcr(formula = Surv(time, status) ~ nodeII + nodeIII + nodeIV - 1 |
#'                 sex + trt + thickness + age, data = e1690)
#'
#'   summary(fit)
#'
#'   ## Cox-Snell residuals
#'   par(mfrow = c(1, 2))
#'   plot(fit, ask = FALSE)
#'   par(mfrow = c(1, 1))
#'
#'   ## Latent variables
#'   plot(fit$latent$M, ylab = "Initial competing causes", pch = 16, cex = 0.8)
#'   plot(fit$latent$D, ylab = "Remaining competing causes", pch = 16, cex = 0.8)
#'
#'   ## Uncorrelated destructive fit (alpha = 0)
#'   fit0 <- cdnbcr(formula = Surv(time, status) ~ nodeII + nodeIII + nodeIV - 1 |
#'                  sex + trt + thickness + age, alpha = FALSE, data = e1690)
#'
#'   summary(fit0)
#'  }
#' @export
#' @author Diego I. Gallardo \email{diego.gallardo.mateluna@gmail.com}
#' @author Rodrigo M. R. de Medeiros \email{rodrigo.matheus@ufrn.br}
#'
cdnbcr <- function(formula, data, subset, na.action, theta.link = "log", p.link = "logit",
                   alpha = TRUE, control = control_EM(...), y = FALSE, x = FALSE, ...)
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

  ## Cox-Snell residuals
  Spop <- pcdnbcr(time, theta, phi, p, alpha, mu, sigma, lower.tail = FALSE)
  residuals <- -log(Spop)

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
              residuals = residuals,
              vcov = vcov,
              logLik = logLik,
              nobs = n,
              df.null = n - sum(c(1 + alpha_id + phi_id)),
              df.residual = n - r1 - r2 - 2 - sum(c(alpha_id + phi_id)),
              convergence = convergence,
              inits = inits,
              control = control,
              EM_iterations = opt$iterations,
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
