# General methods -----------------------------------------------------------------------------
#' @name cdnbcr-methods
#'
#' @title Extract Information From a CDNBCR Fit
#'
#' @description
#'  Methods for \code{"cdnbcr"} objects.
#'
#' @param object An object of class \code{"cdnbcr"}.
#' @param formula A model \code{\link[Formula]{Formula}} or \code{\link[stats]{terms}} object or an
#'     \code{"cdnbcr"} object.
#' @param k Numeric, the penalty per parameter to be used; the default
#'     \code{k = 2} is the classical AIC. See \code{\link[stats]{AIC}}.
#' @param parm A character indicating which regression structure should be used. It can be
#'     \code{"theta"} for the expected initial competing causes regression structure,
#'     \code{"p"} for the activation probability of an initial competing cause regression
#'     submodel, or \code{"full"} for both regression structures.
#' @param ... Additional arguments passed to or from other methods.
#'
#' @author
#' Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
#'
#' @returns
#' \itemize{
#' \item \code{model.frame} returns a \code{data.frame} containing the variables required
#'     by \code{formula} and any additional arguments provided via \code{...}.
#' \item \code{model.matrix} returns the design matrix used in the regression structure,
#'     as specified by the \code{parm} argument.
#' \item \code{coef} returns a numeric vector of estimated regression coefficients, based
#'     on the \code{parm} argument. If \code{parm = "full"}, it returns a list with the
#'     components \code{"theta"} and \code{"p"}, each containing the corresponding
#'     coefficient estimates.
#' \item \code{vcov} returns the asymptotic covariance matrix of the regression coefficients,
#'     based on the \code{parm} argument.
#' \item \code{residuals} returnas a \code{"\link[survival]{Surv}"} object with the Cox-Snell residuals
#'     (Cox and Snell, 1968). If the model is well fitted to the data, the Cox-Snell residuals are
#'      expected to be distributed as a censored random sample from the exponential distribution
#'       with mean 1.
#' \item \code{logLik} returns the log-likelihood value of the fitted model.
#' \item \code{AIC} returns a numeric value representing the Akaike Information Criterion
#'     (AIC), Bayesian Information Criterion, or another criterion, depending on \code{k}.
#' }
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
#'
#' # Model frame
#' mf <- model.frame(fit)
#' mf
#'
#' # Model matrices
#' model.matrix(fit, parm = "theta")
#' model.matrix(fit, parm = "p")
#'
#' # Coef
#' coef(fit)
#' coef(fit, parm = "theta")
#' coef(fit, parm = "p")
#'
#' # vcov
#' vcov(fit)
#' vcov(fit, parm = "theta")
#' vcov(fit, parm = "p")
#'
#' # residuals
#' residuals(fit)
#'
#' # Log-likelihood value
#' logLik(fit)
#'
#' # AIC and BIC
#' AIC(fit)
#' AIC(fit, k = log(fit$nobs))}
NULL

## Model frame
#' @export
#' @rdname cdnbcr-methods
model.frame.cdnbcr <- function(formula, ...) {
  formula$terms <- formula$terms$full
  formula$call$formula <- formula$formula <- formula(formula$terms)
  NextMethod("model.frame")
}


## Model matrix
#' @export
#' @rdname cdnbcr-methods
model.matrix.cdnbcr <- function(object, parm = c("full", "theta", "p"), ...) {
  parm <- match.arg(parm)
  if (!is.null(object$x)){
    rval <- object$x
  } else {
    rval <- list(theta = stats::model.matrix(object$terms[["theta"]], stats::model.frame(object)),
                 p = stats::model.matrix(object$terms[["p"]], stats::model.frame(object)))
  }

  switch(parm,
         "full" = rval,
         "theta" = rval$theta,
         "p" = rval$p)
}

# Parameter estimates
#' @rdname cdnbcr-methods
#' @export
coef.cdnbcr <- function(object,
                        parm = c("full", "theta", "p"), ...) {

  parm <- match.arg(parm, c("full", "theta", "p"))
  switch(parm,
         "full" = object$coefficients,
         "theta" = object$coefficients$theta,
         "p" = object$coefficients$p)
}

#  Variance-covariance matrix
#' @rdname cdnbcr-methods
#' @export
vcov.cdnbcr <- function(object,
                        parm = c("full", "theta", "p"), ...) {

  parm <- match.arg(parm, c("full", "theta", "p"))
  covm <- object$vcov

  r1 <- length(object$coefficients$theta)
  r2 <- length(object$coefficients$p)
  r <- object$nobs - object$df.residual

  param_id <- list(beta1 = 1:r1,
                   phi = r1 + 1,
                   beta2 = 1:r2 + r1 + 1,
                   alpha = r1 + r2 + 1 + length(object$alpha),
                   mu = r1 + r2 + 2 + length(object$alpha),
                   sigma = r)

  switch(parm,
         "theta" = {
           covm[param_id$beta1, param_id$beta1]
         },
         "p" = {
           covm[param_id$beta2, param_id$beta2]
         },
         "full" = {
           covm
         })

}

# Residuals
#' @rdname cdnbcr-methods
#' @export
residuals.cdnbcr <- function(object, ...){

  ## Model specifications
  mf <- stats::model.frame(object)
  time <- stats::model.response(mf)
  theta <- object$fitted$theta
  phi <- object$phi
  p <- object$fitted$p
  alpha <- if (!is.null(object$alpha)) object$alpha else 0
  mu <- object$mu
  sigma <- object$sigma

  ## Cox-Snell residuals
  res <- -log(pcdnbcr(time[, 1], theta, phi, p, alpha, mu, sigma, lower.tail = FALSE))
  survival::Surv(res, event = time[, 2])

}

# Log-likelihood
#' @rdname cdnbcr-methods
#' @export
logLik.cdnbcr <- function(object, ...) {
  structure(object$logLik,
            df = object$nobs - object$df.residual,
            class = "logLik")
}


# AIC
#' @export
#' @rdname cdnbcr-methods
AIC.cdnbcr <- function(object, ..., k = 2) {

  r <- object$nobs - object$df.residual
  AIC <- - 2 * object$logLik + k * r

  class(AIC) <- "AIC"
  return(AIC)
}

# Summary method ------------------------------------------------------------------------------
#' @name summary.cdnbcr
#'
#' @title Summarizing a CDNBCR Fit
#'
#' @description \code{summary} method for class \code{"cdnbcr"}.
#'
#' @param object An object of class \code{"cdnbcr"}, a result of a call to \code{\link{cdnbcr}}.
#' @param x An object of class \code{"summary.cdnbcr"}, a result of a call to \code{summary.cdnbcr}.
#' @param digits A non-null value for digits specifies the minimum number of significant digits to
#'     be printed in values.
#' @param ... Further arguments passed to or from other methods.
#'
#' @returns The function \code{summary.cdnbcr} returns an object of class \code{"summary.cdnbcr"},
#'   a list with the following components:
#' \describe{
#'   \item{call}{The original function call.}
#'   \item{theta}{Summary of the regression coefficients for the \code{theta} submodel.}
#'   \item{p}{Summary of the regression coefficients for the \code{p} submodel.}
#'   \item{par}{Estimates and standard errors for the additional model parameters: \code{phi}
#'     (dispersion of the initial number of competing causes), \code{alpha} (dependence parameter,
#'     when estimated), and the baseline parameters \code{mu} and \code{sigma} (reparameterized Weibull).}
#'   \item{links}{Named list with the link functions used in the \code{theta} and \code{p}
#'     regression submodels.}
#'   \item{residuals}{A \code{"\link[survival]{Surv}"} object with the Cox-Snell residuals
#'     (Cox and Snell, 1968). If the model is well fitted to the data, the Cox-Snell residuals are
#'      expected to be distributed as a censored random sample from the exponential distribution
#'       with mean 1.}
#'   \item{iterations}{Number of EM iterations.}
#'   \item{logLik}{Log-likelihood of the fitted model.}
#'   \item{AIC, BIC}{Akaike and Bayesian information criteria.}
#' }
#' @export
#'
#' @author Diego I. Gallardo \email{diego.gallardo.mateluna@gmail.com}
#' @author Rodrigo M. R. de Medeiros \email{rodrigo.matheus@ufrn.br}
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
#' out <- summary(fit)
#' out
#'
#' ### Summary table for the regression coefficients
#' out$theta
#' out$p
#'
#' ### Summary table for the additional parameters
#' out$par
#'
#' ### Cox-Snell residuals
#' out$residuals
#' class(out$residuals)
#' plot(out$residuals)
#' curve(pexp(x, 1, lower.tail = FALSE), add = TRUE, col = "blue")
#'
#' ## Uncorrelated destructive fit
#' fit0 <- cdnbcr(formula = Surv(time, status) ~ nodeII + nodeIII + nodeIV - 1 |
#'                  sex + trt + thickness + age, data = e1690, alpha = FALSE)
#' summary(fit0)}
summary.cdnbcr <- function(object, ...)
{

  r1 <- length(object$coefficients$theta)
  r2 <- length(object$coefficients$p)
  r <- object$nobs - object$df.residual

  param_id <- list(beta1 = 1:r1,
                   phi = r1 + 1,
                   beta2 = 1:r2 + r1 + 1,
                   alpha = r1 + r2 + 1 + length(object$alpha),
                   mu = r1 + r2 + 2 + length(object$alpha),
                   sigma = r)

  # Summary for beta1
  est.beta1 <- object$coefficients$theta
  se.beta1 <- sqrt(diag(object$vcov)[param_id$beta1])
  tval.beta1 <- est.beta1 / se.beta1
  pval.beta1 <- 2 * stats::pnorm(abs(tval.beta1), lower.tail = FALSE)

  beta1 <- cbind(Estimate = est.beta1,
                 `Std. error` = se.beta1,
                 `t value` = tval.beta1,
                 `Pr(>|t|)` = pval.beta1)

  # Summary for beta2
  est.beta2 <- object$coefficients$p
  se.beta2 <- sqrt(diag(object$vcov)[param_id$beta2])
  tval.beta2 <- est.beta2 / se.beta2
  pval.beta2 <- 2 * stats::pnorm(abs(tval.beta2), lower.tail = FALSE)

  beta2 <- cbind(Estimate = est.beta2,
                 `Std. error` = se.beta2,
                 `t value` = tval.beta2,
                 `Pr(>|t|)` = pval.beta2)

  # Additional model parameters
  se <- sqrt(diag(object$vcov))
  if (is.null(object$alpha)) {
    par <- cbind(phi = c(object$phi, se[param_id$phi]),
                 mu = c(object$mu, se[param_id$mu]),
                 sigma = c(object$sigma, se[param_id$sigma]))
  } else {
    par <- cbind(phi = c(object$phi, se[param_id$phi]),
                 alpha = c(object$alpha, se[param_id$alpha]),
                 mu = c(object$mu, se[param_id$mu]),
                 sigma = c(object$sigma, se[param_id$sigma]))
  }

  rownames(par) <- c("Estimate", "Std. error")

  out <- list(call = object$call,
              theta = beta1,
              p = beta2,
              par = par,
              links = object$links,
              residuals = stats::residuals(object),
              iterations = object$EM_iterations,
              logLik = object$logLik,
              AIC = 2 * (r - object$logLik),
              BIC = log(object$nobs) * r - 2 * object$logLik)


  class(out) <- "summary.cdnbcr"
  out
}

# Print summary
#' @rdname summary.cdnbcr
#' @export
print.summary.cdnbcr <- function(x, digits = getOption("digits"), ...)
{

  cat("Call:\n")
  print(x$call)

  cat("\nSummary for residuals:\n")
  res <- x$residuals
  print(structure(round(stats::quantile(as.numeric(res), probs = c(0, 0.25, 0.5, 0.75, 1)),
                        digits = digits),
                  .Names = c("Min", "1Q", "Median", "3Q", "Max")
  ))

  cat("\nRegression model for theta (", x$links$theta.link, " link):\n", sep = "")
  stats::printCoefmat(round(x$theta, digits))

  cat("\nRegression model for p (", x$links$p.link, " link):\n", sep = "")
  stats::printCoefmat(round(x$p, digits))

  cat("\nDispersion, dependence, and baseline parameters:\n")
  print(round(x$par, digits))

  cat("\n---\nLog-lik value:", round(x$logLik, digits),
      "\nAIC:", round(x$AIC, digits),
      "and BIC:", round(x$BIC, digits),
      "\nEM iterations:", x$it, "\n")

  invisible(x)
}

# Plot method ---------------------------------------------------------------------------------
#' Diagnostic Plots for the CDNBCR Model
#'
#' Produces two diagnostic plots to assess the goodness-of-fit of a Correlated Destructive Negative
#'     Binomial Cure Rate model fit based on the Cox-Snell residuals. Available plots include the
#'     Kaplan-Meier estimate and the cumulative hazard plot.
#'
#' @param x An object of class \code{"cdnbcr"}, a result of a call to \code{\link{cdnbcr}}.
#' @param which Numeric; if a subset of the plots is required, specify a subset
#'     of the numbers \code{1:2}.
#' @param ask Logical; if \code{TRUE}, the user is asked before each plot.
#' @param col.lines A vector with dimension two with the color for empirical and expected lines.
#' @param pch,cex,lwd,... Graphical parameters (see \code{\link[graphics]{par}}).
#'
#' @author Diego I. Gallardo <\email{diego.gallardo.mateluna@gmail.com}>
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
#'
#' @details
#' The function produces two diagnostic plots for assessing model fit:
#' \itemize{
#'    \item{The Kaplan-Meier estimate for Cox-Snell residuals, compared to the
#'     expected survival function of an exponential distribution with mean 1;}
#'    \item{The cumulative hazard function against the Cox-Snell residuals, which
#'     should align approximately with the identity line if the model is well-fitted.}
#' }
#'
#' @export
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
#' @returns \code{plot} method for \code{"cdnbcr"} objects returns two types of
#'     diagnostic plots based on the Cox-Snell residuals.
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
#'               sex + trt + thickness + age, data = e1690)
#'
#' ## Plot of the Cox-Snell residuals
#' oldpar <- par(mfrow = c(1, 2))
#' plot(fit, ask = FALSE)
#' par(oldpar)}
plot.cdnbcr <- function(x, which = 1:2,
                        ask = prod(graphics::par("mfcol")) < length(which) &&
                          grDevices::dev.interactive(),
                        col.lines = c("black", "#56B4E9"),
                        lwd = 2, pch = 16, cex = 0.8,
                        ...)
{

  dots <- list(...)

  if(!is.numeric(which) || any(which < 1) || any(which > 2))
    stop("`which' must be in 1:2")

  ## Reading
  res <- stats::residuals(x)
  KM <- survival::survfit(res ~ 1, res)

  ## Graphical parameters setting
  if (ask) {
    op <- graphics::par(ask = TRUE)
    on.exit(graphics::par(op))
  }

  ## Plots to shown
  show <- rep(FALSE, 2)
  show[which] <- TRUE

  ## Residuals versus Fitted values
  if (show[1]){

    xlab <- if (is.null(dots$xlab)) "Cox-Snell residuals" else dots$xlab
    ylab <- if (is.null(dots$ylab)) "Survival function" else dots$ylab

    dots$xlab <- NULL
    dots$ylab <- NULL

    do.call(plot, c(list(x = KM, xlab = xlab, ylab = ylab, col = col.lines[1]), dots))
    graphics::curve(stats::pexp(x, lower.tail = FALSE), col = col.lines[2], add = TRUE, lwd = 2)
    graphics::legend("topright", c("KM", "Exp(1)"), lty = 1, col = col.lines, bty = "n", lwd = 2)

  }

  ## Residuals versus index observation
  if (show[2]){

    xlab <- if (is.null(dots$xlab)) "Cox-Snell residuals" else dots$xlab
    ylab <- if (is.null(dots$ylab)) "Cumulative hazard" else dots$ylab

    dots$xlab <- NULL
    dots$ylab <- NULL

    do.call(plot, c(list(x = sort(res[, 1]), y = KM$cumhaz, xlab = xlab, ylab = ylab,
                         pch = pch, cex = cex, col = col.lines[1]), dots))
    graphics::abline(0, 1, col = col.lines[2], lwd = lwd)

  }

  invisible(x)

}

# Predict method ------------------------------------------------------------------------------
#' Prediction Method for CDNBCR Models
#'
#' Computes fitted values and predictions from a correlated destructive negative
#' binomial cure rate (CDNBCR) model fitted with \code{\link{cdnbcr}}. Predictions
#' may include the survival function, cure rate, expected number of initial
#' competing causes, or the activation probability of competing causes.
#'
#' @param object An object of class \code{"cdnbcr"}, as returned by \code{\link{cdnbcr}}.
#' @param newdata An optional data frame containing covariate values for which
#'     predictions are required. If omitted, predictions are computed for the
#'     data used in the model fit.
#' @param type Character string indicating the type of prediction. Possible values are:
#'     \describe{
#'       \item{\code{"survival"}}{Predicted survival function evaluated at the values supplied in \code{time}.}
#'       \item{\code{"cure"}}{Predicted cure fraction (probability of being cured).}
#'       \item{\code{"theta"}}{Predicted expected number of initial competing causes.}
#'       \item{\code{"p"}}{Predicted activation probability of an initial competing cause.}
#'     }
#' @param time Numeric vector of time points at which the survival function is
#'     evaluated. Only used when \code{type = "survival"}.
#' @param na.action A function specifying how missing values in \code{newdata}
#'     should be handled. The default is to return \code{NA} for predictions
#'     involving incomplete cases.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return
#' A numeric vector or matrix of predicted values. The structure of the output
#' depends on the selected \code{type}:
#' \describe{
#'   \item{\code{"survival"}}{A matrix of survival probabilities with rows
#'     corresponding to observations in \code{newdata} (or the original data)
#'     and columns corresponding to the values in \code{time}.}
#'   \item{\code{"cure"}}{A numeric vector with the predicted cure fractions.}
#'   \item{\code{"theta"}}{A numeric vector with the predicted expected numbers of
#'     initial competing causes.}
#'   \item{\code{"p"}}{A numeric vector with the predicted activation probabilities.}
#' }
#'
#' @export
#'
#' @author Diego I. Gallardo <\email{diego.gallardo.mateluna@gmail.com}>
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
#'
#' @references
#' Cox, D. R., and Snell, E. J. (1968). A general definition of residuals.
#'     \emph{Journal of the Royal Statistical Society, Series B}, \bold{30}, 248--265.
#'
#' De Medeiros, R. M. R., Bourguignon, M., Gómez, Y. M., and Gallardo, D. I. (2026).
#'     A correlated approach to cancer cell counting in cure rate models.
#'
#' Rodrigues, J., De Castro, M., Balakrishnan, N., and Cancho, V. G. (2011).
#'     Destructive weighted Poisson cure rate models. \emph{Lifetime Data Analysis},
#'     \bold{17}, 333--346.
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
#'
#' ## New data for predictions
#' newdata <- data.frame(trt = c("Control", "Control", "Chemotherapy", "Chemotherapy"),
#'                       age = median(e1690$age),
#'                       sex = c("Male", "Female", "Male", "Female"),
#'                       thickness = median(e1690$thickness),
#'                       nodeII = c(0, 0, 0, 0),
#'                       nodeIII = c(0, 0, 0, 0),
#'                       nodeIV = c(1, 1, 1, 1))
#' newdata
#'
#' ## Fitted survival curves
#' pred <- predict(fit, newdata)
#'
#' plot(pred[1, ], type = "l", ylim = c(0, 1), xlab = "Time", ylab = "Survival")
#' lines(pred[2, ], col = 2, lty = 2)
#' lines(pred[3, ], col = 3, lty = 3)
#' lines(pred[4, ], col = 4, lty = 4)
#' legend("topright", legend = c("trt: Control, sex: Male",
#'                               "trt: Control, sex: Female",
#'                               "trt: Chemotherapy, sex: Male",
#'                               "trt: Chemotherapy, sex: Female"),
#'        col = 1:4, lty = 1:4)
#'
#' ## Predicted cure rates
#' predict(fit, newdata, type = "cure")
#'
#' ## Predicted expected number of initial competing causes
#' predict(fit, newdata, type = "theta")
#'
#' ## Predicted activation probability of an initial competing cause
#' predict(fit, newdata, type = "p")}
predict.cdnbcr <- function(object, newdata = NULL,
                           type = c("survival", "cure", "theta", "p"),
                           time, na.action = stats::na.pass, ...)
{

  mf <- stats::model.frame(object)
  mt <- attr(mf, "terms")
  xlevels <- stats::.getXlevels(mt, mf)

  phi <- object$phi
  alpha <- if (is.null(object$alpha)) 0L else object$alpha
  mu <- object$mu
  sigma <- object$sigma

  theta.link <- object$links$theta.link
  p.link <- object$links$p.link

  if(is.null(newdata)) {

    n <- object$nobs
    theta <- object$fitted$theta
    p <- object$fitted$p

  } else {

    f1 <- stats::formula(object$formula, rhs = 1)
    f2 <- stats::formula(object$formula, rhs = 2)

    rhs1 <- stats::as.formula(call("~", f1[[3]]))
    rhs2 <- stats::as.formula(call("~", f2[[3]]))

    newdata1 <- stats::model.frame(rhs1, newdata, na.action = na.action, xlev = xlevels)
    newdata2 <- stats::model.frame(rhs2, newdata, na.action = na.action, xlev = xlevels)

    Z1 <- stats::model.matrix(rhs1, newdata1, xlev = xlevels)
    Z2 <- stats::model.matrix(rhs2, newdata2, xlev = xlevels)

    n <- nrow(newdata)
    theta <- c(stats::make.link(theta.link)$linkinv(Z1%*%object$coefficients$theta))
    p <- c(stats::make.link(p.link)$linkinv(Z2%*%object$coefficients$p))

  }

  type <- match.arg(type, c("survival", "cure", "theta", "p"))

  rval <- switch(type,

                 "theta" = {
                   theta
                 },

                 "p" = {
                   p
                 },

                 "cure" = {
                   c(cure_rate(theta, phi, p, alpha))
                 },

                 "survival" = {

                   if (missing(time)) {
                     y <- stats::model.response(stats::model.frame(object))
                     time <- seq(0, max(y[,1], na.rm = TRUE), length.out = 100)
                   }


                   pred <- matrix(NA, n, length(time))
                   for (i in 1:n){
                      pred[i, ] <- pcdnbcr(time, theta[i], phi, p[i], alpha, mu, sigma, lower.tail = FALSE)
                   }

                   pred

                 })

  rval

}
