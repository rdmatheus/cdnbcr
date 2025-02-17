#' @name cdnbcr-methods
#'
#' @title Methods for "cdnbcr" objects
#'
#' @param x,object An object of class \code{"cdnbcr"}.
#' @param formula A model \code{\link[Formula]{Formula}} or \code{\link[stats]{terms}} object or an
#'     \code{"cdnbcr"} object.
#' @param k Numeric, the penalty per parameter to be used; the default
#'     \code{k = 2} is the classical AIC. See \code{\link[stats]{AIC}}.
#' @param submodel Specifies which part of the model should be used. Options:
#'     \code{"full"} (default) for all parameters, \code{"theta"} for the regression submodel of
#'     \code{theta}, and \code{"p"} for the regression submodel of\code{p}.
#' @param digits Minimal number of significant digits, see \code{\link[base]{print.default}}.
#' @param ... Additional arguments passed to or from other methods.
#'
#' @author
#' Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
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
model.matrix.cdnbcr <- function(object, submodel = c("theta", "p"), ...) {
  submodel <- match.arg(submodel)
  rval <- if(!is.null(object$x[[submodel]])) object$x[[submodel]]
  else stats::model.matrix(object$terms[[submodel]], stats::model.frame(object))
  return(rval)
}

# Print
#' @rdname cdnbcr-methods
#' @export
print.cdnbcr <- function(x, digits = getOption("digits"), ...)
{

  r <- x$nobs - x$df.residual

  cat("Call:\n")
  print(x$call)

  cat("\ntheta submodel with", x$links$theta.link, "link function:\n")
  print(round(x$coefficients$theta, digits))

  cat("\np submodel with", x$links$p.link, "link function:\n")
  print(round(x$coefficients$p, digits))

  cat("\nAdditional model parameter estimates:\n")
  if (is.null(x$alpha)) {
    par <- cbind(phi = x$phi, mu = x$mu, sigma = x$sigma)
  } else {
    par <- cbind(phi = x$phi, alpha = x$alpha, mu = x$mu, sigma = x$sigma)
  }

  rownames(par) <- " "
  print(round(par, digits))

  cat("\n---\nLog-lik value:", x$logLik,
      "\nAIC:", 2 * (r - x$logLik),
      "and BIC:", log(x$nobs) * r - 2 * x$logLik,
      "\nEM iterations:", x$EM_iterations, "\n")

  invisible(x)
}

# Summary
#' @rdname cdnbcr-methods
#' @export
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

  # Summary for residuals
  rq <- object$residuals
  residuals <- cbind(mean(rq),
                     stats::sd(rq),
                     mean((rq - mean(rq))^3) / (stats::sd(rq)^3),
                     mean((rq - mean(rq))^4) / (stats::sd(rq)^4))
  colnames(residuals) <- c("Mean", "Std. dev.", "Skewness", "Kurtosis")
  rownames(residuals) <- " "

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
              residuals = residuals,
              beta1 = beta1,
              beta2 = beta2,
              par = par,
              links = object$links,
              it = object$EM_iterations,
              logLik = object$logLik,
              AIC = 2 * (r - object$logLik),
              BIC = log(object$nobs) * r - 2 * object$logLik)


  class(out) <- "summary.cdnbcr"
  out
}

# Print summary
#' @rdname cdnbcr-methods
#' @export
print.summary.cdnbcr <- function(x, digits = getOption("digits"), ...)
{

  cat("Call:\n")
  print(x$call)

  cat("\nSummary for residuals:\n")
  print(round(x$residuals, digits))

  cat("\ntheta submodel with", x$links$theta.link, "link function:\n")
  stats::printCoefmat(round(x$beta1, digits))

  cat("\np submodel with", x$links$p.link, "link function:\n")
  stats::printCoefmat(round(x$beta2, digits))

  cat("\nAdditional model parameters:\n")
  print(round(x$par, digits))

  cat("\n---\nLog-lik value:", x$logLik,
      "\nAIC:", x$AIC,
      "and BIC:", x$BIC,
      "\nEM iterations:", x$it, "\n")

  invisible(x)
}

# Parameter estimates
#' @rdname cdnbcr-methods
#' @export
coef.cdnbcr <- function(object,
                        submodel = c("full", "theta", "p"), ...) {

  submodel <- match.arg(submodel, c("full", "theta", "p"))
  switch(submodel,
          "full" = object$coefficients,
         "theta" = object$coefficients$theta,
             "p" = object$coefficients$p)
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

#  Variance-covariance matrix
#' @rdname cdnbcr-methods
#' @export
vcov.cdnbcr <- function(object,
                        submodel = c("full", "theta", "p"), ...) {

  submodel <- match.arg(submodel, c("full", "theta", "p"))
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

  switch(submodel,
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

# Plot
#' Diagnostic Plots for the Correlated Destructive Negative Binomial Cure Rate
#'     Model
#'
#' Produces diagnostic plots to assess the goodness-of-fit of a Correlated Destructive Negative
#'     Binomial Cure Rate model fit based on the Cox-Snell residuals. Available plots include the
#'     Kaplan-Meier estimate and the cumulative hazard plot.
#'
#' @param x An object of class \code{cdnbcr}.
#' @param which numeric; if a subset of the plots is required, specify a subset
#'     of the numbers \code{1:2}.
#' @param ask logical; if \code{TRUE}, the user is asked before each plot.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' The function produces two diagnostic plots for assessing model fit:
#' (1) The Kaplan-Meier estimate for Cox-Snell residuals, compared to the expected
#'     survival function of an exponential distribution with mean 1; and
#' (2) The cumulative hazard function against the Cox-Snell residuals, which
#'     should align approximately with the identity line if the model is well-fitted.
#'
#' @export
#'
plot.cdnbcr <- function(x, which = 1:2,
                        ask = prod(graphics::par("mfcol")) < length(which) &&
                          grDevices::dev.interactive(),
                        ...)
{

  if(!is.numeric(which) || any(which < 1) || any(which > 2))
    stop("`which' must be in 1:2")

  ## Reading
  res <- data.frame(res = x$residuals,
                    status = stats::model.response(stats::model.frame(x))[, 2])

  KM <- survival::survfit(survival::Surv(res, status) ~ 1, res)

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
    graphics::plot(KM, ylab = "Survival function", xlab = "Cox-Snell residuals")
    graphics::curve(stats::pexp(x, lower.tail = FALSE), col = "#56B4E9", add = TRUE, lwd = 2)
    graphics::legend("topright", c("KM", "Exp(1)"), lty = 1, col = c("black", "#56B4E9"), bty = "n", lwd = 2)
  }

  ## Residuals versus index observation
  if (show[2]){
    graphics::plot(sort(res$res), KM$cumhaz, pch = 16, ylab = "Cumulative hazard",
                   xlab = "Cox-Snell residuals", cex = 0.8)
    graphics::abline(0, 1, col = "#56B4E9", lwd = 3)
  }

  invisible(x)

}


#' Predict Method for the Correlated Destructive Negative Binomial Cure
#'     Rate Regression
#'
#' Obtains predictions from a fitted \code{"cdnbcr"} object.
#'
#' @param object A \code{"cdnbcr"} object.
#' @param newdata Optionally, a data frame in which to look for variables
#'     with which to predict. If omitted, the fitted linear predictors are
#'     used.
#' @param type The type of prediction required. The default is the
#'     fitted survival function (\code{"survival"}). The alternative
#'     (\code{"cure"}) provides the fitted cure rate function.
#' @param na.action function determining what should be done with missing
#'     values in \code{newdata}. The default is to predict \code{NA}.
#' @param ...  arguments passed to or from other methods.
#'
#' @return A vector of predictions.
#' @export
#'
predict.cdnbcr <- function(object, newdata = NULL,
                           type = c("survival", "cure"),
                           na.action = stats::na.pass, ...)
{

  mf <- stats::model.frame(object)
  mt <- attr(mf, "terms")
  xlevels <- stats::.getXlevels(mt, mf)

  type <- match.arg(type, c("survival", "cure"))

  phi <- object$phi
  alpha <- if (is.null(object$alpha)) 0L else object$alpha
  mu <- object$mu
  sigma <- object$sigma

  theta.link <- object$links$theta.link
  p.link <- object$links$p.link

  if(is.null(newdata)) {

    time <- stats::model.response(stats::model.frame(object))[, 1]
    theta <- object$fitted$theta
    p <- object$fitted$p

  } else {

    mf <- stats::model.frame(object$terms[["full"]], newdata,
                             na.action = na.action, xlev = xlevels)

    newdata <- newdata[rownames(mf), , drop = FALSE]

    time <- c(stats::model.response(mf)[, 1])
    Z1 <- stats::model.matrix(stats::delete.response(object$terms$theta), mf)
    Z2 <- stats::model.matrix(stats::delete.response(object$terms$p), mf)

    theta <- c(stats::make.link(theta.link)$linkinv(Z1%*%object$coefficients$theta))
    p <- c(stats::make.link(p.link)$linkinv(Z2%*%object$coefficients$p))

  }

  rval <- switch(type,
                 "cure" = {
                   c(cure_rate(theta, phi, p, alpha))
                 },
                 "survival" = {
                   c(pcdnbcr(time, theta, phi, p, alpha, mu, sigma, lower.tail = FALSE))
                 })

  rval

}
