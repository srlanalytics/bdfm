
#' Extractor Functions for Dynamic Factor Models
#'
#' `predict` extracts the predicted values of model. `adjusted` returns the
#' original series with the predicted values substituted if missing. `factors`
#' return the factor(s) of the model.
#'
#' @param object object of class `"dfm"`
#' @param ... unused, to comply with generic
#' @export
#' @examples
#' \dontrun{
#' # predict vs. adjusted
#' library(tsbox)
#' dta0 <- ts_seas(cbind(mdeaths, fdeaths))  # seasonally adjust
#' dta <- dta0
#' dta[1:10, 2] <- NA
#'
#' m <- dfm(dta)
#' ts_plot(predict(m)[, 'fdeaths'], dta0[, 'fdeaths'])
#' ts_plot(adjusted(m)[, 'fdeaths'], dta0[, 'fdeaths'])
#' }
factors <- function(object) {
  stopifnot(inherits(object, "dfm"))
  object$factors
}


#' @export
#' @name factors
adjusted <- function(object) {
  stopifnot(inherits(object, "dfm"))
  object$adjusted
}

#' @name factors
#' @export
#' @method predict dfm
predict.dfm <- function(object, ...) {
  object$values
}

#' @export
#' @method print dfm
print.dfm <- function(x, ...) {

  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n",
      sep = "")
  map_model <- c(
    "bayesian" = "Bayesian dynamic factor model",
    "ml" = "Dynamic factor model (maximum Likelihood)",
    "pc" = "Dynamic factor model (two step)"
  )

  cat("\n")
  cat(unname(map_model[x$method]), "with",
    nrow(x$B), "factor(s) and", ncol(x$B) / nrow(x$B), "lag(s)"
  )
  cat("\n\n")
  invisible(x)
}



#' @export
#' @method summary dfm
#' @importFrom stats printCoefmat
summary.dfm <- function(object, digits = max(3, getOption("digits") - 3), ...) {

  cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"), "\n",
      sep = "")

  map_model <- c(
    "bayesian" = "Bayesian dynamic factor model",
    "ml" = "Dynamic factor model (maximum Likelihood)",
    "pc" = "Dynamic factor model (two step)"
  )

  cat("\n")
  cat(unname(map_model[object$method]), "with",
    nrow(object$B), "factor(s) and", ncol(object$B) / nrow(object$B), "lag(s)"
  )
  cat("\n")

  obs_var <- as.matrix(1 -(object$R / 10000))
  colnames(obs_var) <- "Sample Fit"
  obs_loadings <- object$H
  colnames(obs_loadings) <- paste0("Load (", seq(ncol(obs_loadings)), ")")
  obs <- cbind(obs_var, obs_loadings)

  cat("\n")
  cat("Observation Equation:\n")
  printCoefmat(obs, digits = digits, na.print = "NA")

  cat("\n")

  cat("Log Likelihood:", object$Lik)
  cat("  BIC:", object$BIC)
  cat("\n\n")
  invisible(object)


  # cat(
  #   "Call: \n Bayesian dynamic factor model with", nrow(object$B),
  #   "factor(s) and", ncol(object$B) / nrow(object$B), "lag(s)."
  # )
  # cat("\n \n")
  # cat("Log Likelihood:", object$Lik)
  # cat("\n \n")
  # cat("BIC:", object$BIC)
  # cat("\n \n")
  # cat("Posterior medians for transition equation: \n")
  # cat("\n Coefficients B: \n")
  # print(object$B)
  # cat("\n Covariance Q: \n")
  # print(object$q)
  # cat("\n \n")
  # cat("Posterior medians for observation equation: \n")
  # cat("\n Coefficients H: \n")
  # # H <- data.frame(object$H)
  # # row.names(H) <- colnames(object$values)
  # # colnames(H) <- as.character(seq(1, ncol(H)))
  # print(object$H)
  # cat("\n Shocks R: \n")
  # # r <- data.frame(diag(object$R))
  # # row.names(r) <- colnames(object$values)
  # # colnames(r) <- "Variance of Shocks"
  # print(object$R)
}


