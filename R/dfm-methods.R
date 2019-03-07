
#' extractor function for factors
#' @param x obeject of class `"dfm"`
#' @export
factors <- function(x) {
  stopifnot(inherits(x, "dfm"))
  x$factors
}

# methods
#' @export
#' @method predict dfm
predict.dfm <- function(object, ...) {
  object$values
}

#' @export
#' @method print dfm
print.dfm <- function(x, ...) {

  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
  cat("Bayesian dynamic factor model with",
    nrow(x$B), "factor(s) and", ncol(x$B) / nrow(x$B), "lag(s)."
  )
  cat("\n")
  cat("Log Likelihood:", x$Lik, " ")
  cat("BIC:", x$BIC)
  cat("\n \n")
  cat("Goodness of in-sample fit: \n")
  print(x$R2)
  invisible(x)
}

#' @export
#' @method summary dfm
summary.dfm <- function(object, ...) {
  cat(
    "Call: \n Bayesian dynamic factor model with", nrow(object$B),
    "factor(s) and", ncol(object$B) / nrow(object$B), "lag(s)."
  )
  cat("\n \n")
  cat("Log Likelihood:", object$Lik)
  cat("\n \n")
  cat("BIC:", object$BIC)
  cat("\n \n")
  cat("Posterior medians for transition equation: \n")
  cat("\n Coefficients B: \n")
  print(object$B)
  cat("\n Covariance Q: \n")
  print(object$q)
  cat("\n \n")
  cat("Posterior medians for observation equation: \n")
  cat("\n Coefficients H: \n")
  # H <- data.frame(object$H)
  # row.names(H) <- colnames(object$values)
  # colnames(H) <- as.character(seq(1, ncol(H)))
  print(object$H)
  cat("\n Shocks R: \n")
  # r <- data.frame(diag(object$R))
  # row.names(r) <- colnames(object$values)
  # colnames(r) <- "Variance of Shocks"
  print(object$R)
}


