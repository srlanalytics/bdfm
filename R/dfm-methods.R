
#' Extractor Functions for Dynamic Factor Models
#'
#' @param x obeject of class `"dfm"`
#' @export
factors <- function(x) {
  stopifnot(inherits(x, "dfm"))
  x$factors
}


#' @export
#' @name factors
adjusted <- function(x) {
  stopifnot(inherits(x, "dfm"))
  x$adjusted
}

# methods
#' @export
#' @method predict dfm
predict.dfm <- function(object, newdata = NULL, return_intermediates = FALSE, ...) {
  if(is.null(newdata)){
    return(object$values)
  }else{
    if(NCOL(object$values) != NCOL(newdata)){
      stop("'newdata' must include the same observable series as the original model was fitted with.")
    }

    # logs
    if (!is.null(object$logs)) {
      newdata[, object$logs] <- log(newdata[, object$logs])
    }

    # differences
    if (!is.null(object$diffs)) {
      newdata_lev <- newdata
      newdata[, object$diffs] <- sapply(object$diffs, mf_diff, fq = object$freq, Y = newdata)
    }

    # drop outliers
    newdata[abs(scale(newdata)) > object$outlier_threshold] <- NA

    #scale
    if (object$scale) {
      newdata <- 100*scale(newdata)
      y_scale  <- attr(newdata, "scaled:scale")
      y_center <- attr(newdata, "scaled:center")
    }

    est <- DSmooth(object$B, object$Jb, object$q, object$H, diag(object$R), newdata, object$freq, object$differences)

    if(object$scale){
      est$Ys <- (matrix(1, nrow(est$Ys), 1) %x% t(y_scale)) * (est$Ys / 100) + (matrix(1, nrow(est$Ys), 1) %x% t(y_center))
    }

    # undo differences
    if (!is.null(object$diffs)) {
      est$Ys[,object$diffs] <- sapply(object$diffs, FUN = level, fq = object$freq, Y_lev = newdata_lev, vals = est$Ys)
    }

    # undo logs
    if (!is.null(object$logs)) {
      est$Ys[,object$logs] <- exp(est$Ys[,object$logs])
    }

    #Return intermediate values of low frequency data?
    if (length(unique(object$freq))>1 && !return_intermediates){
      est$Ys[,which(object$freq != 1)] <- do.call(cbind, lapply(X = which(object$freq != 1), FUN = drop_intermediates,
                                                             freq = object$freq, Y_raw = newdata, vals = est$Ys))
    }

    colnames(est$Ys) <- colnames(newdata)

    return(est$Ys)
  }

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
  printCoefmat(obs, digits = digits, signif.stars = signif.stars,
               na.print = "NA")

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


