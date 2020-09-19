# m <- 1
# p <- "auto"
# freq <- "auto"
# method = "bayesian"
# Bp <- NULL
# preD <- 1
# lam_B = 0
# trans_df = 0
# Hp = NULL
# lam_H = 0
# obs_df = NULL
# ID = "pc_long"
# keep_posterior = 1
# reps = 1000
# burn = 500
# verbose = T
# tol = .01
# FC = 3
# # logs = NULL
# # diffs = NULL
# logs = c( 2,  4,  5,  8,  9, 10, 11, 12, 15, 16, 17, 21, 22)
# diffs = c(2, 4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24)
# outlier_threshold <- 4
# scale = TRUE
# orthogonal_shocks = F

dfm_core <- function(Y, m, p, FC = 0, method = "bayesian", scale = TRUE, logs = "auto",
                     outlier_threshold = 4, diffs = "auto", freq = "auto", preD = NULL,
                     Bp = NULL, lam_B = 0, trans_df = 0, Hp = NULL, lam_H = 0, obs_df = NULL, ID = "pc_long",
                     keep_posterior = NULL, reps = 1000, burn = 500, verbose = TRUE,
                     tol = 0.01, interpolate = FALSE, orthogonal_shocks = FALSE) {

  #-------Data processing-------------------------

  data <- Y   # will rename to 'data' later on
  data_orig <- Y

  k <- NCOL(Y) # number of series

  # frequency
  if (freq == "auto") {
    freq <- apply(Y, MARGIN = 2, FUN = get_freq)
  } else if (!is.integer(freq) || length(freq) != ncol(Y)) {
    stop("Argument 'freq' must be 'auto' or integer valued with
         length equal to the number data series")
  }

  if (p == "auto"){
    p <- max(freq)
  }

  # if (FC == "auto"){
  #   FC <- max(freq)
  # }

  # add forecast periods
  if (FC > 0) {
    tmp <- matrix(NA, FC, k)
    Y <- rbind(Y, tmp)
  }


  if((!is.null(logs) && logs == "auto") || (!is.null(diffs) && diffs == "auto")){
    # parse_named_vector <- function(x, name = deparse(substitute(x))) {
    #   z <- paste(paste0("  ", names(x), " = ", x), collapse = ",\n")
    #   paste0(name, " = c(\n", z, "\n)")
    # }
    parse_vector <- function(x, name = deparse(substitute(x))) {
      if (length(x) == 1 && is.null(names(x)) && x == TRUE) return(paste0(name, " = 1"))
      x.char <- names(x)[x]
      if (length(x.char) == 0) return(paste0(name, " = NULL"))
      z <- paste(paste0("  \"", x.char, "\""), collapse = ",\n")
      paste0(name, " = c(\n", z, "\n)")
    }
    do_log_diff <- should_log_diff(Y)
    if (verbose) message("auto log/diff detection, with:", do_log_diff)
    # if(logs == "auto"){
    #   logs <- do_log_diff[1,,drop = FALSE]
    #   print(logs)
    #   if (verbose) message(paste("logs:", names(logs)[logs]), if (diffs == "auto") ",")
    # }
    # if(diffs == "auto"){
    #   diffs <- do_log_diff[2,,drop = FALSE]
    #   if (verbose) message(paste("diffs:", names(diffs)[diffs]))
    # }
  }



  if(is.logical(logs)){
    if(!any(logs)) logs <- NULL
  }

  if(is.logical(diffs)){
    if(!any(diffs)) diffs <- NULL
  }


  # logs
  if (!is.null(logs)) {
    logs <- standardize_index(logs, Y)
    Y[, logs] <- log(Y[, logs])
  }

  
  print(diffs)
  
  
  # differences
  if (!is.null(diffs)) {
    Y_lev <- Y
    diffs <- standardize_index(diffs, Y)
    Y[, diffs] <- sapply(diffs, mf_diff, fq = freq, Y = Y)
  }

  # if (!is.null(trans_df)) {
  #   trans_df <- standardize_index(trans_df, Y)
  # }

  if (!is.null(obs_df)) {
    obs_df <- standardize_numeric(obs_df, Y)
  }

  # specify which series are differenced for mixed frequency estimation

  LD <- rep(0, NCOL(Y))
  if (!is.null(preD)) {
    preD <- standardize_index(preD, Y)
  }
  LD[unique(c(preD, diffs))] <- 1 # in bdfm 1 indicates differenced data, 0 level data

  # drop outliers
  Y[abs(scale(Y)) > outlier_threshold] <- NA

  if (scale) {
    Y <- 100*scale(Y)
    y_scale  <- attr(Y, "scaled:scale")
    y_center <- attr(Y, "scaled:center")
  }

  if (!is.null(keep_posterior)){
    if(length(keep_posterior)>1){
      stop("Length of 'keep_posterior' cannot be greater than 1")
    }
    keep_posterior <- standardize_index(keep_posterior, Y)
  }

  if(all(!ID%in%c("pc_wide", "pc_long", "name"))){
    ID <- standardize_index(ID, Y)
  }

  if (length(unique(freq)) != 1 && method != "bayesian") {
    stop("Mixed freqeuncy models are only supported for Bayesian estimation")
  }

  if (method == "bayesian") {
    est <- bdfm(
      Y = Y, m = m, p = p, Bp = Bp,
      lam_B = lam_B, Hp = Hp, lam_H = lam_H, nu_q = trans_df, nu_r = obs_df,
      ID = ID, keep_posterior = keep_posterior, freq = freq, LD = LD, reps = reps,
      burn = burn, verbose = verbose, orthogonal_shocks = orthogonal_shocks
    )
  } else if (method == "ml") {
    est <- MLdfm(
      Y = Y, m = m, p = p, tol = tol,
      verbose = verbose, orthogonal_shocks = orthogonal_shocks
    )
  } else if (method == "pc") {
    est <- PCdfm(
      Y, m = m, p = p, Bp = Bp,
      lam_B = lam_B, Hp = Hp, lam_H = lam_H, nu_q = trans_df, nu_r = obs_df,
      ID = ID, reps = reps, burn = burn, orthogonal_shocks = orthogonal_shocks
    )
  }

  # any reason why est$Kstore is in such a strange form? why not simple lists?
  # SL: These objects are arma::field<mat> in the Rcpp code, which works much better
  # than a list internally. As I understand it, in R it is in fact already a list, just
  # one in which all the objects are matrices.
  # k_list <- lapply(seq(NROW(est$Kstore)), function(i) est$Kstore[i, 1, drop = FALSE][[1]])
  # pe_list <- lapply(seq(NROW(est$PEstore)), function(i) est$PEstore[i, 1, drop = FALSE][[1]])
  # gain_list <- lapply(k_list, function(e) t(e[1:m, , drop = FALSE]))

  names_list <- lapply(seq(NROW(Y)), function(i) names(Y[i, ])[is.finite(Y[i, ])])

  factor_update <- Map(
    function(g, pe, nm) {
      x <- g * (matrix(1, NROW(g), 1) %x% t(pe))
      colnames(x) <- nm
      x
    },
    g = est$Kstore,
    pe = est$PEstore,
    nm = names_list
  )

  est$Kstore  <- NULL # this is huge and no longer needed, so drop it
  est$PEstore <- NULL

  # get updates to keep_posterior if specified
  if(!is.null(keep_posterior)){
    idx_loading <- est$H[keep_posterior,,drop=FALSE]%*%J_MF(freq[keep_posterior], m = m, ld = LD[keep_posterior], sA = NCOL(est$Jb))
    idx_scale <- if (scale) y_scale[keep_posterior]/100 else 1
    idx_update <- lapply(factor_update, function(x) as.matrix(idx_scale * (idx_loading %*% x)) )
    # same structure as data: missing values as NA
    idx_update <- lapply(idx_update, function(e){
      tmp <- setNames(rep(NA, k), colnames(Y))
      tmp[colnames(e)] <- e
      return(tmp)
    })
    est$idx_update <- do.call(rbind, idx_update)
  }

  est$factor_update <- lapply(factor_update, function(e) e[1:m,]) #return this instead of gain and prediction error far more useful!

  if(!is.null(keep_posterior) && method == "bayesian"){ #fudge untill we allow keep_posterior to be a vector
    est$Ystore <- est$Ystore[,,1]
    est$Ymedian <- est$Ymedian[,1]
  }
  
  # undo scaling
  if(scale){
    est$values <- (matrix(1, nrow(est$values), 1) %x% t(y_scale)) * (est$values / 100) + (matrix(1, nrow(est$values), 1) %x% t(y_center))
    est$R2     <- 1 - est$R/10000
    if(!is.null(keep_posterior) && method == "bayesian"){
      est$Ystore <- est$Ystore*(y_scale[keep_posterior]/100) + y_center[keep_posterior]
      est$Ymedian <- est$Ymedian*(y_scale[keep_posterior]/100) + y_center[keep_posterior]
    }
  }else{
    est$R2 <- 1 - est$R/apply(X = Y, MARGIN = 2, FUN = var, na.rm = TRUE)
  }

  # undo differences
  if (!is.null(diffs)) {
    est$values[,diffs] <- sapply(diffs, FUN = level, fq = freq, Y_lev = Y_lev, vals = est$values)
    if(!is.null(keep_posterior) && method == "bayesian" && keep_posterior%in%diffs){
      est$Ymedian <- level_simple(est$Ymedian, y_lev = Y_lev[,keep_posterior], fq = freq[keep_posterior])
      est$Ystore  <- apply(est$Ystore, MARGIN = 2, FUN = level_simple, y_lev = Y_lev[,keep_posterior], fq = freq[keep_posterior])
    }
  }

  # undo logs
  if (!is.null(logs)) {
    est$values[,logs] <- exp(est$values[,logs])
    if(!is.null(keep_posterior) && method == "bayesian" && keep_posterior%in%logs){
      est$Ymedian <- exp(est$Ymedian)
      est$Ystore  <- exp(est$Ystore)
    }
  }

  #Return intermediate values of low frequency data?
  if (length(unique(freq))>1 && !interpolate){
    est$values[,which(freq != 1)] <- do.call(cbind, lapply(X = which(freq != 1), FUN = drop_intermediates,
                                                           freq = freq, Y_raw = Y, vals = est$values))
  }

  est$freq  <- freq
  est$logs  <- logs
  est$diffs <- diffs
  est$scale <- scale
  est$outlier_threshold <- outlier_threshold
  est$differences <- LD

  colnames(est$values) <- colnames(data)

  # adjusted series: align 'values' with original series
  # browser()
  est$adjusted <- substitute_in_benchmark(est$values, data_orig)

  return(est)
}


substitute_in_benchmark <- function(x, benchmark) {
  nfct <- NROW(x) - NROW(benchmark)
  if (nfct > 0) {
    benchmark <- rbind(benchmark, matrix(NA_real_, ncol = NCOL(benchmark), nrow = nfct))
  }
  stopifnot(identical(NROW(x), NROW(benchmark)))

  benchmark[is.na(benchmark)] <- x[is.na(benchmark)]
  benchmark
}

