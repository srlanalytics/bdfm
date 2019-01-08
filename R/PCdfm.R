PCdfm <- function(Y, m, p, FC = 0, Bp = NULL, lam_B = 0, Hp = NULL, lam_H = 0, nu_q = 0, nu_r = NULL, ID = "PC_full", ITC = T, reps = 1000, burn = 500) {

  # ----------- Preliminaries -----------------
  Y <- as.matrix(Y)
  k <- ncol(Y)
  r <- nrow(Y)
  n_obs <- sum(is.finite(Y))
  if (ITC) {
    itc <- colMeans(Y, na.rm = T)
    Y <- Y - matrix(1, r, 1) %x% t(itc) # De-mean data. Data is not automatically stadardized --- that is left up to the user and should be done before estimation if desired.
  } else {
    itc <- rep(0, k)
  }

  #Estimate principal components

  if (ID == "PC_sub") {
    Ysub <- Y_sub(Y) # submatrix of Y with complete data, i.e. no missing values
    PC   <- PrinComp(Ysub$Ysub, m)
    X    <- matrix(NA, r, m)
    X[Ysub$ind, ] <- PC$components
  } else if (ID == "PC_full") {
    PC <- PrinComp(Y, m)
    X  <- PC$components
    if (!any(!is.na(PC$components))) {
      stop("Every period contains missing data. Try setting ID to PC_sub.")
    }
  }

  # ----------- Format Priors ------------------
  #enter priors multiplicatively so that 0 is a weak prior and 1 is a strong prior (additive        priors are relative to the number of observations)
  lam_B <- r*lam_B + 1
  nu_q  <- r*nu_q  + 1
  lam_H <- r*lam_H + 1
  if(is.null(nu_r)){
    nu_r = rep(1,k)
  }else{
    if(length(nu_r) != k){stop("Length of nu_r must equal the number of observed series")}
    nu_r = r*nu_r + rep(1,k)
  }
  if(is.null(Hp)){
    Hp <- matrix(0,k,m)
  }
  if(is.null(Bp)){
    Bp <- matrix(0,m,m*p)
  }
  # ---------------------------------------------

  #Estimate parameters of the observation equation
  Hest <- BReg_diag(X, Y, Int = F, Bp = Hp, lam = lam_H, nu = nu_r, reps = reps, burn = burn)
  H    <- Hest$B
  R    <- diag(c(Hest$q))

  #Estimate parameters of the transition equation (Bayesian VAR)
  Z    <- stack_obs(X, p = p)
  Z    <- as.matrix(Z[-nrow(Z),])
  xx   <- as.matrix(X[-(1:p),])
  indZ <- which(apply(Z, MARGIN = 1, FUN = AnyNA))
  indX <- which(apply(xx, MARGIN = 1, FUN = AnyNA))
  ind  <- unique(c(indZ, indX)) #index of rows with missing values in Z and xx
  Best <- BReg(Z[-ind,,drop = FALSE], xx[-ind, , drop = FALSE], Int = FALSE, Bp = Bp, lam = lam_B, nu = nu_q, reps = reps, burn = burn)
  B    <- Best$B
  q    <- Best$q

  if(FC>0){
    tmp <- matrix(NA,FC,k)
    Y   <- rbind(Y, tmp)
    r   <- r + FC
  }
  
    Jb <- Matrix::Diagonal(m*p)

    Est <- DSmooth(B = B, Jb = Jb, q = q, H = H, R = R, 
            Y = Y, freq = rep(1,k), LD = rep(0,k))

    BIC <- log(n_obs)*(m*p + m^2 + k*m + k) - 2*Est$Lik

    Out <- list(
      B = B,
      q = q,
      H = H,
      R = R,
      values  = Est$Ys,
      factors = Est$Z[,1:m],
      Qstore  = Best$Qstore, # lets us look at full distribution
      Bstore  = Best$Bstore,
      Rstore  = Hest$Rstore,
      Hstore  = Hest$Hstore,
      Kstore  = Est$Kstr,
      PEstore = Est$PEstr,
      Lik     = Est$Lik,
      BIC     = BIC
    )
  return(Out)
}