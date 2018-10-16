# Example of Bayesian estimation by simulation and comparison with ML estimation via WE(1983)
library(BDFM)

Data <- SimDFM(300)

#-------- Simulation Smoothing vs. Disturbance Smoothing -------------

B <- Data$B
q <- Data$q
H <- Data$H
R <- Data$r

# Durbin Koopman Disturbance Smoother
EstDK <- DSmooth(B, q, H, R, Data$Y_obs)
# plot results

gr(
  ts.plot(cbind(EstDK$Ys[, 1], Data$Y_true[, 1], Data$Y_obs[, 1]), col = c("red", "black", "blue"), lty = c(1, 2, 1))
)
# Put into companion form for Ksmoother
A <- rbind(B, cbind(diag(1, 6, 6), matrix(0, 6, 3)))
Q <- matrix(0, 9, 9)
Q[1:3, 1:3] <- q
HJ <- cbind(H, matrix(0, 10, 6))
# Make sparse for Ksmoother
A <- Matrix(A, sparse = T)
Q <- Matrix(Q, sparse = T)
HJ <- Matrix(HJ, sparse = T)
EstKS <- Ksmoother(A, Q, HJ, R, Data$Y_obs)

gr(
  ts.plot(cbind(EstKS$Ys[, 1], Data$Y_true[, 1], Data$Y_obs[, 1]), col = c("red", "black", "blue"), lty = c(1, 2, 1))
)

Zero <- EstDK$Ys - EstKS$Ys # compare the two --- the difference should be zero!

#---------------- Bayesian Estimation -------------------------
Bayes <- BDFM(Data$Y_obs, m = 3, p = 2, ID = "PC_full")
id <- 1

gr(
  ts.plot(cbind(Bayes$Ys[, id], Data$Y_true[, id], Data$Y_obs[, id]), col = c("red", "black", "blue"), lty = c(1, 2, 1))
)

gr(
  hist(Bayes$Bstore[1, 1, ]) # can look at full posterior distributions
)

#--------------- ML Estimation --------------------------------
ML <- mlDFM(Data$Y_obs, m = 3, p = 2, Loud = T)
ts.plot(cbind(ML$Ys[, 1], Data$Y_true[, 1], Data$Y_obs[, 1]), col = c("red", "black", "blue"), lty = c(1, 2, 1))

#------------------ Simulations -------------------------------

# MLE v. Bayes
r <- 50
reps <- 100
MSE_Bayes <- matrix(0, reps, 10)
MSE_ML <- matrix(0, reps, 10)
set.seed(3030)
for (j in 1:reps) {
  Data <- SimDFM(r)
  Bayes <- BDFM(Data$Y_obs, m = 3, p = 2, lam_B = 25)
  ML <- mlDFM(Data$Y_obs, m = 3, p = 2)
  MSE_Bayes[j, ] <- colMeans((Data$Y_true - Bayes$Ys)^2)
  MSE_ML[j, ] <- colMeans((Data$Y_true - ML$Ys)^2)
  print(j)
}

Bayes_ByCol <- colMeans(MSE_Bayes)
ML_ByCol <- colMeans(MSE_ML)
Bayes_total <- mean(MSE_Bayes)
ML_total <- mean(MSE_ML)
