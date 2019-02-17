
# drop outliers (optional but gets rid of some wierd stuff)
econ_us[abs(scale(econ_us)) > 4] <- NA

logs <- c(
  "W068RCQ027SBEA",
  "PCEDG",
  "PCEND",
  "JTSJOL",
  "INDPRO",
  "CSUSHPINSA",
  "HSN1F",
  "TSIFRGHT",
  "IPG2211S",
  "DGORDER",
  "AMTMNO",
  "CPILFESL",
  "ICSA"
)

diffs <- setdiff(colnames(econ_us), c("A191RL1Q225SBEA", 'W068RCQ027SBEA', "USSLIND"))

m <- dfm(econ_us, factors = 3, pre_differenced = "A191RL1Q225SBEA", logs = logs, diffs = diffs)

# warning: eig_gen(): decomposition failed

# error: max(): object has no elements
# Error in EstDFM(B = B_in, Bp = Bp, Jb = Jb, lam_B = lam_B, q = q, nu_q = nu_q,  :
#   max(): object has no elements
# In addition: Warning message:
# In bdfm(Y = Y, m = m, p = p, Bp = Bp, lam_B = lam_B, Hp = Hp, lam_H = lam_H,  :
#   PC_full not a valid identification string or index vector, defaulting to pc_full
# >

# - [ ] What is wrong?
# - [ ] Better Error that tells me whats wrong

# Are we drawing from a stationary distribution?
# ts.plot(m$Qstore[1,1,])
# ts.plot(m$Hstore[1,1,])

summary(m)


