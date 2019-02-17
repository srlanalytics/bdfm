library(bdfm)

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

# Error in EstDFM(B = B_in, Bp = Bp, Jb = Jb, lam_B = lam_B, q = q, nu_q = nu_q,  :
#   c++ exception (unknown reason)

# - [ ] What is wrong?
# - [ ] Better Error that tells me whats wrong

# Are we drawing from a stationary distribution?
# ts.plot(m$Qstore[1,1,])
# ts.plot(m$Hstore[1,1,])

# summary(m)


