# ============================================================
# Unequal Allocation Sensitivity (Normal Case)
# Addresses Reviewer 2, Comment 7
#
# UMPCU can be used with unequal stage-2 allocation: the Sampson-Sill
# derivation separates nA (stage-1 per arm), nB (stage-2 on selected
# arm M), and nB0 (stage-2 on control). The conditional density
# f_Q(W | Delta_M, X*, T) uses n0 = nA + nB0 and nB separately.
#
# Stage 1 keeps equal allocation across K+1 arms.
# Stage 2 varied:1:1(baseline), 2:1 (treat-heavy), 1:2 (ctrl-heavy)
# with total N_B = 200.
# ============================================================

suppressPackageStartupMessages(library(mvtnorm))
set.seed(2025)

restrictR <- function(numArm, n) sample(rep(0:numArm, n), n*(numArm+1), replace=FALSE)
restrictR_twoArm <- function(n_treat, n_ctrl) sample(c(rep(1,n_treat), rep(0,n_ctrl)))

simRsp <- function(n, numArm, mu, sigma) {
  ssn <- n*(numArm+1); out <- matrix(NA, ssn, numArm+1)
  for (i in 0:numArm) out[,i+1] <- rnorm(ssn, mu[i+1], sigma); out
}

simRspStage2 <- function(n_treat, n_ctrl, mu_treat, mu_ctrl, sigma) {
  ssn <- n_treat + n_ctrl
  out <- matrix(NA, ssn, 2)
  out[, 1] <- rnorm(ssn, mu_ctrl,  sigma)
  out[, 2] <- rnorm(ssn, mu_treat, sigma)
  out
}

armMeans <- function(alloc, rsp, numArm) {
  m <- numeric(numArm+1)
  for (i in 0:numArm) m[i+1] <- mean(rsp[alloc==i, i+1]); m
}

rawPValues <- function(alloc, rsp, numArm) {
  x <- numeric(numArm+1); v <- numeric(numArm+1); nk <- integer(numArm+1)
  for (i in 0:numArm) { idx <- which(alloc==i)
  x[i+1] <- mean(rsp[idx,i+1]); v[i+1] <- var(rsp[idx,i+1]); nk[i+1] <- length(idx) }
  p <- numeric(numArm)
  for (j in 1:numArm) {
    se <- sqrt(v[j+1]/nk[j+1] + v[1]/nk[1])
    t <- (x[j+1] - x[1])/se
    df <- se^4/(v[j+1]^2/(nk[j+1]^2*(nk[j+1]-1)) + v[1]^2/(nk[1]^2*(nk[1]-1)))
    p[j] <- pt(t, df, lower.tail=FALSE)
  }; p
}

simesP <- function(pVec) { r <- length(pVec); min(r*sort(pVec)/seq_len(r)) }

dunnettP <- function(pVec, df = Inf) {
  r <- length(pVec); if (r == 1) return(pVec[1])
  pc <- pmin(pmax(pVec, 1e-15), 1-1e-15)
  z_max <- max(qnorm(1 - pc))
  corr <- matrix(0.5, r, r); diag(corr) <- 1
  if (is.infinite(df)) 1 - pmvnorm(upper=rep(z_max,r), corr=corr)[1]
  else                1 - pmvt(upper=rep(z_max,r), corr=corr, df=df)[1]
}

closureAdj <- function(p_raw, M, numArm, intersect_fn, ...) {
  others <- setdiff(1:numArm, M); n_oth <- length(others)
  if (n_oth == 0) return(p_raw[M])
  p_max <- 0
  for (mask in 0:(2^n_oth - 1)) {
    bits <- as.integer(intToBits(mask))[1:n_oth]
    I <- c(M, others[bits == 1])
    p_I <- intersect_fn(p_raw[I], ...)
    if (p_I > p_max) p_max <- p_I
  }; p_max
}

combFisher <- function(p1,p2,alpha=0.025)
  as.integer(pchisq(-2*log(p1*p2), df=4, lower.tail=FALSE) < alpha)
combInvN <- function(p1,p2,w1,w2,alpha=0.025)
  as.integer((1 - pnorm(w1*qnorm(1-p1) + w2*qnorm(1-p2))) < alpha)

# UMPCU with unequal allocation
dstyQ_u <- function(w, d1, nA, nB, nB0, sigma, Y0bar, Zbar, X2) {
  n0 <- nA + nB0
  mu_w <- (n0 * (nA + nB) * d1) / ((n0 + nA + nB) * sigma^2)
  s2<- (n0 * (nA + nB))/ ((n0 + nA + nB) * sigma^2)
  Tv   <- (n0 * Y0bar + (nA + nB) * Zbar) / (n0 + nA + nB)
  v    <- (sqrt(nA * (nA + nB)) * (sigma^2 * w + (nA + nB) * (Tv - X2))) /(sqrt(nB) * (nA + nB) * sigma)
  sqrt(2*pi*s2) * dnorm(w, mu_w, sqrt(s2)) * pnorm(v)
}

umpcuTest_u <- function(W, nA, nB, nB0, sigma, Y0bar, Zbar, X2, alpha=0.025) {
  CN <- tryCatch(1/integrate(dstyQ_u, 0, nA, nB, nB0, sigma, Y0bar, Zbar, X2,
                             lower=-Inf, upper=Inf)$value,
                 error=function(e) NA_real_)
  if (is.na(CN)) return(NA_integer_)
  pval <- tryCatch(integrate(dstyQ_u, 0, nA, nB, nB0, sigma, Y0bar, Zbar, X2,
                             lower=W, upper=Inf)$value*CN,
                   error=function(e) NA_real_)
  if (is.na(pval)) return(NA_integer_)
  as.integer(pval < alpha)
}

compFunc_unequal <- function(numArm, nA, nB_treat, nB_ctrl, mu, sigma, alpha=0.025) {
  rspA <- simRsp(nA, numArm, mu, sigma); allocA <- restrictR(numArm, nA)
  X <- armMeans(allocA, rspA, numArm)
  ord <- order(X[-1], decreasing=TRUE); M <- ord[1]; X2 <- X[ord[2]+1]
  p_raw <- rawPValues(allocA, rspA, numArm); df1 <- 2*nA - 2
  
  rspB <- simRspStage2(nB_treat, nB_ctrl, mu[M+1], mu[1], sigma)
  allocB <- restrictR_twoArm(nB_treat, nB_ctrl)
  YM <- mean(rspB[allocB==1, 2]); Y0 <- mean(rspB[allocB==0, 1])
  v1 <- var(rspB[allocB==1, 2]); v0 <- var(rspB[allocB==0, 1])
  se2 <- sqrt(v1/nB_treat + v0/nB_ctrl)
  df2 <- se2^4 / (v1^2/(nB_treat^2*(nB_treat-1)) + v0^2/(nB_ctrl^2*(nB_ctrl-1)))
  p2 <- pt((YM - Y0)/se2, df2, lower.tail=FALSE)
  
  p1S <- closureAdj(p_raw, M, numArm, simesP)
  p1D <- closureAdj(p_raw, M, numArm, dunnettP, df=df1)
  eps <- 1e-15
  p1S <- min(max(p1S,eps),1-eps); p1D <- min(max(p1D,eps),1-eps)
  p2c <- min(max(p2, eps),1-eps)
  
  w1 <- sqrt(nA / (nA + min(nB_treat, nB_ctrl)))
  w2 <- sqrt(min(nB_treat, nB_ctrl) / (nA + min(nB_treat, nB_ctrl)))
  
  Zbar<- (nA*X[M+1] + nB_treat*YM) / (nA + nB_treat)
  Y0bar <- (nA*X[1]+ nB_ctrl*Y0)/ (nA + nB_ctrl)
  n0 <- nA + nB_ctrl
  W <- (n0 * (nA + nB_treat)) / ((n0 + nA + nB_treat) * sigma^2) * (Zbar - Y0bar)
  
  c(UMPCU          = umpcuTest_u(W, nA, nB_treat, nB_ctrl, sigma, Y0bar, Zbar, X2, alpha),
    Fisher_Simes   = combFisher(p1S, p2c, alpha),
    Fisher_Dunnett = combFisher(p1D, p2c, alpha),
    InvN_Simes     = combInvN(p1S, p2c, w1, w2, alpha),
    InvN_Dunnett   = combInvN(p1D, p2c, w1, w2, alpha))
}

# ============================================================
nsim <- 5000
sigma <- sqrt(10); alpha <- 0.025; K <- 5; nA <- 100
mu_large <- c(0.70, 1.80, 1.10, 0.43, 0.82, 0.13)
mu_small <- c(0.70, 1.80, 1.75, 0.43, 0.82, 0.13)

settings <- list(
  `1:1 (baseline)`    = c(100, 100),
  `2:1 (treat-heavy)` = c(133, 67),
  `1:2 (ctrl-heavy)`  = c(67, 133)
)

cat("\n=== T1E: mu=0, K=5 ===\n")
for (nm in names(settings)) {
  nB_t <- settings[[nm]][1]; nB_c <- settings[[nm]][2]
  tests <- matrix(NA, nsim, 5)
  for (i in 1:nsim) tests[i,] <- compFunc_unequal(K, nA, nB_t, nB_c, rep(0,6), sigma, alpha)
  rates <- colMeans(tests, na.rm=TRUE)
  cat(sprintf("%-20s nB_t=%d nB_c=%d: ", nm, nB_t, nB_c))
  cat(sprintf("U=%.4f FS=%.4f FD=%.4f NS=%.4f ND=%.4f\n",
              rates[1], rates[2], rates[3], rates[4], rates[5]))
}

cat("\n=== Power: large mu1-mu2, K=5 ===\n")
for (nm in names(settings)) {
  nB_t <- settings[[nm]][1]; nB_c <- settings[[nm]][2]
  tests <- matrix(NA, nsim, 5)
  for (i in 1:nsim) tests[i,] <- compFunc_unequal(K, nA, nB_t, nB_c, mu_large, sigma, alpha)
  rates <- colMeans(tests, na.rm=TRUE)
  cat(sprintf("%-20s nB_t=%d nB_c=%d: ", nm, nB_t, nB_c))
  cat(sprintf("U=%.4f FS=%.4f FD=%.4f NS=%.4f ND=%.4f\n",
              rates[1], rates[2], rates[3], rates[4], rates[5]))
}

cat("\n=== Power: small mu1-mu2, K=5 ===\n")
for (nm in names(settings)) {
  nB_t <- settings[[nm]][1]; nB_c <- settings[[nm]][2]
  tests <- matrix(NA, nsim, 5)
  for (i in 1:nsim) tests[i,] <- compFunc_unequal(K, nA, nB_t, nB_c, mu_small, sigma, alpha)
  rates <- colMeans(tests, na.rm=TRUE)
  cat(sprintf("%-20s nB_t=%d nB_c=%d: ", nm, nB_t, nB_c))
  cat(sprintf("U=%.4f FS=%.4f FD=%.4f NS=%.4f ND=%.4f\n",
              rates[1], rates[2], rates[3], rates[4], rates[5]))
}

cat("\n=== Unequal allocation simulations complete ===\n")
