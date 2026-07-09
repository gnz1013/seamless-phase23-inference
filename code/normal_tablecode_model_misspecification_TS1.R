# ============================================================
# Misspecification Type I Error Check (R1-C2)
# Three distributions, global null, K=5, n_A=100, nsim=10000
# ============================================================

suppressPackageStartupMessages(library(mvtnorm))

# ---- Distribution generators (all standardized to mean=0, var=sigma^2) ----
gen_normal <- function(n, mean, sd) {
  rnorm(n, mean, sd)
}

gen_t4_scaled <- function(n, mean, sd) {
  # t(df=4) has variance df/(df-2) = 2; rescale to target sd
  z <- rt(n, df = 4) / sqrt(4 / (4 - 2))
  z * sd + mean
}

gen_lnorm_shifted <- function(n, mean, sd, sigma_log = 0.5) {
  # Lognormal: scale & shift to mean=0, sd=1, then scale to target
  z <- rnorm(n)
  x <- exp(sigma_log * z)
  mean_x <- exp(sigma_log^2 / 2)
  var_x <- (exp(sigma_log^2) - 1) * exp(sigma_log^2)
  x_std <- (x - mean_x) / sqrt(var_x)
  x_std * sd + mean
}

# ---- Utilities (paper's original, with dist_fn parameterization) ----
restrictR <- function(numArm, n) sample(rep(0:numArm, n), n*(numArm+1), replace=FALSE)

simRsp <- function(n, numArm, mu, sigma, dist_fn) {
  ssn <- n*(numArm+1); out <- matrix(NA, ssn, numArm+1)
  for (i in 0:numArm) out[,i+1] <- dist_fn(ssn, mu[i+1], sigma)
  out
}

armMeans <- function(alloc, rsp, numArm) {
  m <- numeric(numArm+1)
  for (i in 0:numArm) m[i+1] <- mean(rsp[alloc==i, i+1])
  m
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
  }
  p
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
  }
  p_max
}

combFisher <- function(p1, p2, alpha=0.025)
  as.integer(pchisq(-2*log(p1*p2), df=4, lower.tail=FALSE) < alpha)

combInvN <- function(p1, p2, w1, w2, alpha=0.025)
  as.integer((1 - pnorm(w1*qnorm(1-p1) + w2*qnorm(1-p2))) < alpha)

dstyQ <- function(w, d1, nA, nB, sigma, Y0bar, Zbar, X2) {
  n0 <- nA + nB
  mu_w <- (n0*(nA+nB)*d1)/((n0+nA+nB)*sigma^2)
  s2 <- (n0*(nA+nB))/((n0+nA+nB)*sigma^2)
  Tv <- (n0*Y0bar + (nA+nB)*Zbar)/(n0+nA+nB)
  v <- (sqrt(nA*(nA+nB))*(sigma^2*w + (nA+nB)*(Tv-X2)))/(sqrt(nB)*(nA+nB)*sigma)
  sqrt(2*pi*s2)*dnorm(w, mu_w, sqrt(s2))*pnorm(v)
}

umpcuTest <- function(W, nA, nB, sigma, Y0bar, Zbar, X2, alpha=0.025) {
  CN <- tryCatch(1/integrate(dstyQ, 0, nA, nB, sigma, Y0bar, Zbar, X2,
                             lower=-Inf, upper=Inf)$value, error=function(e) NA_real_)
  if (is.na(CN)) return(NA_integer_)
  pval <- tryCatch(integrate(dstyQ, 0, nA, nB, sigma, Y0bar, Zbar, X2,
                             lower=W, upper=Inf)$value*CN,
                   error=function(e) NA_real_)
  if (is.na(pval)) return(NA_integer_)
  as.integer(pval < alpha)
}

# ---- Stallard-Todd ----
f_max_equicorr <- function(m, K, rho = 0.5) {
  integrate(function(w) {
    dnorm(w) * K * dnorm((m - sqrt(rho)*w)/sqrt(1-rho)) / sqrt(1-rho) *
      pnorm((m - sqrt(rho)*w)/sqrt(1-rho))^(K-1)
  }, -Inf, Inf, rel.tol = 1e-6)$value
}

P_upper <- function(u, K, r) {
  integrand <- function(m) {
    sapply(m, function(mi) {
      pnorm((u - mi)/sqrt(r), lower.tail = FALSE) * f_max_equicorr(mi, K, rho = 0.5)
    })
  }
  integrate(integrand, u - 5*sqrt(r), 10, rel.tol = 1e-5)$value
}

find_u_alpha <- function(alpha, K, r, u_lo = 0.5, u_hi = 6) {
  f <- function(u) P_upper(u, K, r) - alpha
  uniroot(f, lower = u_lo, upper = u_hi, tol = 1e-5)$root
}

score_normal <- function(xbar_k, xbar_0, n_k, n_0, sigma) {
  V <- (n_k * n_0) / ((n_k + n_0) * sigma^2)
  Z <- V * (xbar_k - xbar_0)
  list(Z = Z, V = V)
}

# ---- 6-method simulation (with distribution parameter) ----
compFuncAll <- function(numArm, nA, nB, mu, sigma, dist_fn, u_alpha_ST, alpha=0.025) {
  rspA <- simRsp(nA, numArm, mu, sigma, dist_fn); allocA <- restrictR(numArm, nA)
  X <- armMeans(allocA, rspA, numArm)
  ord <- order(X[-1], decreasing=TRUE); M <- ord[1]; X2 <- X[ord[2]+1]
  p_raw <- rawPValues(allocA, rspA, numArm); df1 <- 2*nA-2
  
  rspB <- simRsp(nB, 1, mu=c(mu[1], mu[M+1]), sigma, dist_fn); allocB <- restrictR(1, nB)
  YM <- mean(rspB[allocB==1, 2]); Y0 <- mean(rspB[allocB==0, 1])
  v1 <- var(rspB[allocB==1, 2]); v0 <- var(rspB[allocB==0, 1])
  n1 <- sum(allocB==1); n0s <- sum(allocB==0)
  se2 <- sqrt(v1/n1 + v0/n0s)
  df2 <- se2^4/(v1^2/(n1^2*(n1-1)) + v0^2/(n0s^2*(n0s-1)))
  p2 <- pt((YM-Y0)/se2, df2, lower.tail=FALSE)
  
  p1S <- closureAdj(p_raw, M, numArm, simesP)
  p1D <- closureAdj(p_raw, M, numArm, dunnettP, df=df1)
  eps <- 1e-15
  p1S <- min(max(p1S, eps), 1-eps); p1D <- min(max(p1D, eps), 1-eps)
  p2c <- min(max(p2, eps), 1-eps)
  w1 <- sqrt(nA/(nA+nB)); w2 <- sqrt(nB/(nA+nB))
  Zbar <- (nA*X[M+1] + nB*YM)/(nA+nB); Y0bar <- (nA*X[1] + nB*Y0)/(nA+nB)
  W <- (nA+nB)/(2*sigma^2)*(Zbar - Y0bar)
  
  # Stallard-Todd
  n_M_tot <- nA + n1; n_0_tot <- nA + n0s
  xbar_M_comb <- (nA * X[M+1] + n1 * YM) / n_M_tot
  xbar_0_comb <- (nA * X[1] + n0s * Y0) / n_0_tot
  out2 <- score_normal(xbar_M_comb, xbar_0_comb, n_M_tot, n_0_tot, sigma)
  V1_ref <- nA / (2 * sigma^2)
  Y_std <- out2$Z / sqrt(V1_ref)
  rej_ST <- as.integer(Y_std > u_alpha_ST)
  
  c(UMPCU          = umpcuTest(W, nA, nB, sigma, Y0bar, Zbar, X2, alpha),
    Fisher_Simes   = combFisher(p1S, p2c, alpha),
    Fisher_Dunnett = combFisher(p1D, p2c, alpha),
    InvN_Simes     = combInvN(p1S, p2c, w1, w2, alpha),
    InvN_Dunnett   = combInvN(p1D, p2c, w1, w2, alpha),
    StallardTodd   = rej_ST)
}

reportSim <- function(tests, nsim) {
  rates <- colMeans(tests, na.rm=TRUE)
  nvalid <- apply(tests, 2, function(x) sum(!is.na(x)))
  se <- sqrt(rates*(1-rates)/nvalid)
  out <- rbind(Rate = round(rates, 4), SE = round(se, 4))
  colnames(out) <- c("UMPCU","F-Simes","F-Dunnett","InvN-S","InvN-D","ST")
  out
}

# ============================================================
# SETTINGS
# ============================================================
set.seed(2023)
nsim <- 10000
nA <- 100; nB <- 100; sigma <- sqrt(10); alpha <- 0.025
K_arms <- 5

cat("nsim =", nsim, "\n")
cat("nA=nB=", nA, ", sigma^2=10, K=5 exp arms\n")
cat("All distributions standardized to mean=0, var=10 (per arm)\n\n")

# ---- Pre-compute ST critical value (K=5, r=1) ----
cat("Pre-computing Stallard-Todd critical value for K=5...\n")
u_ST <- find_u_alpha(alpha, K = K_arms, r = 1)
cat(sprintf("  u_alpha = %.4f\n\n", u_ST))

# ============================================================
# MISSPECIFICATION CHECK: T1E under global null
# ============================================================
mu <- rep(0, K_arms + 1)

distributions <- list(
  "Normal (baseline)"       = gen_normal,
  "t(df=4) heavy-tailed"    = gen_t4_scaled,
  "Lognormal (σ_log=0.5)"   = gen_lnorm_shifted
)

cat("=== Misspecification T1E check ===\n")
cat("Global null: mu_0 = mu_1 = ... = mu_K = 0\n\n")

results <- list()
for (dist_name in names(distributions)) {
  dist_fn <- distributions[[dist_name]]
  cat(sprintf("--- %s ---\n", dist_name))
  
  tests <- matrix(NA, nsim, 6)
  for (i in 1:nsim) {
    tests[i,] <- compFuncAll(K_arms, nA, nB, mu, sigma, dist_fn, u_ST, alpha)
  }
  
  r <- reportSim(tests, nsim)
  results[[dist_name]] <- r
  print(r); cat("\n")
}

# ---- Summary table ----
cat("\n=== SUMMARY: T1E across distributions ===\n")
summary_mat <- sapply(results, function(r) r["Rate",])
print(round(summary_mat, 4))
cat("\n")

cat("MC SE at alpha=0.025:", round(sqrt(0.025*0.975/nsim), 4), "\n")
cat("=== Complete ===\n")
