# ============================================================
# CP&CT + Stallard-Todd + UMPCU — BINOMIAL CASE (6 methods, z-test version)
# stage-1/stage-2 p-values use pooled two-proportion
# z-test (standard Phase III practice) instead of Fisher exact.
# ============================================================

suppressPackageStartupMessages(library(mvtnorm))
set.seed(2023)

# ---- Utilities ----
restrictR <- function(numArm, n) {
  sample(rep(0:numArm, n), n * (numArm + 1), replace = FALSE)
}

simRspBinom <- function(n, numArm, pVec) {
  ssn <- n * (numArm + 1)
  out <- matrix(NA, nrow = ssn, ncol = numArm + 1)
  for (i in 0:numArm) out[, i+1] <- rbinom(ssn, 1, pVec[i+1])
  out
}

armCounts <- function(alloc, rsp, numArm) {
  cc <- integer(numArm + 1)
  for (i in 0:numArm) cc[i+1] <- sum(rsp[alloc == i, i+1])
  cc
}

# ---- CHANGED: z-test raw p-values (pooled under H_0k) ----
rawPValuesBinom <- function(alloc, rsp, numArm) {
  idx0 <- which(alloc == 0); n0 <- length(idx0); s0 <- sum(rsp[idx0, 1])
  p <- numeric(numArm)
  for (k in 1:numArm) {
    idx_k <- which(alloc == k); n_k <- length(idx_k); s_k <- sum(rsp[idx_k, k+1])
    p_k_hat <- s_k / n_k; p_0_hat <- s0 / n0
    p_pool <- (s_k + s0) / (n_k + n0)
    se <- sqrt(p_pool * (1 - p_pool) * (1/n_k + 1/n0))
    if (se <= 0) {
      p[k] <- 0.5   # degenerate: both arms all-0 or all-1
    } else {
      z <- (p_k_hat - p_0_hat) / se
      p[k] <- pnorm(z, lower.tail = FALSE)   # one-sided: greater
    }
  }
  p
}

# ---- Intersection tests ----
simesP <- function(pVec) {
  r <- length(pVec); ps <- sort(pVec)
  min(r * ps / seq_len(r))
}

dunnettP_z <- function(pVec) {
  r <- length(pVec)
  if (r == 1) return(pVec[1])
  pc <- pmin(pmax(pVec, 1e-15), 1 - 1e-15)
  z <- qnorm(1 - pc); z_max <- max(z)
  corr <- matrix(0.5, r, r); diag(corr) <- 1
  1 - pmvnorm(upper = rep(z_max, r), corr = corr)[1]
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

combFisher <- function(p1, p2, alpha = 0.025)
  as.integer(pchisq(-2 * log(p1 * p2), df = 4, lower.tail = FALSE) < alpha)

combInvN <- function(p1, p2, w1, w2, alpha = 0.025) {
  z <- w1 * qnorm(1 - p1) + w2 * qnorm(1 - p2)
  as.integer((1 - pnorm(z)) < alpha)
}

# ---- UMPCU (binomial)  ----
tempFuncBinom <- function(z, nA, nB, n0, X2, Tv, d1) {
  x1_range <- max(X2 + 1, z - nB):min(nA, z)
  if (length(x1_range) == 0 || x1_range[1] > x1_range[length(x1_range)]) return(0)
  inner_sum <- sum(choose(nA, x1_range) * choose(nB, z - x1_range))
  choose(n0, Tv - z) * exp(z * d1) * inner_sum
}

calcCN_binom <- function(nA, nB, X2, obsZ, Y0, d1) {
  n0 <- nA + nB; Tv <- Y0 + obsZ
  z_sup <- max(X2 + 1, Tv - n0):min(nA + nB, Tv)
  if (length(z_sup) == 0) return(NA_real_)
  denom <- sum(sapply(z_sup, tempFuncBinom, nA=nA, nB=nB, n0=n0, X2=X2, Tv=Tv, d1=d1))
  if (denom == 0) return(NA_real_)
  1 / denom
}

densityZ_binom <- function(z, nA, nB, X2, obsZ, Y0, d1, CN) {
  CN * tempFuncBinom(z, nA, nB, nA + nB, X2, Y0 + obsZ, d1)
}

umpcuTestBinom <- function(obsZ, nA, nB, X2, Y0, alpha = 0.025) {
  CN <- calcCN_binom(nA, nB, X2, obsZ, Y0, d1 = 0)
  if (is.na(CN)) return(NA_integer_)
  Tv <- Y0 + obsZ; z_up <- min(nA + nB, Tv)
  pval <- sum(sapply(obsZ:z_up, densityZ_binom,
                     nA = nA, nB = nB, X2 = X2,
                     obsZ = obsZ, Y0 = Y0, d1 = 0, CN = CN))
  as.integer(pval < alpha)
}

# ---- Stallard-Todd utilities ----
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

# ---- One-trial simulation (6 methods) ----
compFuncAllBinom <- function(numArm, nA, nB, pVec, u_alpha_ST, alpha = 0.025) {
  rspA <- simRspBinom(nA, numArm, pVec)
  allocA <- restrictR(numArm, nA)
  X <- armCounts(allocA, rspA, numArm)
  ord <- order(X[-1], decreasing = TRUE)
  M <- ord[1]; M2 <- ord[2]; X2 <- X[M2 + 1]
  
  rspB <- simRspBinom(nB, 1, pVec = c(pVec[1], pVec[M + 1]))
  allocB <- restrictR(1, nB)
  
  # CP&CT p-values (z-test)
  p_raw <- rawPValuesBinom(allocA, rspA, numArm)
  p_raw <- pmin(pmax(p_raw, 1e-15), 1 - 1e-15)
  
  p1_simes   <- closureAdj(p_raw, M, numArm, simesP)
  p1_dunnett <- closureAdj(p_raw, M, numArm, dunnettP_z)
  p1_simes   <- min(max(p1_simes,   1e-15), 1 - 1e-15)
  p1_dunnett <- min(max(p1_dunnett, 1e-15), 1 - 1e-15)
  
  p2 <- rawPValuesBinom(allocB, rspB, 1)[1]
  p2 <- min(max(p2, 1e-15), 1 - 1e-15)
  
  w1 <- sqrt(nA / (nA + nB)); w2 <- sqrt(nB / (nA + nB))
  rej_FS <- combFisher(p1_simes,   p2, alpha)
  rej_FD <- combFisher(p1_dunnett, p2, alpha)
  rej_NS <- combInvN(p1_simes,   p2, w1, w2, alpha)
  rej_ND <- combInvN(p1_dunnett, p2, w1, w2, alpha)
  
  # UMPCU (with tie check)
  sortX <- sort(X[-1], decreasing = TRUE)
  obsZ <- X[M + 1] + sum(rspB[which(allocB == 1), 2])
  Y0   <- X[1]     + sum(rspB[which(allocB == 0), 1])
  if (abs(sortX[1] - sortX[2]) <= 1) {
    rej_U <- NA_integer_
  } else {
    rej_U <- umpcuTestBinom(obsZ, nA, nB, X2, Y0, alpha)
  }
  
  # Stallard-Todd
  X_total_s1 <- sum(X); n_total_s1 <- (numArm + 1) * nA
  phat_pool_s1 <- X_total_s1 / n_total_s1
  V1_ref <- (nA / 2) * phat_pool_s1 * (1 - phat_pool_s1)
  Z_M_2 <- 0.5 * (obsZ - Y0)
  if (V1_ref <= 0) {
    rej_ST <- NA_integer_
  } else {
    rej_ST <- as.integer(Z_M_2 / sqrt(V1_ref) > u_alpha_ST)
  }
  
  c(UMPCU          = rej_U,
    Fisher_Simes   = rej_FS,
    Fisher_Dunnett = rej_FD,
    InvN_Simes     = rej_NS,
    InvN_Dunnett   = rej_ND,
    StallardTodd   = rej_ST)
}

# ---- Reporter ----
reportSim6 <- function(tests, nsim) {
  u <- tests[, 1]; u_valid <- u[!is.na(u)]
  rate_u <- mean(u_valid); se_u <- sqrt(rate_u*(1-rate_u)/length(u_valid))
  pct_na <- mean(is.na(u))
  
  others <- colMeans(tests[, 2:6], na.rm = TRUE)
  nvalid <- apply(tests[, 2:6], 2, function(x) sum(!is.na(x)))
  se_o <- sqrt(others*(1-others)/nvalid)
  
  out <- rbind(Rate   = round(c(rate_u, others), 4),
               SE     = round(c(se_u, se_o), 4),
               NA_pct = round(c(pct_na, rep(NA, 5)), 4))
  colnames(out) <- c("UMPCU","F-Simes","F-Dunnett","InvN-S","InvN-D","ST")
  out
}

# ============================================================
# PRE-COMPUTE ST CRITICAL VALUES
# ============================================================
nA <- 100
alpha <- 0.025
nsim <- 10000

nB_vec <- c(300, 100)
K_vec_total <- 3:7
numArm_vec <- K_vec_total - 1

cat("Pre-computing Stallard-Todd critical values...\n")
ST_cache <- list()
for (nB_val in nB_vec) {
  r_val <- nB_val / nA
  for (numArm in numArm_vec) {
    key <- sprintf("K=%d,nB=%d", numArm, nB_val)
    t0 <- Sys.time()
    u <- find_u_alpha(alpha, K = numArm, r = r_val)
    elapsed <- as.numeric(Sys.time() - t0, units = "secs")
    ST_cache[[key]] <- u
    cat(sprintf("  numArm=%d, nB=%d (r=%.2f): u=%.4f (%.1fs)\n",
                numArm, nB_val, r_val, u, elapsed))
  }
}
cat("\n")

run_sim <- function(numArm, nA, nB_val, pVec, nsim) {
  key <- sprintf("K=%d,nB=%d", numArm, nB_val)
  u_alpha_ST <- ST_cache[[key]]
  tests <- matrix(NA, nsim, 6)
  for (i in 1:nsim) {
    tests[i, ] <- compFuncAllBinom(numArm, nA, nB_val, pVec, u_alpha_ST, alpha)
  }
  tests
}

cat("MC SE at alpha=0.025:", round(sqrt(0.025*0.975/nsim), 4), "\n")
cat("MC SE at power=0.70: ", round(sqrt(0.70*0.30/nsim), 4), "\n\n")

# ============================================================
# TABLE 8: Type I Error vs K
# ============================================================
cat("=== TABLE 8: Type I Error (binomial, z-test CP&CT), pi=0.55 ===\n")
for (nB_val in nB_vec) {
  cat(sprintf("\n--- nA=%d, nB=%d ---\n", nA, nB_val))
  for (K in K_vec_total) {
    pVec <- rep(0.55, K + 1)
    tests <- run_sim(K-1, nA, nB_val, pVec, nsim)
    cat("K =", K, "arms:\n"); print(reportSim6(tests, nsim)); cat("\n")
  }
}

# ============================================================
# TABLE 9: Power vs K
# ============================================================
cat("=== TABLE 9: Power vs K (binomial, z-test CP&CT) ===\n")
pi_full <- c(0.55, 0.65, 0.50, 0.32, 0.51, 0.42, 0.21)
for (nB_val in nB_vec) {
  cat(sprintf("\n--- nA=%d, nB=%d ---\n", nA, nB_val))
  for (K in K_vec_total) {
    pVec <- pi_full[1:(K+1)]
    tests <- run_sim(K-1, nA, nB_val, pVec, nsim)
    cat("K =", K, "arms:\n"); print(reportSim6(tests, nsim)); cat("\n")
  }
}

# ============================================================
# TABLE 10: Power vs pi1-pi2
# ============================================================
cat("=== TABLE 10: Power vs pi1-pi2 (binomial, z-test CP&CT) ===\n")
pi0 <- 0.55; pi1 <- 0.65
pi2_vec <- c(0.60, 0.55, 0.50, 0.45, 0.40)
for (nB_val in nB_vec) {
  cat(sprintf("\n--- nA=%d, nB=%d ---\n", nA, nB_val))
  for (p2v in pi2_vec) {
    pVec <- c(pi0, pi1, p2v)
    tests <- run_sim(2, nA, nB_val, pVec, nsim)
    cat("pi2 =", p2v, ":\n"); print(reportSim6(tests, nsim)); cat("\n")
  }
}

cat("=== Binomial 6-method simulations complete ===\n")
