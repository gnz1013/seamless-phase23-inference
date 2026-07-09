# ============================================================
# Heatmap: UMPCU vs InvN-Dunnett Power Difference -- BINOMIAL
# K=5, aligned with R1-C3 / R2-C1: pooled two-proportion z-test
# ============================================================

suppressPackageStartupMessages({
  library(mvtnorm)
  library(ggplot2)
  library(scales)
})
set.seed(2027)

# ---- Utilities ----
restrictR <- function(numArm, n) sample(rep(0:numArm, n), n*(numArm+1), replace=FALSE)

simRspBinom <- function(n, numArm, pVec) {
  ssn <- n*(numArm+1); out <- matrix(NA, ssn, numArm+1)
  for (i in 0:numArm) out[,i+1] <- rbinom(ssn, 1, pVec[i+1]); out
}

armCounts <- function(alloc, rsp, numArm) {
  cc <- integer(numArm+1)
  for (i in 0:numArm) cc[i+1] <- sum(rsp[alloc==i, i+1]); cc
}

# Pooled two-proportion z-test, one-sided (H1: p_k > p_0) -- replaces Fisher
rawPValuesBinom <- function(alloc, rsp, numArm) {
  idx0 <- which(alloc==0); n0 <- length(idx0); s0 <- sum(rsp[idx0,1])
  p <- numeric(numArm)
  for (k in 1:numArm) {
    idx_k <- which(alloc==k); n_k <- length(idx_k); s_k <- sum(rsp[idx_k,k+1])
    p_pool <- (s_k + s0) / (n_k + n0)
    if (p_pool <= 0 || p_pool >= 1) { p[k] <- 1; next }
    se <- sqrt(p_pool * (1-p_pool) * (1/n_k + 1/n0))
    z  <- (s_k/n_k - s0/n0) / se
    p[k] <- pnorm(z, lower.tail=FALSE)
  }
  p
}

dunnettP_z <- function(pVec) {
  r <- length(pVec); if (r==1) return(pVec[1])
  pc <- pmin(pmax(pVec, 1e-15), 1-1e-15)
  z_max <- max(qnorm(1-pc))
  corr <- matrix(0.5, r, r); diag(corr) <- 1
  1 - pmvnorm(upper=rep(z_max, r), corr=corr)[1]
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

combInvN <- function(p1, p2, w1, w2, alpha=0.025)
  as.integer((1 - pnorm(w1*qnorm(1-p1) + w2*qnorm(1-p2))) < alpha)

# ---- UMPCU binomial (Sill 2009) ----
tempFuncBinom <- function(z, nA, nB, n0, X2, Tv, d1) {
  x1_range <- max(X2+1, z-nB):min(nA, z)
  if (length(x1_range)==0 || x1_range[1] > x1_range[length(x1_range)]) return(0)
  choose(n0, Tv-z) * exp(z*d1) *
    sum(choose(nA, x1_range) * choose(nB, z-x1_range))
}

calcCN_binom <- function(nA, nB, X2, obsZ, Y0, d1) {
  n0 <- nA+nB; Tv <- Y0+obsZ
  z_sup <- max(X2+1, Tv-n0):min(nA+nB, Tv)
  if (length(z_sup)==0) return(NA_real_)
  denom <- sum(sapply(z_sup, tempFuncBinom, nA=nA, nB=nB, n0=n0, X2=X2, Tv=Tv, d1=d1))
  if (denom==0) return(NA_real_)
  1/denom
}

umpcuTestBinom <- function(obsZ, nA, nB, X2, Y0, alpha=0.025) {
  CN <- calcCN_binom(nA, nB, X2, obsZ, Y0, d1=0)
  if (is.na(CN)) return(NA_integer_)
  Tv <- Y0+obsZ; z_up <- min(nA+nB, Tv)
  pval <- sum(sapply(obsZ:z_up, function(z) {
    CN * tempFuncBinom(z, nA, nB, nA+nB, X2, Tv, d1=0)
  }))
  as.integer(pval < alpha)
}

# ---- Single trial: returns UMPCU and InvN-Dunnett rejection ----
compFunc2_binom <- function(numArm, nA, nB, pVec, alpha=0.025) {
  rspA <- simRspBinom(nA, numArm, pVec); allocA <- restrictR(numArm, nA)
  X <- armCounts(allocA, rspA, numArm)
  ord <- order(X[-1], decreasing=TRUE); M <- ord[1]; M2 <- ord[2]; X2 <- X[M2+1]

  rspB <- simRspBinom(nB, 1, pVec=c(pVec[1], pVec[M+1]))
  allocB <- restrictR(1, nB)

  p_raw <- rawPValuesBinom(allocA, rspA, numArm)
  p_raw <- pmin(pmax(p_raw, 1e-15), 1-1e-15)
  p1D   <- closureAdj(p_raw, M, numArm, dunnettP_z)
  p1D   <- min(max(p1D, 1e-15), 1-1e-15)
  p2    <- rawPValuesBinom(allocB, rspB, 1)[1]
  p2    <- min(max(p2, 1e-15), 1-1e-15)

  w1 <- sqrt(nA/(nA+nB)); w2 <- sqrt(nB/(nA+nB))
  rej_ND <- combInvN(p1D, p2, w1, w2, alpha)

  sortX <- sort(X[-1], decreasing=TRUE)
  if (abs(sortX[1] - sortX[2]) <= 1) {
    rej_U <- NA_integer_
  } else {
    obsZ <- X[M+1] + sum(rspB[which(allocB==1), 2])
    Y0   <- X[1]   + sum(rspB[which(allocB==0), 1])
    rej_U <- umpcuTestBinom(obsZ, nA, nB, X2, Y0, alpha)
  }
  c(UMPCU=rej_U, InvN_D=rej_ND)
}

# ============================================================
# GRID: K=5, pi0=0.40 fixed, pi1=pi0+delta, pi2=pi1-gap
# other 3 arms at pi0 level minus small decrements
# ============================================================
nsim <- 2000
nA <- 100; nB <- 100; alpha <- 0.025; K <- 5
pi0 <- 0.40

gap_vec   <- c(0.00, 0.05, 0.10, 0.15, 0.20, 0.25)
delta_vec <- c(0.05, 0.10, 0.15, 0.20, 0.25)

results <- expand.grid(gap = gap_vec, delta = delta_vec)
results$UMPCU  <- NA_real_
results$InvN_D <- NA_real_
results$NA_pct <- NA_real_
results$diff   <- NA_real_

t0 <- Sys.time()
for (r in seq_len(nrow(results))) {
  gap <- results$gap[r]; delta <- results$delta[r]
  pi1 <- pi0 + delta
  pi2 <- pi1 - gap
  pi2 <- max(min(pi2, 0.99), 0.01)
  pi_others <- c(pi0 - 0.05, pi0 - 0.10, pi0 - 0.15)
  pi_others <- pmax(pi_others, 0.05)
  pVec <- c(pi0, pi1, pi2, pi_others)

  tests <- matrix(NA, nsim, 2)
  for (i in 1:nsim) tests[i,] <- compFunc2_binom(K, nA, nB, pVec, alpha)

  u_valid <- tests[!is.na(tests[,1]), 1]
  results$UMPCU[r]  <- if (length(u_valid)) mean(u_valid) else NA_real_
  results$InvN_D[r] <- mean(tests[,2])
  results$NA_pct[r] <- mean(is.na(tests[,1]))
  results$diff[r]   <- results$UMPCU[r] - results$InvN_D[r]

  elapsed <- as.numeric(Sys.time() - t0, units="mins")
  cat(sprintf("[%d/%d] %.1f min | gap=%.2f delta=%.2f | U=%.3f ND=%.3f diff=%+.3f NA=%.1f%%\n",
              r, nrow(results), elapsed, gap, delta,
              results$UMPCU[r], results$InvN_D[r], results$diff[r],
              results$NA_pct[r]*100))
}

write.csv(results, "~/Desktop/heatmap_binomial_K5_grid.csv", row.names=FALSE)

# ============================================================
# PLOT
# ============================================================
results$label <- sprintf("%+.2f\n(NA %.0f%%)", results$diff, results$NA_pct*100)
results$unreliable <- results$NA_pct > 0.30

p <- ggplot(results, aes(x = gap, y = delta)) +
  geom_tile(aes(fill = diff), color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 2.8, color = "black") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-0.15, 0.20),
    oob = scales::squish,
    name = "UMPCU -\nInvN-Dunnett"
  ) +
  geom_tile(data = subset(results, unreliable),
            aes(x = gap, y = delta),
            fill = NA, color = "black", linewidth = 1.2, linetype = "dashed") +
  scale_x_continuous(breaks = gap_vec) +
  scale_y_continuous(breaks = delta_vec) +
  labs(
    title = "Power difference: UMPCU vs CP&CT (InvN-Dunnett) -- Binomial case",
    subtitle = sprintf("K=5 arms, nA=nB=100, pi0=0.40, alpha=0.025, nsim=%d. Dashed: UMPCU inapplicable in >30%% of trials.", nsim),
    x = expression(pi[1] - pi[2] ~ "(gap between top two arms)"),
    y = expression(Delta[M] == pi[1] - pi[0] ~ "(treatment effect)")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 9)
  )

ggsave("~/Desktop/heatmap_binomial_K5.pdf", p, width = 10, height = 7)
ggsave("~/Desktop/heatmap_binomial_K5.png", p, width = 10, height = 7, dpi = 200)

cat("\n=== Binomial K=5 heatmap complete ===\n")
cat("Files: ~/Desktop/heatmap_binomial_K5_grid.csv, heatmap_binomial_K5.{pdf,png}\n")
