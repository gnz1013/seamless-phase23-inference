#============================================================
# Heatmap: UMPCU vs InvN-Dunnett Power Difference
#
# Visualizes the power-advantage landscape across:
#   X-axis: mu1 - mu2 (gap between top two arms, 0 to 1.0)
#   Y-axis: Delta_M = mu1 - mu0 (treatment effect, 0.3 to 1.2)
#   Color : UMPCU_power - InvN_Dunnett_power
#red   -> UMPCU wins
#     blue  -> CP&CT (InvN-Dunnett) wins
#     white -> tie
#
# Setting: K=5, nA=nB=100, sigma^2=10, alpha=0.025
# Grid: 11x 10 = 110 cells, nsim=2000 per cell (~20 min on M1/M2 Mac)
# ============================================================

suppressPackageStartupMessages({
  library(mvtnorm)
  library(ggplot2)
})
set.seed(2026)

# ---- Copy all helper functions from cpct_variants_normal_FULL.R ----
restrictR <- function(numArm, n) sample(rep(0:numArm, n), n*(numArm+1), replace=FALSE)

simRsp <- function(n, numArm, mu, sigma) {
  ssn <- n*(numArm+1); out <- matrix(NA,ssn, numArm+1)
  for (i in 0:numArm) out[,i+1] <- rnorm(ssn, mu[i+1], sigma); out
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
  z_max <- max(qnorm(1- pc))
  corr <- matrix(0.5, r, r); diag(corr) <- 1
  if (is.infinite(df)) 1 - pmvnorm(upper=rep(z_max,r), corr=corr)[1]
  else                1 - pmvt(upper=rep(z_max,r), corr=corr, df=df)[1]
}

closureAdj <- function(p_raw, M, numArm, intersect_fn,...) {
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

combInvN <- function(p1,p2,w1,w2,alpha=0.025)
  as.integer((1 - pnorm(w1*qnorm(1-p1) + w2*qnorm(1-p2))) < alpha)

dstyQ <- function(w,d1,nA,nB,sigma,Y0bar,Zbar,X2) {
  n0 <- nA+nB
  mu_w <- (n0*(nA+nB)*d1)/((n0+nA+nB)*sigma^2)
  s2 <- (n0*(nA+nB))/((n0+nA+nB)*sigma^2)
  Tv <- (n0*Y0bar + (nA+nB)*Zbar)/(n0+nA+nB)
  v <- (sqrt(nA*(nA+nB))*(sigma^2*w + (nA+nB)*(Tv-X2)))/(sqrt(nB)*(nA+nB)*sigma)
  sqrt(2*pi*s2)*dnorm(w,mu_w,sqrt(s2))*pnorm(v)
}

umpcuTest <- function(W,nA,nB,sigma,Y0bar,Zbar,X2,alpha=0.025) {
  CN <- tryCatch(1/integrate(dstyQ,0,nA,nB,sigma,Y0bar,Zbar,X2,lower=-Inf,upper=Inf)$value,
                 error=function(e) NA_real_)
  if (is.na(CN)) return(NA_integer_)
  pval <- tryCatch(integrate(dstyQ,0,nA,nB,sigma,Y0bar,Zbar,X2,lower=W,upper=Inf)$value*CN,
                   error=function(e) NA_real_)
  if (is.na(pval)) return(NA_integer_)
  as.integer(pval < alpha)
}

# ---- Only compute UMPCU and InvN-Dunnett (faster) ----
compFunc2 <- function(numArm, nA, nB, mu, sigma, alpha=0.025) {
  rspA <- simRsp(nA,numArm,mu,sigma); allocA <- restrictR(numArm,nA)
  X<- armMeans(allocA,rspA,numArm)
  ord <- order(X[-1], decreasing=TRUE); M <- ord[1]; X2 <- X[ord[2]+1]
  p_raw <- rawPValues(allocA,rspA,numArm); df1 <- 2*nA-2
  
  rspB <- simRsp(nB,1, mu=c(mu[1],mu[M+1]), sigma); allocB <- restrictR(1,nB)
  YM <- mean(rspB[allocB==1,2]); Y0 <- mean(rspB[allocB==0,1])
  v1 <- var(rspB[allocB==1,2]); v0 <- var(rspB[allocB==0,1])
  n1 <- sum(allocB==1); n0s <- sum(allocB==0)
  se2 <- sqrt(v1/n1 + v0/n0s)
  df2 <- se2^4/(v1^2/(n1^2*(n1-1)) + v0^2/(n0s^2*(n0s-1)))
  p2 <- pt((YM-Y0)/se2, df2, lower.tail=FALSE)
  
  p1D <- closureAdj(p_raw,M,numArm,dunnettP,df=df1)
  eps <- 1e-15
  p1D <- min(max(p1D,eps),1-eps); p2c <- min(max(p2, eps),1-eps)
  w1 <- sqrt(nA/(nA+nB)); w2 <- sqrt(nB/(nA+nB))
  
  Zbar <- (nA*X[M+1] + nB*YM)/(nA+nB); Y0bar <- (nA*X[1]+nB*Y0)/(nA+nB)
  W <- (nA+nB)/(2*sigma^2)*(Zbar-Y0bar)
  
  c(UMPCU   = umpcuTest(W,nA,nB,sigma,Y0bar,Zbar,X2,alpha),
    InvN_D  = combInvN(p1D,p2c,w1,w2,alpha))
}

# ============================================================
# GRID SIMULATION
# ============================================================
nsim <- 2000
nA <- 100; nB <- 100; sigma <- sqrt(10); alpha <- 0.025; K <- 5

# Grids
gap_vec<- seq(0.0, 1.0, by = 0.1)   # mu1 - mu2  (11 values)
delta_vec   <- seq(0.3, 1.2, by = 0.1)   # mu1 - mu0  (10 values)

# Fix mu0 = 0.5, mu1 = mu0 + delta; mu2 = mu1 - gap
# Remaining K-2=3 arms: set to small uniform spread well below mu2
# so they don't interfere with top-2 ranking
mu0 <- 0.5

results <- expand.grid(gap = gap_vec, delta = delta_vec)
results$UMPCU <- NA_real_
results$InvN_D <- NA_real_
results$diff <- NA_real_

t0 <- Sys.time()
for (r in seq_len(nrow(results))) {
  gap <- results$gap[r]; delta <- results$delta[r]
  mu1 <- mu0 + delta
  mu2 <- mu1 - gap
  # Other3 arms scattered below mu2 with small effects
  mu_others <- c(0.10, 0.25, 0.40)
  mu <- c(mu0, mu1, mu2, mu_others)# length K+1= 6
  
  tests <- matrix(NA, nsim, 2)
  for (i in 1:nsim) tests[i,] <- compFunc2(K, nA, nB, mu, sigma, alpha)
  rates <- colMeans(tests, na.rm = TRUE)
  results$UMPCU[r]<- rates[1]
  results$InvN_D[r] <- rates[2]
  results$diff[r]   <- rates[1] - rates[2]
  
  if (r %% 10 == 0) {
    elapsed <- as.numeric(Sys.time() - t0, units = "mins")
    cat(sprintf("[%d/%d] elapsed %.1f min | gap=%.1f delta=%.1fU=%.3f ND=%.3f diff=%+.3f\n",
                r, nrow(results), elapsed, gap, delta,
                rates[1], rates[2], rates[1]-rates[2]))
  }
}

# Save raw grid for reproducibility
write.csv(results, "heatmap_grid.csv", row.names = FALSE)

# ============================================================
# PLOT
# ============================================================
p <- ggplot(results, aes(x = gap, y = delta, fill = diff)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%+.2f", diff)), size = 2.5, color = "black") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0,
    limits = c(-0.15, 0.15),
    oob = scales::squish,
    name = "UMPCU -\nInvN-Dunnett"
  ) +
  scale_x_continuous(breaks = gap_vec) +
  scale_y_continuous(breaks = delta_vec) +
  labs(
    title = "Power difference: UMPCU vs CP&CT (InvN-Dunnett)",
    subtitle = sprintf("K=5 arms, nA=nB=100, sigma^2=10, alpha=0.025, nsim=%d", nsim),
    x = expression(mu[1] - mu[2]~"(gap between top two arms)"),
    y = expression(Delta[M] == mu[1] - mu[0]~"(treatment effect)")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave("heatmap_UMPCU_vs_InvN_Dunnett.pdf", p, width = 9, height = 7)
ggsave("heatmap_UMPCU_vs_InvN_Dunnett.png", p, width = 9, height = 7, dpi = 200)
print(p)

cat("\n=== Heatmap complete ===\n")
cat("Files written: heatmap_grid.csv, heatmap_UMPCU_vs_InvN_Dunnett.{pdf,png}\n")
