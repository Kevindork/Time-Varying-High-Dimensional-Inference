source("initialize.R")

pb <- progress_bar$new(
  format = " Simulating [:bar] :percent in :elapsed. ETA :eta",
  total = n.sims * n, clear = F, width = 60)

set.seed(1)
beta <- create.beta(n, p, s) * 2.5
x <- mvrnorm(n, rep(0, p), sigma.x)
errors.array <- mvrnorm(n.sims, rep(0, n), sigma.errors)

t.all <- seq(0, 1, length.out = n)
i.bw.all <- which(t.all < (1 - bw) & t.all > bw)
n.bw.all <- length(i.bw.all)

gamma <- array(0, dim = c(p, n, n.sims))
Omega <- array(0, dim = c(p, n, n.sims))
bias.lasso <- array(0, dim = c(p, n, n.sims))
beta.lasso <- array(0, dim = c(p, n, n.sims))
beta.ridge <- array(0, dim = c(p, n, n.sims))
beta.tv <- array(0, dim = c(p, n, n.sims))
p.values <- array(0, dim = c(p, n, n.sims))
p.values.adj <- array(0, dim = c(p, n, n.sims))


for(sim in 1:n.sims){
  if(n.sims == 1)  {
    errors <- errors.array
  } else  {
    errors <- errors.array[sim, ]
  }
  y <- diag(x %*% beta) + errors

  for(i in 1:n)  {
    i.bw <- which(abs(t.all - t.all[i]) <= bw)
    n.bw <- length(i.bw)
    x.bw <- x[i.bw, ]
    y.bw <- y[i.bw]

    L1.tv <- sqrt(2 * log(p) / n.bw)
    L2.tv <- 1 / n.bw
    cov.x.bw <- t(x.bw) %*% x.bw / n.bw
    A <- solve(cov.x.bw + diag(L2.tv, p))

    # Model fitting
    fit.lasso <- scalreg(x.bw, y.bw, lam0 = L1.tv)
    beta.lasso[, i, sim] <- coef(fit.lasso)
    beta.ridge[, i, sim] <- A %*% t(x.bw) %*% y.bw / n.bw
    sigma.est <- fit.lasso$hsigma
    Omega.full <- sigma.est ^ 2 * A %*% cov.x.bw %*% A / n.bw
    Omega[, i, sim] <- sqrt(diag(Omega.full))

    # Projections
    v <- svd(x.bw)$v
    p.x <- v %*% t(v)
    diag(p.x) <- 0
    bias.lasso[, i, sim] <- p.x %*% beta.lasso[, i, sim]
    for(j in 1:p)  {
      gamma[j, i, sim] <- max(abs(p.x[j, ])) * (L1 / 2) ^ (1 / 2 - xi)
    }

    # Bias correction and inference
    beta.tv[, i, sim] <- beta.ridge[, i, sim] - bias.lasso[, i, sim]
    p.values[, i, sim] <- 2 * (1 - pnorm(abs(abs(beta.tv[, i, sim]) -
                                            gamma[, i, sim]) / Omega[, i, sim]))

    z <- abs(mvrnorm(n.norm, rep(0, p), Omega.full))
    cutoff <- apply(z, 1, function(x) max(x * Omega[, i, sim]))
    p.values.adj[, i, sim] <- sapply(p.values[,i,sim],
                                     function(x) sum(x > cutoff)) / n.norm
    pb$tick()
  }
}

n.FN.tv <- rep(NA, n.sims)
n.FP.tv <- rep(NA, n.sims)
FWER.tv <- rep(NA, n.sims)
RMSE.tv <- rep(NA, n.sims)
for(sim in 1:n.sims)  {
  n.FN.tv[sim] <- sum(p.values.adj[1:s, i.bw.all, sim] > alpha)
  n.FP.tv[sim] <- sum(p.values.adj[-(1:s), i.bw.all, sim] < alpha)
  FWER.tv[sim] <- mean(colSums(p.values.adj[-(1:s), , sim] < alpha) > 0)
  RMSE.tv <- sqrt(mean((beta.tv[, i.bw.all, sim] - beta[, i.bw.all]) ^ 2))
}

FPR.tv <- mean(n.FP.tv) / ((p - s) * n.bw.all * n.sims)
FNR.tv <- mean(n.FN.tv) / (s * n.bw.all * n.sims)

sd.FWER.tv <- sd(FWER.tv) / sqrt(n.sims)
sd.FNR.tv <- sd(n.FN.tv) / sqrt(n.sims) / (s * n.bw.all * n.sims)
sd.FPR.tv <- sd(n.FP.tv) / sqrt(n.sims) / ((p - s) * n.bw.all * n.sims)

results.tv <- data.frame(n.FP = mean(n.FP.tv), FPR = FPR.tv, sd.FPR = sd.FPR.tv,
                         n.FN = mean(n.FN.tv), FNR = FNR.tv, sd.FNR = sd.FNR.tv,
                         FWER = mean(FWER.tv), sd.FWER = sd.FWER.tv,
                         RMSE = mean(RMSE.tv))
