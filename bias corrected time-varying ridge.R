pb <- progress_bar$new(
  format = " Simulating [:bar] :percent in :elapsed. ETA :eta",
  total = n.sims * n, clear = F, width = 60)

set.seed(1)
beta <- create.beta(n, p, s) * b
x <- mvrnorm(n, rep(0, p), sigma.x)
if(errors.type == "t")
{
  errors.array <- array(rt(n * n.sims, 3) / sqrt(3), dim = c(n.sims, n))
} else  {
  errors.array <- mvrnorm(n.sims, rep(0, n), sigma.errors)
}

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
      gamma[j, i, sim] <- max(abs(p.x[j, ])) * (L1.tv / 2) ^ (1 / 2 - xi)
    }

    # Bias correction and inference
    beta.tv[, i, sim] <- beta.ridge[, i, sim] - bias.lasso[, i, sim]
    p.values[, i, sim] <- 2 * (1 - pnorm(abs(abs(beta.tv[, i, sim]) -
                                             gamma[, i, sim]) /
                                         Omega[, i, sim]))

    z <- abs(mvrnorm(n.norm, rep(0, p), Omega.full))
    cutoff <- apply(z, 1, function(x) max(x * Omega[, i, sim]))
    p.values.adj[, i, sim] <- sapply(p.values[, i, sim],
                                     function(x) sum(x > cutoff)) / n.norm
    pb$tick()
  }
}

n.FN.tv <- rep(NA, n.sims)
n.FP.tv <- rep(NA, n.sims)
FWER.tv <- rep(NA, n.sims)

n.FN.tv.adj <- rep(NA, n.sims)
n.FP.tv.adj <- rep(NA, n.sims)
FWER.tv.adj <- rep(NA, n.sims)

RMSE.tv <- rep(NA, n.sims)
for(sim in 1:n.sims)  {
  n.FN.tv[sim] <- sum(p.values[1:s, i.bw.all, sim] > alpha)
  n.FP.tv[sim] <- sum(p.values[-(1:s), i.bw.all, sim] < alpha)
  FWER.tv[sim] <- mean(colSums(p.values[-(1:s), i.bw.all, sim] < alpha) > 0)

  n.FN.tv.adj[sim] <- sum(p.values.adj[1:s, i.bw.all, sim] > alpha)
  n.FP.tv.adj[sim] <- sum(p.values.adj[-(1:s), i.bw.all, sim] < alpha)
  FWER.tv.adj[sim] <- mean(colSums(p.values.adj[-(1:s), i.bw.all, sim] <
                                     alpha) > 0)

  RMSE.tv[sim] <- sqrt(mean((beta.tv[, i.bw.all, sim] - beta[, i.bw.all]) ^ 2))
}

FPR.tv <- mean(n.FP.tv) / ((p - s) * n.bw.all)
FNR.tv <- mean(n.FN.tv) / (s * n.bw.all)

FPR.tv.adj <- mean(n.FP.tv.adj) / ((p - s) * n.bw.all)
FNR.tv.adj <- mean(n.FN.tv.adj) / (s * n.bw.all)

sd.FNR.tv <- sd(n.FN.tv) / sqrt(n.sims) / (s * n.bw.all * n.sims)
sd.FPR.tv <- sd(n.FP.tv) / sqrt(n.sims) / ((p - s) * n.bw.all * n.sims)
sd.FWER.tv <- sd(FWER.tv) / sqrt(n.sims)

sd.FNR.tv.adj <- sd(n.FN.tv.adj) / sqrt(n.sims) / (s * n.bw.all * n.sims)
sd.FPR.tv.adj <- sd(n.FP.tv.adj) / sqrt(n.sims) / ((p - s) * n.bw.all * n.sims)
sd.FWER.tv.adj <- sd(FWER.tv.adj) / sqrt(n.sims)

sd.RMSE.tv <- sd(RMSE.tv) / sqrt(n.sims)

results.tv <- data.frame(n = n, p = p, s = s, b = b,
                         n.FP = mean(n.FP.tv), FPR = FPR.tv, sd.FPR = sd.FPR.tv,
                         n.FN = mean(n.FN.tv), FNR = FNR.tv, sd.FNR = sd.FNR.tv,
                         FWER = mean(FWER.tv), sd.FWER = sd.FWER.tv,
                         RMSE = mean(RMSE.tv), sd.RMSE = sd.RMSE.tv,
                         error = errors.type)

results.tv.adj <- data.frame(n = n, p = p, s = s, b = b,
                             n.FP = mean(n.FP.tv.adj), FPR = FPR.tv.adj,
                             sd.FPR = sd.FPR.tv.adj,
                             n.FN = mean(n.FN.tv.adj), FNR = FNR.tv.adj,
                             sd.FNR = sd.FNR.tv.adj,
                             FWER = mean(FWER.tv.adj), sd.FWER = sd.FWER.tv.adj,
                             RMSE = mean(RMSE.tv), sd.RMSE = sd.RMSE.tv,
                             error = errors.type)

if(print == T)  {
  write.table(results.tv, file = "TV-Ridge.csv", append = T, quote = F, sep = ",",
              eol = "\n", na = "NA", dec = ".", row.names = F,
              col.names = F, qmethod = c("escape", "double"),
              fileEncoding = "")

  write.table(results.tv.adj, file = "TV-Ridge, adjusted.csv", append = T, quote = F,
              sep = ",", eol = "\n", na = "NA", dec = ".", row.names = F,
              col.names = F, qmethod = c("escape", "double"),
              fileEncoding = "")
}
