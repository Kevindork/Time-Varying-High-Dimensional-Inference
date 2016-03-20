source("initialize.R")

pb <- progress_bar$new(
  format = " Simulating [:bar] :percent in :elapsed. ETA :eta",
  total = n.sims, clear = F, width = 60)

set.seed(1)
beta <- create.beta(n, p, s) * b
x <- mvrnorm(n, rep(0, p), sigma.x)
errors.array <- mvrnorm(n.sims, rep(0, n), sigma.errors)

t.all <- seq(0, 1, length.out = n)
i.bw.all <- which(t.all < (1 - bw) & t.all > bw)
n.bw.all <- length(i.bw.all)

gamma <- array(0, dim = c(p, n.sims))
Omega <- array(0, dim = c(p, n.sims))
bias.lasso <- array(0, dim = c(p, n.sims))
beta.lasso <- array(0, dim = c(p, n.sims))
beta.ridge <- array(0, dim = c(p, n.sims))
beta.buhl <- array(0, dim = c(p, n.sims))
p.values <- array(0, dim = c(p, n.sims))
p.values.adj <- array(0, dim = c(p, n.sims))

for(sim in 1:n.sims){
  if(n.sims == 1)  {
    errors <- errors.array
  } else  {
    errors <- errors.array[sim, ]
  }
  y <- diag(x %*% beta) + errors

  L1.buhl <- sqrt(2 * log(p) / n)
  L2.buhl <- 1 / n
  cov.x <- t(x) %*% x / n
  A <- solve(cov.x + diag(L2.buhl, p))

  # Model fitting
  fit.lasso <- scalreg(x, y, lam0 = L1.buhl)
  beta.lasso[, sim] <- coef(fit.lasso)
  beta.ridge[, sim] <- A %*% t(x) %*% y / n
  sigma.est <- fit.lasso$hsigma
  Omega.full <- sigma.est ^ 2 * A %*% cov.x %*% A / n
  Omega[, sim] <- sqrt(diag(Omega.full))

  # Projections
  v <- svd(x)$v
  p.x <- v %*% t(v)
  diag(p.x) <- 0
  bias.lasso[, sim] <- p.x %*% beta.lasso[, sim]
  for(j in 1:p)  {
    gamma[j, sim] <- max(abs(p.x[j, ])) * (L1.buhl / 2) ^ (1 / 2 - xi)
  }

  # Bias correction and inference
  beta.buhl[, sim] <- beta.ridge[, sim] - bias.lasso[, sim]
  p.values[, sim] <- 2 * (1 - pnorm(abs(abs(beta.buhl[, sim]) - gamma[, sim]) /
                                    Omega[, sim]))

  z <- abs(mvrnorm(n.norm, rep(0, p), Omega.full))
  cutoff <- apply(z, 1, function(x) max(x * Omega[, sim]))
  p.values.adj[, sim] <- sapply(p.values[, sim],
                                   function(x) sum(x > cutoff)) / n.norm
  pb$tick()
}

n.FN.buhl <- rep(NA, n.sims)
n.FP.buhl <- rep(NA, n.sims)
FWER.buhl <- rep(NA, n.sims)

n.FN.buhl.adj <- rep(NA, n.sims)
n.FP.buhl.adj <- rep(NA, n.sims)
FWER.buhl.adj <- rep(NA, n.sims)

RMSE.buhl <- rep(NA, n.sims)
for(sim in 1:n.sims)  {
  n.FN.buhl[sim] <- sum(p.values[1:s, sim] > alpha) * n.bw.all
  n.FP.buhl[sim] <- sum(p.values[-(1:s), sim] < alpha) * n.bw.all
  FWER.buhl[sim] <- sum(p.values[-(1:s), sim] < alpha) > 0

  n.FN.buhl.adj[sim] <- sum(p.values.adj[1:s, sim] > alpha) * n.bw.all
  n.FP.buhl.adj[sim] <- sum(p.values.adj[-(1:s), sim] < alpha) * n.bw.all
  FWER.buhl.adj[sim] <- sum(p.values.adj[-(1:s), sim] < alpha) > 0

  RMSE.buhl[sim] <- sqrt(mean((beta.buhl[, sim] - beta[, i.bw.all]) ^ 2))
}

FPR.buhl <- mean(n.FP.buhl) / ((p - s) * n.bw.all)
FNR.buhl <- mean(n.FN.buhl) / (s * n.bw.all)

FPR.buhl.adj <- mean(n.FP.buhl.adj) / ((p - s) * n.bw.all)
FNR.buhl.adj <- mean(n.FN.buhl.adj) / (s * n.bw.all)

sd.FNR.buhl <- sd(n.FN.buhl) / sqrt(n.sims) / (s * n.bw.all)
sd.FPR.buhl <- sd(n.FP.buhl) / sqrt(n.sims) / ((p - s) * n.bw.all)
sd.FWER.buhl <- sd(FWER.buhl) / sqrt(n.sims)

sd.FNR.buhl.adj <- sd(n.FN.buhl.adj) / sqrt(n.sims) / (s * n.bw.all)
sd.FPR.buhl.adj <- sd(n.FP.buhl.adj) / sqrt(n.sims) / ((p - s) * n.bw.all)
sd.FWER.buhl.adj <- sd(FWER.buhl.adj) / sqrt(n.sims)

sd.RMSE.buhl <- sd(RMSE.buhl) / sqrt(n.sims)

results.buhl <- data.frame(n = n, p = p, s = s, b = b,
                           n.FP = mean(n.FP.buhl), FPR = FPR.buhl,
                           sd.FPR = sd.FPR.buhl,
                           n.FN = mean(n.FN.buhl), FNR = FNR.buhl,
                           sd.FNR = sd.FNR.buhl,
                           FWER = mean(FWER.buhl), sd.FWER = sd.FWER.buhl,
                           RMSE = mean(RMSE.buhl), sd.RMSE = sd.RMSE.buhl,
                           error = errors.type)

results.buhl.adj <- data.frame(n = n, p = p, s = s, b = b,
                               n.FP = mean(n.FP.buhl.adj), FPR = FPR.buhl.adj,
                               sd.FPR = sd.FPR.buhl.adj,
                               n.FN = mean(n.FN.buhl.adj), FNR = FNR.buhl.adj,
                               sd.FNR = sd.FNR.buhl.adj,
                               FWER = mean(FWER.buhl.adj),
                               sd.FWER = sd.FWER.buhl.adj,
                               RMSE = mean(RMSE.buhl), sd.RMSE = sd.RMSE.buhl,
                               error = errors.type)


write.table(results.buhl, file = "Ridge.csv", append = T, quote = F, sep = ",",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = F, qmethod = c("escape", "double"),
            fileEncoding = "")

write.table(results.buhl.adj, file = "Ridge, adjusted.csv", append = T, quote = F,
            sep = ",", eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = F, qmethod = c("escape", "double"),
            fileEncoding = "")
