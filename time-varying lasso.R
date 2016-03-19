source("initialize.R")

pb <- progress_bar$new(
  format = " Simulating [:bar] :percent in :elapsed. ETA :eta",
  total = n.sims * n, clear = F, width = 60)

set.seed(1)
beta <- create.beta(n, p, s) * b
x <- mvrnorm(n, rep(0, p), sigma.x)
errors.array <- mvrnorm(n.sims, rep(0, n), sigma.errors)

t.all <- seq(0, 1, length.out = n)
i.bw.all <- which(t.all < (1 - bw) & t.all > bw)
n.bw.all <- length(i.bw.all)

beta.lasso <- array(0, dim = c(p, n, n.sims))

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

    # Model fitting
    fit.lasso <- glmnet(x.bw, y.bw, lambda = L1.lasso, alpha = 1, intercept=F)
    beta.lasso[, i, sim] <- coef(fit.lasso)[-1]
    pb$tick()
  }
}

n.FN.lasso <- rep(NA, n.sims)
n.FP.lasso <- rep(NA, n.sims)
FWER.lasso <- rep(NA, n.sims)
RMSE.lasso <- rep(NA, n.sims)
for(sim in 1:n.sims)  {
  n.FN.lasso[sim] <- sum(beta.lasso[1:s, i.bw.all, sim] == 0)
  n.FP.lasso[sim] <- sum(beta.lasso[-(1:s), i.bw.all, sim] != 0)
  FWER.lasso[sim] <- mean(colSums(beta.lasso[-(1:s), i.bw.all, sim]) != 0)
  RMSE.lasso <- sqrt(mean((beta.lasso[, i.bw.all, sim] - beta[, i.bw.all]) ^ 2))
}

FPR.lasso <- mean(n.FP.lasso) / ((p - s) * n.bw.all)
FNR.lasso <- mean(n.FN.lasso) / (s * n.bw.all)

sd.FWER.lasso <- sd(FWER.lasso) / sqrt(n.sims)
sd.FNR.lasso <- sd(n.FN.lasso) / sqrt(n.sims) / (s * n.bw.all * n.sims)
sd.FPR.lasso <- sd(n.FP.lasso) / sqrt(n.sims) / ((p - s) * n.bw.all * n.sims)

results.lasso <- data.frame(n.FP = mean(n.FP.lasso), FPR = FPR.lasso, sd.FPR = sd.FPR.lasso,
                         n.FN = mean(n.FN.lasso), FNR = FNR.lasso, sd.FNR = sd.FNR.lasso,
                         FWER = mean(FWER.lasso), sd.FWER = sd.FWER.lasso,
                         RMSE = mean(RMSE.lasso))
