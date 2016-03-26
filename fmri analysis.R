library(R.matlab)
library(stringr)
library(network)
library(sna)
library(scalreg)
library(GGally)
library(reshape2)
library(progress)
patient.list.N <- paste("N", str_pad(c(3:5, 7:8, 10:12, 14:15), 3, pad = "0"),
                        sep = "")
patient.list.P <- paste("P", str_pad(c(1:9, 11, 13, 16), 3, pad = "0"),
                        sep = "")
patient.list <- sort(c(patient.list.N,
                       paste(patient.list.P, "A", sep = ""),
                       paste(patient.list.P, "B", sep = "")))

bw <- 0.05
xi <- 0.05
for(id.patient in 1:length(patient.list))  {
  patient <- patient.list[id.patient]
  print(paste("Working on patient", patient))

  file <- paste("Data/BCT_REST_", patient, "_ROIs.mat", sep = "")
  dat <- readMat(file)$data
  n <- dim(dat)[1]
  p <- dim(dat)[2] - 1


  gamma <- array(0, dim = c(p, n, p + 1))
  Omega <- array(0, dim = c(p, n, p + 1))
  bias.lasso <- array(0, dim = c(p, n, p + 1))
  beta.lasso <- array(0, dim = c(p, n, p + 1))
  beta.ridge <- array(0, dim = c(p, n, p + 1))
  beta.tv <- array(0, dim = c(p, n, p + 1))
  p.values <- array(0, dim = c(p, n, p + 1))
  p.values.adj <- array(0, dim = c(p + 1, n, p + 1))


  pb <- progress_bar$new(
    format = " Simulating [:bar] :percent in :elapsed. ETA :eta",
    total = (p + 1) * n, clear = F, width = 60)

  set.seed(1)
  for(node in 1:(p + 1))  {
    x <- dat[, -node]
    y <- dat[, node]

    t.all <- seq(0, 1, length.out = n)
    i.bw.all <- which(t.all < (1 - bw) & t.all > bw)
    n.bw.all <- length(i.bw.all)

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
      beta.lasso[, i, node] <- coef(fit.lasso)
      beta.ridge[, i, node] <- A %*% t(x.bw) %*% y.bw / n.bw
      sigma.est <- fit.lasso$hsigma
      Omega.full <- sigma.est ^ 2 * A %*% cov.x.bw %*% A / n.bw
      Omega[, i, node] <- sqrt(diag(Omega.full))

      # Projections
      v <- svd(x.bw)$v
      p.x <- v %*% t(v)
      diag(p.x) <- 0
      bias.lasso[, i, node] <- p.x %*% beta.lasso[, i, node]
      for(j in 1:p)  {
        gamma[j, i, node] <- max(abs(p.x[j, ])) * (L1.tv / 2) ^ (1 / 2 - xi)
      }

      # Bias correction and inference
      beta.tv[, i, node] <- beta.ridge[, i, node] - bias.lasso[, i, node]
      p.values[, i, node] <- 2 * (1 - pnorm(abs(abs(beta.tv[, i, node]) -
                                                  gamma[, i, node]) / Omega[, i, node]))

      z <- abs(mvrnorm(n.norm, rep(0, p), Omega.full))
      cutoff <- apply(z, 1, function(x) max(x * Omega[, i, node]))
      p.values.adj[-node, i, node] <- sapply(p.values[, i, node],
                                             function(x) sum(x + 0.05 > cutoff)) / n.norm
      pb$tick()
    }
  }

  pb <- progress_bar$new(
    format = " Saving Plots [:bar] :percent in :elapsed",
    total = n.bw.all, clear = FALSE, width = 60)

  colors <- rep("#00BDC4", 52)
  colors[c(1:4, 6, 8, 14, 16, 40)] <- "#F8766D"
  for(t in 1:n.bw.all)  {
    adj <- p.values.adj[, t, ] < 0.05
    net = network(adj, directed = T)
    net %v% "color" = colors
    title <- paste("Patient ", patient, " at Frame ", t, sep = "")
    file <- paste("Graphs/", patient, "/", patient, ",",
                  str_pad(t, 3, pad = "0"), ".png", sep = "")

    CairoPNG(filename = file, width = 600, height = 600,
             pointsize = 12, bg = "transparent")
    plot <- ggnet2(net, mode = "circle", size = 12, color = "color", alpha = 0.5,
                   edge.size = 1, edge.color = "#AAAAAA") +
      geom_point(aes(color = color), size = 9, alpha = 0.8) +
      geom_text(aes(label = 1:52), color = "white", family = "Century Gothic",
                fontface = "bold", size = 5) +
      ggtitle(title)
    print(plot)
    dev.off()

    pb$tick()
  }
}
