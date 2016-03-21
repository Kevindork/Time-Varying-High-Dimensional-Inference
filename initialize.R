library(glmnet)
library(scalreg)
library(progress)
library(MASS)

alpha <- 0.05
n.sims <- 20
b <- 2.5
bw <- 0.1
n <- 300
p <- 400
s <- 3
xi <- 0.05
n.norm <- 1000
errors.type <- "Diag"

# sigma.errors <- diag(n)
sigma.errors <- toeplitz(0.9 ^ c(0 : (n - 1)))
sigma.x <- diag(p)

# Time varying bias corrected ridge lambdas dynamic and unspecified here
L1.lasso <- sqrt(2 * log(p) / n)

if(exists("create.beta") == F)  {
  create.beta <- function(n, p, s, beta.min = 0.25,
                          n.knots.min = 3, n.knots.max = 6)  {
    beta <- array(0, dim = c(p, n))
    n.knots <- sample(n.knots.min : n.knots.max, s, rep = T)
    for(i in 1:s)  {
      i.knot <- round(seq(1, n, length.out = n.knots[i]))
      beta.knot <- runif(n.knots[i], beta.min, 1) *
                   sample(c(-1, 1), n.knots[i], replace = T)
      beta[i,] <- spline(i.knot, beta.knot, xout = 1:n)$y
    }
    return(beta)
  }

  colMax <- function(data)  {
    sapply(data, max, na.rm = TRUE)
  }

  plot.beta <- function(beta, index = "nonzero")  {
    if(index[1] == "nonzero")  {
      index <- which(rowSums(beta) != 0)
    }

    t <- seq(0, 1, length.out = dim(beta)[2])
    plot.dat <- cbind(t, as.data.frame(t(beta[index,])))
    names(plot.dat) <- c("t", paste("Beta", index, sep = " "))
    plot.dat <- melt(plot.dat, id.vars = t)
    names(plot.dat) <- c("t", "Variable", "Beta")
    ggplot(plot.dat, aes(x = t, y = Beta, col = Variable)) +
      geom_line(size = 1.2) +
      geom_hline(aes(yintercept = 0)) +
      scale_colour_discrete(name = "Variable") +
      theme_bw() +
      theme(panel.border = element_rect(colour = "#FFFFFF")) +
      theme(panel.grid.major = element_line(colour = "#D0D0D0", size = .75)) +
      theme(axis.ticks=element_blank()) +
      ggtitle("Beta vs Time") +
      theme(plot.title = element_text(face = "bold", hjust = -.05, vjust = 2,
                                      colour = "#3C3C3C", size = 20,
                                      margin = margin(0, 20, 10, 0))) +
      ylab("Beta") +
      xlab("Time") +
      theme(axis.text.x = element_text(size = 11, colour = "#535353",
                                       face = "bold")) +
      theme(axis.text.y = element_text(size = 11, colour="#535353",
                                       face = "bold")) +
      theme(axis.title.y = element_text(size = 13, colour="#535353",
                                        face = "bold",
                                        margin=margin(0, 15, 0, 0))) +
      theme(axis.title.x = element_text(size = 13, colour = "#535353",
                                        face = "bold",
                                        margin = margin(15, 0, 0, 0))) +
      theme(plot.margin = unit(c(1, 1, .5, .7), "cm"))
  }
}
