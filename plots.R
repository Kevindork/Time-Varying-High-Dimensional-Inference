dat.lasso.magic <- cbind(read.csv("Lasso, magic.csv", header=T),
                         Method = "Lasso (Magic)")[, -5]
dat.lasso.univ <- cbind(read.csv("Lasso, univ.csv", header=T),
                        Method = "Lasso (Univ)")[, -5]
dat.ridge.buhl <- cbind(read.csv("Ridge, adjusted.csv", header=T),
                        Method = "Debiased Ridge")
dat.ridge.tv <- cbind(read.csv("TV-Ridge, adjusted.csv", header=T),
                      Method = "Proposed")

# FPR vs Signal-to-Noise Ratio
dat <- rbind(dat.ridge.tv, dat.ridge.buhl, dat.lasso.univ)
dat$SNR <- dat$b * base::norm(diag(x %*% (beta / b)), "2") / sqrt(n)
ggplot(dat, aes(SNR, FPR, colour = Method)) +
  geom_line(size = 1) +
  geom_point(size = 5) +
  geom_point(size=3, colour="#FFFFFF") +
  scale_y_log10(breaks = 10^(0:-6), limits = c(1e-6, 1)) +
  theme_bw() +
  ggtitle("Type I Error Rate vs Signal to Noise Ratio") +
  theme(plot.title=element_text(size = 15)) +
  ylab("Type I Error Rate") +
  xlab("Signal to Noise Ratio")

# FWER vs Signal-to-Noise Ratio
dat <- rbind(dat.ridge.tv, dat.ridge.buhl, dat.lasso.univ)
dat$SNR <- dat$b * base::norm(diag(x %*% (beta / b)), "2") / sqrt(n)
ggplot(dat, aes(SNR, FWER, colour = Method)) +
  geom_line(size = 1) +
  geom_hline(aes(yintercept = 0.05), linetype = 2) +
  geom_point(size = 5) +
  geom_point(size=3, colour="#FFFFFF") +
  scale_y_log10(breaks = 10^(0:-3), limits = c(1e-3, 1)) +
  theme_bw() +
  ggtitle("Familywise Error Rate vs Signal to Noise Ratio") +
  theme(plot.title=element_text(size = 15)) +
  ylab("FWER") +
  xlab("Signal to Noise Ratio")

# Power vs Signal-to-Noise Ratio
dat <- rbind(dat.ridge.tv, dat.ridge.buhl, dat.lasso.univ, dat.lasso.magic)
dat$SNR <- dat$b * base::norm(diag(x %*% (beta / b)), "2") / sqrt(n)
dat$Power <- 1 - dat$FNR
ggplot(dat, aes(SNR, Power, colour = Method)) +
  geom_line(size = 1) +
  geom_point(size = 5) +
  geom_point(size=3, colour="#FFFFFF") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  ggtitle("Power vs Signal to Noise Ratio") +
  theme(plot.title=element_text(size = 15)) +
  ylab("Power") +
  xlab("Signal to Noise Ratio")
