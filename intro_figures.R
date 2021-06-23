#######

library(MASS)
library(ggplot2)
set.seed(123)

# Generate data
N <- 1000 
rho <- 0.5
mu1 <- 0; s1 <- 1
mu2 <- 0; s2 <- 1
mu <- c(mu1, mu2) # Mean
sigma <- matrix(c(s1^2, s1*s2*rho, s1*s2*rho, s2^2), 2)
bvn <- data.frame(mvrnorm(N, mu = mu, Sigma = sigma))
names(bvn)[1] <- "y1"
names(bvn)[2] <- "y2"

# Plot
q1 <- quantile(bvn[,"y1"], 0.20)
q2 <- quantile(bvn[,"y2"], 0.20)

#FIGURE 2: The means-ends distinction
ggplot(bvn, aes(x=y1, y=y2)) + geom_point() + 
  labs(x = expression(y[2]), y = expression(y[1])) +
  geom_vline(xintercept=q1) + geom_hline(yintercept=q2) +
  annotate("text", x=-2, y=3, label= "IE") +
  annotate("text", x=2, y=3, label= "CE") +
  annotate("text", x=-2, y=-3, label= "CI") +
  annotate("text", x=2, y=-3, label= "EE") +
  theme_bw()
