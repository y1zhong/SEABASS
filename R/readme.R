J = 100
N = 1000
tps = 10
p = 3


nc_J = 10

beta_true <- rnorm(J)
beta_nc <- rnorm(nc_J, 0.5, 0.5)
beta_bias <- beta_true + mean(beta_nc)
beta <- c(beta_bias, beta_nc)

alpha <- matrix(rnorm((J + nc_J)*p), nrow=p)

