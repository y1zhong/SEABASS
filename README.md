# SeaBASS

```R
##########data generation settings###########
set.seed(0)
J = 100
N = 1000
tps = 10
p = 3
nc_J = 10
J_all = J + nc_J

##########generate true effects###########
beta_true <- rnorm(J)
beta_nc <- rnorm(nc_J, 0.5, 0.5)
beta_bias <- beta_true + mean(beta_nc)
beta <- c(beta_bias, beta_nc)
alpha <- matrix(rnorm((J + nc_J)*p), nrow=p)
nc_idx <- c(101:110)

##########generate data###########
seq_data <- gen_data(beta=beta, alpha=alpha, N=N, tps=tps, vax_p=0.5)

##########fit SEABASSmodels###########
seabass_fit <- SEABASS(seq_data[[10]], nc_idx = nc_idx, signal_cutoff = 0, threshold = 0.95)

