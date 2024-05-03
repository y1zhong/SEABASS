#'@title Fit SeaBASS model
#'@param data list with AE count (A), confounders (X), vaccine indication (V), number of each categories
#'@param initial list include the parameter initial values, default is from glm fit results
#' @export
VB = function(data, 
                   initial = NULL,
                   horseshoe = TRUE,
                   cutoff_beta = 4,
                   cutoff_alpha = 12,
                   mu_alpha = NULL,
                   mu_beta = NULL,
                   a_alpha = 0.001,
                   b_alpha  = 0.001,
                   init_inv_sigma2_beta = NULL,
                   init_inv_sigma2_alpha = NULL,
                   init_inv_a_beta = NULL,
                   init_inv_a_tau = NULL,
                   init_inv_tau2 = NULL,
                   nu_0 = 1,
                   eta2_0 = 1,
                   nu_1 = 1,
                   eta2_1 = 1,
                   include_intercept=FALSE,
                   ELBO_stop = 1,
                   max_iter=5000,
                   ELBO_diff_tol=1e-6,
                   paras_diff_tol = 1e-6,
                   verbose=0, save_profile = 1){
  
  
  AEs = colnames(data$A)
  if (include_intercept) {
    data$X = cbind(1, data$X)
  }
  
  J = ncol(data$A)
  p = ncol(data$X)
  
  ### set inital values
  if(is.null(mu_alpha)){
    mu_alpha = matrix(0, ncol = J, nrow=p)
    
  }
  if(is.null(mu_beta)){
    mu_beta = rep(0, J)
  }
  
  if(is.null(init_inv_sigma2_beta)) init_inv_sigma2_beta = rep(1, J)
  if(is.null(init_inv_sigma2_alpha)) init_inv_sigma2_alpha = matrix(1, ncol = J, nrow = p)
  
  if(is.null(init_inv_a_beta)) init_inv_a_beta =rep(0.5 , J)
  if(is.null(init_inv_tau2)) init_inv_tau2 = 1
  if(is.null(init_inv_a_tau)) init_inv_a_tau = 0.5
  
  if(is.null(initial)) {
    initials <- get_glm_coefs(data)
    initials <- coef_cutoff(initials, cutoff_beta, cutoff_alpha)
    init_beta = initials$beta
    init_alpha = initials$alpha
  } else {
    init_beta = initial$beta
    init_alpha = initial$alpha
  }
  
  seabass_fit <- SeaBASS_cpp(data$A, data$X, data$V,data$nn,
                        horseshoe=horseshoe,
                        initial_beta = init_beta,
                        initial_alpha = init_alpha,
                        initial_inv_sigma2_beta = init_inv_sigma2_beta,
                        initial_inv_sigma2_alpha =init_inv_sigma2_alpha,
                        initial_inv_a_beta = init_inv_a_beta,
                        initial_inv_tau2 = init_inv_tau2,
                        initial_inv_a_tau = init_inv_a_tau,
                        mu_alpha = mu_alpha,
                        mu_beta = mu_beta,
                        a_alpha = a_alpha,
                        b_alpha  = b_alpha,
                        nu_0 = nu_0,
                        eta2_0 = eta2_0,
                        nu_1 = nu_1,
                        eta2_1 = eta2_1,
                        ELBO_stop = ELBO_stop,
                        ELBO_diff_tol = ELBO_diff_tol,
                        paras_diff_tol = paras_diff_tol,
                        max_iter = max_iter,
                        verbose = verbose,
                        save_profile = save_profile)
  
  colnames(seabass_fit$data$A) = AEs
  return(seabass_fit)
}


#' @export
Gibbs = function(data, initial, hyper,
                      mcmc_sample = 500,
                      burnin = 5000,
                      thinning = 10,
                      verbose = 500,
                      save_profile = 1){
  

   
  init_beta = initial$beta
  init_alpha = initial$alpha
  initial_inv_sigma2_beta = initial$inv_sigma2_beta
  initial_inv_sigma2_alpha = initial$inv_sigma2_alpha
  initial_inv_a_beta = initial$inv_a_beta
  initial_inv_tau2 = initial$inv_tau2
  initial_inv_a_tau = initial$inv_a_tau
  
  a_alpha = hyper$a_alpha
  b_alpha  = hyper$b_alpha
  nu_0 = hyper$nu_0
  eta2_0 = hyper$eta2_0
  nu_1 = hyper$nu_1
  eta2_1 = hyper$eta2_1
  

  gibbs_fit <- Gibbs_cpp(data$A, data$X, data$V,data$nn,
                              initial_beta = init_beta,
                              initial_alpha = init_alpha,
                              initial_inv_sigma2_beta = initial_inv_sigma2_beta,
                              initial_inv_sigma2_alpha = initial_inv_sigma2_alpha,
                              initial_inv_a_beta = initial_inv_a_beta,
                              initial_inv_tau2 = initial_inv_tau2,
                              initial_inv_a_tau = initial_inv_a_tau,
                              nu_0 = nu_0,
                              eta2_0 = eta2_0,
                              nu_1 = nu_1,
                              eta2_1 = eta2_1,
                              mcmc_sample = mcmc_sample,
                              burnin = burnin,
                              thinning = thinning,
                              verbose = verbose,
                              save_profile = save_profile)
  gibbs_fit[["initial"]] = list(init_beta = init_beta,
                                init_alpha = init_alpha)

  return(gibbs_fit)
}


SEABASS <- function(data, nc_idx, signal_cutoff, threshold=0.95,
                    max_iter=5000,
                    ELBO_diff_tol=1e-8,
                    paras_diff_tol = 1e-6,
                    verbose=0, save_profile = 1,
                    mcmc_sample = 500,
                    burnin = 100,thinning = 1){

  VB_fit <- VB(data,
               initial = NULL,
               horseshoe = TRUE,
               include_intercept=FALSE,
               ELBO_stop = 1,
               max_iter=max_iter,
               ELBO_diff_tol=ELBO_diff_tol,
               paras_diff_tol = paras_diff_tol,
               verbose=verbose, save_profile = save_profile)
  
  init_list <- list(beta=VB_fit$post_mean$beta,
                    alpha=VB_fit$post_mean$alpha,
                    inv_sigma2_beta=VB_fit$post_mean$inv_sigma2_beta,
                    inv_sigma2_alpha=VB_fit$post_mean$inv_sigma2_alpha,
                    inv_a_beta=VB_fit$post_mean$inv_a_beta,
                    inv_tau2=VB_fit$post_mean$inv_tau2,
                    inv_a_tau=VB_fit$post_mean$inv_a_tau)
  init_list <- VB_fit$post_mean
  hyper_list  <- VB_fit$hyper
  Gibbs_fit <- Gibbs(data, initial = init_list,
                     hyper = hyper_list,
                     mcmc_sample = mcmc_sample,
                     burnin = burnin,
                     thinning = thinning,
                     verbose = 0,
                     save_profile = 1)
  
  vb_beta = VB_fit$post_mean$beta
  mcmc_sd = apply(Gibbs_fit$mcmc$beta, 1, sd)
  mcmc_resample = sapply(1:J_all, function(j) rnorm(1000, vb_beta[j], mcmc_sd[j])) 
  vb_nc_mcmc =  rowMeans(mcmc_resample[, nc_idx])
  vb_adjust = colMeans(apply(Gibbs_fit$mcmc$beta, 1, function(x) x - vb_nc_mcmc))  
  
  
  vb_prob = colMeans(apply(Gibbs_fit$mcmc$beta, 1, function(x) x - vb_nc_mcmc)  > cutoff)
  vb_signal = as.numeric(vb_prob >=  threshhold)[-nc_idx]
  
  return(list(signal = vb_signal,
              probability = vb_prob,
              vb_fit = VB_fit,
              gibbs_fit = Gibbs_fit,
              data = data))
}
