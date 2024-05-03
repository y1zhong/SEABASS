#include <RcppArmadillo.h>
#include "rcpp_pgdraw.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
Environment pkg = Environment::namespace_env("stats");
Function rg = pkg["rgamma"];

Environment pgdraw_pkg = Environment::namespace_env("pgdraw");
Function pg_draw = pgdraw_pkg["pgdraw"];

Environment bayelogit_pkg = Environment::namespace_env("BayesLogit");
Function rpg = bayelogit_pkg["rpg"];

class Gibbs_VAX{
private:
  int hs; // 0 for base model, 1 for horseshoe prior
  struct LogitReg_Gibbs{
    int num_reports;
    int num_AEs;
    int num_confounders;
    
    arma::mat A;
    arma::mat X;
    arma::vec V;
    arma::vec nn;
    
    NumericVector nn_rep;
  } dat;
  
  struct HyperParas_Gibbs{
    double a_alpha;
    double b_alpha;
    double nu_0;
    double eta2_0;
    
    double nu_1;
    double eta2_1;
  } hyperparas;
  
  struct LogitRegParas_Gibbs{
    arma::vec beta;
    arma::vec beta_sq;
    arma::vec var_beta;
    
    arma::vec inv_sigma2_beta;
    arma::vec inv_a_beta;
    
    double inv_tau2;
    double inv_a_tau;
    
    arma::mat alpha;
    arma::mat alpha_sq;
    arma::mat inv_sigma2_alpha;
    arma::mat pre_mat_j;
    arma::mat cov_mat_j;
    arma::vec mu_vec_j;
    arma::vec b;
    
    
    arma::mat omega;
    
    arma::mat temp;
    double loglik;
    double logpost;
  } paras;
  
  struct MCMCSample{
    arma::mat beta;
    arma::mat inv_sigma2_beta;
    arma::mat inv_a_beta;
    
    arma::vec inv_tau2;
    arma::vec inv_a_tau;
    
    arma::cube alpha;
    arma::cube inv_sigma2_alpha;
  } paras_sample;
  
  
  struct GibbsSamplerProfile{
    vec loglik;
    vec logpost;
  } gibbs_profile;
  
  struct GibbsSamplerControl{
    int total_iter;
    int burnin;
    int mcmc_sample;
    int thinning;
    int verbose;
    int save_profile;
    int total_profile;
  } gibbs_control;
  
  int iter;
  
public:
  void set_hs(int in_hs){
    if(in_hs==0){
      hs = 0;
    } else if(in_hs==1){
      hs = 1;
    }
  };
  
  int get_hs(){
    return hs;
  }
  
  void load_data_Gibbs(const arma::mat& in_A, const arma::mat& in_X, const arma::vec& in_V, const arma::vec& in_nn){
    dat.A = in_A;
    dat.X = in_X;
    dat.V = in_V;
    dat.nn = in_nn;
    dat.nn_rep = wrap(dat.nn);
    
    
    dat.num_AEs = dat.A.n_cols;
    dat.num_confounders = dat.X.n_cols;
    dat.num_reports = dat.X.n_rows;
    
    dat.nn_rep =rep(dat.nn_rep, dat.num_AEs);
  };
  
  void set_hyperparas_Gibbs(const double& in_a_alpha, const double& in_b_alpha,
                            const double& in_nu_0,const double& in_eta2_0,
                            const double& in_nu_1,const double& in_eta2_1){
    hyperparas.a_alpha = in_a_alpha;
    hyperparas.b_alpha = in_b_alpha;
    
    hyperparas.nu_0 = in_nu_0;
    hyperparas.eta2_0 = in_eta2_0;
    
    hyperparas.nu_1 = in_nu_1;
    hyperparas.eta2_1 = in_eta2_1;
  };
  
  void set_paras_initial_values_Gibbs(const arma::vec& in_beta,
                                      const arma::mat& in_alpha,
                                      const arma::vec& in_inv_sigma2_beta,
                                      const arma::mat& in_inv_sigma2_alpha,
                                      const arma::vec& in_inv_a_beta,
                                      const double& in_inv_tau2,
                                      const double& in_inv_a_tau){
    paras.beta = in_beta;
    paras.beta_sq = in_beta%in_beta;
    paras.inv_sigma2_beta = in_inv_sigma2_beta;
    paras.var_beta = 1/paras.inv_sigma2_beta;
    paras.inv_a_beta = in_inv_a_beta;
    
    paras.inv_tau2 = in_inv_tau2;
    paras.inv_a_tau = in_inv_a_tau;
    
    paras.alpha=in_alpha;
    paras.alpha_sq=in_alpha%in_alpha;
    paras.inv_sigma2_alpha = in_inv_sigma2_alpha;
    paras.omega.zeros(dat.num_reports, dat.num_AEs);
    update_omega();
  };
  
  void set_gibbs_control(int in_mcmc_sample, int in_burnin, int in_thinning,
                         int in_verbose, int in_save_profile){
    gibbs_control.mcmc_sample = in_mcmc_sample;
    gibbs_control.burnin = in_burnin;
    gibbs_control.thinning = in_thinning;
    gibbs_control.total_iter = gibbs_control.burnin;
    gibbs_control.total_iter += gibbs_control.mcmc_sample*gibbs_control.thinning;
    gibbs_control.verbose = in_verbose;
    gibbs_control.save_profile = in_save_profile;
    if(gibbs_control.save_profile > 0){
      gibbs_control.total_profile = gibbs_control.total_iter/gibbs_control.save_profile;
    } else{
      gibbs_control.total_profile = 0;
    }
  };
  
  void update_beta(){
    //std::cout << "update_beta" << std::endl;
    paras.var_beta = 1/(sum(paras.omega.each_col()%dat.V,0).t() + paras.inv_sigma2_beta * paras.inv_tau2);
    paras.temp = (dat.A.each_col() - dat.nn/2) - paras.omega % (dat.X*paras.alpha);
    paras.temp = paras.temp.each_col() % dat.V;
    paras.beta = conv_to< colvec >::from(sum(paras.temp, 0));
    paras.beta %= paras.var_beta;
    paras.beta = randn(dat.num_AEs) % sqrt(paras.var_beta) + paras.beta;
    paras.beta_sq = paras.beta%paras.beta;
  };
  
  void update_inv_sigma2_beta(){
    //std::cout << "update_inv_sigma2_beta" << std::endl;
    // paras.inv_sigma2_beta = as<arma::vec>(rg(dat.num_AEs, (hyperparas.nu+1)*0.5, (hyperparas.nu*paras.inv_a_beta + 0.5*paras.beta_sq)));
    for(int j=0; j<dat.num_AEs; j++){
      paras.inv_sigma2_beta(j) = randg(distr_param((hyperparas.nu_0+1)*0.5, 1.0/(hyperparas.nu_0*paras.inv_a_beta(j) + 0.5*paras.inv_tau2*paras.beta_sq(j))));
      // std::cout << (hyperparas.nu*paras.inv_a_beta + 0.5*paras.beta_sq(j))<< std::endl;
      // std::cout << paras.inv_sigma2_beta(j)<< std::endl;
    }
    
  };
  void update_inv_a_beta(){
    //std::cout << "update_inv_a_beta" << std::endl;
    // paras.inv_a_beta = as<double>(rg(1, (dat.num_AEs*hyperparas.nu+1)/2, (hyperparas.nu*accu(paras.inv_sigma2_beta) + 1/hyperparas.eta2)));
    for(int j=0; j<dat.num_AEs; j++){
      paras.inv_a_beta(j) = randg(distr_param((hyperparas.nu_0+1)*0.5, 1.0/(hyperparas.nu_0*paras.inv_sigma2_beta(j) + 1/hyperparas.eta2_0)));
    }
    //paras.inv_a_beta = randg(distr_param((hyperparas.nu_0+1)/2, 1.0/(hyperparas.nu_0*paras.inv_sigma2_beta + 1/hyperparas.eta2_0)));
  };
  void update_inv_sigma2_alpha(){
    //std::cout << "update_inv_sigma2_alpha" << std::endl;
    // paras.inv_sigma2_alpha = reshape(as<arma::vec>(rg(dat.num_AEs*dat.num_confounders, hyperparas.a_alpha+0.5, (hyperparas.b_alpha + 0.5*vectorise(paras.alpha_sq)))), dat.num_confounders, dat.num_AEs);
    for(int j=0; j<dat.num_AEs; j++){
      for(int i=0; i<dat.num_confounders; i++){
        paras.inv_sigma2_alpha(i,j) = randg(distr_param(hyperparas.a_alpha+0.5, 1.0/(hyperparas.b_alpha + 0.5*paras.alpha_sq(i,j))));
      }
    }
  };
  
  
  void update_alpha(){
    //std::cout << "update_alpha" << std::endl;
    //arma::vec ind = linspace(0, dat.num_AEs-1, dat.num_AEs);
    for(int j=0; j<dat.num_AEs; j++){
      paras.pre_mat_j = dat.X.each_col() % paras.omega.col(j);
      paras.pre_mat_j =  paras.pre_mat_j.t() * dat.X;
      paras.pre_mat_j.diag() += paras.inv_sigma2_alpha.col(j);
      
      paras.mu_vec_j =  dat.X.t() * (dat.A.col(j) - dat.nn/2 - paras.omega.col(j) % dat.V * paras.beta(j));
      
      // paras.cov_mat_j = inv(paras.pre_mat_j);
      // paras.mu_vec_j = paras.cov_mat_j*paras.mu_vec_j;
      // paras.alpha.col(j) = mvnrnd(paras.mu_vec_j, paras.cov_mat_j, 1);
      
      paras.pre_mat_j = chol(paras.pre_mat_j);
      paras.b = solve(paras.pre_mat_j.t(), paras.mu_vec_j);
      paras.alpha.col(j) = randn(dat.num_confounders);
      paras.alpha.col(j) = solve(paras.pre_mat_j, paras.alpha.col(j) + paras.b);
    }
    paras.alpha_sq = paras.alpha%paras.alpha;
  };
  
  
  
  
  
  void update_omega(){
    //std::cout << "update_omega" << std::endl;
    paras.temp = dat.X*paras.alpha + dat.V*paras.beta.t();
    
    paras.omega = reshape(as<arma::vec>(rpg(dat.nn_rep.size(), dat.nn_rep, wrap(vectorise(paras.temp)))), dat.num_reports,dat.num_AEs);
    
    
  };
  
  void update_inv_tau2(){
    paras.inv_tau2 = randg(distr_param((dat.num_AEs + hyperparas.nu_1)*0.5, 1.0/(0.5*accu(paras.beta_sq % paras.inv_sigma2_beta) + hyperparas.nu_1*paras.inv_a_tau)));
  };
  
  void update_inv_a_tau(){
    paras.inv_a_tau = randg(distr_param((hyperparas.nu_1+1)/2, 1.0/(hyperparas.nu_1 * paras.inv_tau2 + 1/hyperparas.eta2_1)));
  };
  
  void update_loglik(){
    paras.loglik = 0.5;
  }
  
  void update_logpost(){
    paras.temp = dat.X * paras.alpha + dat.V * paras.beta.t();
    paras.logpost = accu(paras.temp % dat.A);
    paras.temp = paras.temp.for_each( [](mat::elem_type& val) { val = R::log1pexp(val); } );
    paras.logpost -= accu(paras.temp.each_col() % dat.nn);
  }
  
  void initialize_paras_sample(){
    paras_sample.beta.zeros(paras.beta.n_elem,gibbs_control.mcmc_sample);
    paras_sample.alpha.zeros(paras.alpha.n_rows, paras.alpha.n_cols, gibbs_control.mcmc_sample);
    paras_sample.inv_sigma2_alpha.zeros(paras.alpha.n_rows, paras.alpha.n_cols, gibbs_control.mcmc_sample);
    paras_sample.inv_sigma2_beta.zeros(paras.beta.n_elem,gibbs_control.mcmc_sample);
    paras_sample.inv_a_beta.zeros(paras.beta.n_elem, gibbs_control.mcmc_sample);
    paras_sample.inv_tau2.zeros(gibbs_control.mcmc_sample);
    paras_sample.inv_a_tau.zeros(gibbs_control.mcmc_sample);
  }
  
  void save_paras_sample(){
    if(iter >= gibbs_control.burnin){
      if((iter - gibbs_control.burnin)%gibbs_control.thinning==0){
        int mcmc_iter = (iter - gibbs_control.burnin)/gibbs_control.thinning;
        //std::cout << "paras_sample.beta" << std::endl;
        paras_sample.beta.col(mcmc_iter) = paras.beta;
        paras_sample.alpha.slice(mcmc_iter) = paras.alpha;
        //std::cout << "sigma2_beta" << std::endl;
        paras_sample.inv_sigma2_beta.col(mcmc_iter) = paras.inv_sigma2_beta;
        paras_sample.inv_sigma2_alpha.slice(mcmc_iter) = paras.inv_sigma2_alpha;
        //std::cout << "a_beta.beta" << std::endl;
        paras_sample.inv_a_beta.col(mcmc_iter) = paras.inv_a_beta;
        //std::cout << "tau2" << std::endl;
        paras_sample.inv_tau2(mcmc_iter) = paras.inv_tau2;
        paras_sample.inv_a_tau(mcmc_iter) = paras.inv_a_tau;
        //std::cout << "a_tau" << std::endl;
      }
    }
  };
  
  void initialize_gibbs_profile(){
    //std::cout << "total_profile: " << gibbs_control.total_profile << std::endl;
    if(gibbs_control.save_profile>0){
      gibbs_profile.loglik.zeros(gibbs_control.total_profile);
      gibbs_profile.logpost.zeros(gibbs_control.total_profile);
    }
  }
  
  void save_gibbs_profile(){
    if(gibbs_control.save_profile > 0){
      if(iter%gibbs_control.save_profile==0){
        int profile_iter = iter/gibbs_control.save_profile;
        // update_loglik();
        update_logpost();
        // gibbs_profile.loglik(profile_iter) = paras.loglik;
        gibbs_profile.logpost(profile_iter) = paras.logpost;
      }
    }
  }
  
  void monitor_gibbs(){
    if(gibbs_control.verbose > 0){
      if(iter%gibbs_control.verbose==0){
        std::cout << "iter: " << iter <<  " logpost: "<< paras.logpost << std::endl;
      }
    }
  }
  
  void run_gibbs(){
    initialize_paras_sample();
    initialize_gibbs_profile();
    //std::cout << "total iter:" << gibbs_control.total_iter << std::endl;
    if (get_hs() == 0) {
      for(iter=0; iter<gibbs_control.total_iter;iter++){
        
        //std::cout << "iter:" << iter << std::endl;
        
        //std::cout << "update_inv_sigma2_alpha:"<<std::endl;
        update_inv_sigma2_alpha();
        //std::cout << "update_alpha:"<<std::endl;
        update_alpha();
        
        //std::cout << "update_inv_a_beta:"<<std::endl;
        update_inv_a_beta();
        
        //std::cout << "update_inv_sigma2_beta:"<<std::endl;
        update_inv_sigma2_beta();
        
        //std::cout << "update_beta:"<<std::endl;
        update_beta();
        
        //std::cout << "update_omega:"<<std::endl;
        update_omega();
        
        //std::cout << "save save_paras_sample" << std::endl;
        save_paras_sample();
        
        //std::cout << "save save_gibbs_profile" << std::endl;
        save_gibbs_profile();
        
        //std::cout << "save monitor_gibbs" << std::endl;
        monitor_gibbs();
      }
    } else {
      for(iter=0; iter<gibbs_control.total_iter;iter++){
        //std::cout << "iter:" << iter << std::endl;
        //std::cout << "update_inv_sigma2_alpha:"<<std::endl;
        update_inv_sigma2_alpha();
        //std::cout << "update_alpha:"<<std::endl;
        update_alpha();
        //std::cout << "update_inv_a_tau:"<<std::endl;
        update_inv_a_tau();
        //std::cout << "update_inv_tau2:"<<std::endl;
        update_inv_tau2();
        //std::cout << "update_inv_a_beta:"<<std::endl;
        update_inv_a_beta();
        //std::cout << "update_inv_sigma2_beta:"<<std::endl;
        update_inv_sigma2_beta();
        //std::cout << "update_beta:"<<std::endl;
        update_beta();
        
        //std::cout << "update_omega:"<<std::endl;
        update_omega();
        
        //std::cout << "save save_paras_sample" << std::endl;
        save_paras_sample();
        
        //std::cout << "save save_gibbs_profile" << std::endl;
        save_gibbs_profile();
        
        //std::cout << "save monitor_gibbs" << std::endl;
        monitor_gibbs();
      }
    }
    
  };
  
  
  
  //output for R
  List get_gibbs_post_mean(){
    //std::cout << "get_gibbs_post_mean" << std::endl;
    arma::vec beta = mean(paras_sample.beta,1);
    arma::cube mean_alpha = mean(paras_sample.alpha,2);
    arma::cube mean_inv_sigma2_alpha = mean(paras_sample.inv_sigma2_alpha,2);
    return List::create(Named("beta") = beta,
                        Named("alpha") = mean_alpha.slice(0),
                        Named("inv_sigma2_beta") = mean(paras_sample.inv_sigma2_beta,1),
                        Named("inv_sigma2_alpha") = mean_inv_sigma2_alpha.slice(0),
                        Named("inv_a_beta") = mean(paras_sample.inv_a_beta,1),
                        Named("inv_tau2") = mean(paras_sample.inv_tau2),
                        Named("inv_a_tau") = mean(paras_sample.inv_a_tau));
  };
  
  List get_gibbs_sample(){
    //std::cout << "get_gibbs_sample" << std::endl;
    return List::create(Named("beta") = paras_sample.beta);
    // Named("alpha") = paras_sample.alpha,
    // Named("inv_sigma2_beta") = paras_sample.inv_sigma2_beta,
    // Named("inv_sigma2_alpha") = paras_sample.inv_sigma2_alpha,
    // Named("inv_a_beta") = paras_sample.inv_a_beta,
    // Named("inv_tau2") = paras_sample.inv_tau2,
    // Named("inv_a_tau") = paras_sample.inv_a_tau);
  };
  
  List get_gibbs_trace(){
    //std::cout << "get_gibbs_trace" << std::endl;
    uvec iters = linspace<uvec>(1,gibbs_control.total_iter,gibbs_control.total_profile);
    return List::create(Named("iters") = iters,
                        // Named("loglik") = gibbs_profile.loglik,
                        Named("logpost") = gibbs_profile.logpost);
  }
  
  List get_gibbs_control(){
    //std::cout << "get_gibbs_control" << std::endl;
    return List::create(Named("total_iter") = gibbs_control.total_iter,
                        Named("burnin") = gibbs_control.burnin,
                        Named("mcmc_sample") = gibbs_control.mcmc_sample,
                        Named("thinning") = gibbs_control.thinning,
                        Named("verbose") = gibbs_control.verbose,
                        Named("save_profile") = gibbs_control.save_profile,
                        Named("total_profile") = gibbs_control.total_profile);
  }
  
  List get_gibbs_data(){
    return List::create(Named("X") = dat.X,
                        Named("V") = dat.V,
                        Named("A") = dat.A,
                        Named("nn") = dat.nn);
  };
  
  int get_iter(){
    return iter;
  };
};


// [[Rcpp::export]]
List Gibbs_cpp(arma::mat& A, arma::mat& X, arma::vec& V,  arma::vec& nn,
                    arma::vec initial_beta,
                    arma::mat initial_alpha,
                    arma::vec initial_inv_sigma2_beta,
                    arma::mat initial_inv_sigma2_alpha,
                    arma::vec initial_inv_a_beta,
                    double initial_inv_tau2,
                    double initial_inv_a_tau,
                    double a_alpha = 0.1,
                    double b_alpha = 0.1,
                    double nu_0 = 1,
                    double eta2_0 = 1,
                    double nu_1 = 1,
                    double eta2_1 = 1,
                    bool horseshoe = true,
                    int mcmc_sample = 500,
                    int burnin = 5000,
                    int thinning = 10,
                    double paras_diff_tol = 1e-6,
                    int verbose = 5000,
                    int save_profile = 1){
  
  wall_clock timer;
  timer.tic();
  Gibbs_VAX model;
  
  
  model.load_data_Gibbs(A,X,V, nn);
  model.set_hyperparas_Gibbs(a_alpha, b_alpha, nu_0, eta2_0, nu_1, eta2_1);
  
  
  
  model.set_gibbs_control(mcmc_sample,
                          burnin,
                          thinning,
                          verbose,
                          save_profile);
  model.set_hs(horseshoe);
  
  // std::cout << "set control done" << std::endl;
  
  model.set_paras_initial_values_Gibbs(initial_beta, initial_alpha,
                                       initial_inv_sigma2_beta,
                                       initial_inv_sigma2_alpha,
                                       initial_inv_a_beta,
                                       initial_inv_tau2,
                                       initial_inv_a_tau);
  
  // std::cout << "set initial values" << std::endl;
  
  model.run_gibbs();
  // std::cout << "run finish" << std::endl;
  double elapsed = timer.toc();
  
  List output;
  
  
  output = List::create(Named("post_mean") = model.get_gibbs_post_mean(),
                        Named("mcmc") = model.get_gibbs_sample(),
                        Named("trace") = model.get_gibbs_trace(),
                        Named("data") = model.get_gibbs_data(),
                        Named("mcmc_control") = model.get_gibbs_control(),
                        Named("elapsed") = elapsed);
  //std::cout << "output finish" << std::endl;
  
  return output;
  
}


