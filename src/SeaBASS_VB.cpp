#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

Environment qgam_base = Environment::namespace_env("qgam");
Function log1pexp_base = qgam_base["log1pexp"];

class seabass_VB{
private:
  int hs; // 0 for base model, 1 for horseshoe prior
  struct seabassData{
    int num_reports;
    int num_AEs;
    int num_confounders;
  
    arma::mat A;
    arma::mat X;
    arma::vec V;
    arma::vec nn;
  } dat;

  struct HyperParas{
    double a_alpha;
    double b_alpha;
  
    arma::vec mu_beta;
    arma::mat mu_alpha;
    arma::vec mu_beta_sq;
    arma::mat mu_alpha_sq;
  
    double nu_0;
    double eta2_0;
  
    double nu_1;
    double eta2_1;
  } hyperparas;

  
  struct VBParas{
    arma::vec E_beta;
    arma::vec Var_beta;
    arma::vec E_beta_sq;
    arma::vec beta_E_minus_mu_sq;
  
    arma::vec E_inv_sigma2_beta;
    arma::vec E_inv_a_beta;
  
    double E_inv_tau2;
    double E_inv_a_tau;
  
    arma::mat E_alpha;
    arma::mat Var_alpha;
    arma::mat E_alpha_sq;
    arma::mat cov_mat_j;
    arma::cube cov_mat;
    arma::mat E_inv_sigma2_alpha;
    arma::mat alpha_E_minus_mu_sq;
  
    arma::mat E_omega;
    arma::mat E_phi;
  
    arma::vec log_det_cov_alpha;
    arma::mat log_cosh_b;
    double ELBO;
    double loglik;
  
    arma::mat temp;
  } vb_paras;

  
  struct VBProfile{
    arma::vec ELBO;
    arma::vec loglik;
  } vb_profile;
  
  struct VBControl{
    int max_iter;
    double para_diff_tol;
    int ELBO_stop;
    double ELBO_diff_tol;
    int verbose;
    int save_profile;
    int total_profile;
  } vb_control;
  
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

  
  void load_data(const arma::mat& in_A, const arma::mat& in_X, const arma::vec& in_V, const arma::vec& in_nn){
    dat.A = in_A;
    dat.X = in_X;
    dat.V = in_V;
    dat.nn = in_nn;
  
    dat.num_AEs = dat.A.n_cols;
    dat.num_confounders = dat.X.n_cols;
    dat.num_reports = dat.X.n_rows;
  };
  
  void set_hyperparas(const double& in_a_alpha, const double& in_b_alpha,
                      const arma::vec& in_mu_beta, const arma::mat& in_mu_alpha,
                      const double& in_nu_0,const double& in_eta2_0,
                      const double& in_nu_1,const double& in_eta2_1){
    hyperparas.a_alpha = in_a_alpha;
    hyperparas.b_alpha = in_b_alpha;
    
    hyperparas.mu_beta = in_mu_beta;
    hyperparas.mu_beta_sq = in_mu_beta %in_mu_beta;
    hyperparas.mu_alpha = in_mu_alpha;
    hyperparas.mu_alpha_sq = in_mu_alpha % in_mu_alpha;
  
    hyperparas.nu_0 = in_nu_0;
    hyperparas.eta2_0 = in_eta2_0;
    
    hyperparas.nu_1 = in_nu_1;
    hyperparas.eta2_1 = in_eta2_1;
    
  };

  void set_paras_initial_values(const arma::vec& in_beta,
                                const arma::mat& in_alpha,
                                const arma::vec& in_inv_sigma2_beta,
                                const arma::mat& in_inv_sigma2_alpha,
                                const arma::vec& in_inv_a_beta,
                                const double& in_inv_tau2,
                                const double& in_inv_a_tau){
  
    
    vb_paras.E_inv_sigma2_beta = in_inv_sigma2_beta;
    vb_paras.E_beta = in_beta;
    vb_paras.E_inv_a_beta = in_inv_a_beta;
    
  
    vb_paras.E_inv_tau2 = in_inv_tau2;
    vb_paras.E_inv_a_tau = in_inv_a_tau;
  
    vb_paras.Var_beta = 1 / (vb_paras.E_inv_sigma2_beta * vb_paras.E_inv_tau2);
    vb_paras.E_beta_sq = vb_paras.E_beta % vb_paras.E_beta + vb_paras.Var_beta;
    vb_paras.beta_E_minus_mu_sq = vb_paras.E_beta_sq + hyperparas.mu_beta_sq - 2 * vb_paras.E_beta % hyperparas.mu_beta;
  
    vb_paras.E_inv_sigma2_alpha = in_inv_sigma2_alpha;
    vb_paras.Var_alpha = 1 / vb_paras.E_inv_sigma2_alpha;
    vb_paras.cov_mat.zeros(dat.num_confounders, dat.num_confounders, dat.num_AEs);
    for(int j=0; j<dat.num_AEs; j++){
      vb_paras.cov_mat.slice(j).diag() = vb_paras.Var_alpha.col(j);
    }
  
    vb_paras.E_alpha = in_alpha;
    vb_paras.E_alpha_sq = vb_paras.Var_alpha + vb_paras.E_alpha % vb_paras.E_alpha;
    vb_paras.alpha_E_minus_mu_sq = vb_paras.E_alpha_sq + hyperparas.mu_alpha_sq - 2 * vb_paras.E_alpha % hyperparas.mu_alpha;
  
    vb_paras.E_phi.zeros(dat.num_reports, dat.num_AEs);
    vb_paras.E_omega.zeros(dat.num_reports, dat.num_AEs);
  
    update_E_phi();
    update_E_omega();
    
    vb_paras.log_det_cov_alpha.zeros(dat.num_AEs);
  };
  
  void set_vb_control(int in_max_iter,
                      double in_para_diff_tol,
                      int in_ELBO_stop,
                      double in_ELBO_diff_tol,
                      int in_verbose,
                      int in_save_profile){
    vb_control.max_iter = in_max_iter;
    vb_control.para_diff_tol = in_para_diff_tol;
    vb_control.ELBO_stop = in_ELBO_stop;
    vb_control.ELBO_diff_tol = in_ELBO_diff_tol;
    vb_control.verbose = in_verbose;
    vb_control.save_profile = in_save_profile;
    if(vb_control.save_profile > 0){
      vb_control.total_profile = vb_control.max_iter/vb_control.save_profile;
    } else{
      vb_control.total_profile = 0;
    }
  };

  
  void update_E_beta(){
    for(int j=0; j<dat.num_AEs; j++){
      vb_paras.Var_beta(j) = 1 / (accu(dat.V % vb_paras.E_omega.col(j)) + vb_paras.E_inv_sigma2_beta(j) * vb_paras.E_inv_tau2);
      vb_paras.E_beta(j) = accu((dat.A.col(j) - dat.nn/2 - vb_paras.E_omega.col(j) % (dat.X * vb_paras.E_alpha.col(j))) % dat.V);
      vb_paras.E_beta(j) += hyperparas.mu_beta(j) * vb_paras.E_inv_sigma2_beta(j) * vb_paras.E_inv_tau2;
      vb_paras.E_beta(j) *= vb_paras.Var_beta(j);
      vb_paras.E_beta_sq(j) = vb_paras.E_beta(j) * vb_paras.E_beta(j) + vb_paras.Var_beta(j);
   }
    
    vb_paras.beta_E_minus_mu_sq = vb_paras.E_beta_sq + hyperparas.mu_beta_sq - 2 * vb_paras.E_beta % hyperparas.mu_beta;
  };

  
  void update_E_inv_sigma2_beta(){
    vb_paras.E_inv_sigma2_beta = (0.5 * (1 + hyperparas.nu_0)) / (0.5 * vb_paras.beta_E_minus_mu_sq * vb_paras.E_inv_tau2 + hyperparas.nu_0 * vb_paras.E_inv_a_beta);
  };

  void update_E_inv_sigma2_alpha(){
    vb_paras.E_inv_sigma2_alpha = (hyperparas.a_alpha + 0.5) / (0.5 * vb_paras.alpha_E_minus_mu_sq + hyperparas.b_alpha);
  };

  void update_E_alpha(){
    for(int j=0; j<dat.num_AEs; j++){
      vb_paras.cov_mat_j = dat.X.each_col() % vb_paras.E_omega.col(j);
      vb_paras.cov_mat_j =  vb_paras.cov_mat_j.t() * dat.X;
      vb_paras.cov_mat_j.diag() += vb_paras.E_inv_sigma2_alpha.col(j);
    
      vb_paras.cov_mat_j = inv(vb_paras.cov_mat_j);
      vb_paras.cov_mat.slice(j) = vb_paras.cov_mat_j;
      vb_paras.Var_alpha.col(j) = vb_paras.cov_mat_j.diag();
      vb_paras.E_alpha.col(j) = vb_paras.cov_mat_j * (dat.X.t() * (dat.A.col(j) - dat.nn/2 - vb_paras.E_omega.col(j) % dat.V*vb_paras.E_beta(j)) + vb_paras.E_inv_sigma2_alpha.col(j) % hyperparas.mu_alpha.col(j));
    }
    vb_paras.E_alpha_sq = vb_paras.Var_alpha + vb_paras.E_alpha % vb_paras.E_alpha;
    vb_paras.alpha_E_minus_mu_sq = vb_paras.E_alpha_sq + hyperparas.mu_alpha_sq - 2 * vb_paras.E_alpha % hyperparas.mu_alpha;
  };
  
  
  
  void update_E_inv_a_beta(){
    vb_paras.E_inv_a_beta = ((hyperparas.nu_0 + 1) / 2) / (hyperparas.nu_0 * vb_paras.E_inv_sigma2_beta + 1 / hyperparas.eta2_0);
  };
  
  void update_E_inv_tau2(){
    vb_paras.E_inv_tau2 = 0.5 * (dat.num_AEs + hyperparas.nu_1) / (0.5 * accu(vb_paras.beta_E_minus_mu_sq % vb_paras.E_inv_sigma2_beta) + hyperparas.nu_1 * vb_paras.E_inv_a_tau);
    
  };
  
  void update_E_inv_a_tau(){
    vb_paras.E_inv_a_tau = ((hyperparas.nu_1 + 1)/2) / (hyperparas.nu_1 * vb_paras.E_inv_tau2 + 1 / hyperparas.eta2_1);
  };
  
  void update_E_phi(){
    for(int i=0; i<dat.num_reports; i++){
      for(int j=0; j<dat.num_AEs; j++){
        vb_paras.temp = (dat.X.row(i).as_col() * dat.X.row(i)) % (vb_paras.cov_mat.slice(j) + vb_paras.E_alpha.col(j) * vb_paras.E_alpha.col(j).t());
        vb_paras.E_phi(i,j) = accu(vb_paras.temp) + dat.V(i)*vb_paras.E_beta_sq(j) + 2 * accu(dat.X.row(i)*vb_paras.E_alpha.col(j)) * dat.V(i) * vb_paras.E_beta(j);
        vb_paras.E_phi(i,j) = sqrt(vb_paras.E_phi(i,j));
      }
    }
  }
  
  void update_E_omega(){
    vb_paras.temp = tanh(vb_paras.E_phi/2) / (2*vb_paras.E_phi);
    vb_paras.E_omega = vb_paras.temp.each_col() % dat.nn;
  };
  
  void update_log_cosh_b(){
    vb_paras.temp = log(cosh(vb_paras.E_phi / 2));
    vb_paras.log_cosh_b = vb_paras.temp.each_col() % dat.nn;
  };
  
  
  void update_log_det_cov_mat(){
    for(int j=0; j<dat.num_AEs; j++){
      vb_paras.cov_mat_j = dat.X.each_col() % vb_paras.E_omega.col(j);
      vb_paras.cov_mat_j =  vb_paras.cov_mat_j.t() * dat.X;
      vb_paras.cov_mat_j.diag() += vb_paras.E_inv_sigma2_alpha.col(j);
      vb_paras.log_det_cov_alpha(j) = log_det_sympd(vb_paras.cov_mat_j);
    }
  }
  
  
  void update_ELBO(){
    
    vb_paras.temp = dat.X * vb_paras.E_alpha + dat.V * vb_paras.E_beta.t();
    vb_paras.ELBO = accu((dat.A.each_col() - dat.nn/2) % vb_paras.temp);
    
    vb_paras.ELBO += accu(log(vb_paras.Var_beta)) / 2;
    vb_paras.ELBO -= ((1+hyperparas.nu_0)/2) * accu(log(0.5 * vb_paras.beta_E_minus_mu_sq * vb_paras.E_inv_tau2 + hyperparas.nu_0 * vb_paras.E_inv_a_beta));
    vb_paras.ELBO += hyperparas.nu_0 * accu(vb_paras.E_inv_a_beta % vb_paras.E_inv_sigma2_beta);
    
    update_log_det_cov_mat();
    vb_paras.ELBO -= accu(vb_paras.log_det_cov_alpha) / 2;
    vb_paras.ELBO -= accu(log(hyperparas.b_alpha + 0.5 * vb_paras.alpha_E_minus_mu_sq)) * (hyperparas.a_alpha + 0.5);
    vb_paras.ELBO -= (hyperparas.nu_0 * dat.num_AEs + 1) / 2 * log(hyperparas.nu_0 * accu(vb_paras.E_inv_sigma2_beta) + 1 / hyperparas.eta2_0);
    
    update_log_cosh_b();
    vb_paras.ELBO -= accu(vb_paras.log_cosh_b);
    
    if(hs == 1){
      vb_paras.ELBO -= ((dat.num_AEs + hyperparas.nu_1)/2) * log(0.5 * accu(vb_paras.beta_E_minus_mu_sq % vb_paras.E_inv_sigma2_beta) + hyperparas.nu_1*vb_paras.E_inv_a_tau);
      vb_paras.ELBO += 0.5 * vb_paras.E_inv_tau2 * accu(vb_paras.beta_E_minus_mu_sq % vb_paras.E_inv_sigma2_beta);
      vb_paras.ELBO -= (hyperparas.nu_1 + 1) / 2 * log(hyperparas.nu_1 * vb_paras.E_inv_tau2 + 1 / hyperparas.eta2_1);
      vb_paras.ELBO +=  vb_paras.E_inv_a_tau * vb_paras.E_inv_tau2 * hyperparas.nu_1;
    }
    vb_paras.loglik = accu(dat.A % vb_paras.E_phi) - accu(as<arma::mat>(log1pexp_base(vb_paras.E_phi)));
  };
  
  
  bool compute_paras_diff(double& para_diff, 
                          arma::vec& beta, arma::vec& beta_prev,
                          arma::mat& alpha, arma::mat& alpha_prev){
    double beta_temp = accu((beta - beta_prev) > para_diff);
    double alpha_temp = accu((alpha - alpha_prev) > para_diff);
    return(beta_temp == 0 && alpha_temp == 0);
  };
  
  
  
  void initialize_vb_profile(){
    if(vb_control.save_profile>0){
      vb_profile.ELBO.zeros(vb_control.total_profile);
      vb_profile.loglik.zeros(vb_control.total_profile);
    }
  }
  
  void save_vb_profile(){
    if(vb_control.save_profile > 0){
      if(iter%vb_control.save_profile==0){
        int profile_iter = iter/vb_control.save_profile;
        if(vb_control.ELBO_stop==0){
          update_ELBO();
        }
        vb_profile.ELBO(profile_iter) = vb_paras.ELBO;
        vb_profile.loglik(profile_iter) = vb_paras.loglik;
      }
    }
  }
  
  void monitor_vb(){
    if(vb_control.verbose > 0){
      if(iter%vb_control.verbose==0){
        if(vb_control.ELBO_stop==0){
          update_ELBO();
        }
        std::cout << "iter: " << iter <<  " ELBO: "<< vb_paras.ELBO << std::endl;
      }
    }
  }
  
  
  
  void run_mfvb(){
    initialize_vb_profile();
    for(iter=0; iter<vb_control.max_iter; iter++){
      arma::vec E_beta_prev = vb_paras.E_beta;
      arma::mat E_alpha_prev = vb_paras.E_alpha;
      
      update_E_alpha();
      update_E_inv_sigma2_alpha();
      if(hs == 1){
        update_E_inv_a_tau();
        update_E_inv_tau2();
      }
      
      update_E_beta();
      update_E_inv_sigma2_beta();
      update_E_inv_a_beta();
      
      update_E_phi();
      update_E_omega();
      
      
      if(vb_control.ELBO_stop == 0){
        if(compute_paras_diff(vb_control.para_diff_tol, vb_paras.E_beta, E_beta_prev,
                              vb_paras.E_alpha, E_alpha_prev)){
          update_ELBO();
          save_vb_profile();
          monitor_vb();
          break;
        }
      } else {
        double ELBO_prev = vb_paras.ELBO;
        update_ELBO();
        if(abs(vb_paras.ELBO - ELBO_prev) < vb_control.ELBO_diff_tol){
          save_vb_profile();
          monitor_vb();
          break;
        }
      }
      
      save_vb_profile();
      monitor_vb();
    }
  };
  
  List get_vb_data(){
    return List::create(Named("X") = dat.X,
                        Named("V") = dat.V,
                        Named("A") = dat.A,
                        Named("nn") = dat.nn);
  };
  
  List get_vb_hyperparam(){
    return List::create(Named("a_alpha") = hyperparas.a_alpha,
                        Named("b_alpha") = hyperparas.b_alpha,
                        Named("mu_alpha") = hyperparas.mu_alpha,
                        Named("mu_beta") = hyperparas.mu_beta,
                        Named("nu_0") = hyperparas.nu_0,
                        Named("eta2_0") = hyperparas.eta2_0,
                        Named("nu_1") = hyperparas.nu_1,
                        Named("eta2_1") = hyperparas.eta2_1);
  };
  
  List get_vb_post_mean(){
    return List::create(Named("beta") = vb_paras.E_beta,
                        Named("alpha") = vb_paras.E_alpha,
                        Named("inv_sigma2_beta") = vb_paras.E_inv_sigma2_beta,
                        Named("inv_sigma2_alpha") = vb_paras.E_inv_sigma2_alpha,
                        Named("inv_a_beta") = vb_paras.E_inv_a_beta,
                        Named("inv_tau2") = vb_paras.E_inv_tau2,
                        Named("inv_a_tau") = vb_paras.E_inv_a_tau,
                        Named("beta_var") = vb_paras.Var_beta);
  };
  
  
  List get_vb_trace(){
    int actual_profile_iter = 1;
    if(iter == 0){
      iter = 1;
    }
    if(vb_control.save_profile>0){
      actual_profile_iter = iter/vb_control.save_profile;
    }
    arma::uvec iters = linspace<uvec>(1,iter,actual_profile_iter);
    return List::create(Named("iters") = iters,
                        Named("ELBO") = vb_profile.ELBO.rows(0,actual_profile_iter-1),
                        Named("loglik") = vb_profile.loglik.rows(0,actual_profile_iter-1));
    
  }
  
  List get_vb_control(){
    return List::create(Named("max_iter")= vb_control.max_iter,
                        Named("para_diff_tol") = vb_control.para_diff_tol,
                        Named("ELBO_stop") = vb_control.ELBO_stop,
                        Named("ELBO_diff_tol") = vb_control.ELBO_diff_tol,
                        Named("verbose") = vb_control.verbose,
                        Named("save_profile") = vb_control.save_profile,
                        Named("total_profile") = vb_control.total_profile);
  };
  
  int get_iter(){
    return iter;
  };
};


// ' @useDynLib VBVaccine, .registration = TRUE
// ' @importFrom Rcpp evalCpp
// ' @export
// [[Rcpp::export]]
List SeaBASS_cpp(arma::mat& A, arma::mat& X, arma::vec& V,  arma::vec& nn,
                 arma::vec initial_beta,
                 arma::mat initial_alpha,
                 arma::vec initial_inv_sigma2_beta,
                 arma::mat initial_inv_sigma2_alpha,
                 arma::vec initial_inv_a_beta,
                 double initial_inv_tau2,
                 double initial_inv_a_tau,
                 arma::vec mu_beta,
                 arma::mat mu_alpha,
                 bool horseshoe = false,
                 double a_alpha = 0.1,
                 double b_alpha = 0.1,
                 double nu_0 = 1,
                 double eta2_0 = 1,
                 double nu_1 = 1,
                 double eta2_1 = 1,
                 int max_iter = 1000,
                 double paras_diff_tol = 1e-6,
                 int ELBO_stop = 1,
                 double ELBO_diff_tol = 1e-6,
                 int verbose = 0,
                 int save_profile = 1){
  
  wall_clock timer;
  timer.tic();
  seabass_VB model;
  
  model.load_data(A,X,V, nn);
  model.set_hyperparas(a_alpha, b_alpha, mu_beta, mu_alpha, nu_0, eta2_0, nu_1, eta2_1);
  
  model.set_vb_control(max_iter,
                       paras_diff_tol,
                       ELBO_stop,
                       ELBO_diff_tol,
                       verbose,
                       save_profile);
  model.set_hs(horseshoe);
  // std::cout << "set control done" << std::endl;
  
  model.set_paras_initial_values(initial_beta, initial_alpha,
                                 initial_inv_sigma2_beta,
                                 initial_inv_sigma2_alpha,
                                 initial_inv_a_beta,
                                 initial_inv_tau2,
                                 initial_inv_a_tau);
  
  // std::cout << "set initial values" << std::endl;
  
  model.run_mfvb();
  
  
  double elapsed = timer.toc();
  
  List output;
  
  output = List::create(Named("post_mean") = model.get_vb_post_mean(),
                        Named("data") = model.get_vb_data(),
                        Named("hyper") = model.get_vb_hyperparam(),
                        Named("iter") = model.get_iter(),
                        Named("trace") = model.get_vb_trace(),
                        Named("vb_control") = model.get_vb_control(),
                        Named("elapsed") = elapsed);
  
  return output;
  
}


