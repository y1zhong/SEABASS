#'@title get glm estimates
#'@param data list with AE count (A), confounders (X), vaccine indication (V), number of each categories
#'@returns list of glm estimated beta, alpha coefficients with their p-values and adsusted q-values
#'@export
get_glm_coefs = function(data) {
  A=data$A
  X=data$X
  V=data$V
  nn=data$nn
  if(is.null(nn)) nn = rep(1, nrow(X))
  Alpha_glm = matrix(0, ncol(A), ncol(X) + 1)
  glm_coef = matrix(0, ncol(A), ncol(X) + 1)
  p_value = rep(0, ncol(X))
  b_ind = ncol(Alpha_glm)
  se = rep(0, ncol(X))
  suppressWarnings({
    for (s in 1:ncol(A)) {
      As = A[, s]
      model = glm(cbind(As, nn-As)  ~ 0 + X + V,
                  family = binomial())
      alpha_glm = coef(model)
      glm_coef[s,] =  alpha_glm
      se[s] = summary(model)$coefficients[b_ind, 2]
      p_value[s] = coef(summary(model))[, 4][b_ind]
    }
  })
  
  rownames(Alpha_glm) = NULL
  beta_glm = glm_coef[, b_ind]
  alpha_glm = glm_coef[,-b_ind, drop=F]
  colnames(alpha_glm) = colnames(X)
  q_value = p.adjust(p_value,"BH")
  
  return(list(beta= beta_glm,
              alpha = t(alpha_glm),
              p_value = p_value,
              q_value=q_value,
              se = se))
}

#'@title cut coefficients 
#'@param coefs list with coefficients
#'@param beta_cutoff numeric 
#'@param alpha_cutoff numeric
#'@returns list of glm estimated beta, alpha coefficients with their p-values and adsusted q-values
#'@export
coef_cutoff <- function(coefs, beta_cutoff = 4, alpha_cutoff=12){
  beta <- coefs$beta
  alpha <- coefs$alpha
  
  beta <- pmax(beta, -beta_cutoff)
  beta <- pmin(beta, beta_cutoff)
  
  alpha <- pmax(alpha, -alpha_cutoff)
  alpha <- pmin(alpha, alpha_cutoff)
  
  coefs$beta <- beta
  coefs$alpha <- alpha
  
  return (coefs)
}



#'@import dplyr
gen_data <- function(beta, alpha, N, tps, vax_p=0.5
                     ){
  J = length(beta)
  p = nrow(alpha)
  
  seq_data <- lapply(1:tps, function(s){
    X = data.frame(intercept = 1, 
                   Age =  sample(c(0,1,2), N, prob = c(0.3,0.3,0.4),replace=T),
                   Male = rbinom(N, 1, 0.3),
                   V = rbinom(N, 1, vax_p))
    counts <- count(X, intercept,Age,Male, V)
    V = counts$V
    X = as.matrix(counts[, c(1:3)])
    nn = counts$n
    
    PA = X %*% alpha +  V %*% t(beta)
    PA =  plogis(PA)
    
    A = rbinom(length(PA), size = nn, prob = c(PA))
    A = matrix(A, nrow(PA), ncol(PA))
    
    list(A = A,
         X =  X,
         V = V,
         nn = nn)
  })
  
  for (s in 2:tps){
    seq_data[[s]]$A = seq_data[[s]]$A + seq_data[[s-1]]$A
    seq_data[[s]]$nn = seq_data[[s]]$nn + seq_data[[s-1]]$nn
  }
  
  return(seq_data)
}
