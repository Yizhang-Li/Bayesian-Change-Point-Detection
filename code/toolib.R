# tool fun ----
print_loading = function(x){
  #x:string
  print(paste('loading',x,'......'))
}

UT_flag = 0

# 1. recursion filter theta_t | y_t ----
  # y ~ normal N(\theta_t,1)
  # \theta ~ distribution \pi (a_0,\mu_0) = N(\mu_0,1/a_0)

source('UT_data.R')
recursion_pit = function(yt,p, mu,a0 = 1){
  # return dataframe pit
  constant_a_mu = function(a,mu){
    # c(a,mu) = 1/sqrt(2pi/a) * exp{-0.5a * mu^2}
    c = sqrt(a/(2*pi))*exp(-0.5*a*mu*mu) 
    return(c)
  }
  
  n = length(yt) 
  P_df = matrix(0,n,n)
  for (j in 1:n){
    for (i in 1:j){
      if (i == j){
        mu_tt = (a0*mu + yt[i])/(a0+1)
        P_df[i,j] = p * constant_a_mu(a = a0,mu = mu)/constant_a_mu(a0+1,mu_tt)
      }
      else{
        a_ij_1 = a0 + j - i
        mu_ij_1 = (a0*mu + sum(yt[i:(j-1)]))/(a0+j-i)
        mu_ij = (a0*mu + sum(yt[i:j]))/(a0+j-i+1)
        P_df[i,j] = (1-p) * P_df[i,j-1] * constant_a_mu(a_ij_1,mu_ij_1) / constant_a_mu(a_ij_1+1,mu_ij)
      }
    }
    sum_j = sum(P_df[,j])
    P_df[,j] = P_df[,j]/sum_j
  }
  return(P_df)
}

if (UT_flag==1){
  p_df = recursion_pit(yt,p,mu)
}

# 2.1 recursion filter theta_t | y_{t_11,n} ----
recursion_qjt = function(yt,p,mu,a0 = 1){
  # return dataframe pit
  constant_a_mu = function(a,mu){
    c = sqrt(a/(2*pi))*exp(-0.5*a*mu*mu)
    return(c)
  }
  
  n = length(yt) 
  Q_df = matrix(0,n,n)
  for (j in n:1){
    for (i in n:j){
      if (i == j){
        mu_tt = (a0*mu + yt[i])/(a0+1)
        Q_df[i,j] = p * constant_a_mu(a = a0,mu = mu)/constant_a_mu(a0+1,mu_tt)
      }
      else{
        a_j1_i = a0 + i - j
        mu_j1_i = (a0*mu + sum(yt[(j+1):i]))/a_j1_i
        mu_j_i = (a0*mu + sum(yt[j:i]))/(a_j1_i+1)
        Q_df[i,j] = (1-p) * Q_df[i,j+1] * constant_a_mu(a_j1_i,mu_j1_i) / constant_a_mu(a_j1_i+1,mu_j_i)
      }
    }
    sum_j = sum(Q_df[,j])
    Q_df[,j] = Q_df[,j]/sum_j
  }
  return(Q_df)
}

if (UT_flag==1){
  q_df = recursion_qjt(yt,p,mu)
}

# 2.2 recursion filter mu_t | y_{t_1,n} ----
beta_3d = function(yt,p_df,q_df,p,mu,a0 = 1){
  # return list
  # beta_arrary:3d-arrary beta_ijk, normalized_beta, P_list
  result = list()
  
  constant_a_mu = function(a,mu){
    c = sqrt(a/(2*pi))*exp(-0.5*a*mu*mu)
    return(c)
  }
  
  n = length(yt) 
  beta_arrary = array(0,dim=c(n,n,n))
  
  for (i in 1:n){
    for (j in 1:n){
      for (k in 1:n){
        if ((i<=k)&(k==j)){
          beta_arrary[i,j,k] = p*p_df[i,k]
          }else if ((i<=k)&(k<j)) {
          Numerator = constant_a_mu(a0+k-i+1,(a0*mu + sum(yt[i:k]))/(a0+k-i+1)) * constant_a_mu(a0+j-k,(a0*mu + sum(yt[(k+1):j]))/(a0+j-k))
          Denominator = constant_a_mu(a0+j-i+1,(a0*mu + sum(yt[i:j]))/(a0+j-i+1)) * constant_a_mu(a0,mu)
          beta_arrary[i,j,k] = (1-p) * p_df[i,k] * q_df[j,k+1] * Numerator/Denominator
          }else{
          beta_arrary[i,j,k] = 0}
      }
    }
  }
  
  result$beta_arrary = beta_arrary
  
  P_list = rep(0,n)
  for (k in 1:n){
    tmp = 0
    for (i in 1:n){
      for (j in 1:n){
        if ((i<=k) & (k<j)){
          tmp = tmp + beta_arrary[i,j,k]
        }
      }
    }
    tmp = tmp + p
    P_list[k] = tmp
  }
  result$P_list = P_list
  
  normalized_beta = beta_arrary
  for (k in 1:n){
    normalized_beta[,,k] = beta_arrary[,,k] / P_list[k]
  }
  result$normalized_beta = normalized_beta
  
  return(result)
}

if (UT_flag==1){
  beta_result = beta_3d(yt,p_df,q_df,p,mu)
  normalized_beta = beta_result$normalized_beta
}


# 2.3 expectation \mu_t | y_{1,n}
mu_t_n = function(yt,normalized_beta,mu,a0=1){
  # return mu_1,...,mu_n
  n = length(yt)
  mu_list = rep(0,n)
  
  y_bar_ij = function(yt,mu,a0,i,j){
    y_bar = (a0*mu + sum(yt[i:j])) / (a0 + j -i + 1)
    return(y_bar)
  }

  for (k in 1:n){
    tmp = 0
    for (i in 1:n){
      for (j in 1:n){
        if((i<=k) & (k<=j)){
          tmp = tmp + normalized_beta[i,j,k] * y_bar_ij(yt,mu,a0,i,j)
        }
      }
    }
    mu_list[k] = tmp
  }
  return(mu_list)
}

if (UT_flag==1){
  mu_list = mu_t_n(yt,normalized_beta,mu)
  par(mfrow=c(1,2))
  plot(yt,main = 'y')
  plot(mu_list,main = 'mu_est')
}

# 3. p_argmax ----
p_argmax = function(yt,mu,j0,j1,a0=1){
  #return list
  result = list()
  t = length(yt)
  if (j0 <= j1){
    j_list = j0:j1
  }else{
    j_list = j1:j0
  }
  log_likelihood_list = rep(0,length(j_list))
  for (j_id in 1:length(j_list)){
    j = j_list[j_id]
    p = min(1,2^j/t)
    p_df = recursion_pit(yt,p,mu,a0 = a0)
    l_p = 0
    for (k in 1:t){
      l_p_k = 0 
      for (i in 1:k){
        l_p_k = l_p_k + p_df[i,k]
      }
      l_p_k = log(l_p_k)
      l_p = l_p + l_p_k
    }
    log_likelihood_list[j_id] = l_p
  }
  #print(log_likelihood_list)  
  best_j = j_list[which.max(log_likelihood_list)]
  best_p = min(1,2^best_j/t)
  result$best_j = best_j
  result$best_p = best_p
  result$log_likelihood = max(log_likelihood_list)
  return(result)
}

if (UT_flag==1){
  p_best = p_argmax(yt,mu,-5,5)$best_p
}

# 4. bootstrap_no_cpts
bootsrap_test = function(yt,B,mu,alpha=0.05){
  # yt ~ N(\theta,1)
  result = list()
  L0_result = p_argmax(yt,mu,-5,5)
  L0 = L0_result$log_likelihood
  p = L0_result$best_p
  
  L0 = L0 - sum(log(pnorm(yt,mean = mu,sd = 1)))
  
  L_list = rep(0,B)
  for (b in 1:B){
    sim_yt = rnorm(length(yt),mean = mu,sd = 1)
    L_result = p_argmax(sim_yt,mu,-5,5)
    L = L_result$log_likelihood
    L = L0 - sum(log(pnorm(sim_yt,mean = mu,sd = 1)))
    
    L_list[b] = L
  }
  alpha_hat = length(L_list[L_list>L0])/B
  result$alpha_hat = alpha_hat
  if(alpha_hat < alpha){
    print('Reject H0:no cpts')
  }
  return(result)
}

#bootsrap_np_cpts = bootsrap_test(yt,20,mu)

# 5. Bayesian_cp
Bayesian_kcp = function(yt,best_p,cp_num = 10,beta_c = 1,plot_flag = FALSE,d=10){
  # return:result --cp_list
  result = list()
  
  # b(p) = \beta m(p) = \beta abs(log(p))^(1+eps)
  eps = 0.001
  bp = max(round(beta_c * abs(log(best_p))^(1+eps),0),1)
  
  # mu_list
  mu = mean(yt)
  p_df = recursion_pit(yt,best_p,mu)
  q_df = recursion_qjt(yt,best_p,mu)
  beta_result = beta_3d(yt,p_df,q_df,best_p,mu)
  normalized_beta = beta_result$normalized_beta
  mu_list = mu_t_n(yt,normalized_beta,mu)
  result$mu_list = mu_list
  if (plot_flag == TRUE){
    png('filter_y.png')
    plot(mu_list,main = 'filter_y')
    dev.off()
  }
  
  # location
  t = length(yt)
  valid_location = 1:t
  valid_location = valid_location[(bp+1):(t-bp-1)]
  
  cp_location_list = rep(1,cp_num)
  for (k in 1:cp_num){
    delta_list = rep(0,length(valid_location))
    for (ind in 1:length(valid_location)){
      location = valid_location[ind]
      start_id = location-bp
      end_id = location+bp
      delta_location = (mu_list[end_id] - mu_list[start_id])^2
      delta_list[ind] = delta_location
    }
    #print(valid_location)
    #print(which.max(delta_list))
    if (length(valid_location)==0){
      cp_location_list = cp_location_list[1:(k-1)]
      break
    }else{
      cp_location = valid_location[which.max(delta_list)]
      cp_location_list[k] = cp_location
      
      tmp = array(0,c(length(valid_location),k))
      for (col in 1:k){
        tmp[,col] = abs(valid_location-cp_location_list[col])
      }
      
      tmp = apply(tmp,1,min)
      
      valid_flag = ifelse(tmp>=bp,1,0)
      valid_location = valid_location[valid_flag==1]
    }
  }
  
  cp_location_list = sort(cp_location_list)
  result$cp_location = cp_location_list
  
  # loglikelihood - penalty
  pdf_y = function(x,mu){
    pdf_t_t = 1/sqrt(2*pi)*exp(-(x-mu)^2/2)
    return(pdf_t_t)
  }#y_t ~ N(\theta_t,1)
  
  lambda_k = -(cp_num+1)*log(t)/(2*d)
  seg_num = length(cp_location_list)+1
  
  start_ind = 1
  end_ind = cp_location_list[1]-1
  
  for (seg in 1:seg_num){
    
    theta_seg = mean(yt[start_ind:end_ind])
    lambda_k = lambda_k + sum(pdf_y(yt[start_ind:end_ind],theta_seg))
    
    if (seg==seg_num){
      break
    }else{
      start_ind = end_ind+1
      if (seg<length(cp_location_list)){
        end_ind = cp_location_list[seg+1]-1
        #print(paste('seg:',seg,'end_ind:',end_ind,sep = ""))
      }else{
        end_ind = t
      }
    }
  }
  result$lambda_k = lambda_k
  return(result)
}

library(doParallel)
if (UT_flag==1){
  max_num = 10
  
  cl = makeCluster(detectCores())
  registerDoParallel(cl)
  res_list = foreach(kcp=1:max_num) %dopar% Bayesian_kcp(yt,p_best,cp_num = kcp)
  stopCluster(cl)
  rm(cl)
  res_vec = unlist(res_list)
  res_vec = res_vec[names(res_vec)=='lambda_k']
  #print(res_vec)
  best_k = which.max(res_vec)
  best_cp = Bayesian_kcp(yt,p_best,cp_num = best_k)
}

  
 




