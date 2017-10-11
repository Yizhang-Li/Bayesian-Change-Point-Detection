source('toolib.R')
library(readr)

# 1. loading data ----
Treasury_monthly = read_csv("treasury_rates.csv",col_types = cols(DATE = col_date(format = "%Y/%m/%d")))
Treasury_monthly = Treasury_monthly[Treasury_monthly$DATE <= as.Date('2010-06-01','%Y-%m-%d'),]
Treasury_monthly = Treasury_monthly[Treasury_monthly$DATE >= as.Date('1984-01-01','%Y-%m-%d'),]
Treasury_monthly$index = 1:dim(Treasury_monthly)[1]

# 2. 
for (col in 2:9){
  print(colnames(Treasury_monthly)[col])
  # 2.1 smooth
  Treasury_df = Treasury_monthly[,c(1,col,10)]
  colnames(Treasury_df) = c('DATE','rate','index')
  Treasury_df$avg_rate = 0
  for (i in 1:dim(Treasury_df)[1]){
    if (Treasury_df$index[i] %% 4 == 3){
      Treasury_df$avg_rate[i] = mean(Treasury_df$rate[(i-2):i])
    }
  }
  rm(i)
  Treasury_df = Treasury_df[Treasury_df$index %% 4 == 3,]
  Treasury_df$index = (1:dim(Treasury_df)[1] - 1) %% 4 + 1
  Treasury_df$index = paste(substr(Treasury_df$DATE,1,4),'Q',Treasury_df$index,sep = "")
  Treasury_df$rate = Treasury_df$avg_rate/sd(Treasury_df$avg_rate)
  
  # 2.2 cp_location
  yt = Treasury_df$rate
  mu = mean(yt)
  
  best_p_result = p_argmax(yt,mu,0,5)
  best_p = best_p_result$best_p
  best_j = best_p_result$best_j

  library(doParallel)
  max_num = 10  
  
  cl = makeCluster(detectCores())
  registerDoParallel(cl)
  res_list = foreach(kcp=1:max_num) %dopar% Bayesian_kcp(yt,best_p,cp_num = kcp,beta_c=1)
  stopCluster(cl)
  rm(cl)
  res_vec = unlist(res_list)
  res_vec = res_vec[names(res_vec)=='lambda_k']
  best_k = which.max(res_vec)
  best_cp = Bayesian_kcp(yt,best_p,cp_num = best_k,beta_c=2,plot_flag = FALSE)
  mu_list = best_cp$mu_list
  cp_location = best_cp$cp_location
  cp_date = Treasury_df$index[cp_location]
  #print(cp_date)
  
  png(paste(colnames(Treasury_monthly)[col],'_filter_y.png',sep = ""))
  plot(mu_list,main = colnames(Treasury_monthly)[col])
  dev.off()

  png(paste(colnames(Treasury_monthly)[col],'_y.png'))
  plot(Treasury_df$DATE,Treasury_df$rate,main = colnames(Treasury_monthly)[col])
  dev.off()
  
  write.csv(cp_date,paste('cp_date_',colnames(Treasury_monthly)[col],'.csv',sep = ""))

}

