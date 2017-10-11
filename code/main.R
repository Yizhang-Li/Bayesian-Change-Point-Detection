source('toolib.R')
library(readr)

# 1. loading data ----
Tbill_3mth_monthly = read_csv('Tbill_3mth_monthly.csv')
colnames(Tbill_3mth_monthly) = c('date','rate')
Tbill_3mth_monthly = Tbill_3mth_monthly[Tbill_3mth_monthly$date <= as.Date('2010-06-01','%Y-%m-%d'),]
Tbill_3mth_monthly = Tbill_3mth_monthly[Tbill_3mth_monthly$date >= as.Date('1984-01-01','%Y-%m-%d'),]
Tbill_3mth_monthly$index = 1:dim(Tbill_3mth_monthly)[1]
Tbill_3mth_monthly$avg_rate = 0
for (i in 1:dim(Tbill_3mth_monthly)[1]){
  if (Tbill_3mth_monthly$index[i] %% 4 == 3){
    Tbill_3mth_monthly$avg_rate[i] = mean(Tbill_3mth_monthly$rate[(i-2):i])
  }
}
rm(i)
Tbill_3mth_monthly = Tbill_3mth_monthly[Tbill_3mth_monthly$index %% 4 == 3,]
Tbill_3mth_monthly$index = (1:dim(Tbill_3mth_monthly)[1] - 1) %% 4 + 1
Tbill_3mth_monthly$index = paste(substr(Tbill_3mth_monthly$date,1,4),'Q',Tbill_3mth_monthly$index,sep = "")
Tbill_3mth_monthly$rate = Tbill_3mth_monthly$avg_rate/sd(Tbill_3mth_monthly$avg_rate)

# 2. cp_location
yt = Tbill_3mth_monthly$rate
mu = mean(yt)

best_p = p_argmax(yt,mu,-5,5)$best_p
#bootsrap_result = bootsrap_test(yt,20,mu)
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
best_cp = Bayesian_kcp(yt,best_p,cp_num = best_k,beta_c=2,plot_flag = TRUE)
cp_location = best_cp$cp_location
cp_date = Tbill_3mth_monthly$index[cp_location]
print(cp_date)

# 3. plot
plot(Tbill_3mth_monthly$date,Tbill_3mth_monthly$rate)

write.csv(cp_date,'cp_date_3mth.csv')
