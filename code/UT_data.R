# UT data ----
yt = c(rnorm(n = 25,mean = 1,sd=1),rnorm(n = 25,mean = -1,sd = 1))
j = 2
p = 2^j/length(yt)
mu = mean(yt)
rm(j)