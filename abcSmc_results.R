abcSmc_results <- readRDS(paste0(dir,"IIDIPUS_Results/abcsmc_2022-08-24_002627"))

dim(abcSmc_results$W)
dim(abcSmc_results$omegastore) # 14 parameters in Omega, 5 particles, 20 steps

abcSmc_results$omegastore[,,3]==abcSmc_results$omegastore[,,4]


abcSmc_results <- readRDS(paste0(dir,"IIDIPUS_Results/abcsmc_2022-08-23_140027")) # 14 parameters in Omega, 100 particles, 20 steps
# particles are just a collection of N random samples from a sequence of probability distributions, i.e. they are approximations of these
# distributions

abcSmc_results$omegastore[,,1]==abcSmc_results$omegastore[,,3]

abcSmc_results$omegastore[,,2]

plot(x=1:100, y=(abcSmc_results$omegastore[,,1][,1]))

W <- abcSmc_results$W
d <- abcSmc_results$d
s <- 1
Npart <- 5
Npart <- AlgoParams$smc_Npart
i_zeroweight <- which(W==0,  arr.ind=TRUE)
plot(rep(1, Npart), apply(d[,,1], c(1), min), xlim=c(1, s), ylim=c(min(d), max(d[which(d<Inf)])))
for (i in 2:s){
  i_s_zeroweight <- i_zeroweight[which(i_zeroweight[,2]==i),1]
  if (length(i_s_zeroweight)==0){
    points(rep(i, Npart), apply(d[,,i], c(1), min))
  }
  else{
    points(rep(i, Npart-length(i_s_zeroweight)), apply(d[-i_s_zeroweight,,i], c(1), min))
  }
}
