################                WARNING!!! THIS SCRIPT TAKES A FEW WEEKS/DAYS TO COMPLETE AND NEEDS LOTS OF MEMEORY!!!!
###############               WARNING!!! THIS SCRIPT TAKES A FEW WEEKS/DAYS TO COMPLETE AND NEEDS LOTS OF MEMEORY!!!!
#############                WARNING!!! THIS SCRIPT TAKES A FEW WEEKS/DAYS TO COMPLETE AND NEEDS LOTS OF MEMEORY!!!!
###########                 WARNING!!! THIS SCRIPT TAKES A FEW WEEKS/DAYS TO COMPLETE AND NEEDS LOTS OF MEMEORY!!!!



## first, I create 15 random dataframes that each represent a cluster when printed next to each other
## for this I am using the clusterGeneration 
## this first step runs kluster_sim on 2-30 clusters data with 0.1 to 0.7 separation values 
# for 100 sampling iterations simulated 10 times. Cool stuff! 
simit <- function (n,smpl, smpl_iter, iter_simu){
sepvalue = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7) 
output = list()
for (j in 1:length(sepvalue)){
  sepValu <- sepvalue[j]
  dat <- genRandomClust(numClust = n,
                          clustszind=2,
                          clustSizeEq=100,
                          sepVal= sepValu,
                          numReplicate = 1,
                          rangeN=c(50,200))$datList[1]
    
    da <- data.frame(dat)
    rm(dat)
    colnames(da) <- c("x", "y")
    run = stringi::stri_rand_strings(1, 7)
    out <- kluster_sim(data = da, clusters = n, iter_sim = iter_simu, iter_klust = smpl_iter, smpl = smpl)$sim
    out$sepVal = sepValu
    out$run = run
    out$sample = smpl
    out$sampling_iteration = smpl_iter
    out$simulation_iteration = iter_simu    
    jpeg(paste0(getwd(),"/plots/",run,".jpg"))
    plot(da)
    dev.off()
  output[[j]] <- out
  rm(da,run,out)
}
result <- do.call(rbind, lapply(output, data.frame, stringsAsFactors=FALSE))
write.csv(result, paste0(getwd(),"/simulation_results/result",n,smpl, smpl_iter, iter_simu,".csv"))
}





for (t in 3:15){
  simit(n=t, smpl = 100,smpl_iter = 100,iter_simu = 10)
}


#####testing sensitivity to sample size

## first, I create 15 random dataframes that each represent a cluster when printed next to each other
## for this I am using the clusterGeneration 
## this first step runs kluster_sim on 2-30 clusters data with 0.1 to 0.7 separation values 
# for 100 sampling iterations simulated 10 times. Cool stuff! 
simit2 <- function (n, smpl_iter, iter_simu){
  sepvalue = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7) 
  output2 = list()
  for (j in 1:length(sepvalue)){
    sepValu <- sepvalue[j]
    dat <- genRandomClust(numClust=n, 
                          sepVal=sepValu, 
                          numNonNoisy=2, 
                          numNoisy=0, 
                          numOutlier=0, 
                          numReplicate=1, 
                          fileName="test",  
                          clustszind=1, # equal cluster size
                          clustSizeEq=n*10000, # total data points is nClust^10
                          rangeN=c(50,200), 
                          clustSizes=NULL, 
                          covMethod=c("eigen"), 
                          rangeVar=c(1, 10), 
                          lambdaLow=1, 
                          ratioLambda=10,  
                          alphad=1,
                          eta=1,
                          rotateind=TRUE, 
                          iniProjDirMethod=c("SL"), 
                          projDirMethod=c("newton"), 
                          alpha=0.05, 
                          ITMAX=20, 
                          eps=1.0e-10, 
                          quiet=TRUE, 
                          outputDatFlag=FALSE, 
                          outputLogFlag=FALSE, 
                          outputEmpirical=FALSE, 
                          outputInfo=FALSE)$datList[1]
    
    da <- data.frame(dat)
    rm(dat)
    colnames(da) <- c("x", "y")
    run = stringi::stri_rand_strings(1, 7)
    out100 <- kluster(data = da, clusters = n, iter_sim = iter_simu, iter_klust = smpl_iter, smpl = 100)$sim
    out100$sample = 100
    out500 <- kluster(data = da, clusters = n, iter_sim = iter_simu, iter_klust = smpl_iter, smpl = 500)$sim
    out500$sample = 500
    out1000 <- kluster(data = da, clusters = n, iter_sim = iter_simu, iter_klust = smpl_iter, smpl = 1000)$sim
    out1000$sample = 1000
    out2 = rbind(out100,out500,out1000)
    out2$sepVal = sepValu
    out2$run = run
    out2$sampling_iteration = smpl_iter
    out2$simulation_iteration = iter_simu    
    jpeg(paste0(getwd(),"/plots/",run,".jpg"))
    plot(da)
    dev.off()
    output2[[j]] <- out2
    rm(da,run,out2)
  }
  result2 <- do.call(rbind, lapply(output2, data.frame, stringsAsFactors=FALSE))
  write.csv(result2, paste0(getwd(),"/simulation_results2/result_kluster_",n, smpl_iter, iter_simu,".csv"))
}


for (t in 3:15){
  simit2(n=t,smpl_iter = 100,iter_simu = 1)
}


