
kluster_sim <- function(data,
                        clusters, #number of clusters we know
                        iter_sim, # number of simulation iterations
                        iter_klust, #number of iterations for clustering with sample size x
                        smpl #size of the sample to be taken with replacement out of data
) {
  
  size <- nrow(data)

##starting to store results from different algorithms
tic()
##now let's compute optimal ks with BIC
BIC.best <- dim(Mclust(as.matrix(data), G=3:15)$z)[2]
tBIC <- toc()
tBIC <- as.numeric(tBIC$toc - tBIC$tic)

tic()
pamk.best <- pamk(data)$nc
tpamk <- toc()
tpamk <- as.numeric(tpamk$toc - tpamk$tic)

tic()
calinski.best <- as.numeric(which.max(cascadeKM(data, 1, 15, iter = 1000)$results[2,]))
tcalinski <- toc()
tcalinski <- as.numeric(tcalinski$toc - tcalinski$tic)

tic()
apclus.best <- length(apcluster(negDistMat(r=2), data)@clusters)
tap <- toc()
tap <- as.numeric(tap$toc - tap$tic)


###kluster procedure
##bic
kbicssum <- list()
kbics2 <- list()
t2 <- list()
tBIC_kluster <- list()
for (j in 1:iter_sim) {
  for (i in 1:iter_klust) {
    tic()
    dat2 <- as.matrix(sample(data, smpl, replace = T))
    kbics2[[i]] <- dim(Mclust(dat2, G=3:15)$z)[2]
    rm(dat2)
    t <- toc()
    t2[[i]] <- t$toc - t$tic
  }
  kbicssum[[j]] <- mean(as.numeric(kbics2))
}
tBIC_kluster <- mean(unlist(t2))
mean_BIC_kluster <- round(mean(as.numeric(kbics2)),0)

table(kbicssum)


###pam method
kpamsum <- list()
kpam2 <- list()
t2 <- list()
tpam_kluster <- list()
for (j in 1:iter_sim) {
  for (i in 1:iter_klust) {
    tic()
    dat2 <- sample(data, smpl, replace = T)
    kpam2[[i]] <- pamk(dat2)$nc
    rm(dat2)
    t <- toc()
    t2[[i]] <- t$toc - t$tic
  }
  kpamsum[[j]] <- mean(as.numeric(kpam2))
}
tpam_kluster <- mean(unlist(t2))
mean_pam_kluster <- round(mean(as.numeric(kpam2)),0)

table(kpamsum)



kcalsum <- list()
kcal2 <- list()
t2 <- list()
tcal_kluster <- list()
for (j in 1:iter_sim) {
  for (i in 1:iter_klust) {
    tic()
    dat2 <- sample(data, smpl, replace = T)
    kcal2[[i]] <- as.numeric(which.max(cascadeKM(dat2, 1, 15, iter = 1000)$results[2,]))
    rm(dat2)
    t <- toc()
    t2[[i]] <- t$toc - t$tic
  }
  kcalsum[[j]] <- mean(as.numeric(kcal2))
}
tcal_kluster <- mean(unlist(t2))
mean_cal_kluster <- round(mean(as.numeric(kcal2)),0)

table(kcalsum)

###AP
kapsum <- list()
kap2 <- list()
t2 <- list()
tap_kluster <- list()
for (j in 1:iter_sim) {
  for (i in 1:iter_klust) {
    tic()
    dat2 <- sample(data, smpl, replace = T)
    kap2[[i]] <- length(apcluster(negDistMat(r=2), dat2)@clusters)
    rm(dat2)
    t <- toc()
    t2[[i]] <- t$toc - t$tic
  }
  kapsum[[j]] <- mean(as.numeric(kap2))
}
tap_kluster <- mean(unlist(t2))
mean_ap_kluster <- round(mean(as.numeric(kap2)),0)

table(kapsum)




algorithm <- c("BIC.best","pamk.best","calinski.best","apclus.best",
               "BIC_kluster",
               "cal_kluster",
               "pam_kluster",
               "ap_kluster")
ptime <- c(tBIC,tpamk,tcalinski,tap,
           tBIC_kluster,tcal_kluster,
           tpam_kluster,tap_kluster
)
k_num <- c(BIC.best,pamk.best,calinski.best,apclus.best,
           mean_BIC_kluster,mean_cal_kluster,
           mean_pam_kluster,
           mean_ap_kluster)


sim <- data.frame(algorithm,k_num,ptime)
sim$k_orig <- clusters
sim$e <- k_num - clusters
sim$n <- size


# 
# 
# 
# boxplot(kcalsum)
# boxplot(round(as.numeric(kapsum),0))
# which.max(table(round(as.numeric(kbicssum),0)))



####cluster and visualize the performance
km1 <- hkmeans(data, k = BIC.best,iter.max = 300) 
km2 <- hkmeans(data, k = pamk.best,iter.max = 300) 
km3 <- hkmeans(data, k = calinski.best,iter.max = 300) 
km4 <- hkmeans(data, k = apclus.best,iter.max = 300)
km5 <- hkmeans(data, k = mean_BIC_kluster,iter.max = 300) 
km6 <- hkmeans(data, k = mean_pam_kluster,iter.max = 300) 
km7 <- hkmeans(data, k = mean_cal_kluster,iter.max = 300) 
km8 <- hkmeans(data, k = mean_ap_kluster,iter.max = 300)



par(mfrow=c(2,4))
plot(data, col = km1$cluster, pch = 19, frame = FALSE,
     main = paste0("HK-means with k = ",BIC.best,"")) 
points(km1$centers, col = 1:BIC.best, pch = 8, cex = 3)

plot(data, col = km2$cluster, pch = 19, frame = FALSE,
     main = paste0("HK-means with k = ",pamk.best,"")) 
points(km2$centers, col = 1:pamk.best, pch = 8, cex = 3)

plot(data, col = km3$cluster, pch = 19, frame = FALSE,
     main = paste0("HK-means with k = ",calinski.best,"")) 
points(km3$centers, col = 1:calinski.best, pch = 8, cex = 3)

plot(data, col = km4$cluster, pch = 19, frame = FALSE,
     main = paste0("HK-means with k = ",apclus.best,"")) 
points(km4$centers, col = 1:apclus.best, pch = 8, cex = 3)

plot(data, col = km5$cluster, pch = 19, frame = FALSE,
     main = paste0("HK-means with kluster BIC process, k = ",mean_BIC_kluster,"")) 
points(km5$centers, col = 1:mean_BIC_kluster, pch = 8, cex = 3)

plot(data, col = km6$cluster, pch = 19, frame = FALSE,
     main = paste0("HK-means with kluster PAM process, k = ",mean_pam_kluster,"")) 
points(km6$centers, col = 1:mean_pam_kluster, pch = 8, cex = 3)

plot(data, col = km7$cluster, pch = 19, frame = FALSE,
     main = paste0("HK-means with kluster Calinski process, k = ",mean_cal_kluster,"")) 
points(km7$centers, col = 1:mean_cal_kluster, pch = 8, cex = 3)

plot(data, col = km8$cluster, pch = 19, frame = FALSE,
     main = paste0("HK-means with kluster AP process, k = ",mean_ap_kluster,"")) 
points(km8$centers, col = 1:mean_ap_kluster, pch = 8, cex = 3)


par(mfrow=c(2,2))
boxplot(round(as.numeric(kbicssum),0),
        main = paste0("kluster optimum cluster number from BIC w/ resampling.
                      Mean resampling estimate = ",mean_BIC_kluster,"
                      ",
                      " and ordinary BIC suggested ",BIC.best," clusters.")) 

boxplot(round(as.numeric(kpamsum),0),
        main = paste0("kluster optimum cluster number from PAM w/ resampling.
                      Mean resampling estimate = ",mean_pam_kluster,"
                      ",
                      " and ordinary PAM suggested ",pamk.best," clusters.")) 
boxplot(round(as.numeric(kcalsum),0),
        main = paste0("kluster optimum cluster number from Calinski w/ resampling.
                      Mean resampling estimate = ",mean_cal_kluster,"
                      and ordinary Calinski suggested ",calinski.best," clusters.")) 

boxplot(round(as.numeric(kapsum),0),
        main = paste0("kluster optimum cluster number from AP w/ resampling.
                      Mean resampling estimate = ",mean_ap_kluster,"
                      and ordinary AP suggested ",apclus.best," clusters.")) 


return(sim)

rm(size,clusters,tBIC,tpamk,tcalinski,tap,
   tBIC_kluster,tcal_kluster,
   tpam_kluster,tap_kluster,BIC.best,pamk.best,calinski.best,apclus.best,
   mean_BIC_kluster,mean_cal_kluster,
   mean_pam_kluster,
   mean_ap_kluster,df)
}


###############################################################
###############################################################
###############################################################



kluster_sim_sole <- function(data,
                        clusters, #number of clusters we know
                        iter_sim, # number of simulation iterations
                        iter_klust, #number of iterations for clustering with sample size x
                        smpl, #size of the sample to be taken with replacement out of data
                        algorithm #select analysis algorithm from BIC, PAMK, CAL, and AP
) {
  
  size <- nrow(data)
  if (algorithm == 'BIC') {
  ##starting to store results from different algorithms
  tic()
  ##now let's compute optimal ks with BIC
  BIC.best <- dim(Mclust(as.matrix(data), G=3:15)$z)[2]
  tBIC <- toc()
  tBIC <- as.numeric(tBIC$toc - tBIC$tic)
  
  ###kluster procedure
  ##bic
  kbicssum <- list()
  kbics2 <- list()
  t2 <- list()
  tBIC_kluster <- list()
  for (j in 1:iter_sim) {
    for (i in 1:iter_klust) {
      tic()
      dat2 <- as.matrix(sample(data, smpl, replace = T))
      kbics2[[i]] <- dim(Mclust(dat2, G=3:15)$z)[2]
      rm(dat2)
      t <- toc()
      t2[[i]] <- t$toc - t$tic
    }
    kbicssum[[j]] <- mean(as.numeric(kbics2))
  }
  tBIC_kluster <- mean(unlist(t2))
  mean_BIC_kluster <- round(mean(as.numeric(kbics2)),0)
  
  table(kbicssum)
  
  
  method <- c("BIC.best","BIC_kluster")
  ptime <- c(tBIC,tBIC_kluster)
  k_num <- c(BIC.best,mean_BIC_kluster)
  
  
  sim <- data.frame(method,k_num,ptime)
  sim$k_orig <- clusters
  sim$e <- k_num - clusters
  sim$n <- size
  
  
  # 
  # 
  # 
  # boxplot(kcalsum)
  # boxplot(round(as.numeric(kapsum),0))
  # which.max(table(round(as.numeric(kbicssum),0)))
  
  
  
  ####cluster and visualize the performance
  km1 <- hkmeans(data, k = BIC.best,iter.max = 300) 
  km5 <- hkmeans(data, k = mean_BIC_kluster,iter.max = 300) 
  
  
  
  par(mfrow=c(1,2))
  plot(data, col = km1$cluster, pch = 19, frame = FALSE,
       main = paste0("HK-means with k = ",BIC.best,"")) 
  
  plot(data, col = km5$cluster, pch = 19, frame = FALSE,
       main = paste0("HK-means with kluster BIC process, k = ",mean_BIC_kluster,"")) 
  
  
  
  par(mfrow=c(1,1))
  boxplot(round(as.numeric(kbicssum),0),
          main = paste0("kluster optimum cluster number from BIC w/ resampling.
                      Mean resampling estimate = ",mean_BIC_kluster,"
                      ",
                        " and ordinary BIC suggested ",BIC.best," clusters.")) 
  
  
  return(sim)
  
  rm(size,clusters,tBIC,tBIC_kluster,BIC.best,mean_BIC_kluster,df)
  
  
  
  
  } else
    if (algorithm == "PAMK") {
  
  tic()
  pamk.best <- pamk(data)$nc
  tpamk <- toc()
  tpamk <- as.numeric(tpamk$toc - tpamk$tic)
  
  ###pam method
  kpamsum <- list()
  kpam2 <- list()
  t2 <- list()
  tpam_kluster <- list()
  for (j in 1:iter_sim) {
    for (i in 1:iter_klust) {
      tic()
      dat2 <- sample(data, smpl, replace = T)
      kpam2[[i]] <- pamk(dat2)$nc
      rm(dat2)
      t <- toc()
      t2[[i]] <- t$toc - t$tic
    }
    kpamsum[[j]] <- mean(as.numeric(kpam2))
  }
  tpam_kluster <- mean(unlist(t2))
  mean_pam_kluster <- round(mean(as.numeric(kpam2)),0)
  
  table(kpamsum)
  
  
  method <- c("pamk.best","pam_kluster")
  ptime <- c(tpamk,tpam_kluster)
  k_num <- c(pamk.best,mean_pam_kluster)
  
  
  sim <- data.frame(method,k_num,ptime)
  sim$k_orig <- clusters
  sim$e <- k_num - clusters
  sim$n <- size
  
  
  
  ####cluster and visualize the performance
  
  km2 <- hkmeans(data, k = pamk.best,iter.max = 300) 
  km6 <- hkmeans(data, k = mean_pam_kluster,iter.max = 300) 
  
  
  
  par(mfrow=c(1,2))
  
  plot(data, col = km2$cluster, pch = 19, frame = FALSE,
       main = paste0("HK-means with k = ",pamk.best,"")) 
  points(km2$centers, col = 1:pamk.best, pch = 8, cex = 3)
  
  plot(data, col = km6$cluster, pch = 19, frame = FALSE,
       main = paste0("HK-means with kluster PAM process, k = ",mean_pam_kluster,"")) 
  points(km6$centers, col = 1:mean_pam_kluster, pch = 8, cex = 3)
  

  
  par(mfrow=c(1,1))
  
  boxplot(round(as.numeric(kpamsum),0),
          main = paste0("kluster optimum cluster number from PAM w/ resampling.
                      Mean resampling estimate = ",mean_pam_kluster,"
                      ",
                        " and ordinary PAM suggested ",pamk.best," clusters.")) 
  
  
  
  return(sim)
  
  rm(size,clusters,tpamk,mean_pam_kluster,df)
  
    } else 
      if (algorithm == "CAL") {
        tic()
        calinski.best <- as.numeric(which.max(cascadeKM(data, 1, 15, iter = 1000)$results[2,]))
        tcalinski <- toc()
        tcalinski <- as.numeric(tcalinski$toc - tcalinski$tic)    
        
        
        kcalsum <- list()
        kcal2 <- list()
        t2 <- list()
        tcal_kluster <- list()
        for (j in 1:iter_sim) {
          for (i in 1:iter_klust) {
            tic()
            dat2 <- sample(data, smpl, replace = T)
            kcal2[[i]] <- as.numeric(which.max(cascadeKM(dat2, 1, 15, iter = 1000)$results[2,]))
            rm(dat2)
            t <- toc()
            t2[[i]] <- t$toc - t$tic
          }
          kcalsum[[j]] <- mean(as.numeric(kcal2))
        }
        tcal_kluster <- mean(unlist(t2))
        mean_cal_kluster <- round(mean(as.numeric(kcal2)),0)
        
        table(kcalsum)
        
        
        method <- c("calinski.best","cal_kluster")
        ptime <- c(tcalinski,tcal_kluster)
        k_num <- c(calinski.best,mean_cal_kluster)
        
        
        sim <- data.frame(method,k_num,ptime)
        sim$k_orig <- clusters
        sim$e <- k_num - clusters
        sim$n <- size
  
        
        ####cluster and visualize the performance
        km3 <- hkmeans(data, k = calinski.best,iter.max = 300) 
        km7 <- hkmeans(data, k = mean_cal_kluster,iter.max = 300) 

        
        
        par(mfrow=c(1,2))
        
        plot(data, col = km3$cluster, pch = 19, frame = FALSE,
             main = paste0("HK-means with k = ",calinski.best,"")) 
        points(km3$centers, col = 1:calinski.best, pch = 8, cex = 3)
        
        
        plot(data, col = km7$cluster, pch = 19, frame = FALSE,
             main = paste0("HK-means with kluster Calinski process, k = ",mean_cal_kluster,"")) 
        points(km7$centers, col = 1:mean_cal_kluster, pch = 8, cex = 3)
        
        
        
        par(mfrow=c(1,1))

        boxplot(round(as.numeric(kcalsum),0),
                main = paste0("kluster optimum cluster number from Calinski w/ resampling.
                      Mean resampling estimate = ",mean_cal_kluster,"
                      and ordinary Calinski suggested ",calinski.best," clusters.")) 
        
        
        return(sim)
        
        rm(size,clusters,tcal_kluster,calinski.best,mean_cal_kluster,df)
        
        
      } else 
        
        if (algorithm == "AP") {
          
          tic()
          apclus.best <- length(apcluster(negDistMat(r=2), data)@clusters)
          tap <- toc()
          tap <- as.numeric(tap$toc - tap$tic)
          
          
          ###AP
          kapsum <- list()
          kap2 <- list()
          t2 <- list()
          tap_kluster <- list()
          for (j in 1:iter_sim) {
            for (i in 1:iter_klust) {
              tic()
              dat2 <- sample(data, smpl, replace = T)
              kap2[[i]] <- length(apcluster(negDistMat(r=2), dat2)@clusters)
              rm(dat2)
              t <- toc()
              t2[[i]] <- t$toc - t$tic
            }
            kapsum[[j]] <- mean(as.numeric(kap2))
          }
          tap_kluster <- mean(unlist(t2))
          mean_ap_kluster <- round(mean(as.numeric(kap2)),0)
          
          table(kapsum)       
          
          method <- c("apclus.best","ap_kluster")
          ptime <- c(tap,tap_kluster)
          k_num <- c(apclus.best,mean_ap_kluster)
          
          
          sim <- data.frame(method,k_num,ptime)
          sim$k_orig <- clusters
          sim$e <- k_num - clusters
          sim$n <- size
          
      
          ####cluster and visualize the performance
          
          km4 <- hkmeans(data, k = apclus.best,iter.max = 300)
          
          km8 <- hkmeans(data, k = mean_ap_kluster,iter.max = 300)
          
          
          
          par(mfrow=c(1,2))
          
          plot(data, col = km4$cluster, pch = 19, frame = FALSE,
               main = paste0("HK-means with k = ",apclus.best,"")) 
          points(km4$centers, col = 1:apclus.best, pch = 8, cex = 3)
          
          plot(data, col = km8$cluster, pch = 19, frame = FALSE,
               main = paste0("HK-means with kluster AP process, k = ",mean_ap_kluster,"")) 
          points(km8$centers, col = 1:mean_ap_kluster, pch = 8, cex = 3)
          
          
          par(mfrow=c(1,1))
          
          boxplot(round(as.numeric(kapsum),0),
                  main = paste0("kluster optimum cluster number from AP w/ resampling.
                      Mean resampling estimate = ",mean_ap_kluster,"
                      and ordinary AP suggested ",apclus.best," clusters.")) 
          
          
          return(sim)
          
          rm(size,clusters,tap,tap_kluster,apclus.best,mean_ap_kluster,df)
          
          
          
        } 
  
}



