source("load.R") #seed is set in load

####scenario 1D
# generating one-dimensional data with different number of clusters

###################################################################
#######################  sc1d1 --> 4 clusters


df <- rbind(
  cbind(rnorm(121,1),rnorm(121,-5,3)),
  cbind(rnorm(263,1),rnorm(263,-9,1)),
  cbind(rnorm(363,1),rnorm(363,4,2)),
  cbind(rnorm(63,1),rnorm(63,35,3)),
  cbind(rnorm(21,1),rnorm(21,29,2)),
  cbind(rnorm(20,1),rnorm(20,80,3)),
  cbind(rnorm(3,1),rnorm(3,120,2))
)

colnames(df) <- c("x", "y")
df <- data.frame(df)
size <- nrow(df)
df$x <- 1
### see the plot
plot(df)

#enter the actual number of clusters
clusters <- 4

##starting to store results from different algorithms
tic()
  ##now let's compute optimal ks with BIC
BIC.best <- dim(Mclust(as.matrix(df), G=3:15)$z)[2]
tBIC <- toc()
tBIC <- tBIC$toc - tBIC$tic

tic()
pamk.best <- pamk(df)$nc
tpamk <- toc()
tpamk <- tpamk$toc - tpamk$tic

tic()
calinski.best <- as.numeric(which.max(cascadeKM(df, 1, 15, iter = 1000)$results[2,]))
tcalinski <- toc()
tcalinski <- tcalinski$toc - tcalinski$tic


tic()
apclus.best <- length(apcluster(negDistMat(r=2), df)@clusters)
tap <- toc()
tap <- tap$toc - tap$tic






kbicssum <- list()
kbics2 <- list()
t2 <- list()
tBIC_kluster <- list()
for (j in 1:3) {
      for (i in 1:5) {
        tic()
        dat2 <- as.matrix(sample(df, 50, replace = T))
        kbics2[[i]] <- dim(Mclust(dat2, G=3:15)$z)[2]
        t <- toc()
        t2[[i]] <- t$toc - t$tic
      }
    kbicssum[[j]] <- mean(as.numeric(kbics2))
  }
tBIC_kluster <- mean(unlist(t2))
mean_BIC_kluster <- round(mean(as.numeric(kbics2)),0)

table(kbicssum)



cols <- c("BIC.best","tBIC","pamk.best","tpamk","calinski.best","tcalinski","apclus.best","tap",
          "mean_BIC_kluster","tBIC_kluster")
results <- c(BIC.best,tBIC,pamk.best,tpamk,calinski.best,tcalinski,apclus.best,tap,
             mean_BIC_kluster,tBIC_kluster)
sim1 <- data.frame(cols,results)
sim1$n <- clusters



boxplot(kbicssum)
boxplot(round(as.numeric(kbicssum),0))
which.max(table(round(as.numeric(kbicssum),0)))



km1 <- hkmeans(df, k = 4,iter.max = 300) 
km2 <- hkmeans(df, k = 4,iter.max = 300) 


par(mfrow=c(1,2))
plot(df, col = km1$cluster, pch = 19, frame = FALSE,
     main = "K-means with k = 9") #fix the title
points(km9$centers, col = 1:6, pch = 8, cex = 3)


plot(df, col = km2$cluster, pch = 19, frame = FALSE,
     main = "K-means with k = 4") # fix the title
points(km4$centers, col = 1:6, pch = 8, cex = 3)

