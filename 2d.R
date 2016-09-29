source("load.R") #seed is set in load

####scenario 2D
# generating two-dimensional data with different number of clusters

###################################################################
#######################  sc2d1 --> 7 clusters
df <- rbind(
  cbind(rnorm(10601,1,2),rnorm(10601,-5,3)),
  cbind(rnorm(4063,20,1),rnorm(4063,-9,1)),
  cbind(rnorm(8063,36,2),rnorm(8063,4,2)),
  cbind(rnorm(6063,23,2),rnorm(6063,35,3)),
  cbind(rnorm(2061,3,1),rnorm(2061,29,2)),
  cbind(rnorm(200,14,2),rnorm(200,80,3)),
  cbind(rnorm(30,49,1),rnorm(30,120,2))
)

colnames(df) <- c("x", "y")
df <- data.frame(df)
size <- nrow(df)
### see the plot
plot(df)


time1 <- system.time(
  ##now let's compute optimal ks with BIC
  kBIC <- dim(Mclust(as.matrix(df), G=3:15)$z)[2]
)




kbicssum <- list()
kbics2 <- list()
time2 <- list()
time3 <- system.time(
  for (j in 1:20) {
    time2[[j]] <- system.time(
      for (i in 1:10) {
        dat2 <- as.matrix(sample(df, 500, replace = T))
        kbics2[[i]] <- dim(Mclust(dat2, G=3:15)$z)[2]
      }
    )
    kbicssum[[j]] <- mean(as.numeric(kbics2))
  }
)




boxplot(kbicssum)
boxplot(round(as.numeric(kbicssum),0))
table(round(as.numeric(kbicssum),0))
mean(as.numeric(kbics2))



km1 <- clara(df, KBIC, samples=200) 
km2 <- clara(df, "xxxxx", samples=200)


par(mfrow=c(1,2))
plot(df, col = km1$cluster, pch = 19, frame = FALSE,
     main = "K-means with k = 9") ##fix the title
points(km9$centers, col = 1:6, pch = 8, cex = 3)


plot(df, col = km2$cluster, pch = 19, frame = FALSE,
     main = "K-means with k = 4") ##fix the title
points(km4$centers, col = 1:6, pch = 8, cex = 3)




