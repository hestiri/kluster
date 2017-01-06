source("load.R") #seed is set in load
source("kluster.R")


##database 2 with 4 clusters 1000X
db4 <- rbind(
  cbind(rnorm(90000,-50,10),rnorm(90000,0,50))
  ,cbind(rnorm(80000,120,10),rnorm(80000,50,40))
  ,cbind(rnorm(70000,90,10),rnorm(70000,550,25))
  ,cbind(rnorm(50000,-30,5),rnorm(50000,600,30))
)

colnames(db4) <- c("x", "y")
db4 <- data.frame(db4)
# plot(db4)


##database 3 with 5 close clusters 100X
db5 <- rbind(
  cbind(rnorm(9000,-20,3),rnorm(9000,0,50))
  ,cbind(rnorm(8000,5,3),rnorm(8000,90,40))
  ,cbind(rnorm(7000,20,3),rnorm(7000,550,25))
  ,cbind(rnorm(4000,-10,3),rnorm(4000,350,35))
  ,cbind(rnorm(5000,-15,3),rnorm(5000,600,30))
  ,cbind(rnorm(3500,0,3),rnorm(3500,800,30))
  
)

colnames(db5) <- c("x", "y")
db5 <- data.frame(db5)
# plot(db5)


##database 6 with 8 clusters 10X
db8 <- rbind(
  cbind(rnorm(400,0,2),rnorm(400,150,50))
  ,cbind(rnorm(650,30,2),rnorm(650,50,60))
  ,cbind(rnorm(300,40,2),rnorm(300,550,50))
  ,cbind(rnorm(350,-15,3),rnorm(350,600,30))
  ,cbind(rnorm(2000,10,3),rnorm(2000,890,60))
  ,cbind(rnorm(250,5,2),rnorm(250,-300,60))
  ,cbind(rnorm(500,-20,3),rnorm(500,0,60))
  ,cbind(rnorm(150,35,3),rnorm(150,1300,20))
)

colnames(db8) <- c("x", "y")
db8 <- data.frame(db8)
# plot(db8)





sim4 <- kluster(data = db4, clusters = 4, iter_sim = 20, iter_klust = 5, smpl = 100)
sim5 <- kluster(data = db5, clusters = 5, iter_sim = 20, iter_klust = 5, smpl = 100)
sim8 <- kluster(data = db8, clusters = 8, iter_sim = 20, iter_klust = 5, smpl = 100)

###increase iteration from 5 to 1000
sim4_iter1000 <- kluster(data = db4, clusters = 4, iter_sim = 20, iter_klust = 1000, smpl = 100)
sim5_iter1000 <- kluster(data = db5, clusters = 5, iter_sim = 20, iter_klust = 1000, smpl = 100)
sim8_iter1000 <- kluster(data = db8, clusters = 8, iter_sim = 20, iter_klust = 1000, smpl = 100)

###increase iteration from 5 to 100
sim4_iter100 <- kluster(data = db4, clusters = 4, iter_sim = 20, iter_klust = 100, smpl = 100)
sim5_iter100 <- kluster(data = db5, clusters = 5, iter_sim = 20, iter_klust = 100, smpl = 100)
sim8_iter100 <- kluster(data = db8, clusters = 8, iter_sim = 20, iter_klust = 100, smpl = 100)

###increase iteration from 20 to 100 and 100
sim4_iter100100 <- kluster(data = db4, clusters = 4, iter_sim = 100, iter_klust = 100, smpl = 100)
sim5_iter100100 <- kluster(data = db5, clusters = 5, iter_sim = 100, iter_klust = 100, smpl = 100)
sim8_iter100100 <- kluster(data = db8, clusters = 8, iter_sim = 100, iter_klust = 100, smpl = 100)

###increase sample from 100 to 300 
sim4_n300 <- kluster(data = db4, clusters = 4, iter_sim = 100, iter_klust = 100, smpl = 300)
sim5_n300 <- kluster(data = db5, clusters = 5, iter_sim = 100, iter_klust = 100, smpl = 300)
sim8_n300 <- kluster(data = db8, clusters = 8, iter_sim = 100, iter_klust = 100, smpl = 300)


write.csv(sim4$sim, file = "results/kluster_4_20_5_100.csv")
write.csv(sim5$sim, file = "results/kluster_5_20_5_100.csv")
write.csv(sim8$sim, file = "results/kluster_8_20_5_100.csv")
write.csv(sim4_iter1000$sim, file = "results/kluster_4_20_1000_100.csv")
write.csv(sim5_iter1000$sim, file = "results/kluster_5_20_1000_100.csv")
write.csv(sim8_iter1000$sim, file = "results/kluster_8_20_1000_100.csv")
write.csv(sim4_iter100$sim, file = "results/kluster_4_20_100_100.csv")
write.csv(sim5_iter100$sim, file = "results/kluster_5_20_100_100.csv")
write.csv(sim8_iter100$sim, file = "results/kluster_8_20_100_100.csv")
write.csv(sim4_iter100100$sim, file = "results/kluster_4_100_100_100.csv")
write.csv(sim5_iter100100$sim, file = "results/kluster_5_100_100_100.csv")
write.csv(sim8_iter100100$sim, file = "results/kluster_8_100_100_100.csv")
write.csv(sim4_n300$sim, file = "results/kluster_4_100_100_300.csv")
write.csv(sim5_n300$sim, file = "results/kluster_5_100_100_300.csv")
write.csv(sim8_n300$sim, file = "results/kluster_8_100_100_300.csv")