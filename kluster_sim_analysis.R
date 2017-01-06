source("load.R") #seed is set in load
source("kluster.R")


####phase I with samll data sets
##database 1 with 3 clusters
db3 <- rbind(
  cbind(rnorm(263,0,1),rnorm(263,0,5))
  ,cbind(rnorm(363,0,1),rnorm(363,60,5))
  ,cbind(rnorm(63,0,2),rnorm(63,120,5))
)

colnames(db3) <- c("x", "y")
db3 <- data.frame(db3)
# plot(db3)

##database 2 with 4 clusters
db4 <- rbind(
  cbind(rnorm(90,-50,10),rnorm(90,0,50))
  ,cbind(rnorm(80,120,10),rnorm(80,50,40))
  ,cbind(rnorm(70,90,10),rnorm(70,550,25))
  ,cbind(rnorm(50,-30,5),rnorm(50,600,30))
)

colnames(db4) <- c("x", "y")
db4 <- data.frame(db4)
# plot(db4)

##database 3 with 5 close clusters
db5 <- rbind(
  cbind(rnorm(90,-20,3),rnorm(90,0,50))
  ,cbind(rnorm(80,5,3),rnorm(80,90,40))
  ,cbind(rnorm(70,20,3),rnorm(70,550,25))
  ,cbind(rnorm(40,-10,3),rnorm(40,350,35))
  ,cbind(rnorm(50,-15,3),rnorm(50,600,30))
  ,cbind(rnorm(35,0,3),rnorm(35,800,30))
  
)

colnames(db5) <- c("x", "y")
db5 <- data.frame(db5)
# plot(db5)

##database 4 with 6 clusters
db6 <- rbind(
  cbind(rnorm(90,-50,10),rnorm(90,150,50))
  ,cbind(rnorm(80,60,10),rnorm(80,50,50))
  ,cbind(rnorm(70,70,10),rnorm(70,550,40))
  ,cbind(rnorm(50,-30,5),rnorm(50,600,30))
  ,cbind(rnorm(20,10,5),rnorm(20,-200,40))
  ,cbind(rnorm(80,10,10),rnorm(80,890,50))
  
)

colnames(db6) <- c("x", "y")
db6 <- data.frame(db6)
# plot(db6)


##database 5 with 7 clusters
db7 <- rbind(
  cbind(rnorm(90,0,5),rnorm(90,150,50))
  ,cbind(rnorm(80,60,5),rnorm(80,50,50))
  ,cbind(rnorm(70,40,5),rnorm(70,550,40))
  ,cbind(rnorm(50,-30,5),rnorm(50,600,30))
  ,cbind(rnorm(20,10,5),rnorm(20,890,40))
  ,cbind(rnorm(80,10,5),rnorm(80,-300,30))
  ,cbind(rnorm(50,-40,5),rnorm(50,0,50))
  
)

colnames(db7) <- c("x", "y")
db7 <- data.frame(db7)
# plot(db7)


##database 6 with 8 clusters
db8 <- rbind(
  cbind(rnorm(40,0,2),rnorm(40,150,50))
  ,cbind(rnorm(65,30,2),rnorm(65,50,60))
  ,cbind(rnorm(30,40,2),rnorm(30,550,50))
  ,cbind(rnorm(35,-15,3),rnorm(35,600,30))
  ,cbind(rnorm(200,10,3),rnorm(200,890,60))
  ,cbind(rnorm(25,5,2),rnorm(25,-300,60))
  ,cbind(rnorm(50,-20,3),rnorm(50,0,60))
  ,cbind(rnorm(15,35,3),rnorm(15,1300,20))
)

colnames(db8) <- c("x", "y")
db8 <- data.frame(db8)
# plot(db8)


## database 9 with 9 clusters
db9 <- rbind(
  cbind(rnorm(100,-20,0.5),rnorm(100,-50,30))
  ,cbind(rnorm(130,-10,1),rnorm(130,100,35))
  ,cbind(rnorm(200,20,2),rnorm(200,330,50))
  ,cbind(rnorm(80,1,1),rnorm(80,650,50))
  ,cbind(rnorm(90,5,1),rnorm(90,-250,30))
  ,cbind(rnorm(29,5,1),rnorm(29,170,50))
  ,cbind(rnorm(45,12,2),rnorm(45,1000,60))
  ,cbind(rnorm(9,-7,1),rnorm(9,1350,50))
  ,cbind(rnorm(35,0,1),rnorm(35,-750,30))
  
)

colnames(db9) <- c("x", "y")
db9 <- data.frame(db9)
# plot(db9)


sim3 <- kluster_sim(data = db3, clusters = 3, iter_sim = 20, iter_klust = 5, smpl = 100)
sim4 <- kluster_sim(data = db4, clusters = 4, iter_sim = 20, iter_klust = 5, smpl = 100)
sim5 <- kluster_sim(data = db5, clusters = 5, iter_sim = 20, iter_klust = 5, smpl = 100)
sim6 <- kluster_sim(data = db6, clusters = 6, iter_sim = 20, iter_klust = 5, smpl = 100)
sim7 <- kluster_sim(data = db7, clusters = 7, iter_sim = 20, iter_klust = 5, smpl = 100)
sim8 <- kluster_sim(data = db8, clusters = 8, iter_sim = 20, iter_klust = 5, smpl = 100)
sim9 <- kluster_sim(data = db9, clusters = 9, iter_sim = 20, iter_klust = 5, smpl = 100)

write.csv(sim3$sim, file = paste("kluster_sim_3_20_5_100_",nrow(db3),".csv", sep=""))
write.csv(sim4$sim, file = paste("kluster_sim_4_20_5_100_",nrow(db4),".csv", sep=""))
write.csv(sim5$sim, file = paste("kluster_sim_5_20_5_100_",nrow(db5),".csv", sep=""))
write.csv(sim6$sim, file = paste("kluster_sim_6_20_5_100_",nrow(db6),".csv", sep=""))
write.csv(sim7$sim, file = paste("kluster_sim_7_20_5_100_",nrow(db7),".csv", sep=""))
write.csv(sim8$sim, file = paste("kluster_sim_8_20_5_100_",nrow(db8),".csv", sep=""))
write.csv(sim9$sim, file = paste("kluster_sim_9_20_5_100_",nrow(db9),".csv", sep=""))


####phase II with 10x data sets
##database 1 with 3 clusters
db30 <- rbind(
  cbind(rnorm(2630,0,1),rnorm(2630,0,5))
  ,cbind(rnorm(3630,0,1),rnorm(3630,60,5))
  ,cbind(rnorm(630,0,2),rnorm(630,120,5))
)

colnames(db30) <- c("x", "y")
db30 <- data.frame(db30)
# plot(db3)

##database 2 with 4 clusters
db40 <- rbind(
  cbind(rnorm(900,-50,10),rnorm(900,0,50))
  ,cbind(rnorm(800,120,10),rnorm(800,50,40))
  ,cbind(rnorm(700,90,10),rnorm(700,550,25))
  ,cbind(rnorm(500,-30,5),rnorm(500,600,30))
)

colnames(db40) <- c("x", "y")
db40 <- data.frame(db40)
# plot(db4)

##database 3 with 5 close clusters
db50 <- rbind(
  cbind(rnorm(900,-20,3),rnorm(900,0,50))
  ,cbind(rnorm(800,5,3),rnorm(800,90,40))
  ,cbind(rnorm(700,20,3),rnorm(700,550,25))
  ,cbind(rnorm(400,-10,3),rnorm(400,350,35))
  ,cbind(rnorm(500,-15,3),rnorm(500,600,30))
  ,cbind(rnorm(350,0,3),rnorm(350,800,30))
  
)

colnames(db50) <- c("x", "y")
db50 <- data.frame(db50)
# plot(db5)

##database 4 with 6 clusters
db60 <- rbind(
  cbind(rnorm(900,-50,10),rnorm(900,150,50))
  ,cbind(rnorm(800,60,10),rnorm(800,50,50))
  ,cbind(rnorm(700,70,10),rnorm(700,550,40))
  ,cbind(rnorm(500,-30,5),rnorm(500,600,30))
  ,cbind(rnorm(200,10,5),rnorm(200,-200,40))
  ,cbind(rnorm(800,10,10),rnorm(800,890,50))
  
)

colnames(db60) <- c("x", "y")
db60 <- data.frame(db60)
# plot(db6)


##database 5 with 7 clusters
db70 <- rbind(
  cbind(rnorm(90,0,5),rnorm(90,150,50))
  ,cbind(rnorm(80,60,5),rnorm(80,50,50))
  ,cbind(rnorm(70,40,5),rnorm(70,550,40))
  ,cbind(rnorm(50,-30,5),rnorm(50,600,30))
  ,cbind(rnorm(20,10,5),rnorm(20,890,40))
  ,cbind(rnorm(80,10,5),rnorm(80,-300,30))
  ,cbind(rnorm(50,-40,5),rnorm(50,0,50))
  
)

colnames(db70) <- c("x", "y")
db70 <- data.frame(db70)
# plot(db7)


##database 6 with 8 clusters
db80 <- rbind(
  cbind(rnorm(400,0,2),rnorm(400,150,50))
  ,cbind(rnorm(650,30,2),rnorm(650,50,60))
  ,cbind(rnorm(300,40,2),rnorm(300,550,50))
  ,cbind(rnorm(350,-15,3),rnorm(350,600,30))
  ,cbind(rnorm(2000,10,3),rnorm(2000,890,60))
  ,cbind(rnorm(250,5,2),rnorm(250,-300,60))
  ,cbind(rnorm(500,-20,3),rnorm(500,0,60))
  ,cbind(rnorm(150,35,3),rnorm(150,1300,20))
)

colnames(db80) <- c("x", "y")
db80 <- data.frame(db80)
# plot(db8)


## database 9 with 9 clusters
db90 <- rbind(
  cbind(rnorm(1000,-20,0.5),rnorm(1000,-50,30))
  ,cbind(rnorm(1300,-10,1),rnorm(1300,100,35))
  ,cbind(rnorm(2000,20,2),rnorm(2000,330,50))
  ,cbind(rnorm(800,1,1),rnorm(800,650,50))
  ,cbind(rnorm(900,5,1),rnorm(900,-250,30))
  ,cbind(rnorm(290,5,1),rnorm(290,170,50))
  ,cbind(rnorm(450,12,2),rnorm(450,1000,60))
  ,cbind(rnorm(90,-7,1),rnorm(90,1350,50))
  ,cbind(rnorm(350,0,1),rnorm(350,-750,30))
  
)

colnames(db90) <- c("x", "y")
db90 <- data.frame(db90)
# plot(db9)


sim30 <- kluster_sim(data = db30, clusters = 3, iter_sim = 20, iter_klust = 5, smpl = 100)
sim40 <- kluster_sim(data = db40, clusters = 4, iter_sim = 20, iter_klust = 5, smpl = 100)
sim50 <- kluster_sim(data = db50, clusters = 5, iter_sim = 20, iter_klust = 5, smpl = 100)
sim60 <- kluster_sim(data = db60, clusters = 6, iter_sim = 20, iter_klust = 5, smpl = 100)
sim70 <- kluster_sim(data = db70, clusters = 7, iter_sim = 20, iter_klust = 5, smpl = 100)
sim80 <- kluster_sim(data = db80, clusters = 8, iter_sim = 20, iter_klust = 5, smpl = 100)
sim90 <- kluster_sim(data = db90, clusters = 9, iter_sim = 20, iter_klust = 5, smpl = 100)

write.csv(sim30$sim, file = paste("kluster_sim_3_20_5_100_",nrow(db30),".csv", sep=""))
write.csv(sim40$sim, file = paste("kluster_sim_4_20_5_100_",nrow(db40),".csv", sep=""))
write.csv(sim50$sim, file = paste("kluster_sim_5_20_5_100_",nrow(db50),".csv", sep=""))
write.csv(sim60$sim, file = paste("kluster_sim_6_20_5_100_",nrow(db60),".csv", sep=""))
write.csv(sim70$sim, file = paste("kluster_sim_7_20_5_100_",nrow(db70),".csv", sep=""))
write.csv(sim80$sim, file = paste("kluster_sim_8_20_5_100_",nrow(db80),".csv", sep=""))
write.csv(sim90$sim, file = paste("kluster_sim_9_20_5_100_",nrow(db90),".csv", sep=""))




####phase III with 100x data sets
##database 1 with 3 clusters
db300 <- rbind(
  cbind(rnorm(26300,0,1),rnorm(26300,0,5))
  ,cbind(rnorm(36300,0,1),rnorm(36300,60,5))
  ,cbind(rnorm(6300,0,2),rnorm(6300,120,5))
)

colnames(db300) <- c("x", "y")
db300 <- data.frame(db300)
# plot(db3)

##database 2 with 4 clusters
db400 <- rbind(
  cbind(rnorm(9000,-50,10),rnorm(9000,0,50))
  ,cbind(rnorm(8000,120,10),rnorm(8000,50,40))
  ,cbind(rnorm(7000,90,10),rnorm(7000,550,25))
  ,cbind(rnorm(5000,-30,5),rnorm(5000,600,30))
)

colnames(db400) <- c("x", "y")
db400 <- data.frame(db400)
# plot(db4)

##database 3 with 5 close clusters
db500 <- rbind(
  cbind(rnorm(9000,-20,3),rnorm(9000,0,50))
  ,cbind(rnorm(8000,5,3),rnorm(8000,90,40))
  ,cbind(rnorm(7000,20,3),rnorm(7000,550,25))
  ,cbind(rnorm(4000,-10,3),rnorm(4000,350,35))
  ,cbind(rnorm(5000,-15,3),rnorm(5000,600,30))
  ,cbind(rnorm(3500,0,3),rnorm(3500,800,30))
  
)

colnames(db500) <- c("x", "y")
db500 <- data.frame(db500)
# plot(db5)

##database 4 with 6 clusters
db600 <- rbind(
  cbind(rnorm(9000,-50,10),rnorm(9000,150,50))
  ,cbind(rnorm(8000,60,10),rnorm(8000,50,50))
  ,cbind(rnorm(7000,70,10),rnorm(7000,550,40))
  ,cbind(rnorm(5000,-30,5),rnorm(5000,600,30))
  ,cbind(rnorm(2000,10,5),rnorm(2000,-200,40))
  ,cbind(rnorm(8000,10,10),rnorm(8000,890,50))
  
)

colnames(db600) <- c("x", "y")
db600 <- data.frame(db600)
# plot(db6)


##database 5 with 7 clusters
db700 <- rbind(
  cbind(rnorm(9000,0,5),rnorm(9000,150,50))
  ,cbind(rnorm(8000,60,5),rnorm(8000,50,50))
  ,cbind(rnorm(7000,40,5),rnorm(7000,550,40))
  ,cbind(rnorm(5000,-30,5),rnorm(5000,600,30))
  ,cbind(rnorm(2000,10,5),rnorm(2000,890,40))
  ,cbind(rnorm(8000,10,5),rnorm(8000,-300,30))
  ,cbind(rnorm(5000,-40,5),rnorm(5000,0,50))
  
)

colnames(db700) <- c("x", "y")
db700 <- data.frame(db700)
# plot(db7)


##database 6 with 8 clusters
db800 <- rbind(
  cbind(rnorm(4000,0,2),rnorm(4000,150,50))
  ,cbind(rnorm(6500,30,2),rnorm(6500,50,60))
  ,cbind(rnorm(3000,40,2),rnorm(3000,550,50))
  ,cbind(rnorm(3500,-15,3),rnorm(3500,600,30))
  ,cbind(rnorm(20000,10,3),rnorm(20000,890,60))
  ,cbind(rnorm(2500,5,2),rnorm(2500,-300,60))
  ,cbind(rnorm(5000,-20,3),rnorm(5000,0,60))
  ,cbind(rnorm(1500,35,3),rnorm(1500,1300,20))
)

colnames(db800) <- c("x", "y")
db800 <- data.frame(db800)
# plot(db8)


## database 9 with 9 clusters
db900 <- rbind(
  cbind(rnorm(10000,-20,0.5),rnorm(10000,-50,30))
  ,cbind(rnorm(13000,-10,1),rnorm(13000,100,35))
  ,cbind(rnorm(20000,20,2),rnorm(20000,330,50))
  ,cbind(rnorm(8000,1,1),rnorm(8000,650,50))
  ,cbind(rnorm(9000,5,1),rnorm(9000,-250,30))
  ,cbind(rnorm(2900,5,1),rnorm(2900,170,50))
  ,cbind(rnorm(4500,12,2),rnorm(4500,1000,60))
  ,cbind(rnorm(900,-7,1),rnorm(900,1350,50))
  ,cbind(rnorm(3500,0,1),rnorm(3500,-750,30))
  
)

colnames(db900) <- c("x", "y")
db900 <- data.frame(db900)
# plot(db9)


sim300 <- kluster_sim(data = db300, clusters = 3, iter_sim = 20, iter_klust = 5, smpl = 100)
sim400 <- kluster_sim(data = db400, clusters = 4, iter_sim = 20, iter_klust = 5, smpl = 100)
sim500 <- kluster_sim(data = db500, clusters = 5, iter_sim = 20, iter_klust = 5, smpl = 100)
sim600 <- kluster_sim(data = db600, clusters = 6, iter_sim = 20, iter_klust = 5, smpl = 100)
sim700 <- kluster_sim(data = db700, clusters = 7, iter_sim = 20, iter_klust = 5, smpl = 100)
sim800 <- kluster_sim(data = db800, clusters = 8, iter_sim = 20, iter_klust = 5, smpl = 100)
sim900 <- kluster_sim(data = db900, clusters = 9, iter_sim = 20, iter_klust = 5, smpl = 100)

write.csv(sim300$sim, file = paste("kluster_sim_3_20_5_100_",nrow(db300),".csv", sep=""))
write.csv(sim400$sim, file = paste("kluster_sim_4_20_5_100_",nrow(db400),".csv", sep=""))
write.csv(sim500$sim, file = paste("kluster_sim_5_20_5_100_",nrow(db500),".csv", sep=""))
write.csv(sim600$sim, file = paste("kluster_sim_6_20_5_100_",nrow(db600),".csv", sep=""))
write.csv(sim700$sim, file = paste("kluster_sim_7_20_5_100_",nrow(db700),".csv", sep=""))
write.csv(sim800$sim, file = paste("kluster_sim_8_20_5_100_",nrow(db800),".csv", sep=""))
write.csv(sim900$sim, file = paste("kluster_sim_9_20_5_100_",nrow(db900),".csv", sep=""))





