# kluster on public data


# Breast Cancer Wisconsin (Diagnostic) Data Se
dat = read.csv(paste0(getwd(),"/Wisconsin/data.csv"))



plot = ggplot(dat)+
  geom_point(aes(x=area_mean, y= texture_mean,colour=diagnosis))+
  ggtitle("") +
  theme_bw() +
  theme(panel.grid.major.y = element_line(colour = "gray"),
        panel.grid.minor.y = element_blank(),
        axis.line = element_line(size=1, colour = "black"),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma", face = "bold"),
        axis.text.x=element_text(colour="black", size = 10, face="plain"),
        axis.text.y=element_text(colour="black", size = 10, face="plain"),
        legend.position="bottom")

ggsave(filename="publicdata/Wisconsin.PNG", 
       plot=plot, dpi = 300, width = 12, height = 8)

k = klust(data = dat[,c("area_mean","texture_mean")],iter_klust = 100, smpl=100, algorithm = "BIC")$f_BIC_k
BreastCancerWisconsinDB = kluster_sim(data = dat[,c("area_mean","texture_mean")], clusters = 2, iter_sim = 1, iter_klust = 100, smpl = 100)$sim
BreastCancerWisconsinDB$data = "Breast Cancer Wisconsin (Diagnostic) Data Set"

write.csv(BreastCancerWisconsinDB, paste0(getwd(),"/Wisconsin/kluster_BreastCancerWisconsinDB.csv"))


