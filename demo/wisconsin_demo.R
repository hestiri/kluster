if (!require("easypackages")) install.packages('easypackages', repos = "http://cran.rstudio.com/")
packages("data.table","devtools","dplyr","DT","ggplot2","gridExtra","factoextra",
         "rmarkdown","reshape2","visNetwork","rmdformats","ggrepel",
         prompt = F)

devtools::install_github("hestiri/kluster")



dat = read.csv("data/Breast_Cancer_Wisconsin.csv")
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
plot


BreastCancerWisconsinDB = kluster_sim(data = dat[,c("area_mean","texture_mean")], clusters = 2, iter_sim = 1, iter_klust = 100, smpl = 100)$sim
BreastCancerWisconsinDB$data = "Breast Cancer Wisconsin (Diagnostic) Data Set"


k = kluster(data = dat[,c("area_mean","texture_mean")],iter_klust = 100, smpl=100, algorithm = "BIC")$f_BIC_k
