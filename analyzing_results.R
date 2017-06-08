##### binding all result tables. results on 100 samples, 100 sampling iterations, 10 iterations all

path = paste0(getwd(),"/simulation_results")
names <- list.files(path)
names <- names[-1] #dropping the first name, which is plots folder
n <- length(names)

output <- list()
N <- length(names)
for (n in 1:N) {
  output[[n]] = data.frame(read.csv(paste0(path,"/",names[n],sep="")))
}
allresults <- do.call(rbind, lapply(output, data.frame, stringsAsFactors=FALSE))

##creating a percent error to compare overall performance
allresults$err <- round((allresults$e/allresults$k_orig),3)


## recalculating the actual time to iteratte sampling for 100 times and iterate the whole thing 10 times
# processing time is for a single run with a single iteration. so to reproduce the actual time for
# 100 iterations, we need to multiply the time by 100
# it was mistakenly divided by number of iterations, so multiply by 10
# it was iterated 10 times, so multiply by 10
allresults$ptime2 = ifelse(allresults$method %in% c("BIC.best","pamk.best","calinski.best","apclus.best"),
                           allresults$ptime, allresults$ptime*10000)

##removing data with more than 15 clusters -- because the algorithms were designated to consider 
# 3-15 clusters
allresults = subset(allresults, allresults$k_orig < 16)

### let's look at time of processing for original algorithms
origs = subset(allresults,allresults$method %in% c("BIC.best","pamk.best","calinski.best","apclus.best"))

p0 = ggplot(origs, aes(x=n, y=ptime))+
  # geom_smooth(aes(group=as.factor(method)),alpha = 0.7)+
  # geom_hline(yintercept=0.05, size = 1, colour = "#FF3721",linetype = "dashed",alpha = 0.5)+
  # geom_hline(yintercept=-0.05, size = 1, colour = "#FF3721",linetype = "dashed",alpha = 0.5)+
  # stat_smooth(aes(color=as.factor(method0)),method="lm",
  #             formula = y ~ x + poly(x, 2),fullrange=TRUE, alpha=0.2)+
  geom_line(aes(color=as.factor(method0)),alpha = 0.7) +
  geom_point(size = 0.6, alpha = 0.8, color = "gray") +
  # geom_vline(xintercept=0, size = 2000, colour = "#FF3721",linetype = "dashed",alpha = 0.5)+
  # ?geom_smooth
  # geom_boxplot(alpha = 0.7)+
  # geom_point(aes(x=n, y=scale(ptime2)), color = "blue", size = 1, alpha = 0.4 ) +
  scale_y_continuous(limits=c(0,50))+
  scale_colour_hue(h=c(0, 360),guide = F) + 
  theme_bw() +
  theme(panel.grid.major.y = element_line(colour = "gray"),
        panel.grid.minor.y = element_blank(),
        axis.line = element_line(size=1, colour = "black"),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma", face = "bold"),
        axis.text.x=element_text(colour="black", size = 10, face="plain"),
        axis.text.y=element_text(colour="black", size = 10, face="plain"),
        legend.position="bottom") +
  xlab("") +
  ylab("processing time (seconds)")


## now project the processing time for data up to 2 million rows
options(scipen=999)
p00 = ggplot(origs, aes(x=n, y=ptime))+
  stat_smooth(aes(color=as.factor(method0)),method="glm",
              formula = y ~ poly(x, 1),fullrange=TRUE, alpha=0.2)+
  # geom_point(aes(x=n, y=scale(ptime2)), color = "blue", size = 1, alpha = 0.4 ) +
  scale_x_continuous(limits=c(0,2000000))+
  scale_colour_hue(h=c(0, 360),guide = guide_legend(title = "Algorithm", nrow= 1)) + 
  theme_bw() +
  theme(panel.grid.major.y = element_line(colour = "gray"),
        panel.grid.minor.y = element_blank(),
        axis.line = element_line(size=1, colour = "black"),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma", face = "bold"),
        axis.text.x=element_text(colour="black", size = 10, face="plain"),
        axis.text.y=element_text(colour="black", size = 10, face="plain", angle = 90),
        legend.position="bottom") +
  xlab("Size of database (number of rows)") +
  ylab("processing time (seconds)")

p000 = grid.arrange(p0,p00)

ggsave(filename=paste0(getwd(),"/figures/proc_time_algorithm.PNG"), 
       plot=p000, dpi = 300, width = 12, height = 8)



## looking at good and bad results by sepVal and k_orig

aggregate(err ~ method+sepVal, data = allresults, length)
agg1 <- aggregate(err ~ method+sepVal, data = allresults, mean)
ggplot(agg1,aes(x=as.factor(sepVal),y=err, group=method,color=method)) + geom_path() + geom_point()
agg2 <- aggregate(err ~ method+k_orig, data = allresults, mean)

goods = subset(allresults, err <= 0.01 & err >= -0.01)
sepVal_good = aggregate(err ~ method+sepVal, data = goods, length)
k_orig_good = aggregate(err ~ method+k_orig, data = goods, length)
ggplot(k_orig_good, aes(x= as.numeric(reorder(k_orig,err)), y= err)) + geom_point(aes(color = method)) + scale_x_continuous(limits=c(0,30))

## looking at good performance by 
  ### number of clusters
# k_orig_good2 = as.data.frame.matrix(table(goods$method,goods$k_orig))
k_orig_good = data.frame(table(goods$method,goods$k_orig))
k_orig_good$ratio = round(as.numeric(k_orig_good$Freq/7),4)
p1 = ggplot(k_orig_good, aes(Var2, reorder(Var1, ratio))) + 
  geom_tile(aes(fill = ratio),colour = "gray") + 
  scale_fill_gradient(low = "white", high = "darkblue") + 
  labs(fill='Prediction ratio with 90%+ accuracy') + #+ geom_label(label = k_orig_good$ratio)
  xlab("Number of Clusters") +
  ylab("") +
  theme(legend.position="bottom",
        axis.text=element_text(face="bold", size = 12)) +
  ggtitle("")


ggsave(filename=paste0(getwd(),"/figures/goodsclusternumber.PNG"), 
       plot=p1, dpi = 300, width = 7.5, height = 7)


  ### separation value
# table(goods$method,goods$sepVal)
sepVal_good = data.frame(table(goods$method,goods$sepVal))
sepVal_good$ratio = round(as.numeric(sepVal_good$Freq/13),4)
p2 = ggplot(sepVal_good, aes(Var2, reorder(Var1, ratio))) + 
  geom_tile(aes(fill = ratio),colour = "gray") + 
  scale_fill_gradient(limits= c(0,1), low = "white", high = "darkblue") + 
  labs(fill='Prediction ratio with 90%+ accuracy') + #+ geom_label(label = k_orig_good$ratio)
  xlab("Separation Value") +
  ylab("") +
  theme(legend.position="bottom",
        axis.text=element_text(face="bold", size = 12)) +
  ggtitle("")

ggsave(filename=paste0(getwd(),"/figures/goodssepval.PNG"), 
       plot=p2, dpi = 300, width = 5.5, height = 7)

                                      # ## looking at bad performance by 
                                      # ### number of clusters
                                      # bads = subset(allresults, err >= 0.2 | err <= -0.2)
                                      # # table(bads$method,bads$k_orig)
                                      # k_orig_bads = data.frame(table(bads$method,bads$k_orig))
                                      # k_orig_bads$ratio = round(as.numeric(k_orig_bads$Freq/7),4)
                                      # ggplot(k_orig_bads, aes(Var2, reorder(Var1, ratio))) + geom_tile(aes(fill = ratio),colour = "white") + scale_fill_gradient(low = "white", high = "orange")
                                      # 
                                      # ### separation value
                                      # # table(bads$method,bads$sepVal)
                                      # sepVal_bads = data.frame(table(bads$method,bads$sepVal))
                                      # sepVal_bads$ratio = round(as.numeric(sepVal_bads$Freq/13),4)
                                      # ggplot(sepVal_bads, aes(Var2, reorder(Var1, ratio))) + geom_tile(aes(fill = ratio),colour = "white") + scale_fill_gradient(low = "white", high = "orange")
                                      # 
                                      # 
                                      # ## average performance by method + sepVal
                                      # ggplot(agg1, aes(x=reorder(method,err), y=err))+
                                      #   geom_point(aes(color=as.factor(agg1$sepVal)), size = 4) +
                                      #   geom_hline(yintercept=0.1)+
                                      #   geom_hline(yintercept=-0.1)+
                                      #   scale_y_continuous(limits=c(-1,1))
                                      # 
                                      # ## average performance by method+k_orig
                                      # ggplot(agg2, aes(x=reorder(method,err), y=err))+
                                      #   geom_point(aes(color=as.factor(agg2$k_orig)), size = 4) +
                                      #   geom_hline(yintercept=0.1)+
                                      #   geom_hline(yintercept=-0.1)+
                                      #   scale_y_continuous(limits=c(-1,1))
                                      # 
                                      # 
                                      # ## performance by method+k_orig+sepVal
                                      # ggplot(allresults, aes(x=reorder(method,err), y=err))+
                                      #   geom_point(aes(color=allresults$k_orig, size = as.factor(allresults$sepVal))) +
                                      #   geom_hline(yintercept=0.1)+
                                      #   geom_hline(yintercept=-0.1)+
                                      #   scale_y_continuous(limits=c(-2,2))
                                      # 
                                      #   
                                      # #performance by method
                                      # ggplot(allresults, aes(x=err))+
                                      # #   geom_point(aes(color=as.factor(sepVal)), size = 2, alpha= 0.5) +
                                      #   stat_density(aes(ymax = ..density..,  ymin = -..density..),
                                      #                geom = "ribbon", alpha = 0.6,
                                      #                position = "identity") +
                                      #   
                                      #   geom_vline(xintercept=0.1, size = 1, colour = "#FF3721",linetype = "dashed")+
                                      #   geom_vline(xintercept=-0.1, size = 1, colour = "#FF3721",linetype = "dashed")+
                                      #   # scale_y_continuous(limits=c(-6,6))+
                                      #   coord_cartesian(ylim=c(0, 3))+
                                      #   facet_grid(~method) + coord_flip()



##density plot by method
fill <- "#4271AE"
lines <- "#1F3552"
p3 = ggplot(allresults, aes(x =err)) +
  geom_density(colour = lines, fill = fill, 
               size = 1,alpha = 0.9) + 
  scale_x_continuous(name = "Ratio of error to actual cluster number",
                     breaks = seq(-1, 1, .5),
                     limits=c(-1, 1)) +
  # scale_y_continuous(name = "Density") +
  coord_cartesian(ylim=c(0, 8))+
  # ggtitle("Density plot of cluster number identification error") +
  geom_vline(xintercept = 0.1, size = 1, colour = "#FF3721",
             linetype = "dashed") +
  geom_vline(xintercept = -0.1, size = 1, colour = "#FF3721",
             linetype = "dashed") +
  theme_bw() +
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma", face = "bold"),
        axis.text.x=element_text(colour="black", size = 9),
        axis.text.y=element_text(colour="black", size = 9)) +
  facet_wrap(~ method , ncol = 3)

ggsave(filename=paste0(getwd(),"/figures/densityplot.PNG"), 
       plot=p3, dpi = 300, width = 14, height = 7)



                      # 
                      # # performance by method over number of clusters
                      # ggplot(allresults, aes(x=as.factor(k_orig), y=err))+
                      #   geom_point(aes(color=as.factor(sepVal)), size = 2, alpha = 0.4 ) +
                      #   geom_boxplot(alpha = 0.5)+
                      #   geom_hline(yintercept=0.1, size = 1, colour = "#FF3721",linetype = "dashed")+
                      #   geom_hline(yintercept=-0.1, size = 1, colour = "#FF3721",linetype = "dashed")+
                      #   scale_y_continuous(limits=c(-1,1))+
                      #   facet_wrap(~method, ncol = 3)


# performance by method over number of clusters
p4 = ggplot(allresults, aes(x=k_orig, y=err))+
  geom_hline(yintercept=0.1, size = 1, colour = "#FF3721",linetype = "dashed",alpha = 0.5)+
  geom_hline(yintercept=-0.1, size = 1, colour = "#FF3721",linetype = "dashed",alpha = 0.5)+
  geom_line(aes(color=as.factor(sepVal)),alpha = 0.7) +
  geom_point(aes(color=as.factor(sepVal)),shape = 21, fill = "white",size = 2, alpha = 0.7 ) +
  scale_y_continuous(limits=c(-2,2))+
  scale_x_continuous(breaks = seq(3,30,3))+
  scale_colour_hue(h=c(180, 270),guide = guide_legend(title = "Separation value", nrow= 1)) + 
  theme_bw() +
  theme(panel.grid.major.y = element_line(colour = "gray"),
        panel.grid.minor.y = element_blank(),
        axis.line = element_line(size=1, colour = "black"),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma", face = "bold"),
        axis.text.x=element_text(colour="black", size = 7, face="plain"),
        axis.text.y=element_text(colour="black", size = 7, face="plain"),
        legend.position="bottom") +
    xlab("Actual number of clusters") +
    ylab("Ratio of error to actual cluster number") +
  facet_wrap(~method, ncol = 3) 

ggsave(filename=paste0(getwd(),"/figures/resultsbymethod.PNG"), 
       plot=p4, dpi = 300, width = 7, height = 7)
  


### adding a new variable to be able to combine the four models
allresults$method0 = "BIC"
allresults$method0 = ifelse(allresults$method %in% c("pamk.best","pam_kluster_frq","pam_kluster_mean"), "PAM", allresults$method0)
allresults$method0 = ifelse(allresults$method %in% c("calinski.best","cal_kluster_mean","cal_kluster_frq"), "CAL", allresults$method0)
allresults$method0 = ifelse(allresults$method %in% c("apclus.best","ap_kluster_mean","ap_kluster_frq"), "AP", allresults$method0)

##adding a new variable to distinguish time for kluster method
allresults$algorithm = "original algorithm"
allresults$algorithm = ifelse(allresults$method %in% c("ap_kluster_mean","ap_kluster_frq","cal_kluster_mean","cal_kluster_frq",
                                                       "pam_kluster_frq","pam_kluster_mean", "BOC_kluster_frq", "BIC_kluster_mean"), 
                              "kluster proc",
                              allresults$algorithm)

                                        # # performance by method over data size
                                        # ggplot(allresults, aes(x=as.factor(n), y=err))+
                                        #   geom_point(aes(color=as.factor(sepVal)), size = 2, alpha = 0.4 ) +
                                        #   geom_boxplot(alpha = 0.7)+
                                        #   geom_hline(yintercept=0.05, size = 1, colour = "#FF3721",linetype = "dashed")+
                                        #   geom_hline(yintercept=-0.05, size = 1, colour = "#FF3721",linetype = "dashed")+
                                        #   scale_y_continuous(limits=c(-1,1))+
                                        #   facet_wrap(~method, ncol = 3)

### processing time over size of data
p5 = ggplot(allresults, aes(x=n, y=ptime2))+
  # geom_hline(yintercept=0.05, size = 1, colour = "#FF3721",linetype = "dashed",alpha = 0.5)+
  # geom_hline(yintercept=-0.05, size = 1, colour = "#FF3721",linetype = "dashed",alpha = 0.5)+
  geom_line(aes(color=as.factor(algorithm)),alpha = 0.7) +
  geom_point(size = 1, alpha = 0.4 ) +
  # geom_boxplot(alpha = 0.7)+
  # geom_point(aes(x=n, y=scale(ptime2)), color = "blue", size = 1, alpha = 0.4 ) +
  # scale_y_continuous(limits=c(-2,2))+
  scale_colour_hue(h=c(0, 360),guide = guide_legend(title = "Procedure", nrow= 1)) + 
  theme_bw() +
  theme(panel.grid.major.y = element_line(colour = "gray"),
        panel.grid.minor.y = element_blank(),
        axis.line = element_line(size=1, colour = "black"),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma", face = "bold"),
        axis.text.x=element_text(colour="black", size = 7, face="plain"),
        axis.text.y=element_text(colour="black", size = 7, face="plain"),
        legend.position="bottom") +
  xlab("Size of database (number of rows)") +
  ylab("processing time (seconds)") +
  facet_wrap(~method0, ncol = 2)

ggsave(filename=paste0(getwd(),"/figures/timeeoversize.PNG"), 
       plot=p5, dpi = 300, width = 7, height = 7)


                                            # # performance by method over sepVal
                                            # ggplot(allresults, aes(x=as.factor(sepVal), y=err))+
                                            #   geom_point(aes(color=as.factor(sepVal)), size = 2, alpha = 0.4 ) +
                                            #   geom_boxplot(alpha = 0.7)+
                                            #   geom_hline(yintercept=0.05, size = 1, colour = "#FF3721",linetype = "dashed")+
                                            #   geom_hline(yintercept=-0.05, size = 1, colour = "#FF3721",linetype = "dashed")+
                                            #   scale_y_continuous(limits=c(-1,1))+
                                            #   facet_wrap(~method, ncol = 3)
                                            # 
                                            # # speed by method over size
                                            # ggplot(allresults, aes(x=as.factor(n), y=ptime))+
                                            #   geom_point(aes(color=as.factor(sepVal)), size = 2, alpha = 0.4 ) +
                                            # #   geom_boxplot(alpha = 0.7)+
                                            # #   geom_hline(yintercept=0.05)+
                                            # #   geom_hline(yintercept=-0.05) +
                                            # #   scale_y_continuous(limits=c(-1,1))+
                                            #   facet_wrap(~method, ncol = 3)



                                            # ggplot(data = allresults) +
                                            #   geom_quasirandom(aes(method,as.numeric(err),color=sepVal), alpha=0.7) +
                                            #   geom_hline(yintercept=0.1, size = 1, colour = "#FF3721",linetype = "dashed")+
                                            #   geom_hline(yintercept=-0.1, size = 1, colour = "#FF3721",linetype = "dashed") +
                                            #   scale_y_continuous(limits=c(-1,1)) +
                                            #   facet_wrap(~method, ncol = 3, scales = "free_x") 
                                            #   
                                            #   
                                            #   xlab("Observation") +
                                            #   ylab("Observation Value") +
                                            #   ggtitle(paste('Unsupervised Learning w/ Bayesian Information Criterion
                                            #                   Implausible Values in "', observ[i], '"
                                            #                   Hopkins statistic = ',mean(observation[(observation$observation_source_value == observ[i]),"hopkins"]),'', sep=''))



# # examine results on clinical-like data
# allresults_clin <- subset(allresults, allresults$k_orig < 6 & allresults$sepVal < 0.3)
# ##density plot by method
# fill <- "#4271AE"
# lines <- "#1F3552"
# p6 = ggplot(allresults_clin, aes(x = err)) +
#   geom_density(aes(fill = as.factor(method)), colour = lines, 
#                size = 1,alpha = 0.7) +
#   scale_x_continuous(name = "Ratio of error to actual cluster number",
#                      breaks = seq(-2, 2, .5),
#                      limits=c(-2, 2)) +
#   scale_y_continuous(name = "Density") +
#   ggtitle("") +
#   geom_vline(xintercept = 0.1, size = 1, colour = "#FF3721",
#              linetype = "dashed") +
#   geom_vline(xintercept = -0.1, size = 1, colour = "#FF3721",
#              linetype = "dashed") +
#   theme_bw() +
#   theme(axis.line = element_line(size=1, colour = "black"),
#         panel.grid.major = element_line(colour = "#d3d3d3"),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(), panel.background = element_blank(),
#         plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
#         text=element_text(family="Tahoma", face = "bold"),
#         axis.text.x=element_text(colour="black", size = 9),
#         axis.text.y=element_text(colour="black", size = 9),
#         legend.position="bottom") +
#   scale_fill_discrete(name = "Method") +
#   # scale_colour_hue(h=c(0, 360),guide = guide_legend(title = "Method", nrow= 3)) +
#   facet_wrap(~ method0 , ncol = 2)
# 
# 
# ggsave(filename=paste0(getwd(),"/figures/densityforclinical.PNG"), 
#        plot=p6, dpi = 300, width = 7, height =7)



###########################                                                                ###########################
#########################################                                    #########################################
#########################################################    #########################################################
#                                                results on lab data
#########################################################    #########################################################
#########################################                                    #########################################
###########################                                                                ###########################


# results from the following simulations are from 100 smaples and 100 times iterations -- i have simulations with 10 and 1 iterations too.
LOINCS = fread(paste0(getwd(),"/cl_results/result_LOINC_klust.csv"))
### adding a new variable to be able to combine the four models
LOINCS$method0 = "BIC_kluster"
LOINCS$method0 = ifelse(LOINCS$method %in% c("pamk.best","pam_kluster_frq","pam_kluster_mean"), "pam_kluster", LOINCS$method0)
LOINCS$method0 = ifelse(LOINCS$method %in% c("calinski.best","cal_kluster_mean","cal_kluster_frq"), "cal_kluster", LOINCS$method0)
LOINCS$method0 = ifelse(LOINCS$method %in% c("apclus.best","ap_kluster_mean","ap_kluster_frq"), "ap_kluster", LOINCS$method0)



##density plot by method
fill <- "#4271AE"
lines <- "#1F3552"
p7 = ggplot(LOINCS, aes(x =k_num)) +
  geom_density(colour = lines, fill = fill, 
               size = 1,alpha = 0.9) + 
  scale_x_continuous(name = "Estimated number of clusters",
                     breaks = seq(2, 22, 10),
                     limits=c(1, 22)) +
  # scale_y_continuous(name = "Density") +
  coord_cartesian(ylim=c(0, 1))+
  # ggtitle("Density plot of cluster number identification error") +
  geom_vline(xintercept = 2, size = 1, colour = "#FF3721",
             linetype = "dashed") +
  geom_vline(xintercept = 5, size = 1, colour = "#FF3721",
             linetype = "dashed") +
  theme_bw() +
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma", face = "bold"),
        axis.text.x=element_text(colour="black", size = 9),
        axis.text.y=element_text(colour="black", size = 9)) +
  facet_wrap(~ method , ncol = 2)

ggsave(filename=paste0(getwd(),"/figures/densityplot_kluster_loinc.PNG"), 
       plot=p7, dpi = 300, width = 10, height = 7)



# speed over data size
## 100 iterations for simulation and kluster
p8 = ggplot(LOINCS, aes(x=n, y=ptime*100))+
  # geom_hline(yintercept=0.05, size = 1, colour = "#FF3721",linetype = "dashed",alpha = 0.5)+
  # geom_hline(yintercept=-0.05, size = 1, colour = "#FF3721",linetype = "dashed",alpha = 0.5)+
  geom_line(aes(color=as.factor(method0)),alpha = 0.7) +
  geom_point(aes(color=as.factor(method0)),size = 2, alpha = 0.4 ) +
  # geom_boxplot(alpha = 0.7)+
  # geom_point(aes(x=n, y=scale(ptime2)), color = "blue", size = 1, alpha = 0.4 ) +
  # scale_y_continuous(limits=c(-2,2))+
  scale_colour_hue(h=c(0, 360),guide = guide_legend(title = "Kluster Method", nrow= 1)) + 
  theme_bw() +
  theme(panel.grid.major.y = element_line(colour = "gray"),
        panel.grid.minor.y = element_blank(),
        axis.line = element_line(size=1, colour = "black"),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma", face = "bold"),
        axis.text.x=element_text(colour="black", size = 7, face="plain"),
        axis.text.y=element_text(colour="black", size = 7, face="plain"),
        legend.position="bottom") +
  xlab("Size of database (number of rows)") +
  ylab("processing time (seconds)") 

ggsave(filename=paste0(getwd(),"/figures/speed_kluster_loinc.PNG"), 
       plot=p8, dpi = 300, width = 7, height = 7)

##100 sampling iterations
LOINCS$ptime = LOINCS$ptime*100
aggregate(ptime ~ method0, data = LOINCS, mean)

p81 = ggplot(LOINCS) +
  geom_boxplot(aes(x= method0, y= ptime)) +
  theme_bw() +
  theme(panel.grid.major.y = element_line(colour = "gray"),
        panel.grid.minor.y = element_blank(),
        axis.line = element_line(size=1, colour = "black"),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma", face = "bold"),
        axis.text.x=element_text(colour="black", size = 7, face="plain"),
        axis.text.y=element_text(colour="black", size = 7, face="plain")) +
  xlab("kluster procedure") +
  ylab("processing time (seconds)") 

ggsave(filename=paste0(getwd(),"/figures/loincs_ptime.PNG"), 
       plot=p81, dpi = 300, width = 7, height = 7)




###########################                                                                ###########################
#########################################                                    #########################################
#########################################################    #########################################################
#                       checking sensitivity with sampling at different sizes on big data
#########################################################    #########################################################
#########################################                                    #########################################
###########################                                                                ###########################
### 

##### binding all result tables.
path2 = paste0(getwd(),"/simulation_results2")
names2 <- list.files(path2)
names2 <- names2[-1] #dropping the first name, which is plots folder
n2 <- length(names2)

output2 <- list()
N2 <- length(names2)
for (n in 1:N2) {
  output2[[n]] = data.frame(read.csv(paste0(path2,"/",names2[n],sep="")))
}
allresults2 <- do.call(rbind, lapply(output2, data.frame, stringsAsFactors=FALSE))


##creating a percent error to compare overall performance
allresults2$err <- round((allresults2$e/allresults2$k_orig),3)


## recalculating the actual time to iteratte sampling for 100 times and iterate the whole thing 1 time
# processing time is for a single run with a single iteration. so to reproduce the actual time for
# 100 iterations, we need to multiply the time by 100
# it was iterated 1 time, so multiply by 1
allresults2$ptime2 = ifelse(allresults2$method %in% c("BIC.best","pamk.best","calinski.best","apclus.best"),
                           allresults2$ptime, allresults2$ptime*100)

##removing data with more than 15 clusters -- because the algorithms were designated to consider 
# 3-15 clusters
allresults2 = subset(allresults2, allresults2$k_orig < 16)

allresults2$method0 = "BIC_kluster"
allresults2$method0 = ifelse(allresults2$method %in% c("pam_kluster_frq","pam_kluster_mean"), "pam_kluster", allresults2$method0)
allresults2$method0 = ifelse(allresults2$method %in% c("cal_kluster_mean","cal_kluster_frq"), "cal_kluster", allresults2$method0)
allresults2$method0 = ifelse(allresults2$method %in% c("ap_kluster_mean","ap_kluster_frq"), "ap_kluster", allresults2$method0)

time2= aggregate(ptime2 ~ method0+n, data = allresults2, mean)
ggplot(time2)+
  geom_line(aes(x=n,y=ptime2, color=method0))
