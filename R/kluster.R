###### if you want to do kluster only for application

kluster <- function(data,
                  iter_klust = 100, #number of iterations for clustering with sample_n size x
                  smpl = 100, #size of the sample_n to be taken with replacement out of data
                  algorithm = "BIC" #select analysis algorithm from BIC, PAMK, CAL, and AP
) {

    smpl = ifelse(nrow(data) >= 3000, 500, smpl)

    if (algorithm == 'BIC') {
        ###kluster procedure
        ##bic
        kbics <- list()
        for (i in 1:iter_klust) {
            dat2 <- as.matrix(data[sample(nrow(data), smpl, replace=TRUE), ])
            kbics[[i]] <- dim(Mclust(dat2, G=1:15)$z)[2]
            rm(dat2)
        }
        m_BIC_k <- round(mean(as.numeric(kbics)),0)
        f_BIC_k <- as.numeric(names(which.max(table(unlist(kbics)))))

        method <- c("BIC_kluster")
        k_mean <- c(m_BIC_k)
        k_freq <- c(f_BIC_k)

        sim <- data.frame(method,k_mean,k_freq)

        return(
            list("sim"=sim,
                 "m_BIC_k"=m_bic_k,
                 "f_BIC_k"=f_bic_k)
        )


    } else
        if (algorithm == "PAMK") {

            ###pam method
            kpam <- list()
            for (i in 1:iter_klust) {
                dat2 <- sample_n(data, smpl, replace = T)
                kpam[[i]] <- pamk(dat2)$nc
                rm(dat2)
            }
            m_pam_k <- round(mean(as.numeric(kpam)),0)
            f_pam_k <- as.numeric(names(which.max(table(unlist(kpam)))))

            method <- c("PAM_kluster")
            k_mean <- c(m_pam_k)
            k_freq <- c(f_pam_k)
            sim <- data.frame(method,k_mean,k_freq)


            return(
                list("sim"=sim,
                     "m_pam_k"=m_pam_k,
                     "f_pam_k"=f_pam_k)
            )


        } else
            if (algorithm == "CAL") {

                kcal <- list()
                for (i in 1:iter_klust) {

                    dat2 <- sample_n(data, smpl, replace = T)
                    kcal[[i]] <- as.numeric(which.max(cascadeKM(dat2, 1, 15, iter = 1000)$results[2,]))
                    rm(dat2)

                }
                m_cal_kluster <- round(mean(as.numeric(kcal)),0)
                f_cal_k <- as.numeric(names(which.max(table(unlist(kcal)))))

                method <- c("CAL_kluster")
                k_mean <- c(m_cal_k)
                k_freq <- c(f_cal_k)

                sim <- data.frame(method,k_mean,k_freq)


                return(
                    list("sim"=sim,
                         "m_cal_k"=m_cal_k,
                         "f_cal_k"=f_cal_k)
                )


            } else

                if (algorithm == "AP") {

                    kap <- list()
                    for (i in 1:iter_klust) {
                        dat2 <- sample_n(data, smpl, replace = T)
                        kap2[[i]] <- length(apcluster(negDistMat(r=2), dat2)@clusters)
                        rm(dat2)
                    }
                    m_ap_k <- round(mean(as.numeric(kap)),0)
                    f_ap_k <- as.numeric(names(which.max(table(unlist(kap)))))


                    method <- c("AP_kluster")
                    k_mean <- c(m_ap_k)
                    k_freq <- c(f_ap_k)

                    sim <- data.frame(method,k_mean,k_freq)

                    return(
                        list("sim"=sim,
                             "m_ap_k"=m_ap_k,
                             "f_ap_k"=f_ap_k)
                    )


                }
    # else
    #   if (algorithm == "Default") {
    #
    #     kbics <- list()
    #     kpam <- list()
    #     kap <- list()
    #       for (i in 1:iter_klust) {
    #         dat2 <- sample_n(data, smpl, replace = T)
    #         kbics[[i]] <- dim(Mclust(dat2, G=1:15)$z)[2]
    #
    #
    #
    #         tic()
    #         kpam2[[i]] <- pamk(dat2)$nc
    #         t2i <- toc()
    #         t2[[i]] <- t2i$toc - t2i$tic
    #
    #         # tic()
    #         # kcal2[[i]] <- as.numeric(which.max(cascadeKM(dat2, 1, 15, iter = 1000)$results[2,]))
    #         # t3i <- toc()
    #         # t3[[i]] <- t3i$toc - t3i$tic
    #
    #         tic()
    #         kap2[[i]] <- length(apcluster(negDistMat(r=2), dat2)@clusters)
    #         t4i <- toc()
    #         t4[[i]] <- t4i$toc - t4i$tic
    #
    #         rm(dat2)
    #
    #       }
    #       kbicssum[[j]] <- mean(as.numeric(kbics2))
    #       kpamsum[[j]] <- mean(as.numeric(kpam2))
    #       # kcalsum[[j]] <- mean(as.numeric(kcal2))
    #       kapsum[[j]] <- mean(as.numeric(kap2))
    #
    #     }
    #     tBIC_kluster <- mean(unlist(t1))
    #     m_BIC_k <- round(mean(as.numeric(kbics2)),0)
    #     tpam_kluster <- mean(unlist(t2))
    #     m_pam_k <- round(mean(as.numeric(kpam2)),0)
    #     # tcal_kluster <- mean(unlist(t3))
    #     # m_cal_kluster <- round(mean(as.numeric(kcal2)),0)
    #     tap_kluster <- mean(unlist(t4))
    #     m_ap_k <- round(mean(as.numeric(kap2)),0)
    #
    #     f_BIC_k <- as.numeric(names(which.max(table(unlist(kbics2)))))
    #     f_pam_k <- as.numeric(names(which.max(table(unlist(kpam2)))))
    #     f_ap_k <- as.numeric(names(which.max(table(unlist(kap2)))))
    #
    #     method <- c("BIC_kluster","pam_kluster","ap_kluster")
    #     ptime <- c(tBIC_kluster/iter_sim,tpam_kluster/iter_sim,tap_kluster/iter_sim)
    #     k_mean <- c(m_BIC_k,m_pam_k,m_ap_k)
    #     k_freq <- c(f_BIC_k,f_pam_k,f_ap_k)
    #
    #
    #
    #     sim <- data.frame(method,k_mean,k_freq,ptime)
    #     sim$k_orig <- clusters
    #     sim$e_mean <- k_mean - clusters
    #     sim$e_freq <- k_freq - clusters
    #
    #     sim$n <- size
    #
    #     #
    #     #
    #     #
    #     # ####cluster and visualize the performance
    #     # km5 <- hkmeans(data, k = mean_BIC_kluster,iter.max = 300)
    #     # km6 <- hkmeans(data, k = mean_pam_kluster,iter.max = 300)
    #     # km7 <- hkmeans(data, k = mean_cal_kluster,iter.max = 300)
    #     # km8 <- hkmeans(data, k = mean_ap_kluster,iter.max = 300)
    #     #
    #     #
    #     #
    #     #
    #     #
    #     # par(mfrow=c(2,2))
    #     # # plot(data, col = km1$cluster, pch = 19, frame = FALSE,
    #     # #      main = paste0("HK-means with k = ",BIC.best,""))
    #     # #
    #     # plot(data, col = km5$cluster, pch = 19, frame = FALSE,
    #     #      main = paste0("HK-means with kluster BIC process, k = ",mean_BIC_kluster,""))
    #     #
    #     # plot(data, col = km6$cluster, pch = 19, frame = FALSE,
    #     #      main = paste0("HK-means with kluster PAM process, k = ",mean_pam_kluster,""))
    #     # points(km6$centers, col = 1:mean_pam_kluster, pch = 8, cex = 3)
    #     #
    #     # plot(data, col = km7$cluster, pch = 19, frame = FALSE,
    #     #      main = paste0("HK-means with kluster Calinski process, k = ",mean_cal_kluster,""))
    #     # points(km7$centers, col = 1:mean_cal_kluster, pch = 8, cex = 3)
    #     #
    #     # plot(data, col = km8$cluster, pch = 19, frame = FALSE,
    #     #      main = paste0("HK-means with kluster AP process, k = ",mean_ap_kluster,""))
    #     # points(km8$centers, col = 1:mean_ap_kluster, pch = 8, cex = 3)
    #     #
    #     #
    #     #
    #     #
    #     # boxplot(round(as.numeric(kbicssum),0),
    #     #         main = paste0("kluster optimum cluster number from BIC w/ resampling.
    #     #                       Mean resampling estimate = ",mean_BIC_kluster,"
    #     #                       ",
    #     #                       " and ordinary BIC suggested ",BIC.best," clusters."))
    #     # boxplot(round(as.numeric(kpamsum),0),
    #     #         main = paste0("kluster optimum cluster number from PAM w/ resampling.
    #     #                       Mean resampling estimate = ",mean_pam_kluster,"
    #     #                       ",
    #     #                       " and ordinary PAM suggested ",pamk.best," clusters."))
    #     # boxplot(round(as.numeric(kcalsum),0),
    #     #         main = paste0("kluster optimum cluster number from Calinski w/ resampling.
    #     #                       Mean resampling estimate = ",mean_cal_kluster,"
    #     #                       and ordinary Calinski suggested ",calinski.best," clusters."))
    #     # boxplot(round(as.numeric(kapsum),0),
    #     #         main = paste0("kluster optimum cluster number from AP w/ resampling.
    #     #                       Mean resampling estimate = ",mean_ap_kluster,"
    #     #                       and ordinary AP suggested ",apclus.best," clusters."))
    #
    #
    #     return(
    #       list("sim"=sim,
    #            "m_BIC_k"=m_BIC_k,
    #            "m_pam_k"=m_pam_k,
    #            "m_ap_k"=m_ap_k,
    #            "f_BIC_k"=f_BIC_k,
    #            "f_pam_k"=f_pam_k,
    #            "f_ap_k"=f_ap_k,
    #            "BICsimk"=kbics2,
    #            "PAMsimk"=kpam2,
    #            "APsimk"=kap2)
    #     )
    #
    #
    #   }


}




