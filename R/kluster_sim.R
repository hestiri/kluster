######## 3 kluster functions



kluster_sim <- function(data,
                        clusters, #number of clusters we know
                        iter_sim, # number of simulation iterations
                        iter_klust, #number of iterations for clustering with sample_n size x
                        smpl, #size of the sample_n to be taken with replacement out of data
                        algorithm = "Default", #select analysis algorithm from BIC, PAMK, CAL, and AP
                        cluster = FALSE # if TURE it'll do clustering which will take a lot longer!
) {

  size <- nrow(data)
  if (algorithm == 'BIC') {
    ##starting to store results from different algorithms
    tic()
    ##now let's compute optimal ks with BIC
    BIC.best <- dim(Mclust(as.matrix(data), G=1:15)$z)[2]
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
        dat2 <- as.matrix(data[sample(nrow(data), smpl, replace=TRUE), ])
        kbics2[[i]] <- dim(Mclust(dat2, G=1:15)$z)[2]
        rm(dat2)
        t <- toc()
        t2[[i]] <- t$toc - t$tic
      }
      kbicssum[[j]] <- mean(as.numeric(kbics2))
    }
    tBIC_kluster <- mean(unlist(t2))
    m_BIC_k <- round(mean(as.numeric(kbics2)),0)
    f_BIC_k <- as.numeric(names(which.max(table(unlist(kbics2)))))



    method <- c("BIC.best","BIC_kluster_mean","BIC_kluster_frq")
    ptime <- c(tBIC,tBIC_kluster/iter_sim,tBIC_kluster/iter_sim)
    k_num <- c(BIC.best,m_BIC_k,f_BIC_k)


    sim <- data.frame(method,k_num,ptime)
    sim$k_orig <- clusters
    sim$e  <- k_num - clusters
    sim$n <- size




    if (cluster == TRUE) {

      ####cluster and visualize the performance
      km1 <- hkmeans(data, k = BIC.best,iter.max = 300)
      km5 <- hkmeans(data, k = m_BIC_k,iter.max = 300)



      par(mfrow=c(1,2))
      plot(data, col = km1$cluster, pch = 19, frame = FALSE,
           main = paste0("HK-means with k = ",BIC.best,""))

      plot(data, col = km5$cluster, pch = 19, frame = FALSE,
           main = paste0("HK-means with kluster BIC process, k = ",m_BIC_k,""))



      par(mfrow=c(1,1))
      boxplot(round(as.numeric(kbicssum),0),
              main = paste0("kluster optimum cluster number from BIC w/ resampling.
                            Mean resampling estimate = ",m_BIC_k,"
                            ",
                            " and ordinary BIC suggested ",BIC.best," clusters."))

    }
    return(
      list("sim"=sim,
           "m_kluster"=m_bic_k,
           "alg_orig"=bic.best,
           "f_kluster"=f_bic_k)
    )




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
          dat2 <- data[sample(nrow(data), smpl, replace=TRUE), ]
          kpam2[[i]] <- pamk(dat2)$nc
          rm(dat2)
          t <- toc()
          t2[[i]] <- t$toc - t$tic
        }
        kpamsum[[j]] <- mean(as.numeric(kpam2))
      }
      tpam_kluster <- mean(unlist(t2))
      m_pam_k <- round(mean(as.numeric(kpam2)),0)
      f_pam_k <- as.numeric(names(which.max(table(unlist(kpam2)))))



      method <- c("pamk.best","pamk_kluster_mean","pamk_kluster_frq")
      ptime <- c(tpamk,tpam_kluster/iter_sim,tpam_kluster/iter_sim)
      k_num <- c(pamk.best,m_pam_k,f_pam_k)


      sim <- data.frame(method,k_num,ptime)
      sim$k_orig <- clusters
      sim$e  <- k_num - clusters
      sim$n <- size


      if (cluster == TRUE) {
        ####cluster and visualize the performance

        km2 <- hkmeans(data, k = pamk.best,iter.max = 300)
        km6 <- hkmeans(data, k = m_pam_k,iter.max = 300)



        par(mfrow=c(1,2))

        plot(data, col = km2$cluster, pch = 19, frame = FALSE,
             main = paste0("HK-means with k = ",pamk.best,""))
        points(km2$centers, col = 1:pamk.best, pch = 8, cex = 3)

        plot(data, col = km6$cluster, pch = 19, frame = FALSE,
             main = paste0("HK-means with kluster PAM process, k = ",m_pam_k,""))
        points(km6$centers, col = 1:m_pam_k, pch = 8, cex = 3)



        par(mfrow=c(1,1))

        boxplot(round(as.numeric(kpamsum),0),
                main = paste0("kluster optimum cluster number from PAM w/ resampling.
                              Mean resampling estimate = ",m_pam_k,"
                              ",
                              " and ordinary PAM suggested ",pamk.best," clusters."))

      }

      return(
        list("sim"=sim,
             "m_kluster"=m_pam_k,
             "alg_orig"=pamk.best,
             "f_kluster"=f_pam_k)
      )


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
            dat2 <- data[sample(nrow(data), smpl, replace=TRUE), ]
            kcal2[[i]] <- as.numeric(which.max(cascadeKM(dat2, 1, 15, iter = 1000)$results[2,]))
            rm(dat2)
            t <- toc()
            t2[[i]] <- t$toc - t$tic
          }
          kcalsum[[j]] <- mean(as.numeric(kcal2))
        }
        tcal_kluster <- mean(unlist(t2))
        m_cal_k <- round(mean(as.numeric(kcal2)),0)
        f_cal_k <- as.numeric(names(which.max(table(unlist(kcal2)))))



        method <- c("calinski.best","cal_kluster_mean","cal_kluster_frq")
        ptime <- c(tcalinski,tcal_kluster/iter_sim,tcal_kluster/iter_sim)
        k_num <- c(calinski.best,m_cal_k,f_cal_k)


        sim <- data.frame(method,k_num,ptime)
        sim$k_orig <- clusters
        sim$e <- k_num - clusters
        sim$n <- size

        if (cluster == TRUE) {
          ####cluster and visualize the performance
          km3 <- hkmeans(data, k = calinski.best,iter.max = 300)
          km7 <- hkmeans(data, k = m_cal_k,iter.max = 300)



          par(mfrow=c(1,2))

          plot(data, col = km3$cluster, pch = 19, frame = FALSE,
               main = paste0("HK-means with k = ",calinski.best,""))
          points(km3$centers, col = 1:calinski.best, pch = 8, cex = 3)


          plot(data, col = km7$cluster, pch = 19, frame = FALSE,
               main = paste0("HK-means with kluster Calinski process, k = ",m_cal_k,""))
          points(km7$centers, col = 1:m_cal_k, pch = 8, cex = 3)



          par(mfrow=c(1,1))

          boxplot(round(as.numeric(kcalsum),0),
                  main = paste0("kluster optimum cluster number from Calinski w/ resampling.
                                Mean resampling estimate = ",m_cal_k,"
                                and ordinary Calinski suggested ",calinski.best," clusters."))

        }
        return(
          list("sim"=sim,
               "m_kluster"=m_cal_k,
               "alg_orig"=calinski.best,
               "f_kluster"=f_cal_k)
        )



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
              dat2 <- data[sample(nrow(data), smpl, replace=TRUE), ]
              kap2[[i]] <- length(apcluster(negDistMat(r=2), dat2)@clusters)
              rm(dat2)
              t <- toc()
              t2[[i]] <- t$toc - t$tic
            }
            kapsum[[j]] <- mean(as.numeric(kap2))
          }
          tap_kluster <- mean(unlist(t2))
          m_ap_k <- round(mean(as.numeric(kap2)),0)
          f_ap_k <- as.numeric(names(which.max(table(unlist(kap2)))))



          method <- c("apclus.best","ap_kluster_mean","ap_kluster_frq")
          ptime <- c(tap,tap_kluster/iter_sim,tap_kluster/iter_sim)
          k_num <- c(apclus.best,m_ap_k,f_ap_k)


          sim <- data.frame(method,k_num,ptime)
          sim$k_orig <- clusters
          sim$e <- k_num - clusters
          sim$n <- size

          if (cluster == TRUE) {
            ####cluster and visualize the performance

            km4 <- hkmeans(data, k = apclus.best,iter.max = 300)

            km8 <- hkmeans(data, k = m_ap_k,iter.max = 300)



            par(mfrow=c(1,2))

            plot(data, col = km4$cluster, pch = 19, frame = FALSE,
                 main = paste0("HK-means with k = ",apclus.best,""))
            points(km4$centers, col = 1:apclus.best, pch = 8, cex = 3)

            plot(data, col = km8$cluster, pch = 19, frame = FALSE,
                 main = paste0("HK-means with kluster AP process, k = ",m_ap_k,""))
            points(km8$centers, col = 1:m_ap_k, pch = 8, cex = 3)


            par(mfrow=c(1,1))

            boxplot(round(as.numeric(kapsum),0),
                    main = paste0("kluster optimum cluster number from AP w/ resampling.
                                  Mean resampling estimate = ",m_ap_k,"
                                  and ordinary AP suggested ",apclus.best," clusters."))

          }
          return(
            list("sim"=sim,
                 "m_kluster"=m_ap_k,
                 "alg_orig"=apclus.best,
                 "f_kluster"=f_ap_k)

          )

        } else
          if (algorithm == "Default") {


            ##starting to store results from different algorithms
            tic()
            ##now let's compute optimal ks with BIC
            BIC.best <- dim(Mclust(as.matrix(data), G=1:15)$z)[2]
            tBIC <- toc()
            tBIC <- as.numeric(tBIC$toc - tBIC$tic)

            tic()
            pamk.best <- pamk(data)$nc
            tpamk <- toc()
            tpamk <- as.numeric(tpamk$toc - tpamk$tic)

            tic()
            calinski.best <- as.numeric(which.max(cascadeKM(data, 1, 15, iter = 500)$results[2,]))
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
                dat2 <- as.matrix(data[sample(nrow(data), smpl, replace=TRUE), ])
                kbics2[[i]] <- dim(Mclust(dat2, G=1:15)$z)[2]
                rm(dat2)
                t <- toc()
                t2[[i]] <- t$toc - t$tic
              }
              kbicssum[[j]] <- mean(as.numeric(kbics2))
            }
            tBIC_kluster <- mean(unlist(t2))
            m_BIC_k <- round(mean(as.numeric(kbics2)),0)
            f_BIC_k <- as.numeric(names(which.max(table(unlist(kbics2)))))



            ###pam method
            kpamsum <- list()
            kpam2 <- list()
            t2 <- list()
            tpam_kluster <- list()
            for (j in 1:iter_sim) {
              for (i in 1:iter_klust) {
                tic()
                dat2 <- data[sample(nrow(data), smpl, replace=TRUE), ]
                kpam2[[i]] <- pamk(dat2)$nc
                rm(dat2)
                t <- toc()
                t2[[i]] <- t$toc - t$tic
              }
              kpamsum[[j]] <- mean(as.numeric(kpam2))
            }
            tpam_kluster <- mean(unlist(t2))
            m_pam_k <- round(mean(as.numeric(kpam2)),0)
            f_pam_k <- as.numeric(names(which.max(table(unlist(kpam2)))))




            kcalsum <- list()
            kcal2 <- list()
            t2 <- list()
            tcal_kluster <- list()
            for (j in 1:iter_sim) {
              for (i in 1:iter_klust) {
                tic()
                dat2 <- data[sample(nrow(data), smpl, replace=TRUE), ]
                kcal2[[i]] <- as.numeric(which.max(cascadeKM(dat2, 1, 15, iter = 1000)$results[2,]))
                rm(dat2)
                t <- toc()
                t2[[i]] <- t$toc - t$tic
              }
              kcalsum[[j]] <- mean(as.numeric(kcal2))
            }
            tcal_kluster <- mean(unlist(t2))
            m_cal_k <- round(mean(as.numeric(kcal2)),0)
            f_cal_k <- as.numeric(names(which.max(table(unlist(kcal2)))))


            ###AP
            kapsum <- list()
            kap2 <- list()
            t2 <- list()
            tap_kluster <- list()
            for (j in 1:iter_sim) {
              for (i in 1:iter_klust) {
                tic()
                dat2 <- data[sample(nrow(data), smpl, replace=TRUE), ]
                kap2[[i]] <- length(apcluster(negDistMat(r=2), dat2)@clusters)
                rm(dat2)
                t <- toc()
                t2[[i]] <- t$toc - t$tic
              }
              kapsum[[j]] <- mean(as.numeric(kap2))
            }
            tap_kluster <- mean(unlist(t2))
            m_ap_k <- round(mean(as.numeric(kap2)),0)
            f_ap_k <- as.numeric(names(which.max(table(unlist(kap2)))))





            method <- c("BIC.best","pamk.best","calinski.best","apclus.best",
                        "BIC_kluster_mean","cal_kluster_mean","pam_kluster_mean","ap_kluster_mean",
                        "BIC_kluster_frq","cal_kluster_frq","pam_kluster_frq","ap_kluster_frq")
            ptime <- c(tBIC,tpamk,tcalinski,tap,
                       tBIC_kluster/iter_sim,tcal_kluster/iter_sim,tpam_kluster/iter_sim,tap_kluster/iter_sim,
                       tBIC_kluster/iter_sim,tcal_kluster/iter_sim,tpam_kluster/iter_sim,tap_kluster/iter_sim
            )
            k_num <- c(BIC.best,pamk.best,calinski.best,apclus.best,
                       m_BIC_k,m_cal_k,m_pam_k,m_ap_k,
                       f_BIC_k,f_cal_k,f_pam_k,f_ap_k)


            sim <- data.frame(method,k_num,ptime)
            sim$k_orig <- clusters
            sim$e <- k_num - clusters
            sim$n <- size





            if (cluster == TRUE) {
              ####cluster and visualize the performance
              km1 <- hkmeans(data, k = BIC.best,iter.max = 300)
              km2 <- hkmeans(data, k = pamk.best,iter.max = 300)
              km3 <- hkmeans(data, k = calinski.best,iter.max = 300)
              km4 <- hkmeans(data, k = apclus.best,iter.max = 300)
              km5 <- hkmeans(data, k = m_BIC_k,iter.max = 300)
              km6 <- hkmeans(data, k = m_pam_k,iter.max = 300)
              km7 <- hkmeans(data, k = m_cal_k,iter.max = 300)
              km8 <- hkmeans(data, k = m_ap_k,iter.max = 300)



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
                   main = paste0("HK-means with kluster BIC process, k = ",m_BIC_k,""))
              points(km5$centers, col = 1:m_BIC_k, pch = 8, cex = 3)

              plot(data, col = km6$cluster, pch = 19, frame = FALSE,
                   main = paste0("HK-means with kluster PAM process, k = ",m_pam_k,""))
              points(km6$centers, col = 1:m_pam_k, pch = 8, cex = 3)

              plot(data, col = km7$cluster, pch = 19, frame = FALSE,
                   main = paste0("HK-means with kluster Calinski process, k = ",m_cal_k,""))
              points(km7$centers, col = 1:m_cal_k, pch = 8, cex = 3)

              plot(data, col = km8$cluster, pch = 19, frame = FALSE,
                   main = paste0("HK-means with kluster AP process, k = ",m_ap_k,""))
              points(km8$centers, col = 1:m_ap_k, pch = 8, cex = 3)


              par(mfrow=c(2,2))
              boxplot(round(as.numeric(kbicssum),0),
                      main = paste0("kluster optimum cluster number from BIC w/ resampling.
                                    Mean resampling estimate = ",m_BIC_k,"
                                    ",
                                    " and ordinary BIC suggested ",BIC.best," clusters."))

              boxplot(round(as.numeric(kpamsum),0),
                      main = paste0("kluster optimum cluster number from PAM w/ resampling.
                                    Mean resampling estimate = ",m_pam_k,"
                                    ",
                                    " and ordinary PAM suggested ",pamk.best," clusters."))
              boxplot(round(as.numeric(kcalsum),0),
                      main = paste0("kluster optimum cluster number from Calinski w/ resampling.
                                    Mean resampling estimate = ",m_cal_k,"
                                    and ordinary Calinski suggested ",calinski.best," clusters."))

              boxplot(round(as.numeric(kapsum),0),
                      main = paste0("kluster optimum cluster number from AP w/ resampling.
                                    Mean resampling estimate = ",m_ap_k,"
                                    and ordinary AP suggested ",apclus.best," clusters."))

            }
            return(
              list("sim"=sim,
                   "m_BIC_k"=m_bic_k,
                   "m_pam_k"=m_pam_k,
                   "m_ap_k"=m_ap_k,
                   "f_BIC_k"=f_bic_k,
                   "f_pam_k"=f_pam_k,
                   "f_ap_k"=f_ap_k)
            )


          }

          }


