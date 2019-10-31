InDisc <- function(SCO, nquad = 30, model = "linear", approp = FALSE, display = TRUE){

  ## checks

  # check if the data is a data.frame

  if (is.data.frame(SCO)){
    SCO <- as.matrix(SCO)
  }

  # check if the data is a matrix

  f1 <- size(SCO)[1]
  f2 <- size(SCO)[2]

  if (f1 == 1 || f2 == 1){
    stop("The SCO argument has to be a matrix")
  }

  # check the number of items

  if (f2 < 20){ #binary 20, graded o continuo 15
    # stop or warning?
    warning("The number of items is probably too small for accurate PDD estimation.")
  }


  ################
  ## begin analysis


  #mean of each item
  mu <- optimbase::transpose(colMeans(SCO))

  if (model == "linear"){

    #minres using correlation
    out_fa <-psych::fa(SCO)

    #minres using covariance
    lam <- psych::fa(SCO, covar=T)$loadings
  }

  if (model == "graded"){

    # when graded, check if min is 0
    if (min(SCO)==0){
      SCO <- SCO + 1
    }

    tmp <- max(SCO) - min(SCO)
    if (tmp ==1){
      #minres using tetrachoric correlation
      out_fa<-psych::fa(SCO, cor = "tet")
    }
    if (tmp >1){
      #minres using polyhcoric correlation
      out_fa<-psych::fa(SCO, cor = "poly")
    }

    #minres using covariance
    lam <- psych::fa(SCO, covar=T)$loadings

  }

  alpha<-out_fa$loadings
  gof_dof <- out_fa$dof # degrees of freedom
  gof_chi <- out_fa$chi # empirical chi square
  gof_RMSR <- out_fa$rms
  gof_TLI <- out_fa$TLI
  gof_RMSEA <- as.numeric(out_fa$RMSEA[1])

  #1st priorgral for estimating residual variance and nodes

  OUT<-priorgral(alpha,nquad)

  # WHEN NON CONVERGENCY, RECOMPUTE USING A MORE STRICT PRIOR ON priorgral
  # OUT <- priorgral(alpha, nquad, 6, 4)

  vres<-OUT$resi
  nodos1<-OUT$nodth
  nodos2<-OUT$nodichi
  var_nodos<-OUT$var_nodos
  EVARI <- OUT$mvari

  # reap

  if (model=="linear"){

    FAC <- reap17c(mu, SCO, lam, vres, nodos1, nodos2, display)

  }
  if (model == "graded"){

    #1st threesholds
    THRES <- t(thresholds(SCO,min(SCO),max(SCO)))
    n <- dim(THRES)[2]
    zeros <- matrix(0,1,n)
    THRES <- rbind(THRES,zeros)

    FAC <- reap17g(THRES, SCO, alpha, vres, nodos1, nodos2, display)

  }

  # WHEN NON CONVERGENCY, FIND WHO AND RECOMPUTE USING A MORE STRICT PRIOR ON priorgral

  if (any(is.na(FAC))){

    is_NaN_FAC <- is.na(FAC)

    f1 <- size(FAC)[1]
    rows_NaN = ((is_NaN_FAC-1) %% f1) + 1

    OUT <- priorgral(alpha, nquad, 6, 4)

    vresb<-OUT$resi
    nodos1b<-OUT$nodth
    nodos2b<-OUT$nodichi

    if (model == "linear"){
      FAC_NaN <- reap17c(mu, SCO[rows_NaN,], lam, vresb, nodos1b, nodos2b, FALSE)
    }

    if (model == "graded"){
      FAC_NaN <- reap17g(THRES, SCO[rows_NaN,], alpha, vresb, nodos1b, nodos2b, FALSE)
    }

    if (any(is.na(FAC_NaN))){

      #still non convergency, more strict prior
      OUT <- priorgral(alpha, nquad, 10, 8)

      vresc<-OUT$resi
      nodos1c<-OUT$nodth
      nodos2c<-OUT$nodichi

      if (model == "linear"){
        FAC_NaN2 <- reap17c(mu, SCO[rows_NaN,], lam, vresc, nodos1c, nodos2c, FALSE)
      }

      if (model == "graded"){
        FAC_NaN2 <- reap17g(THRES, SCO[rows_NaN,], alpha, vresc, nodos1c, nodos2c, FALSE)
      }

      for (i in 1:(size(rows_NaN)[2])){
        FAC[rows_NaN[i],] <- FAC_NaN2[i,]
      }

    }

    else {
      #after the first cut, all values converge, replace the NaN values with the obtained
      for (i in 1:(size(rows_NaN)[2])){
        FAC[rows_NaN[i],] <- FAC_NaN[i,]
      }

    }

  }

  ##

  FAC <- as.matrix(FAC)
  colnames(FAC) <- c("theta","PDD","PSD (theta)","PSD (PDD)","th reli", "PDD reli")

  theta <- FAC[,1]
  person_var <- FAC[,2]
  PSD_theta <- FAC[,3]
  PSD_PDD <- FAC[,4]
  reli_th_i<-FAC[,5]
  reli_PDD_i<-FAC[,6]

  ## reliability

  reli_theta <- 1 - mean(PSD_theta^2)

  reli_PDD <- var_nodos / (var_nodos + mean(PSD_PDD^2))

  if (reli_PDD < 0){ #when the variance between person variability is extremely low
    reli_PDD <- 0
  }

  #average of the individual reliabilities

  aver_r_theta <- mean(reli_th_i)

  aver_r_PDD <- mean(reli_PDD_i)

  ############ appropiateness indices ############

  if (approp == TRUE){

    cvar <- mean(person_var)

    #calculate th using the mean of person variance

     if (model == "linear"){
       th0 <- reapth17c(mu, SCO, lam, vres, nodos1, cvar, display)[,1]
       OUT2 <- rlrtest(SCO, mu, lam, cvar ,vres, th0, theta, person_var)
     }
     else { #graded
       th0 <- reapth17g(THRES, SCO, alpha, vres, nodos1, cvar, display)[,1]
       OUT2 <-rlrtesg(THRES, SCO, alpha, cvar, vres, th0, theta, person_var)
     }

    LR <- OUT2[,1]
    S <- OUT2[,2]

    LR_stat <- mean(LR) # fit, the closer to 0, the better
    Q_Chi_square <- sum(S) # approximate chi square with N degrees of freedom

    OUT<-list("INDIES"=FAC, "degrees_of_freedom"=gof_dof, "Model_Chi_square"=gof_chi, "RMSR"=gof_RMSR, "TLI"=gof_TLI, "RMSEA"=gof_RMSEA, "EVARI"=EVARI, "reli_theta"=reli_theta, "aver_r_theta" = aver_r_theta, "reli_PDD"=reli_PDD, "aver_r_PDD"=aver_r_PDD,"LR_stat"=LR_stat, "Q_Chi_square"=Q_Chi_square)

  }
  else {

    OUT<-list("INDIES"=FAC, "Degrees_of_freedom"=gof_dof, "Model_Chi_square"=gof_chi, "RMSR"=gof_RMSR, "TLI"=gof_TLI, "RMSEA"=gof_RMSEA, "EVARI"=EVARI, "reli_theta"=reli_theta, "aver_r_theta" = aver_r_theta, "reli_PDD"=reli_PDD, "aver_r_PDD"=aver_r_PDD)

  }

  if (display==TRUE){
    return(OUT)
  }
  else {
    invisible(OUT)
  }

}
