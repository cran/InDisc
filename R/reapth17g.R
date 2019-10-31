reapth17g <- function(DIF, RES, alpha, vres, nodos1, vari, display){

  #begin timer
  ptm <- proc.time()

  n<-size(RES)[1]

  for (i in 1:n){

    if (display==TRUE){
      if (i==1){
        #waitbar
        cat('Computing Appropiateness indices: Please wait \n')
        pb <- txtProgressBar(min = 0, max = n-1, style = 3)
      }
    }

    OUT <- eap17gth(DIF, RES[i,], alpha, vres, nodos1, vari)
    tmp1<-OUT$th
    tmp2<-OUT$st

    if (i == 1){
      FAC <- c(tmp1,tmp2)
    }
    else {
      FAC<-rbind(FAC,c(tmp1,tmp2))
    }


    if (display==TRUE){
      Sys.sleep(0.1)
      # update progress bar
      setTxtProgressBar(pb, i-1)
    }
  }

  if (display==TRUE){
    close(pb)

    time.taken <- proc.time() - ptm

    seconds<-time.taken[3]
    hours<-floor(seconds/3600)
    minutes<-floor(seconds/60)
    seconds<-floor(seconds-(minutes*60))
    total_time<-sprintf('%02.0f:%02.0f:%02.0f',hours,minutes,seconds)
    cat('\n')
    cat(sprintf('Computing time: %s \n\n',total_time))
  }


  return(FAC)



}
