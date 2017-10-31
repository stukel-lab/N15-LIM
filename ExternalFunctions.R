
## Function sets the A matrix equations involved in N15 equations
ResetRN15 <- function(Aa, d15N, del15NupNO3, del15NMes, del15NDoc, del15NDet){
    
  ## Set Known Fractionation Values
  del15NNH4 <- d15N[1]
  del15NCya <- d15N[2]
  del15NDtm <- d15N[3]
  del15NHnf <- d15N[4]
  del15NMic <- d15N[5]
  del15NNO3 <- d15N[6]
  
  ## Set fractionation coefficients
  RN2 <- 0.0036765
  RupNO3 <- del15NupNO3 / 1000 * RN2 + RN2
  RNO3 <- del15NNO3 / 1000 * RN2 + RN2
  RNH4 <- del15NNH4 / 1000 * RN2 + RN2
  RDtm <- del15NDtm / 1000 * RN2 + RN2
  RCya <- del15NCya / 1000 * RN2 + RN2
  RHnf <- del15NHnf / 1000 * RN2 + RN2
  RMic <- del15NMic / 1000 * RN2 + RN2
  RMes <- del15NMes / 1000 * RN2 + RN2
  RDet <- del15NDet / 1000 * RN2 + RN2
  RDoc <- del15NDoc / 1000 * RN2 + RN2
  
  ## Set Epsilons
  Eps_NO3up <- -5
  Eps_NH4up <- -10
  Eps_Exc <- -5
  Eps_Eg <- -2
  Eps_Remin <- -1
#  Eps_TL <- 3.5
    
  ## Calculate coefficinets
  Alpha_NO3up <- exp(Eps_NO3up / 1000)
  Alpha_NH4up <- exp(Eps_NH4up / 1000)
  Alpha_Exc <- exp(Eps_Exc / 1000)
  Alpha_Eg <- exp(Eps_Eg / 1000)
  Alpha_Remin <- exp(Eps_Remin / 1000)
#  Alpha_TL <- exp(Eps_TL/1000)
  
  ## Manually apply functions
  Aa[1,1] <- RupNO3
  Aa[1,3] <- -RNO3 * Alpha_NO3up
  Aa[1,10] <- -RNO3 * Alpha_NO3up
  Aa[2,4] <- -RNH4 * Alpha_NH4up
  Aa[2,11] <- -RNH4 * Alpha_NH4up
  Aa[2,18] <- RHnf * Alpha_Exc
  Aa[2,21] <- RMic * Alpha_Exc
  Aa[2,25] <- RMes * Alpha_Exc
  Aa[2,32] <- RDoc * Alpha_Remin
  Aa[3,9] <- RN2
  Aa[3,10] <- RNO3 * Alpha_NO3up
  Aa[3,11] <- RNH4 * Alpha_NH4up
  Aa[3,12] <- -RCya
  Aa[3,13] <- -RCya
  Aa[3,14] <- -RCya
  Aa[3,15] <- -RCya
  Aa[4,2] <- RN2
  Aa[4,3] <- RNO3 * Alpha_NO3up
  Aa[4,4] <- RNH4 * Alpha_NH4up
  Aa[4,5] <- -RDtm
  Aa[4,6] <- -RDtm
  Aa[4,7] <- -RDtm
  Aa[4,8] <- -RDtm
  Aa[5,12] <- RCya
  Aa[5,16] <- -RHnf
  Aa[5,17] <- -RHnf
  Aa[5,18] <- -(RHnf * Alpha_Exc)
  Aa[5,19] <- -(RHnf * Alpha_Eg)
  Aa[5,20] <- -(RHnf * Alpha_Exc)
  Aa[5,28] <- RDet
  Aa[6,5] <- RDtm
  Aa[6,13] <- RCya
  Aa[6,16] <- RHnf
  Aa[6,21] <- -(RMic * Alpha_Exc)
  Aa[6,22] <- -RMic
  Aa[6,23] <- -(RMic * Alpha_Eg)
  Aa[6,24] <- -(RMic * Alpha_Exc)
  Aa[6,29] <- RDet
  Aa[7,6] <- RDtm
  Aa[7,17] <- RHnf
  Aa[7,22] <- RMic
  Aa[7,25] <- -(RMes * Alpha_Exc)
  Aa[7,26] <- -(RMes * Alpha_Eg)
  Aa[7,27] <- -(RMes * Alpha_Exc)
  Aa[7,30] <- RDet
  Aa[7,34] <- -RMes
  Aa[8,7] <- RDtm
  Aa[8,14] <- RCya
  Aa[8,19] <- RHnf * Alpha_Eg
  Aa[8,23] <- RMic * Alpha_Eg
  Aa[8,26] <- RMes * Alpha_Eg
  Aa[8,28] <- -RDet
  Aa[8,29] <- -RDet
  Aa[8,30] <- -RDet
  Aa[8,31] <- -RDet * Alpha_Remin
  Aa[8,33] <- -RDet
  Aa[9,8] <- RDtm
  Aa[9,15] <- RCya
  Aa[9,20] <- RHnf * Alpha_Exc
  Aa[9,24] <- RMic * Alpha_Exc
  Aa[9,27] <- RMes * Alpha_Exc
  Aa[9,31] <- RDet * Alpha_Remin
  Aa[9,32] <- -RDoc * Alpha_Remin
  Aa[10,1] <- -RupNO3
  Aa[10,2] <- -RN2
  Aa[10,9] <- -RN2
  Aa[10,33] <- RDet
  Aa[10,34] <- RMes
  
  
  return(Aa)
}

### Notify Function ###
pb = txtProgressBar(0, 1, style=3) ## Initialized on source (meant to be run from commandline, e.g. rscript, bat file, sh, etc.)
notify = function(i, n) {
    setTxtProgressBar(pb, i/n)
}
