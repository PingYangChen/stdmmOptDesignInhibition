
rOptDesignInhibitionPSO <- function(D_INFO, PSO_INFO) {
  
  # model, dSupp, nSupp, dPara, sUB(vec), sLB(vec), pUB(vec), pLB(vec)
  DESIGNINFO <- c(D_INFO$model, D_INFO$dSupp, 
                  D_INFO$nSupp, D_INFO$dPara, 
                  D_INFO$dsUpper, D_INFO$dsLower, 
                  D_INFO$paraUpper, D_INFO$paraLower)

  # LoopID, Maximize, Terminate by MaxIter, dSwarm, nSwarm, maxIter, Tolerence, c1, c2, vmaxK, w0, w1, wv, chi, numNB
  PSOOPTIONS <- rbind(
    c(1, 1, 1, D_INFO$nSupp*(D_INFO$dSupp + 1) - 1,  
      PSO_INFO$OUTER$nSwarm, PSO_INFO$OUTER$maxIter, 1e-6, 
      PSO_INFO$OUTER$c1, PSO_INFO$OUTER$c2, PSO_INFO$OUTER$vk, 
      PSO_INFO$OUTER$w0, PSO_INFO$OUTER$w1, PSO_INFO$OUTER$wv, 1, 0), # OUTER
    c(2, 0, 1, D_INFO$dPara,  
      PSO_INFO$INNER$nSwarm, PSO_INFO$INNER$maxIter, 1e-6, 
      PSO_INFO$INNER$c1, PSO_INFO$INNER$c2, PSO_INFO$INNER$vk, 
      PSO_INFO$INNER$w0, PSO_INFO$INNER$w1, PSO_INFO$INNER$wv, 1, 0),# INNER
    c(3, 1, 1, D_INFO$nSupp*(D_INFO$dSupp + 1) - 1, 
      PSO_INFO$LOACL$nSwarm, PSO_INFO$LOACL$maxIter, 1e-6, 
      PSO_INFO$LOACL$c1, PSO_INFO$LOACL$c2, PSO_INFO$LOACL$vk, 
      PSO_INFO$LOACL$w0, PSO_INFO$LOACL$w1, PSO_INFO$LOACL$wv, 1, 0),# LOCAL
    c(4, 0, 1, 1, 64, 400, 1e-6, 2, 2, 2, .95, .2, .8, 1, 0)  # ASSIST
  )
  
  pout <- aout <- dout <- NULL
  if (all(D_INFO$paraUpper == D_INFO$paraLower)) {
    # Locally D-optimal Design
    cputime <- system.time({
      pout <- NESTEDPSO(3, PSOOPTIONS, DESIGNINFO, DESIGNINFO[9:(9+D_INFO$dPara-1)])
    })[3]
    dout <- QUICKDISPERSION(pout$fGBest, pout$GBest, DESIGNINFO[9:(9+D_INFO$dPara-1)], DESIGNINFO)
  } else {
    # Standardized maximin Design
    cputime <- system.time({
      pout <- NESTEDPSO(1, PSOOPTIONS, DESIGNINFO, 0)
    })[3]
    aout <- MUSEARCH(pout$fGBest, pout$GBest, PSOOPTIONS, DESIGNINFO)
    dout <- QUICKDISPERSION(pout$fGBest, pout$GBest, c(aout$Mu, aout$Mu_w), DESIGNINFO)
  }
 
  list(DESIGN = matrix(c(pout$GBest, 1 - sum(pout$GBest[(D_INFO$nSupp*2+1):(D_INFO$nSupp*3-1)])), 3, D_INFO$nSupp, byrow = TRUE),
       CRIVAL = pout$fGBest,
       MAXDISP = dout$maxVal,
       PSO_CONV = pout$fGBestHist,
       CPUTIME = cputime,
       PSO_RESULT = pout,
       ANSWERING_SET = aout,
       DISPERSION = dout)
}

##
getD_INFO <- function(model = 1, paraUpper = c(1, 5, 3), paraLower = NULL) {
  modelInfo <- cbind(c("Competitive", "Noncompetitive", "Uncompetitive", "Mixed-type"),
                     c("Vmax, km, kic", "Vmax, km, kic", "Vmax, km, kiu", "Vmax, km, kic, kiu"))
  if (!(model %in% 1:4)) {
    stop("Model Index\n 
          1: Competitive model\n 2: noncompetitive model\n
          3: Uncompeitive model\n4: mixed-type model\n")
  }
  if ((model < 4) & !(length(paraUpper) == 3)) {
    stop(paste0("For ", modelInfo[1,model], "model, the model parameters are ", modelInfo[2,model]))
  }
  if ((model == 4) & !(length(paraUpper) == 4)) {
    stop(paste0("For ", modelInfo[1,model], "model, the model parameters are ", modelInfo[2,model]))
  }
  if (is.null(paraLower)) paraLower <- paraUpper
  
  list(model = model, dSupp = 2, 
       nSupp = ifelse(model < 4, 3, 4), 
       dPara = ifelse(model < 4, 3, 4),  
       dsUpper = c(30, 60), dsLower = c(0, 0), 
       paraUpper = paraUpper, paraLower = paraLower)
}
##
getPSO_INFO <- function(nSwarm = c(64, 32, 64), maxIter = c(100, 50, 100)) {
  PSO_INFO <- lapply(1:3, function(i) {
    list(nSwarm = nSwarm[i], maxIter = maxIter[i], 
         c1 = 2, c2 = 2, vk = 2, w0 = .95, w1 = .2, wv = .8)
  })
  names(PSO_INFO) <- c("OUTER", "INNER", "LOCAL")
  PSO_INFO
}
##
drawEquiv <- function(PSORESULT, D_INFO) {
  
  nSupp <- D_INFO$nSupp
  dLB <- D_INFO$dsLower
  dUB <- D_INFO$dsUpper
  
  pout <- PSORESULT$PSO_RESULT
  dout <- PSORESULT$DISPERSION
  cL <- 25
  colormap <- rgb(c(rep(0,1*cL), seq(0,0.95,length=1*cL), 1,      #R
                    rep(0.95,1*cL), seq(0.95,0.5,length = 1*cL)),   
                  c(rep(0,1*cL), seq(0,0.9,length=1*cL), 1,				#G
                    seq(0.9,0,length=1*cL), rep(0,1*cL)), 
                  c(seq(0.5,0.95,length = 1*cL), rep(0.95,1*cL), 1, #B
                    seq(0.95,0,length = 1*cL), rep(0,1*cL)), 
                  maxColorValue = 1)
  
  cM <- ifelse(max(abs(dout$dispVal)) == 0, 1, max(abs(dout$dispVal)))
  
  designOut <- matrix(c(pout$GBest, 1 - sum(pout$GBest[(nSupp*2+1):(nSupp*3-1)])), 3, 
                      nSupp, byrow = TRUE)
  
  lay <- layout(matrix(1:2, 1, 2), widths = c(13, 2))
  par(mar = c(5, 4, 1, 1) + 0.1, las = 1)
  image(dout$sGrid, dout$iGrid, dout$dispVal, 
        col = colormap, zlim = c(-1,1)*cM, axes = FALSE, xlab = "[s]", ylab = "[i]")
  contour(dout$sGrid, dout$iGrid, dout$dispVal, add = TRUE, col = "gray90", labcex = 1.2)
  points(t(designOut[1:2,]), pch = 19, cex = 1.8, col = "forestgreen")
  points(t(designOut[1:2,]), pch = 19, cex = 1.2, col = "chartreuse1")
  axis(1, at = seq(dLB[1], dUB[1], length = 5), labels = seq(dLB[1], dUB[1], length = 5),
       lty = 0, cex.axis = 1, line = -.5)
  axis(2, at = seq(dLB[2], dUB[2], length = 5), labels = seq(dLB[2], dUB[2], length = 5),
       lty = 0, cex.axis = 1, line = -.5)
  par(mar = c(5, 0, 1, 2) + 0.1, las = 1)
  image(1, seq(-cM, cM, length = 100), t(as.matrix(seq(-cM, cM, length = 100))), 
        col = colormap, axes = F, xlab = "", ylab = "", zlim = c(-cM, cM))
  mtext(round(cM, 1), 3, cex = .8); mtext(-round(cM, 1), 1, cex = .8); mtext(0, 4, cex = .8)
  
}