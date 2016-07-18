library(shiny)
library(Rcpp)
library(inline)
library(RcppArmadillo)

sourceCpp(file.path(getwd(), "kernel/psokernel.cpp"))

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
	
	# Define reactive objects
	appInput <- reactiveValues(designInfo = 0, psoOptions = NULL)
	appMsg <- reactiveValues(errMsg = NULL, result = NULL, applog = NULL, cpustatus = 0)
	appOutput <- reactiveValues(pout = NULL, aout = NULL, dout = NULL, cputime = NULL)
	
	# Collect Input value
	COLLECTINPUT <- reactive({
		modelIdx <- as.numeric(input$inhibition)
		nSupp <- dPara <- ifelse(modelIdx > 3, 4, 3)
		dDesign <- nSupp*3 - 1
		pUB <- {switch(input$inhibition,
												"0" = {0},
												"1" = {c(input$v, input$kmu, input$kcu)},
												"2" = {c(input$v, input$kmu, input$kcu)},
												"3" = {c(input$v, input$kmu, input$kuu)},
												"4" = {c(input$v, input$kmu, input$kcu, input$kuu)}
											)}
		pLB <- {switch(input$inhibition,
												"0" = {0},
												"1" = {c(input$v, input$kml, input$kcl)},
												"2" = {c(input$v, input$kml, input$kcl)},
												"3" = {c(input$v, input$kml, input$kul)},
												"4" = {c(input$v, input$kml, input$kcl, input$kul)}
											)}
		dUB <- c(input$su, input$iu)
		dLB <- c(input$sl, input$il)
		appInput$designInfo <- c(modelIdx, 2, nSupp, dPara, dUB, dLB, pUB, pLB) 
		# LoopID, Maximize, Terminate by MaxIter, dSwarm, nSwarm, maxIter, Tolerence, c1, c2, vmaxK, w0, w1, wv, chi, numNB
		appInput$psoOptions <- rbind(
			c(1, 1, 1, dDesign, input$L1ns, input$L1nt, 1e-6, input$c1, input$c2, input$vk, input$w0, input$w1, input$wv, 1, 0), # OUTER
			c(2, 0, 1, 	 dPara, input$L2ns, input$L2nt, 1e-6, input$c1, input$c2, input$vk, input$w0, input$w1, input$wv, 1, 0), # INNER
			c(3, 1, 1, dDesign, input$L3ns, input$L3nt, 1e-6, input$c1, input$c2, input$vk, input$w0, input$w1, input$wv, 1, 0), # LOCAL
			c(4, 0, 1, 			 1, input$L4ns, input$L4nt, 1e-6, input$c1, input$c2, input$vk, input$w0, input$w1, input$wv, 1, 0)  # ASSIST
		)
	})
	
	# Check Input value
	CHECKINPUT <- reactive({
		dPara <- ifelse(appInput$designInfo[1] > 3, 4, 3)
		dUB <- appInput$designInfo[5:6]
		dLB <- appInput$designInfo[7:8]
		pUB <- appInput$designInfo[9:(9 + dPara - 1)]
		pLB <- appInput$designInfo[(9 + dPara):length(appInput$designInfo)]
		nSwarm <- appInput$psoOptions[,5]
		nIter <- appInput$psoOptions[,6]
		tmpMsg <- NULL
		if (appInput$designInfo[1] == 0) {
			tmpMsg <- c(tmpMsg, "Please choose one inhibition model. \n")
		} else {
			if (!all(pUB >= pLB)) {
				tmpMsg <- c(tmpMsg, "ERROR: The upper bound of the parameter space \n should be larger or equal to its lower bound. \n")
			} 
			if (!all(dUB > dLB)) {
				tmpMsg <- c(tmpMsg, "ERROR: The upper bound of the design space \n should be larger or equal to its lower bound. \n")
			}
			if (!all(nSwarm > 0)) {
				tmpMsg <- c(tmpMsg, "ERROR: The swarm size of PSO should be positive. \n")
			}
			if (!all(nIter >= 0)) {
				tmpMsg <- c(tmpMsg, "ERROR: The iterations of PSO should be zero(initialization) or positive. \n")
			}		
		}
		if (is.null(tmpMsg)) {
			appMsg$errMsg <- NULL
			appMsg$cpustatus <- 1
		} else {
			appMsg$errMsg <- tmpMsg
			appMsg$cpustatus <- 0
		}
	})
	
	# The Main Function
	DOSEARCH <- reactive({
		dPara <- appInput$designInfo[4]
		pUB <- appInput$designInfo[9:(9 + dPara - 1)]
		pLB <- appInput$designInfo[(9 + dPara):length(appInput$designInfo)]
		if (all((pUB - pLB) == 0)) {
			# Locally Design
			appOutput$cputime <- system.time({
				appOutput$pout <- NESTEDPSO(3, appInput$psoOptions, appInput$designInfo, pUB)
			})[3]
			appOutput$dout <- QUICKDISPERSION(appOutput$pout$fGBest, appOutput$pout$GBest, c(pUB, 1), appInput$designInfo)
		} else {
			# Standardized maximin Design
			appOutput$cputime <- system.time({
				appOutput$pout <- NESTEDPSO(1, appInput$psoOptions, appInput$designInfo, pUB)
			})[3]
			appOutput$aout <- MUSEARCH(appOutput$pout$fGBest, appOutput$pout$GBest, appInput$psoOptions, appInput$designInfo)
			appOutput$dout <- QUICKDISPERSION(appOutput$pout$fGBest, appOutput$pout$GBest, c(appOutput$aout$Mu, appOutput$aout$Mu_w), appInput$designInfo)
		}
		appMsg$cpustatus <- 0
	})
	
	# UI Reactive Objects
	observeEvent(input$run, {
		isolate({
			COLLECTINPUT()
			CHECKINPUT()
			if (appInput$designInfo[1] > 0) {
				if (appMsg$cpustatus == 1) {
					DOSEARCH()
				}
			}
		})
	})
		
	output$msg <- renderText({
		ifelse(appMsg$cpustatus == 1, "CPU Status: Busy", "CPU Status: Idle")
	})		
	
	output$result <- renderPrint({
		if ((appInput$designInfo[1] > 0) & (!is.null(appOutput$pout))) {
			nSupp <- appInput$designInfo[3]
			dPara <- appInput$designInfo[4]
			dUB <- appInput$designInfo[5:6]
			dLB <- appInput$designInfo[7:8]
			pUB <- appInput$designInfo[9:(9 + dPara - 1)]
			pLB <- appInput$designInfo[(9 + dPara):length(appInput$designInfo)]			
			inhiInfo <- {switch(appInput$designInfo[1],
													"1" = list(NAME = "Competitive Model", 
																		 PARNAME = c("km", "kic"),
																		 PARSPNAME = c("(km, kic) in ")
																		),
													"2" = list(NAME = "Noncompetitive Model",
																		 PARNAME = c("km", "kic"),
																		 PARSPNAME = c("(km, kic) in ")
																		),
													"3" = list(NAME = "Uncompetitive Model",
																		 PARNAME = c("km", "kiu"),
																		 PARSPNAME = c("(km, kiu) in ")
																		),
													"4" = list(NAME = "Mixed-type Model",
																		 PARNAME = c("km", "kic", "kiu"),
																		 PARSPNAME = c("(km, kic, kiu) in ")
																		)
												)}
			
			if (all((pUB - pLB) == 0)) {
				paraMsg <- paste0("Parameter: ", paste0(inhiInfo$PARNAME, " = ", pUB[-1]))
				designName <- "Locally "
			} else {
				paraMsg <- paste0("Parameter Space: ", 
										paste0(inhiInfo$PARSPNAME, 
											paste0(paste0("[", pLB[-1], ", ", pUB[-1], "]"), collapse = " x ")))
				designName <- "Standarized Maximin "
			}
			
			cat(paste0("Inhibition Model: ", inhiInfo$NAME, "\n"))
			cat(paste0(paraMsg, "\n"))
			cat(paste0("Design Space: (s, i) in [", dLB[1], ", ", dUB[1], "] x [", dLB[2], ", ", dUB[2], "]\n"))
			cat("\n")
			if (is.null(appMsg$errMsg)) {
				designOut <- matrix(c(appOutput$pout$GBest, 1 - sum(appOutput$pout$GBest[(nSupp*2+1):(nSupp*3-1)])), 3, nSupp, byrow = TRUE)
				dimnames(designOut) <- list(c("s", "i", "weight"), paste0("Point ", 1:nSupp))
				designOut <- round(designOut, 4)
				cat(paste0("PSO-generated ", designName, "Design: \n"))
				print(designOut)
				cat("\n")
				cat(paste0("Maximal Dispersion Value: ", round(appOutput$dout$maxVal, 4), "\n"))
				cat(paste0("CPU Time: ", round(appOutput$cputime, 2)," seconds. \n"))
			} else {
				for (i in 1:length(appMsg$errMsg)) {
					cat(appMsg$errMsg[i])
				}
			}
		} else if ((is.null(appOutput$pout)) & (!is.null(appMsg$errMsg))) {
			for (i in 1:length(appMsg$errMsg)) {
				cat(appMsg$errMsg[i])
			}
		}
	})
	
	output$gpath <- renderPlot({
		if ((appInput$designInfo[1] > 0) & (!is.null(appOutput$pout))) {
			dPara <- appInput$designInfo[4]
			dUB <- appInput$designInfo[5:6]
			dLB <- appInput$designInfo[7:8]
			pUB <- appInput$designInfo[9:(9 + dPara - 1)]
			pLB <- appInput$designInfo[(9 + dPara):length(appInput$designInfo)]					
			par(mar = c(5, 4, 1, 2) + 0.1, las = 1)
			if (all((pUB - pLB) == 0)) {
				plot(1:length(appOutput$pout$fGBestHist), appOutput$pout$fGBestHist, 
					type = "l", xlab = "Iteration", ylab = "log-Determinant")
			} else {
				plot(1:length(appOutput$pout$fGBestHist), exp(appOutput$pout$fGBestHist)*100, 
					type = "l", xlab = "Iteration", ylab = "D-eff. (%)")
			}
		}
	})

	output$outcome <- renderPlot({
		if ((appInput$designInfo[1] > 0) & (!is.null(appOutput$dout))) {
			dPara <- appInput$designInfo[4]
			dUB <- appInput$designInfo[5:6]
			dLB <- appInput$designInfo[7:8]
			pUB <- appInput$designInfo[9:(9 + dPara - 1)]
			pLB <- appInput$designInfo[(9 + dPara):length(appInput$designInfo)]			
			
			lay <- layout(matrix(1:2, 1, 2), widths = c(13, 2))
			
			cL <- 25
			colormap <- rgb(c(rep(0,1*cL), seq(0,0.95,length=1*cL), 1,      #R
												rep(0.95,1*cL), seq(0.95,0.5,length=1*cL)),   
											c(rep(0,1*cL), seq(0,0.9,length=1*cL), 1,				#G
												seq(0.9,0,length=1*cL), rep(0,1*cL)), 
											c(seq(0.5,0.95,length = 1*cL), rep(0.95,1*cL), 1, #B
												seq(0.95,0,length = 1*cL), rep(0,1*cL)), 
											maxColorValue = 1)
										
			cM <- ifelse(max(abs(appOutput$dout$dispVal)) == 0, 1, max(abs(appOutput$dout$dispVal)))
			
			par(mar = c(5, 4, 1, 1) + 0.1, las = 1)
			image(appOutput$dout$sGrid, appOutput$dout$iGrid, appOutput$dout$dispVal, 
				col = colormap, zlim = c(-1,1)*cM, axes = FALSE, xlab = "s", ylab = "i")
			axis(1, at = seq(dLB[1], dUB[1], length = 5), labels = seq(dLB[1], dUB[1], length = 5),
					 lty = 0, cex.axis = 1, line = -.5)
			axis(2, at = seq(dLB[2], dUB[2], length = 5), labels = seq(dLB[2], dUB[2], length = 5),
					 lty = 0, cex.axis = 1, line = -.5)
			par(mar = c(5, 0, 1, 2) + 0.1, las = 1)
			image(1, seq(-cM, cM, length = 100), t(as.matrix(seq(-cM, cM, length = 100))), 
						col = colormap, axes = F, xlab = "", ylab = "", zlim = c(-cM, cM))
			mtext(round(cM, 1), 3, cex = .8); mtext(-round(cM, 1), 1, cex = .8); mtext(0, 4, cex = .8)
		}
	})
	
})
