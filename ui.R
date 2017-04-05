library(shiny)

shinyUI(fluidPage(
  
	fluidRow(
		column(12,
			tags$h4('Design Configuration'),
			column(2, 
				selectInput("inhibition", "Inhibition Model", 
										c("Select One" = 0,
											"Competitive" 		= 1,
											"Noncompetitive"	= 2,
											"Uncompetitive" 	= 3,
											"Mixed-type" 			= 4),
										0, FALSE, width = "180px"),	
				conditionalPanel("input.inhibition != '0'",
					numericInput("v",  "Vmax",          1, min = 1e-2, max = NA, width = "180px")
				)
			),
			column(4,
				conditionalPanel("input.inhibition == '1'",
					column(4,
						numericInput("kmu1", "Km Upper",  5, min = 1e-2, max = NA, width = "100px"),
						numericInput("kml1", "Km Lower",  4, min = 1e-2, max = NA, width = "100px")
					),
					column(4,
						numericInput("kcu1", "Kic Upper", 3, min = 1e-2, max = NA, width = "100px"),	
						numericInput("kcl1", "Kic Lower", 2, min = 1e-2, max = NA, width = "100px")
					)
				),
				conditionalPanel("input.inhibition == '2'",
					column(4,
						numericInput("kmu2", "Km Upper",  5, min = 1e-2, max = NA, width = "100px"),
						numericInput("kml2", "Km Lower",  4, min = 1e-2, max = NA, width = "100px")
					),
					column(4,
						numericInput("kcu2", "Kic Upper", 3, min = 1e-2, max = NA, width = "100px"),	
						numericInput("kcl2", "Kic Lower", 2, min = 1e-2, max = NA, width = "100px")
					)
				),
				conditionalPanel("input.inhibition == '3'",
					column(4,
						numericInput("kmu3", "Km Upper",  5, min = 1e-2, max = NA, width = "100px"),
						numericInput("kml3", "Km Lower",  4, min = 1e-2, max = NA, width = "100px")
					),
					column(4,
						numericInput("kuu3", "Kiu Upper", 5, min = 1e-2, max = NA, width = "100px"),
						numericInput("kul3", "Kiu Lower", 4, min = 1e-2, max = NA, width = "100px")
					)
				),	
				conditionalPanel("input.inhibition == '4'",
					column(4,
						numericInput("kmu4", "Km Upper",  5, min = 1e-2, max = NA, width = "100px"),
						numericInput("kml4", "Km Lower",  4, min = 1e-2, max = NA, width = "100px")
					),
					column(4,
						numericInput("kcu4", "Kic Upper", 3, min = 1e-2, max = NA, width = "100px"),	
						numericInput("kcl4", "Kic Lower", 2, min = 1e-2, max = NA, width = "100px")
					),				
					column(4,
						numericInput("kuu4", "Kiu Upper", 5, min = 1e-2, max = NA, width = "100px"),
						numericInput("kul4", "Kiu Lower", 4, min = 1e-2, max = NA, width = "100px")
					)
				)
			),
			column(3, 
				column(6, 
					conditionalPanel("input.inhibition != '0'",
						numericInput("su", "[s] Upper", 30, min = 0, max = NA, width = "100px"),
						numericInput("sl", "[s] Lower",  0, min = 0, max = NA, width = "100px")
					)
				),
				column(6, 
					conditionalPanel("input.inhibition != '0'",
						numericInput("iu", "[i] Upper", 60, min = 0, max = NA, width = "100px"),
						numericInput("il", "[i] Lower",  0, min = 0, max = NA, width = "100px")
					)
				)
			)
		)
	),
	fluidRow(
		column(6, 
			tags$h4('PSO Configuration'),
			column(3, 
				tags$h5('Swarm Size'),
				conditionalPanel("(input.inhibition == '1' & (input.kmu1 != input.kml1 | input.kcu1 != input.kcl1)) |
													(input.inhibition == '2' & (input.kmu2 != input.kml2 | input.kcu2 != input.kcl2)) |
													(input.inhibition == '3' & (input.kmu3 != input.kml3 | input.kuu3 != input.kul3)) |
													(input.inhibition == '4' & (input.kmu4 != input.kml4 | input.kcu4 != input.kcl4 | input.kuu4 != input.kul4))",
					numericInput("L1ns", "Outer Loop",  64, min = 4, max = 256, step = 1, width = "100px"),
					numericInput("L2ns", "Inner Loop",  32, min = 4, max = 256, step = 1, width = "100px")
					#numericInput("L4ns", "Mu Loop",    128, min = 4, max = 256, step = 1, width = "100px")
				),
				conditionalPanel("(input.inhibition == '1' & input.kmu1 == input.kml1 & input.kcu1 == input.kcl1) |
													(input.inhibition == '2' & input.kmu2 == input.kml2 & input.kcu2 == input.kcl2) |
													(input.inhibition == '3' & input.kmu3 == input.kml3 & input.kuu3 == input.kul3) |
													(input.inhibition == '4' & input.kmu4 == input.kml4 & input.kcu4 == input.kcl4 & input.kuu4 == input.kul4)",
					numericInput("L3ns", "Local Loop",  64, min = 4, max = 256, step = 1, width = "100px")
				)
			),
			column(3, 
				tags$h5('Max. Iteration'),
				conditionalPanel("(input.inhibition == '1' & (input.kmu1 != input.kml1 | input.kcu1 != input.kcl1)) |
													(input.inhibition == '2' & (input.kmu2 != input.kml2 | input.kcu2 != input.kcl2)) |
				                  (input.inhibition == '3' & (input.kmu3 != input.kml3 | input.kuu3 != input.kul3)) |
				                  (input.inhibition == '4' & (input.kmu4 != input.kml4 | input.kcu4 != input.kcl4 | input.kuu4 != input.kul4))",
					numericInput("L1nt", "Outer Loop", 100, min = 10, max = 500, step = 5, width = "100px"),
					numericInput("L2nt", "Inner Loop",  50, min = 10, max = 500, step = 5, width = "100px")
					#numericInput("L4nt", "Mu Loop",    100, min = 10, max = 500, step = 5, width = "100px")
				),
				conditionalPanel("(input.inhibition == '1' & input.kmu1 == input.kml1 & input.kcu1 == input.kcl1) |
													(input.inhibition == '2' & input.kmu2 == input.kml2 & input.kcu2 == input.kcl2) |
													(input.inhibition == '3' & input.kmu3 == input.kml3 & input.kuu3 == input.kul3) |
													(input.inhibition == '4' & input.kmu4 == input.kml4 & input.kcu4 == input.kcl4 & input.kuu4 == input.kul4)",
					numericInput("L3nt", "Local Loop", 100, min = 10, max = 500, step = 5, width = "100px")
				)				
			),
			column(6, 
				tags$h5('Miscellaneous'),
				column(6,
					conditionalPanel("input.inhibition != '0'",
						numericInput("c1", "c1",     2, min = 1, max = 5, width = "100px"),
						numericInput("c2", "c2",     2, min = 1, max = 5, width = "100px"),
						numericInput("vk", "Vel. Clamp",  2, min = 2, max = 5, width = "100px")
					)
				),
				column(6, 
					conditionalPanel("input.inhibition != '0'",
						numericInput("w0", "IW(start)",  0.95, min = 0.5, max = 2, width = "100px"),
						numericInput("w1", "IW(end)",  0.20, min = 0.5, max = 1, width = "100px"),
						numericInput("wv", "IW(varfor)",  0.80, min = 0.1, max = 1, width = "100px")
					)				
				)
				

			)
		),		
		column(6, 
     tabsetPanel(
				tabPanel("Instruction",
					#tags$p('To use this app, please see the following instructions.'), 
					tags$ol(
						tags$li('This app will take a few seconds to initiate NestedPSO kernel via Rcpp.'),
						tags$li('Choose the inhibition model in ', tags$strong('Select One'), '.'),
						tags$li('Specify the model parameters:'),
						tags$ul(
							tags$li(
								'The first parameter, ', tags$strong('Vmax'),
							  'is linear to the ', tags$i('D'), '-criterion.  One can verify this by observing results from different',
							  tags$strong('Vmax'), 'values.'),
							tags$li(
								'When all values in ', tags$strong('Lower'), ' equal to that in ', tags$strong('Upper'), 
								', the PSO search for locally ', tags$i('D'), 
								'-optimal design will be performed.'),
							tags$li(
								'When at leat one value in ', tags$strong('Lower'), ' is smaller than that in ', tags$strong('Upper'), 
								', the NestedPSO search for standardized maximin ', tags$i('D'), 
								'-optimal design will be performed.')
						),
						tags$li('Specify the design space for each component ', tags$strong('[s]'), ' and ', tags$strong('[i]'), '.'),
						tags$li('Specify the PSO parameters:'),
						tags$ul(
							tags$li(
								'For finding locally ', tags$i('D'), '-optimal design, ', 
								tags$strong('swarm size'), ' and ', tags$strong('maximal number of iterations'), ' for the PSO (', 
								tags$strong('Local Loop'), ') precedure are needed.'),
							tags$li(
								'For finding the standardized maximin ', tags$i('D'), 
								'-optimal design, the NestedPSO consists of two nested PSO loop, ', 
								tags$strong('Outer'), ' and ', tags$strong('Inner'), 
								'loops.  For each loop, please specify ', tags$strong('swarm size'), ' and ', tags$strong('maximal number of iterations.'),
								tags$strong('Note that, by default seeting, the NestedPSO takes about 40 seconds and 
								it should be enough to find the optimal design.')),
							tags$li(
								'Some miscellaneous options are available for altering, such as the cognitive and the social parameters, 
								the velocity clamping constant and the descending mode of inertia weight.')
						),						
						tags$li('Click on ', tags$strong('Execute'), ' button to start the NestedPSO search.')
					)
				),
        tabPanel("Run & Result", 
					tags$h5('Press "Execute" to Start'),
					actionButton("run", "Execute"),
					textOutput("msg"),
					tags$h4('Result'),
					verbatimTextOutput("result")
				), 
        tabPanel("Directional Derivative Plot", 
					#tags$h5('Dispersion Plot'),
					tags$p('The Directional Derivative Plot is a vital tool to check whether a design is optimal or not.  
								The plot for a (nearly) optimal design should only has non-positive values (blue area) on the design space and
								zero values (white area) occur at the support points (green dots).'),
					imageOutput("outcome", height = 360, width = 360)
				)#,
        #tabPanel("PSO Search Path", 
					#tags$h5('PSO Search Path'),
					#imageOutput("gpath", height = 360, width = 360)
				#) 
      )		
		)
	)
))
