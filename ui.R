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
					numericInput("v",  "a",          1, min = 1e-2, max = NA, width = "180px")
				)
			),
			column(4,
				conditionalPanel("input.inhibition == '1'",
					column(4,
						numericInput("kmu", "Km Upper",  5, min = 1e-2, max = NA, width = "100px"),
						numericInput("kml", "Km Lower",  4, min = 1e-2, max = NA, width = "100px")
					),
					column(4,
						numericInput("kcu", "Kic Upper", 3, min = 1e-2, max = NA, width = "100px"),	
						numericInput("kcl", "Kic Lower", 2, min = 1e-2, max = NA, width = "100px")
					)
				),
				conditionalPanel("input.inhibition == '2'",
					column(4,
						numericInput("kmu", "Km Upper",  5, min = 1e-2, max = NA, width = "100px"),
						numericInput("kml", "Km Lower",  4, min = 1e-2, max = NA, width = "100px")
					),
					column(4,
						numericInput("kcu", "Kic Upper", 3, min = 1e-2, max = NA, width = "100px"),	
						numericInput("kcl", "Kic Lower", 2, min = 1e-2, max = NA, width = "100px")
					)
				),
				conditionalPanel("input.inhibition == '3'",
					column(4,
						numericInput("kmu", "Km Upper",  5, min = 1e-2, max = NA, width = "100px"),
						numericInput("kml", "Km Lower",  4, min = 1e-2, max = NA, width = "100px")
					),
					column(4,
						numericInput("kuu", "Kiu Upper", 5, min = 1e-2, max = NA, width = "100px"),
						numericInput("kul", "Kiu Lower", 4, min = 1e-2, max = NA, width = "100px")
					)
				),	
				conditionalPanel("input.inhibition == '4'",
					column(4,
						numericInput("kmu", "Km Upper",  5, min = 1e-2, max = NA, width = "100px"),
						numericInput("kml", "Km Lower",  4, min = 1e-2, max = NA, width = "100px")
					),
					column(4,
						numericInput("kcu", "Kic Upper", 3, min = 1e-2, max = NA, width = "100px"),	
						numericInput("kcl", "Kic Lower", 2, min = 1e-2, max = NA, width = "100px")
					),				
					column(4,
						numericInput("kuu", "Kiu Upper", 5, min = 1e-2, max = NA, width = "100px"),
						numericInput("kul", "Kiu Lower", 4, min = 1e-2, max = NA, width = "100px")
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
				conditionalPanel("(input.inhibition == '1' & (input.kmu != input.kml || input.kcu != input.kcl)) |
													(input.inhibition == '2' & (input.kmu != input.kml || input.kcu != input.kcl)) |
													(input.inhibition == '3' & (input.kmu != input.kml || input.kuu != input.kul)) |
													(input.inhibition == '4' & (input.kmu != input.kml || input.kcu != input.kcl || input.kuu != input.kul))",
					numericInput("L1ns", "Outer Loop",  64, min = 4, max = 256, step = 1, width = "100px"),
					numericInput("L2ns", "Inner Loop",  32, min = 4, max = 256, step = 1, width = "100px"),
					numericInput("L4ns", "Mu Loop",    128, min = 4, max = 256, step = 1, width = "100px")
				),
				conditionalPanel("(input.inhibition == '1' & input.kmu == input.kml & input.kcu == input.kcl) |
													(input.inhibition == '2' & input.kmu == input.kml & input.kcu == input.kcl) |
													(input.inhibition == '3' & input.kmu == input.kml & input.kuu == input.kul) |
													(input.inhibition == '4' & input.kmu == input.kml & input.kcu == input.kcl & input.kuu == input.kul)",
					numericInput("L3ns", "Local Loop",  64, min = 4, max = 256, step = 1, width = "100px")
				)
			),
			column(3, 
				tags$h5('Max. Iteration'),
				conditionalPanel("(input.inhibition == '1' & (input.kmu != input.kml || input.kcu != input.kcl)) |
													(input.inhibition == '2' & (input.kmu != input.kml || input.kcu != input.kcl)) |
													(input.inhibition == '3' & (input.kmu != input.kml || input.kuu != input.kul)) |
													(input.inhibition == '4' & (input.kmu != input.kml || input.kcu != input.kcl || input.kuu != input.kul))",
					numericInput("L1nt", "Outer Loop", 100, min = 10, max = 500, step = 5, width = "100px"),
					numericInput("L2nt", "Inner Loop",  50, min = 10, max = 500, step = 5, width = "100px"),
					numericInput("L4nt", "Mu Loop",    100, min = 10, max = 500, step = 5, width = "100px")
				),
				conditionalPanel("(input.inhibition == '1' & input.kmu == input.kml & input.kcu == input.kcl) |
													(input.inhibition == '2' & input.kmu == input.kml & input.kcu == input.kcl) |
													(input.inhibition == '3' & input.kmu == input.kml & input.kuu == input.kul) |
													(input.inhibition == '4' & input.kmu == input.kml & input.kcu == input.kcl & input.kuu == input.kul)",
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
				  tags$div(
				    tags$p("To use this app, please see the following instructions."),
				    tags$ol(
				      tags$li("This app will take a few seconds to initiate."),
				      tags$li("Choose the interested inhibition model in 'Select One'."),
				      tags$li("Specify the model parameters.
				              The users can change the value of the first parameter, a, 
                      to verify that it is linear to the D-criterion.
				              For the lower and upper bounds for the remaining parameters, 
				              please find the illustation below."),
				      tags$ul(
				        tags$li("When some values in 'Lower' is smaller than that in 'Upper', 
				                 the app will search the standardized maximin D-optimal design 
				                 among the specified parameter space."),
				        tags$li("When all values in 'Lower' equal to that in 'Upper', 
				                 the app will search the locally D-optimal design at the 
				                 specified parameter point.")
				      ),
				      tags$li("Specify the design space for each component, s and i."),
				      tags$li("Set the PSO options:"),
				      tags$ul(
				        tags$li("For searching the standardized maximin D-optimal design, 
                         there are two PSO procedures.  
				                 First, specify the wanted swarm size and maximal iteration number for the 
				                 'Outer' and 'Inner' loops of design search procedure, the NestedPSO.  
				                 Then, one also need to specify these two settings for the 'Mu Loop', 
				                 which is another PSO procedure for obtaining the assistant design
                         used in the equivalence theorem checking.
				                 <b>Note that, the NestedPSO takes about 40 seconds if the user runs 
                         it with default settings, and, it should be enough to find the optimal 
				                 design.</b>"),
				        tags$li("For searching the locally <i>D</i>-optimal design, the swarm size and maximal iteration number
				                  for the PSO (<strong>'Local Loop'</strong>) precedure is needed."),
				        tags$li("Some miscellaneous options are available for altering, such as the cognitive and the social parameters, 
				                 the velocity clamping constant and the descending mode of inertia weight")
				      ),
				      tags$li("As the settings are done, press <strong>'Execute'</strong> to start the PSO search.") 
				    )
				  )
				),
        tabPanel("Run & Result", 
					tags$h5('Press "Execute" to Start'),
					actionButton("run", "Execute"),
					textOutput("msg"),
					tags$h4('Result'),
					verbatimTextOutput("result")
				), 
        tabPanel("PSO Search Path", 
					#tags$h5('PSO Search Path'),
					plotOutput("gpath", height = 360, width = 360)
				), 
        tabPanel("Dispersion Plot", 
					#tags$h5('Dispersion Plot'),
					plotOutput("outcome", height = 360, width = 360)
				)
      )		
		)
	)
))
