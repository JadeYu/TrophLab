shinyUI(fluidPage(
  titlePanel("TrophLab: create your food web"),
  
  sidebarLayout(
    sidebarPanel(
     selectInput("mac_attr", 
        label = "Macroecological attributes",
        choices = c("Number of species"="nsp","Metabolic distribution type"="tdistr","Local resource"="R","Immigration sequence"="immig","--Please select--"="wait"),
        selected = "wait"),
      
        
      conditionalPanel(
      condition = "input.mac_attr == 'nsp'",
      selectInput("nsp", 
        label = "Please select the number of species in the pool", 
        choices = c(30,50,100),
        selected = 30)
        ),
        
      
       conditionalPanel(
      condition = "input.mac_attr == 'tdistr'",
      selectInput("tdistr",
      	label="Metabolic distribution among species in the pool",
      	choices = c("Normal"="poisson","Geometric"="geometric","Uniform"="uniform"),
        selected = "uniform")
      ),
      
      conditionalPanel(
      condition = "input.mac_attr == 'immig'",
      selectInput("immig",
      	label="Immigration sequence from the species pool to the local community",
      	choices = c("Random"="random","Bigger species first"="descend","Smaller species first"="ascend"),
        selected = "random")
      ),
      
       conditionalPanel(
      condition = "input.mac_attr == 'R'",
      sliderInput("R",
      	label="Total resource in the local community",
      	min=100,max=10000,value=5000
      )
      ),  
      
      selectInput("Dr_attr", 
        label = "Relative individual distinguishability",
        choices = c("Mean","Dispersion","--Please select--"="wait"),
        selected = "wait"),
                   
      conditionalPanel(
      condition = "input.Dr_attr == 'Mean'",
      sliderInput("Dr", 
        label = "Mean value for relative individual distinguishability",
        min=0.1,max=0.9,value=0.5)
      ),
      
      conditionalPanel(
      condition = "input.Dr_attr == 'Dispersion'",
      sliderInput("disperse.Dr", 
        label = "Dispersion for relative individual distinguishability",
        min=0.001,max=1,value=0.01)
      ),

	  selectInput("gamma_attr", 
        label = "Generalization cost",
        choices = c("Mean","Dispersion","--Please select--"="wait"),
        selected = "wait"),
                   
      conditionalPanel(
      condition = "input.gamma_attr == 'Mean'",
      sliderInput("gamma", 
        label = "Mean value for generalization cost",
        min=0.001,max=1,value=0.1)
      ),
      
      conditionalPanel(
      condition = "input.gamma_attr == 'Dispersion'",
      sliderInput("disperse.gamma", 
        label = "Dispersion for generalization cost",
        min=0.001,max=1,value=0.01)
      ),
      
      sliderInput("tau_u", 
        label = "Trophic efficiency (food transferred to biomass)",
        min=0.01,max=0.5,value=0.15),
        
      actionButton("do", "Get my web"),
      
      helpText("Codes can be found on GitHub in:"),
      a(href="https://github.com/JadeYu/TrophLab.git","https://github.com/JadeYu/TrophLab.git")
        ),
    mainPanel(
      conditionalPanel(
      condition = "input.do == 0",
      imageOutput("wait")
      ),
      conditionalPanel(
      condition = "input.do > 0",
     fluidRow(
          column(5,plotOutput("foodweb",width=500,height=450)),
     column(3,
     br(),
     br(),
     tableOutput("table"),
     helpText("*L/S^2 is calculated excluding ungrazed plants."),
     helpText("**Mean metabolic rate = 1 in the species pool."),offset=3)
     
     ),
     
    fluidRow(
    column(4,plotOutput("dynamics",width=350,height=300)),
     column(4,plotOutput("dlink",width=350,height=300),offset=2)
     ),

      h5(helpText("Simulate trophic links among species based on optimal competition outcome derived from maximizing resource allocation microstates (Zhang and Harte, 2015):"),a(href="http://dx.doi.org/10.1016/j.tpb.2015.07.003","http://dx.doi.org/10.1016/j.tpb.2015.07.003"))
    )
    )
  )
))
