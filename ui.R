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
        min=0.1,max=0.9,value=0.7)
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
        min=0.001,max=1,value=0.5)
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
      #actionButton("do", "Get my web"),
      
      helpText("Codes can be found on GitHub in:"),
      a("https://github.com/JadeYu/TrophLab.git")
        ),
    mainPanel(
      plotOutput("foodweb"),
      h5(helpText("Simulate trophic links among species based on optimal competition outcome derived from maximizing resource allocation microstates (Zhang and Harte, 2015):")),a("http://dx.doi.org/10.1016/j.tpb.2015.07.003")
    )
  )
))
