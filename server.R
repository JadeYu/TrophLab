source("FWA_simu.R")

shinyServer(function(input, output) {
	 output$default <- renderImage({
      	
        list(src = filename,
         contentType = 'image/png',
         width = 600,
         height = 350,
         alt = "Image is not ready.")
  }, deleteFile = F)
  
     output$foodweb <- renderPlot({
     		one_food_web(input$nsp,mu=1,input$tdistr,input$immig,e=1,input$tau_u,input$Dr,input$gamma,input$R,disperse=c(input$disperse.Dr,input$disperse.gamma),color="rainbow")
     })

  }
)
