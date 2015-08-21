source("FWA_simu.R")

input <- list(nsp=30,tdistr="uniform",immig="random",tau_u=0.2,Dr=0.5,gamma=0.2,R=1000,disperse.Dr=0.1,disperse.gamma=0.1)

shinyServer(function(input, output) {
  output$wait <- renderImage({
  	    filename = paste("wait",sample(1:6,1),".png",sep="")
  	    
        list(src = filename,
         contentType = 'image/png',
         width = 650,
         height = 450,
         alt = "This category hasn't been done yet.")
  }, deleteFile = F)
  
  observeEvent(input$do,{
     scenario = one_food_web(input$nsp,mu=1,input$tdistr,input$immig,e=1,input$tau_u,input$Dr,input$gamma,input$R,disperse=c(input$disperse.Dr,input$disperse.gamma))
     
     output$foodweb <- renderPlot({
     		show_food_web(scenario$FW,zoom=1.5,color="heat")
     })
     
     output$dynamics <- renderPlot({
     		community_dynamics(scenario$invade,scenario$extinct,zoom=0.8)
     })
     
     output$dlink <- renderPlot({
     		link_distr(scenario$FW,zoom=0.8)
     })
     
    output$table <- renderTable({
    	link_mat <- scenario$FW$link_mat
    	sp_mat <- scenario$FW$sp_mat
    	table <- matrix(nrow=4,ncol=1)
    	rownames(table) <- c("#species","#trophic links","L/S^2*","mean metabolic rate**")
    	colnames(table) <- c("Value")
    	S0 = sum(sp_mat$tl>0)
    		S = sum(sp_mat$tl>=2)+sum(sp_mat$tl==1&rowSums(link_mat)>0)
    	L = sum(link_mat[1:dim(link_mat)[2],]>0)
    		B = round(mean(sp_mat$theta[sp_mat$tl>0]),3)
    	table[,1] <- c(S0,L,round(L/S^2,3),B)
    	table
    })
     })

  }
)
