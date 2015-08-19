one_food_web <- function(nsp,mu,tdistr,immig,e,tau_u,Dr,gamma,R,disperse=c(0.1,0.1),color){
	Y <- try(try_one_food_web(nsp,mu,tdistr,immig,e,tau_u,Dr,gamma,R,disperse,color),silent=T)
	while(class(Y)=="try-error"){
		Y <- try(try_one_food_web(nsp,mu,tdistr,immig,e,tau_u,Dr,gamma,R,disperse,color))	
	}
}

try_one_food_web <- function(nsp,mu,tdistr,immig,e,tau_u,Dr,gamma,R,disperse=c(0.1,0.1),color){
	##calculate all species paramters
	theta_seq <- pool.theta(nsp,mu,tdistr)
	sequence <- immig_seq(theta_seq,1:nsp,immig)
	R_mat <- data.frame(abundance=R*mu,theta=1,Dr=0)
	##here multiple resources can also be used
	theta = rescale(theta_seq,c(-3,3))
	Dr_seq = sigmoid(theta,disperse[1])
	Dr_seq = Dr_seq/mean(Dr_seq)*Dr
	nDr = rescale(1/Dr_seq,c(-3,3))
	gamma_seq <- sigmoid(nDr,disperse[2])
	gamma_seq = gamma_seq/mean(gamma_seq)*gamma
	##simulate the community and food web
	community <- generate.pool(theta_seq,Dr_seq,gamma_seq,R_mat) 
	scenario <-FWA_switch_cost(community,e,tau_u,sequence,omni=F)
##omni specifies whether the species can consume the fundamental resource or not
	show_food_web(scenario$FW,tau_u,lab="",zoom=1.5,color)
}

##functions################
rescale <- function(x_seq,span){
	(x_seq-min(x_seq))/(max(x_seq)-min(x_seq))*(span[2]-span[1])+span[1]
}

sigmoid <- function(z,a){
	1/(1+exp(-a*z))
}

immig_seq <- function(x,sequence,immig){
	if(immig=="random"){
		sequence = sample(sequence)
	}else{
		sequence = sequence[order(x,decreasing=(immig=="descend"))]
	}
	sequence
}

pool.theta <- function(nsp,mu,distr){
	if(distr=="poisson"){
		theta_seq <- rpois(nsp,mu)+1
	}else if(distr=="geometric"){
		theta_seq <- rgeom(nsp,1/(mu+1))+1
	}else if(distr=="uniform"){
		theta_seq <- runif(nsp,0,2*mu)
	}else{
		print("no such bs distribution")
	}
	theta_seq
}

generate.pool <- function(theta_seq,Dr_seq,gamma_seq,R_mat){
	nsp <- length(theta_seq)
	nR <- dim(R_mat)[1]
	sp_pool <- matrix(nrow=nsp+nR,ncol=8)
	sp_pool[,1] <- 1:(nsp+nR)
	sp_pool[,2] <- c(theta_seq,R_mat$theta)
	sp_pool[,3] <- c(Dr_seq,rep(R_mat$Dr,nR))
	sp_pool[,4] <- c(gamma_seq,rep(0,nR))
	sp_pool[,5] <- c(rep(-1,nsp),rep(0,nR))
	sp_pool[,6] <- c(rep(0,nsp),R_mat$abundance)
	sp_pool[,7] <- c(rep(0,nsp),R_mat$abundance)
	sp_pool[,8] <- 1
	colnames(sp_pool) <- c("index","theta","Dr","gamma","tl","BPN","APN","tau_c")##The new sp_mat contains information on the species and resources (basal species)
	link_mat <- matrix(0,nrow=nsp+nR,ncol=nsp)
	colnames(link_mat) <- paste("Sp",1:nsp)
	rownames(link_mat) <- c(paste("Sp",1:nsp),paste("R",1:nR))
	list(sp_mat=data.frame(sp_pool),link_mat=link_mat)
}

FWA_switch_cost <- function(community0,e,tau_u,sequence,omni){
	community <- community0
	invade <- extinct <- 0
	for(index in sequence){
		Y <- add_species(community,index,e,tau_u,omni)
		community <- Y$community
		invade <- c(invade,invade[length(invade)]+Y$invade)
		extinct <- c(extinct,extinct[length(extinct)]+Y$extinct)
	}
	community <- flow_pattern(Y$community,tau_u)
	list(FW=community,invade=invade,extinct=extinct)
}

add_species <- function(community,index,e,tau_u,omni){
	invade <- 0
	extinct <- 0
	new_comm <- invasion(community,index,e,tau_u,omni)
	if(sum(new_comm$sp_mat$tl>0)!=sum(community$sp_mat$tl>0)){
		invade <- 1
		new_comm <- extinction(new_comm,e,tau_u)
		extinct <- sum(community$sp_mat$tl>0)+1 - sum(new_comm$sp_mat$tl>0)
	}
	list(community=new_comm,invade=invade,extinct=extinct)
}

invasion <- function(community,index,e,tau_u,omni){
	sp_mat <- community$sp_mat
	link_mat <- community$link_mat
	Y <- try_add(community,index,tau_u,omni)
	if(Y$SSN>e){
		link_mat[Y$added_prey,index] <- 1
		sp_mat$tl[index] <- Y$tl
		sp_mat$tau_c[index] <- Y$tau_c
		sp_mat$BPN[index] <- sp_mat$APN[index] <- Y$SSN
	}
	list(sp_mat=sp_mat,link_mat=link_mat)
}

extinction <- function(community,e,tau_u){
	new_comm <- update(community,e,tau_u)
	c <- 1
	while(sum(new_comm$sp_mat$tl>0)!=sum(community$sp_mat$tl>0)|sum(new_comm$sp_mat$tl>0)!= sum(colSums(new_comm$link_mat)>0)){##extinction cascade
		#print(paste("cascade",c))
		community <- new_comm
		new_comm <- update(community,e,tau_u)
		c <- c+1
	}
	new_comm
}

update <- function(community,e,tau_u){
	community <- abundance_update(community,tau_u)
	sp_mat <- community$sp_mat
	link_mat <- community$link_mat
	for(extinct in sp_mat$index[sp_mat$tl>0&sp_mat$APN<e]){
		sp_mat$tl[extinct] <- -1
		sp_mat$BPN[extinct] <- sp_mat$APN[extinct] <- 0
		sp_mat$tau_c[extinct] <- 1
		link_mat[extinct,] <- 0
		link_mat[,extinct] <- 0
	}
	community <- list(sp_mat=sp_mat,link_mat=link_mat)
	community <- trophic_update(community,tau_u)
	community
}

abundance_update <- function(community,tau_u){
	sp_mat <- community$sp_mat
	link_mat <- community$link_mat
	ntl <- max(sp_mat$tl)
	if(ntl>0){
		for(tl in 1:ntl){##from bottom to top
			for(sp in sp_mat$index[sp_mat$tl==tl]){
				sp_mat$BPN[sp] <- SSN_before_predation(sp,list(sp_mat=sp_mat,link_mat=link_mat),tau_u)
			}
		}
		for(tl in 0:ntl){##notice that these two loops cannot be merged; resource remaining will also be updated
			for(sp in sp_mat$index[sp_mat$tl==tl]){
				sp_mat$APN[sp] <- SSN_after_predation(sp,list(sp_mat=sp_mat,link_mat=link_mat),tau_u)
			}
		}
	}
	list(sp_mat=sp_mat,link_mat=link_mat)
}

trophic_update <- function(community,tau_u){
	sp_mat <- community$sp_mat
	link_mat <- community$link_mat
	ntl <- max(sp_mat$tl)
	if(ntl>0){
		for(tl in 1:ntl){##from bottom to top
			for(sp in sp_mat$index[sp_mat$tl==tl]){
				if(sum(sp_mat$tl[link_mat[,sp]>0]>=0)){
					sp_mat$tl[sp] <- max(sp_mat$tl[link_mat[,sp]>0])+1
				}else{
					sp_mat$tl[sp] <- -1
				}
			}
		}
	}
	list(sp_mat=sp_mat,link_mat=link_mat)
}

try_add <- function(community,index,tau_u,omni){
	sp_mat <- community$sp_mat
	link_mat <- community$link_mat
	if(omni){
		preys <- c(sp_mat$index[sp_mat$tl>=0])
	}else{
		preys <- c(sp_mat$index[sp_mat$tl>0])
	}
	select <- select_prey(community,preys,index,tau_u) 
	if(!omni){
		basal <- select_prey(community,sp_mat$index[sp_mat$tl==0],index,tau_u)
		if(basal$SSN>select$SSN){
			select <- basal
		}
	}
	list(SSN=select$SSN,added_prey=select$added_prey,tau_c=select$tau_c,tl=max(sp_mat$tl[select$added_prey])+1)
	##assuming new species will not be preyed on when they first enter the community, but can be preyed on by species coming in later
}


select_prey <- function(community,preys,index,tau_u){
	sp_mat <- community$sp_mat
	tau_c <- 1
	increase <- SSN <- 0
	added_prey <- c()
	SSNs <- numeric(length(preys))
	if(length(preys)>0){
		for(p in 1:length(preys)){
			SSNs[p] <- try_prey(preys[p],community,index,tau_c*exp(-sp_mat$gamma[preys[p]]*(sp_mat$Dr[preys[p]])),tau_u)
		}
		preys <- preys[order(SSNs,decreasing=T)]##determine priority 
	}
	while(increase>=0&&length(preys)>0){
		increase <- sum(unlist(lapply(c(added_prey,preys[1]),try_prey,tau_c=tau_c*exp(-sp_mat$gamma[preys[1]]*(sp_mat$Dr[preys[1]])),community=community,index=index,tau_u=tau_u))) - SSN
		if(increase>0){
			added_prey <- c(added_prey,preys[1])
			tau_c <- tau_c*exp(-sp_mat$gamma[preys[1]]*(sp_mat$Dr[preys[1]]))
			preys <- preys[-1]##remove the best prey	
			SSN <- SSN+increase
		}
	}
	list(SSN=SSN,tau_c=tau_c,added_prey=added_prey)
}

try_prey <- function(prey,community,index,tau_c,tau_u){
	sp_mat <- community$sp_mat
	link_mat <- community$link_mat
	predators <- sp_mat$index[1:dim(link_mat)[2]][link_mat[prey,]>0]
	predators <- predators[predators!=index]
	if(length(predators)==0){##no other predators
		SSN <- sp_mat$BPN[prey]*sp_mat$theta[prey]/2*tau_c*tau_u/sp_mat$theta[index]
	}else{
		theta_seq <- c(sp_mat$theta[index]/tau_c,(sp_mat$theta/sp_mat$tau_c)[predators])
		Dr_seq <- sp_mat$Dr[c(index,predators)]
		R_prey <- sp_mat$theta[prey]*sp_mat$BPN[prey]/2*tau_u
		SSN <-solve.analytical(theta_seq,Dr_seq,R_prey,1)[1]
	}
	SSN
}

##caution!! self reference in this and try_prey functions
SSN_before_predation <- function(sp,community,tau_u){##calculate SSN before a species is preyed on
	sp_mat= community$sp_mat
	link_mat = community$link_mat
	if(sp_mat$tl[sp]==0){
		SSN <-sp_mat$BPN[sp]
	}else{
		preys <- sp_mat$index[link_mat[,sp]>0]
		SSN <- 0
		if(length(preys)>0){
			SSN <- sum(unlist(lapply(preys,try_prey,community=community,index=sp,tau_c=sp_mat$tau_c[sp],tau_u=tau_u)))
		}
	}
	SSN
}

SSN_after_predation <- function(sp,community,tau_u){##half the prey population and add the uncaptured prey back to its steady state abundance
	sp_mat = community$sp_mat
	link_mat = community$link_mat
	SSN <- sp_mat$BPN[sp]
	predators <- sp_mat$index[1:dim(link_mat)[2]][link_mat[sp,]>0]
	if(length(predators)>0){
		SSN <- 0.5*SSN
		for(i in 1:length(predators)){
			tau_c <- sp_mat$tau_c[predators[i]]
			SSN <- SSN+try_prey(sp,community,predators[i],tau_c,tau_u)*sp_mat$theta[predators[i]]/sp_mat$theta[sp]/tau_u*(1-tau_c)/tau_c##the uncaptured amount from a given predator
		}
	}
	SSN
}

flow_pattern <- function(community,tau_u){
	sp_mat <- community$sp_mat
	link_mat <- community$link_mat
	for(sp in sp_mat$index[sp_mat$tl>0]){
		preys <- sp_mat$index[link_mat[,sp]>0]
		link_mat[preys,sp] <- unlist(lapply(preys,try_prey,community=community,index=sp,tau_c=sp_mat$tau_c[sp],tau_u=tau_u))*sp_mat$theta[sp]/tau_u
	}
	list(sp_mat=sp_mat,link_mat=link_mat)
}

######functions for graphic demonstration
show_food_web <- function(community,tau_u,lab="",zoom,color="gray"){
	sp_mat <- community$sp_mat
	link_mat <- community$link_mat
	ntl <- max(sp_mat$tl)
	plot(1~1,col=0,xlim=c(0,1),ylim=c(0,ntl+1),xaxt='n',ylab="trophic level",xlab=lab,cex=zoom,cex.lab=zoom,cex.axis=zoom)
	energy <- sp_mat$theta*sp_mat$APN
	col <- rep(0,length(energy))
	mu <- mean(sp_mat$theta)
	if(color=="gray"){
		col[sp_mat$tl>0] <- paste("gray",round(rescale(energy[sp_mat$tl>0],c(80,0))),sep="")
	}else if(color=="heat"){
			col[sp_mat$tl>0] <- heat.colors(50)[round(rescale(energy[sp_mat$tl>0],c(50,1)))]
	}else if(color=="topo"){
			col[sp_mat$tl>0] <- topo.colors(50)[round(rescale(energy[sp_mat$tl>0],c(1,50)))]
	}else if(color=="rainbow"){
			col[sp_mat$tl>0] <- rainbow(50,start=0.2)[round(rescale(energy[sp_mat$tl>0],c(1,50)))]
	}
	for(tl in 1:ntl){
		abline(h=tl,lty=2,lwd=zoom)
		sps <- sp_mat$index[sp_mat$tl==tl]
		xseq <- 1:length(sps)/(length(sps)+1)
		sizes <- sp_mat$theta[sps]
		sizes <- sizes/mu*2
		
		points(xseq,rep(tl,length(sps)),pch=16,col=col[sps],cex=sizes*zoom)
	}
	for(sp in sp_mat$index[sp_mat$tl>1]){
		tl <- sp_mat$tl[sp]
		sps <- sp_mat$index[sp_mat$tl==tl]
		x <- (1:length(sps)/(length(sps)+1))[sp==sps]
		for(prey in sp_mat$index[link_mat[1:dim(link_mat)[2],sp]>0]){
			tl_prey <- sp_mat$tl[prey]
			sps_prey <- sp_mat$index[sp_mat$tl==tl_prey]
			x_prey <- (1:length(sps_prey)/(length(sps_prey)+1))[prey==sps_prey]
			flow <- link_mat[prey,sp]/mean(link_mat[1:dim(link_mat)[2],][link_mat[1:dim(link_mat)[2],]>0])
						arrows(x,tl,x_prey,tl_prey,length=0,lwd=flow*zoom,col="gray40")
		}
	  }				
}

community_dynamics <- function(invasion,extinction,zoom){
	time <- 1:length(invasion)
	plot(invasion~time,ylim=c(0,length(invasion)),xlab="time",pch=16,col=3,ylab="#species",cex=zoom,cex.lab=zoom,cex.axis=zoom)
	points(time,extinction,pch=16,col=2,cex=zoom)
	points(time,(invasion-extinction),pch=16,col=1,cex=zoom)
}

graph_output_FWA <- function(scenarios,tau_u,labs,indice,filename,layout,zoom){
	if(!is.na(filename)){
		png(filename,width=300*layout[2],height=300*layout[1])
	par(mfrow=layout)
	}
	for(i in 1:length(indice)){
		show_food_web(scenarios[[indice[i]]]$FW,tau_u,lab=labs[i],zoom)
	}

	for(i in indice){
		community_dynamics(scenarios[[i]]$invade,scenarios[[i]]$extinct,zoom)
	}
	legend("topleft",c("invasion","extinction","existing"),pch=16,col=3:1,cex=zoom)
	if(!is.na(filename)){
		dev.off()
	}
}

##metrics for the FWA process
get_communities_seq <- function(R_seq,theta_seq,Dr,gamma,e,tau_u,immig_seq,omni){
	communities <- list()
	for(i in 1:length(R_seq)){
		R_mat <- data.frame(abundance=R_seq[i],theta=1,Dr=0)
		community <- generate.pool(theta_seq,Dr,R_mat)
		communities[[i]] <- FWA_switch_cost(community,gamma,e,tau_u,immig_seq,omni)
		print(paste("R level",i))
	}
	communities
}

FW_structure <- function(communities,Dr,gamma,e,tau_u){
	S_seq <- L_seq <- B_seq <- F1_seq <- F2_seq <- numeric()
	for(i in 1:length(communities)){
		scenario <- communities[[i]]$FW
		S_seq[i] <- sum(scenario$sp_mat$tl>=2)+sum(scenario$sp_mat$tl==1&rowSums(scenario$link_mat)>0)
		F1_seq[i] <- sum(scenario$sp_mat$tl==1&rowSums(scenario$link_mat)>0)/sum(scenario$sp_mat$tl==1)
		F2_seq[i] <- 1-sum(scenario$sp_mat$APN[scenario$sp_mat$tl==1&rowSums(scenario$link_mat)>0])/sum(scenario$sp_mat$BPN[scenario$sp_mat$tl==1&rowSums(scenario$link_mat)>0])##fraction of consumed plant species that is not eaten
		L_seq[i] <- sum(scenario$link_mat[1:dim(scenario$link_mat)[2],]>0)##only higher levels
		B_seq[i] <- mean(scenario$sp_mat$theta[scenario$sp_mat$tl>0])
	}
	data.frame(S=S_seq,L=L_seq,mbs=B_seq,F1=F1_seq,F2=F2_seq,Dr=rep(Dr,length(communities)),gamma=rep(gamma,length(communities)),tau_u=rep(tau_u,length(communities)))
	
}

compare_metrics <- function(metrics,metric_name1,metric_name2,metric_pch,metric_col,fit=0,zoom,new,frac=F,graph=NA){	
	if(!is.na(graph)){
		png(graph,width=600,height=300)
		par(mfrow=c(1,2))
	}
	pch <- metrics[,match(metric_pch,names(metrics))]
	index_pch <- match(pch,unique(pch))
	col <- metrics[,match(metric_col,names(metrics))]
	index_col <- match(col,unique(col))+1
	metric1 <- metrics[,match(metric_name1,names(metrics))]
	metric2 <- metrics[,match(metric_name2,names(metrics))]
	if(new){
		if(frac){
			plot(metric2~metric1,xlab=paste(metric_name1,"(",metric_pch,"=",unique(pch),")",sep=""),xlim=c(0.4,1),ylim=c(0.25,0.5),ylab=metric_name2,pch=index_pch,col=index_col,cex=zoom,cex.lab=zoom,cex.axis=zoom)
		}else{
			plot(metric2~metric1,xlab=paste(metric_name1,"(",metric_pch,"=",unique(pch),")",sep=""),ylab=metric_name2,pch=index_pch,col=index_col,cex=zoom,cex.lab=zoom,cex.axis=zoom)
		}
	}
	
	if(fit>0){
		C_seq <- c()
		for(i in unique(index_col)){
			y <- metric2[index_col==i]
			x <- metric1[index_col==i]^fit
			Y <- lm(y~x-1)
			C_seq <- c(C_seq,round(summary(Y)$coefficients[1,1],2))
			lines(metric1[index_col==i][order(x)],predict(Y)[order(x)],col=i,lwd=zoom)
			abline(a=-1,b=1,lty=2)
			lines(metric1[order(metric1)],0.1*metric1[order(metric1)]^2,lty=3)
			
		}
		#plot.new()
		
		legend("topleft",paste(metric_name2,"=",C_seq,metric_name1,"^",fit,sep=""),col=unique(index_col),lty=1,lwd=zoom,cex=zoom)
	}
	if(!is.na(graph)){
		dev.off()
	}
}

show_LSR <- function(metrics,metric_graph,metric_col,graph=NA){
	zoom=1
	if(!is.na(graph)){
		png(graph,width=600,height=600)
		par(mfrow=c(2,2))
		zoom =1.5
	}
	gr <- metrics[,match(metric_graph,names(metrics))]
	
	for(i in unique(gr)){
		compare_metrics(metrics[gr==i,],"S","L",metric_graph,metric_col,fit=2,zoom=zoom,new=T)
		#compare_metrics(metrics[gr==i,],"S","L",metric_graph,metric_col,fit=1,zoom=zoom,new=F)
	}
	plot.new()
	
	col <- metrics[,match(metric_col,names(metrics))]
	index_col <- match(col,unique(col))+1
	legend("bottomleft",c("L=S-1","L=0.1S^2",paste(metric_col,"=",unique(col))),col=c(1,1,unique(index_col)),lty=c(2,3,rep(1,length(unique(index_col)))),cex=zoom)
	if(!is.na(graph)){
		dev.off()
	}
}

show_fraction <- function(metrics,metric_graph,metric_col,graph=NA){
	zoom=1
	if(!is.na(graph)){
		png(graph,width=600,height=600)
		par(mfrow=c(2,2))
		zoom =1.5
	}
	gr <- metrics[,match(metric_graph,names(metrics))]
	for(i in unique(gr)){
		compare_metrics(metrics[gr==i,],"F1","F2",metric_graph,metric_col,fit=0,zoom=zoom,frac=T,new=T)
	}
	plot.new()
	
	col <- metrics[,match(metric_col,names(metrics))]
	index_col <- match(col,unique(col))+1
	legend("bottomleft",c(paste(metric_col,"=",unique(col))),col=c(unique(index_col)),pch=1,cex=zoom)
	if(!is.na(graph)){
		dev.off()
	}
}



##functions for MERA SSN#####################
library(nleqslv)

solve.analytical <- function(theta_seq,Dr_seq,R,r){##given r solve for N+g; r=1 is the special case for equilibrium where N+g=Ne
	if(length(theta_seq)==1){
		sol <- R/theta_seq[1]
	}
	r_seq <- c(r,rep(1,length(theta_seq)-length(r)))
	C <- get.C(theta_seq,Dr_seq,R,r_seq)
	sol <- analytical.sol(C,theta_seq,Dr_seq,r_seq)
	sol
}

analytical.sol <- function(C,theta_seq,Dr_seq,r_seq){
	2/exp(1)*(C*theta_seq*r_seq^(0.5/theta_seq))^(1/(Dr_seq-1))
}

get.C <- function(theta_seq,Dr_seq,R,r_seq){
	init.C <- 1/mean(theta_seq)
	C <- nleqslv(init.C,R_constraint,theta_seq=theta_seq,Dr_seq=Dr_seq,R=R,r_seq=r_seq)
	if(C[[3]]!=1&&C[[3]]!=2){
		print(C[[4]])
	}
	C[[1]]
}

R_constraint <- function(C,theta_seq,Dr_seq,R,r_seq){
	R-sum(theta_seq*analytical.sol(C,theta_seq,Dr_seq,r_seq))
}

steady.state.mR <- function(Sp_mat,R_mat,combine=T){
	SSN_mat <- matrix(ncol=dim(R_mat)[1],nrow=dim(Sp_mat)[1])
	for(i in 1:dim(R_mat)[1]){
		SSN_mat[,i] <- solve.analytical(Sp_mat[,1]/R_mat[i,2],Sp_mat[,2],R_mat[i,1],1)
	}
	rownames(SSN_mat) <- paste("Sp",1:dim(SSN_mat)[1])
	colnames(SSN_mat) <- paste("R",1:dim(SSN_mat)[2])
	if(combine){
		SSN_mat <- rowSums(SSN_mat)
	}
	SSN_mat
}



