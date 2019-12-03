# ----------------------------------------------------------------------------------------------------
# Markov Chains
# Multiple Diagnosis concerning the transition matrix of a Markov chain
# ----------------------------------------------------------------------------------------------------
# written
# Gabriel LEMYRE
# ----------------------------------------------------------------------------------------------------
# Under the supervision of :
# Maciej AUGUSTYNIAK
# ----------------------------------------------------------------------------------------------------
# First version : November 13th, 2019
# Last version  : November 14th, 2019
# ----------------------------------------------------------------------------------------------------

#load numDeriv package for computing numerical derivatives
#enables the use of functions grad and hessian
# install.packages('numDeriv')
library(numDeriv)
#load Rsolnp package for numerical optimization with solnp function
# install.packages('Rsolnp')
library(Rsolnp)
# install.packages('rootSolve')
library(rootSolve)
# install.packages('data.table')
library(data.table) # Permet d'utiliser rbindlist
# install.packages('matrixcalc') # To test weither matrix is positive definite
library(matrixcalc) # Permet d'utiliser rbindlist
# install.packages('magic') #combiner matrice par blocs
library(magic)

# install.packages('gridExtra') # add tables in plots
library(gridExtra)
# install.packages('grid') # add tables in plots
library(grid)
# install.packages('plotrix') # add tables in plots
library(plotrix)
# install.packages('igraph') # basics of network analysis and visualization
library(igraph)


# -------------------------------------------------------
# 
# 
# -------------------------------------------------------
# Sourcing private libraries
# set Sourcing directory
# -------------------------------------------------------

path <- '~/Documents/GitHub/R-Tools'
setwd(path.expand(path)) # Setting Sourcing path
source("LaTeXTable.R")
source("RunTime.R")
source("BasicFunctions.R")

# -------------------------------------------------------
# 
# 
# -------------------------------------------------------
# set working directory
# -------------------------------------------------------

path <- '~/Documents/GitHub/HHMM'
setwd(path.expand(path)) # Setting path

# -------------------------------------------------------
#
# -------------------------------------------------------

Diagnosis <- function(gamma, additional.info=NULL, plotName=NULL){
	custom.Colors.vec <- matrix(c("firebrick3","blue4",NA,NA,
								  "firebrick3","deepskyblue1","blue4",NA,
								  "firebrick3","darksalmon","deepskyblue1","blue4"),ncol=4,byrow=T)
	
	t <- 1 # Initial entry row number
	n <- dim(gamma)[1]
	# print(paste("Gamma dimension :",n,sep=" "))
	
	temp.eigen <- eigen(gamma)
	print(temp.eigen)
	
	eigen.values <- temp.eigen$values
	eigen.vectors.D <- temp.eigen$vectors
	eigen.vectors.G <- solve(eigen.vectors.D)
	
	# print(eigen.vectors.D%*%diag(eigen.values)%*%eigen.vectors.G)
	# print(eigen.vectors.D%*%eigen.vectors.G)
	# print(eigen.vectors.D[,1]*eigen.vectors.G[1,])
	# print(gamma%*%eigen.vectors.D)
	# print(eigen.vectors.D%*%diag(eigen.values))
	
	# TESTING FOR PI
	# print(eigen.vectors.D %*% diag(c(1,rep(0,n-1))) %*% eigen.vectors.G)
	
	# STATIONNARY DISTRIBUTION COMPUTATION
	vec.unitaire <- rep(1,n)
	dist.asympt <- t(vec.unitaire)%*%solve(diag(n)-gamma+vec.unitaire%*%t(vec.unitaire))
	print(dist.asympt)
	
	
	# -------------------------------------------------------
	# GRAPHICAL PRESENTATION
	# Creating the output file and arranging the graphs
	# -------------------------------------------------------
	
	# Graphical options
	title.size <- 6
	cex.lab <- 7.5
	image.width <- 5000
	image.heigth <- 1000*5/3
	
	# Construction du nom du fichier exporté
	if (is.null(plotName)){
		name <- paste("test",Sys.time(),sep="_")
	}else{
		name <- PlotName
	}
	
	jpeg(paste(path,'/2 - Graphiques/',name,".jpg",sep=""), width = image.width, height = image.heigth)
	
	layout(mat = matrix(c(seq(from=1,to=7,by=1),seq(from=15,to=21,by=1),12,10,8,11,9,13,14),ncol=7,nrow = 3,byrow=T),
		   widths = c(5, 0.5, 5, 0.5, 5, 0.5, 5),
		   height = matrix(c(5*rep(1,7),5/3*rep(1,7),5/3*rep(1,7)),ncol=7,nrow = 3,byrow=T))
	
	# GAMMA
	par(mar=c(5,5,5,5))
	color2D.matplot(gamma,cellcolors="#ffffff",
					main="Matrice de transition",
					cex.main=title.size,
					show.values=3,
					vcol=rgb(0,0,0),
					axes=FALSE,
					vcex=cex.lab, xlab="", ylab="")
	
	# EQUAL SIGN
	par(mar = c(0,0,0,0))
	
		# ann - Display Annotoations (set to FALSE)
		# bty - Border Type (none)
		# type - Plot Type (one that produces no points or lines)
		# xaxt - x axis type (none)
		# yaxt - y axis type (none)
		plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
		text(x = 0.5, y = 0.5, paste("="), 
			 cex = 15, col = "black")
		
		
	# RIGHT EIGEN VECTORS
	par(mar=c(5,5,5,5))
		
	cellcol<-matrix(rep("#ffffff",n^2),nrow=n)
	cellcol[,1]<-rep("#c5e1a5",n)
	
	color2D.matplot(eigen.vectors.D,cellcolors=cellcol,
					main="Vecteurs propres droits",
					cex.main=title.size,
					show.values=3,
					vcol=rgb(0,0,0),
					axes=FALSE,
					vcex=cex.lab, xlab="", ylab="")
	
	# MULTIPLY SIGN
	par(mar = c(0,0,0,0))
	
	# ann - Display Annotoations (set to FALSE)
	# bty - Border Type (none)
	# type - Plot Type (one that produces no points or lines)
	# xaxt - x axis type (none)
	# yaxt - y axis type (none)
	plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
	text(x = 0.5, y = 0.5, "X", 
		 cex = 15, col = "black")
	
	
	# EIGEN VALUES ON DIAGONAL
	par(mar=c(5,5,5,5))
	
	cellcol<-matrix(rep("#ffffff",n^2),nrow=n)
	diag(cellcol)<-c("#fff59D","#81d4fa",rep("#ffffff",(n-2)))
	
	color2D.matplot(diag(eigen.values),cellcolors=cellcol,
					main="Valeurs propres",
					cex.main=title.size,
					show.values=3,
					vcol=rgb(0,0,0),
					axes=FALSE,
					vcex=cex.lab, xlab="", ylab="")
	
	# MULTIPLY SIGN
	par(mar = c(0,0,0,0))
	
	# ann - Display Annotoations (set to FALSE)
	# bty - Border Type (none)
	# type - Plot Type (one that produces no points or lines)
	# xaxt - x axis type (none)
	# yaxt - y axis type (none)
	plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
	text(x = 0.5, y = 0.5, "X", 
		 cex = 15, col = "black")
	
	# LEFT EIGEN VECTORS
	par(mar=c(5,5,5,5))
	
	cellcol<-matrix(rep("#ffffff",n^2),nrow=n)
	cellcol[1,]<-rep("#fff59D",n)
	
	color2D.matplot(eigen.vectors.G,
					main="Vecteurs propres gauches",
					cex.main=title.size,
					cellcolors=cellcol,
					show.values=3,
					vcol=rgb(0,0,0),
					axes=FALSE,
					vcex=cex.lab, xlab="", ylab="")
	
	# DISTRIBUTION STATIONNAIRE
	par(mar=c(5,5,5,5))
	
	color2D.matplot(dist.asympt,
					main="Distribution Asymptotique",
					cex.main=title.size,
					cellcolors=rep("#ffffff",n),
					show.values=3,
					vcol=rgb(0,0,0),
					axes=FALSE,
					vcex=cex.lab, xlab="", ylab="")
	
	# TEMPS MOYEN DANS L'ÉTAT
	par(mar=c(5,5,5,5))
	
	color2D.matplot(1/dist.asympt,
					main="E[ Temps de retour i | Départ en i ]",
					cex.main=title.size,
					cellcolors=rep("#ffffff",n),
					show.values=3,
					vcol=rgb(0,0,0),
					axes=FALSE,
					vcex=cex.lab, xlab="", ylab="")
	
	
	# Closing the plot
	dev.off()
	
	
	# layout(mat = matrix(c(1),ncol=5,nrow = 1))
	# 
	# links <- data.frame(matrix(0,nrow=sum(gamma>0),ncol=5))
	# names(links) <- c("Source","Target","weight","edge.color","edge.loop.angle")
	# 
	# # Getting links and their weights
	# for (i in 1:n){
	# 	for (j in 1:n){
	# 		if (gamma[i,j]>0){
	# 			links[t,] <- c(i,j,gamma[i,j], custom.Colors.vec[(n-1),i],
	# 						   if(i==j){-(i-1)*pi/(n/2)}else{0} # Rotation de la boucle de -2pi/n radian
	# 						   )
	# 			t <- t+1
	# 		}
	# 	}
	# }
	# 
	# links[,c(1:3,5)] <- sapply(links[,c(1:3,5)], as.numeric)
	# 
	# # Building the network
	# network <- graph_from_data_frame(d=links, vertices=1:n, directed=T)
	# print(E(network)$edge.color)
	# V(network)$size <- colSums(gamma)/sum(colSums(gamma))*100
	# E(network)$width <- E(network)$weight*10
	# 
	# # ceb <- cluster_edge_betweenness(network) 
	# 
	# # Plotting the network
	# plot(network, 
	# 	 edge.arrow.size=.4,
	# 	 edge.curved=seq(-0.5, 0.5, 
	# 	 				length = ecount(network)), 
	# 	 edge.color=E(network)$edge.color,
	# 	 vertex.color=custom.Colors.vec[(n-1),],
	# 	 vertex.label.color="white",
	# 	 layout=layout.circle, 
	# 	 edge.loop.angle = E(network)$edge.loop.angle
	# 	 )
	
	# if(!is.null(additional.info)){
	# 	add.info.string <- 
	# 	legend(x=-1.5, y=-1.1, add.info.string)
	# }
	
	return(list(eigen.values=eigen.values,
				eigen.vectors.D=eigen.vectors.D,
				eigen.vectors.G=eigen.vectors.G,
				dist.asympt=dist.asympt))
}


# -------------------------------------------------------
#
# -------------------------------------------------------
mu.Test <- c(-0.4882, 0.0862, 0.0170, 0.1967)
sigma.Test <- c(2.9746, 1.3410, 0.9686, 0.5169)

Gamma.Test <- matrix(c(9.498620e-01, 5.013800e-02, 4.657267e-54, 2.447848e-30,
					   8.677874e-34, 9.530261e-01, 4.697389e-02, 4.935129e-11,
					   5.880273e-03, 3.434210e-30, 9.227759e-01, 7.134384e-02,
					   1.556352e-13, 1.348732e-17, 3.212027e-02, 9.678797e-01), nrow=4, byrow=T)

# Gamma.Test <- matrix(c(0.9048, 0.0952, 0.0000000001,
# 					   0.1, 0.8, 0.1,
# 					   0.0000000001, 0.0429, 0.9571),
# 					 byrow=T,nrow=3)

# Gamma.Test <- matrix(c(2, 4, 3,
# 					   -4, -6, -3,
# 					   3, 3, 1),
# 					 byrow=T,nrow=3)

# Gamma.Test=matrix(c(0.4, 0.6,
# 					 0.02, 0.98),
# 					 byrow=T,nrow=2)

test <- Diagnosis(Gamma.Test, additional.info = rbind(mu.Test,sigma.Test))

