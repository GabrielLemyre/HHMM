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

Diagnosis <- function(gamma, additional.info=NULL){
	custom.Colors.vec <- matrix(c("firebrick3","blue4",NA,NA,
								  "firebrick3","deepskyblue1","blue4",NA,
								  "firebrick3","darksalmon","deepskyblue1","blue4"),ncol=4,byrow=T)
	
	links <- data.frame(matrix(0,nrow=sum(gamma>0),ncol=5))
	names(links) <- c("Source","Target","weight","edge.color","edge.loop.angle")
	
	t <- 1 # Initial entry row number
	n <- dim(gamma)[1]
	# print(paste("Gamma dimension :",n,sep=" "))
	
	# Getting links and their weights
	for (i in 1:n){
		for (j in 1:n){
			if (gamma[i,j]>0){
				links[t,] <- c(i,j,gamma[i,j], custom.Colors.vec[(n-1),i],
							   if(i==j){-(i-1)*pi/(n/2)}else{0} # Rotation de la boucle de -2pi/n radian
							   )
				t <- t+1
			}
		}
	}
	
	links[,c(1:3,5)] <- sapply(links[,c(1:3,5)], as.numeric)
	
	# Building the network
	network <- graph_from_data_frame(d=links, vertices=1:n, directed=T)
	print(E(network)$edge.color)
	V(network)$size <- colSums(gamma)/sum(colSums(gamma))*100
	E(network)$width <- E(network)$weight*10
	
	# ceb <- cluster_edge_betweenness(network) 
	
	# Plotting the network
	plot(network, 
		 edge.arrow.size=.4,
		 edge.curved=seq(-0.5, 0.5, 
		 				length = ecount(network)), 
		 edge.color=E(network)$edge.color,
		 vertex.color=custom.Colors.vec[(n-1),],
		 vertex.label.color="white",
		 layout=layout.circle, 
		 edge.loop.angle = E(network)$edge.loop.angle
		 )
	
	# if(!is.null(additional.info)){
	# 	add.info.string <- 
	# 	legend(x=-1.5, y=-1.1, add.info.string)
	# }
	
	return(net2)
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

Gamma.Test <- matrix(c(0.9048, 0.0952, 0.0000000001,
					 0.4122, 0.5485, 0.0393,
					 0.0000000001, 0.0429, 0.9571),
					 byrow=T,nrow=3)

# Gamma.Test=matrix(c(0.98, 0.02,
# 					 0.02, 0.98),
# 					 byrow=T,nrow=2)

test <- Diagnosis(Gamma.Test, additional.info = rbind(mu.Test,sigma.Test))

