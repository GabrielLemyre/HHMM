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

# --------------------------------------------------------
# DIAGNOSTIQUES SUR LA CHAÎNE DE MARKOV
# --------------------------------------------------------
# Fonction pour calculer les moments des v.a. de la chaîne de Markov
#     Calcul de l'espérance, de l'autocovariance et autocorrélation 
#     de chaîne homogène/stationnaire ou non.
# --------------------------------------------------------
Diagnosis <- function(gamma, 
                      additional.info=NULL, 
                      plotName=NULL,
                      special.Weights=NULL,
                      network.print=FALSE,
                      docPath=NULL){
  
  
  custom.Colors.vec <- matrix(c("firebrick3","blue4",NA,NA,
                                "firebrick3","deepskyblue1","blue4",NA,
                                "firebrick3","darksalmon","deepskyblue1","blue4"),ncol=4,byrow=T)
  
  t <- 1 # Initial entry row number
  n <- dim(gamma)[1]
  # print(paste("Gamma dimension :",n,sep=" "))
  
  temp.eigen <- eigen(gamma)
  # print(temp.eigen)
  
  eigen.values <- temp.eigen$values
  eigen.vectors.D <- temp.eigen$vectors
  eigen.vectors.G <- solve(eigen.vectors.D)
  
  # -------------------------------------
  # -------------------------------------
  # -------------------------------------
  # Temps de séjour dans un état TODO
  # -------------------------------------
  # -------------------------------------
  # -------------------------------------
  
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
  # print(dist.asympt)
  
  
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
  
  # Création de l'image en format jpg et spécification de son emplacement
  FilePathName <- paste(docPath,'/2 - Graphiques/',name,".jpg",sep="")
  jpeg(FilePathName, width = image.width, height = image.heigth)
  
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
  
  if (network.print){ # TODO : Add the space to accomodate the network
    
    links <- data.frame(matrix(0,nrow=sum(gamma>0),ncol=5))
    names(links) <- c("Source","Target","weight","edge.color","edge.loop.angle")
    
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
    # print(E(network)$weight)
    
    temp.weights <- gamma-diag(gamma)
    temp.weights.normalized <- temp.weights/rowSums(temp.weights)+diag(gamma)
    
    if (!is.null(special.Weights)){
      print(V(network))
      Vertex.size <- special.Weights/sum(special.Weights)*100
    } else {
      Vertex.size <- colSums(gamma)/sum(colSums(gamma))*100
    }
    
    E(network)$width <- temp.weights.normalized*5
    
    # ceb <- cluster_edge_betweenness(network)
    
    # Plotting the network
    plot(network,
         edge.arrow.size=1,
         edge.curved=seq(0, 0.8,
                         length = ecount(network)),
         edge.color=E(network)$edge.color,
         vertex.color=custom.Colors.vec[(n-1),],
         vertex.label.color="white",
         vertex.label.cex=Vertex.size/10,
         layout=layout.circle,
         edge.loop.angle = E(network)$edge.loop.angle
    )
    
    # if(!is.null(additional.info)){
    # 	add.info.string <-
    # 	legend(x=-1.5, y=-1.1, add.info.string)
    # }
  }
  
  
  # Closing the plot
  dev.off()
  
  cat("Image was saved as ",FilePathName,sep="")
  
  return(list(eigen.values=eigen.values,
              eigen.vectors.D=eigen.vectors.D,
              eigen.vectors.G=eigen.vectors.G,
              dist.asympt=dist.asympt))
}
