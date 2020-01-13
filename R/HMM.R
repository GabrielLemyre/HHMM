# ----------------------------------------------------------------------------------------------------
# Hidden Markov Models
# Basic architecture
# ----------------------------------------------------------------------------------------------------
# written
# Gabriel LEMYRE
# ----------------------------------------------------------------------------------------------------
# Under the supervision of :
# Maciej AUGUSTYNIAK
# ----------------------------------------------------------------------------------------------------
# Last version : april 15th, 2019
# Last version : january 13th, 2020
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
# install.packages('timeDate') # add tables in plots
library(timeDate)


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



# ------------------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////
# TRAINING THE MODEL
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# ------------------------------------------------------------------------------------

# --------------------------------------------------------
# FULL HMM EXPERIENCE : TRAINING FOR THE PARAMETERS
# --------------------------------------------------------
HMM.Train = function(index,
                     data,
                     start.date=NULL,
                     end.date=NULL,
                     frequency=5,
                     mult=1,
                     mu0=NULL,
                     sigma0=NULL,
                     Gamma0=NULL,
                     distribution="Normal",
                     initial.Distribution=NULL,
                     Transition.Type="Homogeneous",
                     nbStepsBack=1,
                     nu=4,
                     nbRegime=NULL,
                     type="HMM",
                     auto.assign=FALSE,
                     data.Origin="MATLAB",
                     PlotName = NULL,
                     special.suffix.name=NULL,
                     path="",
                     nbTicks=20,
                     FormatData=FALSE){

    if (type=="HHMM"){
        nbRegime <- 4
    } else if (is.null(nbRegime)){
        nbRegime <- dim(Gamma0)[1]
    }

    if (Transition.Type!="Homogeneous"){
        nbRegime <- 2
    }

    
    if (FormatData){
        # --------------------------------------------
        # Data manipulation
        # --------------------------------------------
        n <- length(data[,1])
        # Date bounds for training
        (if (is.null(start.date)){start.date <- data[1,1]})
        (if (is.null(start.date)){end.date <- data[n,1]})
    
        if (data.Origin=="MATLAB"){
            # converting the dates from MATLAB's format to R's
            dates <- as.Date(as.POSIXct((data[,1] - 719529)*86400, origin = "1970-01-01", format="%Y-%m-%d"))
            print("MATLAB")
        }else{
            dates <- as.Date(as.character(data[,1]))
        }
    
        st <- which(dates>=start.date)[1]
    
        en <- which(dates<=end.date)
        en <- en[length(en)]
    
        dates.FULL <- dates[st:en]
        sequence.index <- seq(from=frequency, to=length(dates.FULL),by=frequency)
        dates <- dates.FULL[sequence.index]
    
        if (data.Origin=="MATLAB"){
            Pt.FULL <- as.numeric(as.character(data[st:en,index]))
        }else{
            Pt.FULL <- as.numeric(as.character(data[st:en,2]))
        }
    
    
        Pt <- Pt.FULL[sequence.index]
        names(Pt) <- dates
        Pt <- Pt[!is.na(Pt)]
        dates <- names(Pt)
    
        logR <- log( Pt[-1] / Pt[-length(Pt)] )
        logR <- logR * mult
        dates <- dates[-1]
        #the length of logR is one less than Pt
    
        #convert dates into a numeric vector (useful for figures)
        dates.Char <- as.POSIXlt(c(names(logR)), format="%Y-%m-%d")
        n.days <- rep(365, length(dates.Char))
        n.days[((1900+dates.Char$year) %% 4)==0] <- 366 #leap years
        dates.Char  <- 1900 + dates.Char$year + (dates.Char$yday+1)/n.days
    } else {
        logR <- data;
        dates.FULL <- as.Date(names(logR));
        dates.Char <- as.POSIXlt(c(names(logR)), format="%Y-%m-%d")
        n.days <- rep(365, length(dates.Char))
        n.days[((1900+dates.Char$year) %% 4)==0] <- 366 #leap years
        dates.Char  <- 1900 + dates.Char$year + (dates.Char$yday+1)/n.days
    }
    
    # --------------------------------------------
    # NAME CONSTRUCTION FOR PLOT EXPORT
    # --------------------------------------------
    if (is.null(PlotName)){
        dates.Range <- as.character(c(dates.FULL[1],dates.FULL[length(dates.FULL)]))
        dates.split <- matrix(unlist(strsplit(dates.Range,split="\\-")),ncol=3,byrow=T)
        dates.Plot <- paste(dates.split[,3],month.abb[as.numeric(dates.split[,2])],dates.split[,1],sep='-')

        if (frequency==5){
            freq <- "2-Hebdo"
        } else if (frequency==20){
            freq <- "3-Month"
        } else if (frequency==1){
            freq <- "1-Daily"
        } else {
            freq <- "NULL"
        }
        PlotName <- paste(index,
                          "_(logRx",mult,")(",
                          freq,")(dist",distribution,")(",
                          dates.Plot[1],")(",dates.Plot[2],
                          ")_",nbRegime,"_regimes_",
                          type,"_",Transition.Type,"_Transitions",sep="")
        if (!is.null(nbStepsBack)){
            PlotName <- paste(PlotName,"_",nbStepsBack,"_step",sep="")
            if (nbStepsBack>1){
                PlotName <- paste(PlotName,"s",sep="")
            }
        }
        if (!is.null(initial.Distribution)){
            PlotName <- paste(PlotName,"_init_(",paste(initial.Distribution,collapse=","),")",sep="")
        }
        if (!is.null(special.suffix.name)){
            PlotName <- paste(PlotName,"_",special.suffix.name,sep="")
        }
    }


    # --------------------------------------------
    # PARAMETER AUTOMATIC ASSIGNMENT
    # --------------------------------------------
    # Testing for correct specification of the initial parameters
    # and initialization if requested : init.assign=TRUE
    if (auto.assign){
        param0 <- init.Params(logR,nbRegime)
        mu0 <- param0$mu0
        sigma0 <- param0$sigma0
    }
    # Testing for correct specification
    if (is.null(mu0)){stop("You must define mu0, or specify init.assign=TRUE")}
    if (is.null(sigma0) || (sum(sigma0<0)>0)){stop("You must define sigma0 correctly (all positive elements), or specify init.assign=TRUE")}


    # --------------------------------------------
    # TRAINING FOR OPTIMAL PARAMETERS (HAMILTON)
    # --------------------------------------------

    start_time <- as.numeric(Sys.time()) # Start timer
    # start_time <- get.Sys.Time() # Start timer

    HMM.Train <- normal.HMM.mle(data=logR,
                                mu=mu0,
                                sigma=sigma0,
                                Gamma=Gamma0,
                                type=type,
                                distribution=distribution,
                                nu=nu,
                                Transition.Type=Transition.Type,
                                nbStepsBack=nbStepsBack,
                                initial.Distribution=initial.Distribution)

    end_time <- as.numeric(Sys.time()) # End timer
    # end_time <- get.Sys.Time() # End timer

    timeSpent.String <- get.Time.Diff(start_time, end_time)

    cat("Time for full estimation :",timeSpent.String,"\n ----------------------------- \n",sep="")


    print(HMM.Train)

    # If the model is the basic HMM, the parameters are ordered in decreasing
    # order of the volatility values
    if (type=="HMM"){
        # Getting sigma order
        augmented.sigma <- cbind(HMM.Train$sigma,c(1:nbRegime))
        order.matrix <- augmented.sigma[order(augmented.sigma[,1], decreasing = TRUE),]
        order.vector <- order.matrix[,2]

        # Ordering the parameters in decreasing order of volatility
        mu.F <- HMM.Train$mu[order.vector]
        sigma.F <- HMM.Train$sigma[order.vector]
        if (Transition.Type=="Homogeneous"){
            Gamma.F <- HMM.Train$Gamma[order.vector,order.vector]
        } else if (Transition.Type=="Diebold"){
            Gamma.F <- HMM.Train$Gamma[order.vector,]
        } else if (Transition.Type=="Diebold.w.filter"){
            Gamma.F <- HMM.Train$Gamma[order.vector,]
        }

        # Reassignment of the parameters for the returned parameters in
        # HMM.Train to match the ordered values passed in the other functions
        HMM.Train$mu <- mu.F
        HMM.Train$sigma <- sigma.F
        HMM.Train$Gamma <- Gamma.F
    }

    # --------------------------------------------
    # KIM FILTER -> SMOOTHED PROBABILITIES P[C_t|x_1:T]
    # --------------------------------------------
    print("Post Train")
    smooth.Prob <- normal.HMM.KimFilter(data=logR,
                                        type=type,
                                        distribution=distribution,
                                        mu=HMM.Train$mu,
                                        sigma=HMM.Train$sigma,
                                        Gamma=HMM.Train$Gamma,
                                        nu=nu,
                                        Transition.Type=Transition.Type,
                                        nbStepsBack=nbStepsBack,
                                        initial.Distribution=initial.Distribution)

    print("Post Kim")
    
    # --------------------------------------------
    # CREATING NEW FOLDER TO ADD THE GRAPHS
    # --------------------------------------------
    newFolderPlot <- paste(path,PlotName,sep="/")
    
    # Testing if file exists
    if(!file.exists(newFolderPlot)){
        dir.create(newFolderPlot) 
    }
    
    path <- newFolderPlot
    
    
    # --------------------------------------------
    # CREATING THE GRID PLOT
    # --------------------------------------------
    Grid.Plot.HMM.Full(timeseries=logR,
                       dates=dates,
                       state.Probs=smooth.Prob$p.ct.x1T,
                       mu=HMM.Train$mu,
                       sigma=HMM.Train$sigma,
                       nu=nu,
                       distribution=distribution,
                       Gamma=HMM.Train$Gamma,
                       name=PlotName,
                       path=path,
                       nbTicks=nbTicks,
                       Transition.Type=Transition.Type,
                       nbStepsBack=nbStepsBack,
                       type=type,
                       optimResults=HMM.Train,
                       diag.pers.Gamma=smooth.Prob$pers.Gamma,
                       Time.String=paste("Time for full estimation :",timeSpent.String,sep=""))
    
    
    # --------------------------------------------
    # CREATING THE DIAGNOSIS PLOT
    # --------------------------------------------
    DiagnosisPlotName <- paste(PlotName,"Diagnosis",sep="_")
    Diagnosis(gamma=HMM.Train$Gamma, 
              # additional.info = rbind(HMM.Train$mu,HMM.Train$sigma), 
              # special.Weights=HMM.Train$sigma,
              # network.print = FALSE,
              plotName=DiagnosisPlotName,
              docPath=path)
    
    
    # --------------------------------------------
    # RETURNING RESULTS TO THE FUNCTION NAME
    # --------------------------------------------
    return(list(HMM.Train=HMM.Train,
                smooth.Prob=smooth.Prob$p.ct.x1T,
                filtered.Prob=smooth.Prob$p.ct.x1t,
                dates.char=dates.Char[-1],
                dates=dates,
                timeSpent.String=timeSpent.String))

}


# ------------------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////
# DIAGNOSTIQUE ET MANIPULATION CHAINES DE MARKOV
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# ------------------------------------------------------------------------------------

# --------------------------------------------------------
# DIAGNOSTIQUE SUR LA MATRICE DE TRANSITION
# --------------------------------------------------------
# Fonction qui retourne plusieurs quantités d'intérêt
#   Retourne les valeurs et vecteurs propres, la distribution
#   stationnaire ainsi qu'une vsualisation de toutes ces 
#   quantité en format jpg
# --------------------------------------------------------
Diagnosis <- function(gamma,                # Transition matrix
                      plotName=NULL,        # Name to give the plot, if absent it's replace by test and the date/time
                      special.Weights=NULL, # Used to define thickness of lines in network
                      network.print=FALSE,  # TODO : Should the network be printed 
                      docPath=NULL,            # Path to where the output should be saved
                      additional.info=NULL){# TODO : make it possible to add the parameters on the graph       
    
    
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
        name <- plotName
    }
    
    # Création de l'image en format jpg et spécification de son emplacement
    FilePathName <- paste(docPath,'/',name,".jpg",sep="")
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
    }
    
    
    # Closing the plot
    dev.off()
    
    cat("Image was saved as ",FilePathName,sep="")
    
    return(list(eigen.values=eigen.values,
                eigen.vectors.D=eigen.vectors.D,
                eigen.vectors.G=eigen.vectors.G,
                dist.asympt=dist.asympt))
}



# --------------------------------------------------------
# MOMENTS DE LA CHAÎNE DE MARKOV
# --------------------------------------------------------
# Fonction pour calculer les moments des v.a. de la chaîne de Markov
#     Calcul de l'espérance, de l'autocovariance et autocorrélation 
#     de chaîne homogène/stationnaire ou non.
# --------------------------------------------------------
MomentsMC = function(Gamma,                         # Matrice de transition de la chaîne de Markov
                     initial.Distribution=NULL,     # Distribution initiale de la MC. Si absente, utilisation de la dist asymptotique
                     Transition.Type="Homogeneous", # Type de transition
                     Data=NULL,                     # Si transition non-homogène, nécessaire pour construire matrice de transition
                     nbStepsEsp=100,                # Nombre de pas en avant pour calcul de l'espérance dans la MC
                     nbStepsACF=20){                # Nombre de pas pour le calcul de Cov[C_t, C_(t+k)]
    
    K <- dim(Gamma)[1] # Dimension de l'espace d'état
    
    upsilon <- 1:K
    Upsilon <- diag(upsilon)
    
    vecUn <- rep(1,K)
    
    if (Transition.Type=="Homogeneous"){
        distribution.Distribution <- vecUn%*%solve(diag(K)-Gamma+(vecUn%*%t(vecUn))) 
        
        if(is.null(initial.Distribution)){
            
            # Utilisation de la distribution stationnaire comme distribution initiale
            initial.Distribution <- distribution.Distribution
            
            if (nbStepsEsp>1){
                cat("Avertissement : \n",
                    "  Puisque la distribution stationnaire est utilisée comme \n",
                    "  distribution initiale et puisque celle-ci multipliée par \n",
                    "  la matrice de transition donne elle même, l'espérance de \n",
                    "  chaque pas est la même et ne sera donc calculée qu'une fois.\n",
                    "  Ainsi la variable 'nbStepsEsp=",nbStepsEsp,"' est remplacée par 1.\n",sep="")
                
                nbStepsEsp <- 1
            }
        }
    }
    
    Esperance.Ct <- rep(NA,nbStepsEsp) # Initialisation du vecteur des espérances
    ACF.Ct <- rep(NA,nbStepsACF) # Initialisation du vecteur des espérances
    Gamma.m.step <- diag(K) # Matrice identité pour premier pas
    
    for (i in 1:max(nbStepsEsp,nbStepsACF)){
        # E[Ct] = u(1) %*% Gamma^(t-1) %*% c(1, ..., K)
        Gamma.m.step <- Gamma.m.step%*%Gamma
        if (i<=nbStepsEsp){
            Esperance.Ct[i] <- initial.Distribution%*%Gamma.m.step%*%upsilon
        }
        if (i<=nbStepsACF){
            ACF.Ct[i] <- distribution.Distribution%*%Upsilon%*%Gamma.m.step%*%upsilon-(distribution.Distribution%*%upsilon)^2
        }
        
    }
    plot(Esperance.Ct, ylim=c(0.999,K+0.001), type="l")
    plot(ACF.Ct, type="h")
    
    MomentsMC <- list(Esperance.Ct=Esperance.Ct,
                      ACF.Ct=ACF.Ct,
                      distribution.Distribution=distribution.Distribution)
}



# ------------------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////
# TRANSFORMATION DES PARAMÈTRES
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# ------------------------------------------------------------------------------------

# --------------------------------------------------------
# PARAMÈTRES CONTRAINTS -> PARAMÈTRES NON-CONTRAINTS
# --------------------------------------------------------
# Function to convert from natural parameters to their working equivalents
#   The idea is that by changing the constrained parameters to
#   unconstrained versions, the optimization can be done without
#   constraints
# --------------------------------------------------------
normal.HMM.N2W = function(mu, sigma, Gamma, type, Transition.Type=NULL){
    
    if (length(mu)!=length(sigma)){stop("Les dimensions des vecteurs mu et sigma ne concordent pas.")}
    if ((length(mu)!=dim(Gamma)[1])&(type!="HHMM")){stop("Les dimensions de la matrice Gamma et des vecteurs mu et sigma ne concordent pas.")}

    tsigma <- log(sigma)
    tmu <- mu # Aucune transformation nécessaire
    tGamma <- NULL # Construit plus loin

    if (type=='HMM'){
        if ((Transition.Type=="Homogeneous")||is.null(Transition.Type)){
        # --------------
        # CLASSIC HMM
        # --------------
            nbRegime <- length(mu)
            trans.Gamma <- log(Gamma/Gamma[,nbRegime])
            # Retrait de la dernière colonne de notre matrice de transition
            tGamma      <- as.vector(t(trans.Gamma[,-nbRegime]))
        } else if (Transition.Type=="Diebold" || Transition.Type=="Diebold.w.filter"){
        # --------------
        # NON-HOMOGENOUS
        # --------------
            nbRegime <- length(mu)
            tGamma <- as.vector(t(Gamma)) # Because Beta_0 and Beta_1 are unconstrained vectors of parameters
        }
        
    } else if (type=='FHMM'){
        # --------------
        # FACTORIAL HMM
        # --------------
        p1 <- Gamma[1,1]+Gamma[1,2]
        q1 <- Gamma[3,1]+Gamma[3,3]
        p2 <- Gamma[3,3]+Gamma[3,4]
        p2 <- Gamma[2,2]+Gamma[2,4]
        tGamma <- log(c(p1, p2, q1, q2)/(1-c(p1, p2, q1, q2)))

    } else if (type=='HHMM'){
        # --------------
        # HIERARCHICAL HMM
        # --------------
        Gamma.temp <- Gamma[1:4,]
        trans.Gamma <- log(Gamma.temp/Gamma.temp[,3])
        tGamma.temp <- trans.Gamma[,1:2]
        
        for (i in 1:4){
            tGamma=cbind(tGamma, t(tGamma.temp[i,]))
        }
        
        Gamma.temp <- Gamma[5:8,1:2]
        trans.Gamma <- log(Gamma.temp/Gamma.temp[,2])
        tGamma.temp <- matrix(trans.Gamma[,1],ncol=1)
        
        for (i in 1:4){
            tGamma=cbind(tGamma, tGamma.temp[i,])
        }
    }

    return(c(tmu,tsigma,tGamma))
}

# --------------------------------------------------------
# PARAMÈTRES NON-CONTRAINTS -> PARAMÈTRES CONTRAINTS
# --------------------------------------------------------
# Function to convert from Working parameters to their natural equivalents
# --------------------------------------------------------
normal.HMM.W2N = function(parvect, type, Transition.Type=NULL){
    if (type=="HMM"){

        if ((Transition.Type=="Homogeneous")||is.null(Transition.Type)){
        # --------------
        # CLASSIC HMM
        # --------------
            nparvect=length(parvect)
            nbRegime=(-1+sqrt(1+4*nparvect))/2

            nat.param  <- exp(parvect)
            mu         <- log(nat.param[1:nbRegime])
            sigma      <- nat.param[(nbRegime+1):(2*nbRegime)]

            Gamma   <- t(matrix(nat.param[(2*nbRegime+1):(2*nbRegime+nbRegime*(nbRegime-1))],
                                nrow = nbRegime-1, ncol = nbRegime))
            Gamma   <- cbind(Gamma, rep(1, nbRegime))
            Gamma   <- Gamma/apply(Gamma,1,sum)

        } else if (Transition.Type=="Diebold" || Transition.Type=="Diebold.w.filter"){
        # --------------
        # NON-HOMOGENOUS
        # --------------
            nparvect=length(parvect)
            nbRegime=2

            nat.param  <- exp(parvect)
            mu         <- log(nat.param[1:nbRegime])
            sigma      <- nat.param[(nbRegime+1):(2*nbRegime)]

            Gamma   <- log(matrix(nat.param[(2*nbRegime+1):length(nat.param)],
                                  nrow = nbRegime,byrow=T))
        }

    } else if (type=="FHMM"){ # TODO
        # --------------
        # FACTORIAL HMM
        # --------------
        nbRegime <- 4;
        mu <- parvect[1:4];

        nat.param  <- exp(parvect)
        sigma <- nat.param[5:8];

        Gamma <- matrix(nat.param[9:12],ncol=2,byrow=T)
        Gamma <- cbind(Gamma, rep(1, nbRegime))
        Gamma.temp.1 <- Gamma/apply(Gamma,1,sum)
        

    } else if (type=="HHMM"){ 
        # --------------
        # HIERARCHICAL HMM
        # --------------
        nbRegime <- 4;
        mu <- parvect[1:4];

        nat.param  <- exp(parvect)
        sigma <- nat.param[5:8];

        Gamma <- matrix(nat.param[9:16],ncol=2,byrow=T)
        Gamma <- cbind(Gamma, rep(1, nbRegime))
        Gamma.temp.1 <- Gamma/apply(Gamma,1,sum)

        Gamma <- matrix(nat.param[17:20],ncol=1,byrow=T)
        Gamma <- cbind(Gamma, rep(1, nbRegime))
        Gamma.temp.2 <- Gamma/apply(Gamma,1,sum)

        # Padding the second part to fit dimensions of the first one and
        # then stacking them
        Gamma = rbind(Gamma.temp.1, cbind(Gamma.temp.2,rep(0, nbRegime)));
    }

    return(list(mu=mu, sigma=sigma, Gamma=Gamma))
}

# --------------------------------------------------------
#   Construction de la matrice de transition (Gamma)
# --------------------------------------------------------
#   Cette fonction construit la matrice de transition d'un modèle HHMM à 2
#   états parents avec 2 états enfants chaques par le biais des
#   probabilités de transitions horizontales (p1-p2) et transitions
#   verticales (q1-q2)
# -------------------------------
# STRUCTURE POUR LE HIERARCHICAL HMM
# -------------------------------
#     Gamma.i = [a^2.11  +  a^2.12  +  e1     ;  = 1
#                a^2.22  +  a^2.21  +  e2     ;  = 1
#                a^2.33  +  a^2.34  +  e3     ;  = 1
#                a^2.44  +  a^2.43  +  e4     ;  = 1
#                a^1.11  +  a^1.12  +  0      ;  = 1
#                a^1.22  +  a^1.21  +  0      ;  = 1
#                pi.1    +  pi.2    +  0      ;  = 1
#                pi.3    +  pi.4    +  0    ] ;  = 1
# -------------------------------
# STRUCTURE POUR LE FACTORIAL HMM
# -------------------------------
# #          p1          (1-p1)          (1-p2)        p2
# Gamma.1 = [prob.i(1,1) 1-prob.i(1,1) ; 1-prob.i(1,2) prob.i(1,2)];
# #          q1          (1-q1)          (1-q2)        q2
# Gamma.2 = [prob.i(1,3) 1-prob.i(1,3) ; 1-prob.i(1,4) prob.i(1,4)];
# # Produit de Kronecker entre les 2 matrices
# Gamma=kron(Gamma.1,Gamma.2);
# -------------------------
# --------------------------------------------------------
Gamma.Build = function(prob.i, type = NULL,
                       Transition.Type,
                       nbStepsBack=0, 
                       data=NULL){

    if (Transition.Type=="Homogeneous"){
        if (type=="HHMM"){ # Récupération des paramètres pour HHMM
            a2.11 <- prob.i[1,1]
            a2.12 <- prob.i[1,2]
            e1 <- prob.i[1,3]
            a2.22 <- prob.i[2,1]
            a2.21 <- prob.i[2,2]
            e2 <- prob.i[2,3]
            a2.33 <- prob.i[3,1]
            a2.34 <- prob.i[3,2]
            e3 <- prob.i[3,3]
            a2.44 <- prob.i[4,1]
            a2.43 <- prob.i[4,2]
            e4 <- prob.i[4,3]
            a1.11 <- prob.i[5,1]
            a1.12 <- prob.i[5,2]
            a1.22 <- prob.i[6,1]
            a1.21 <- prob.i[6,2]
            pi.1 <- prob.i[7,1]
            pi.2 <- prob.i[7,2]
            pi.3 <- prob.i[8,1]
            pi.4 <- prob.i[8,2]

            # Construction de la matrice de transition à partir de ces paramètres
            Gamma.Build = matrix(c(a2.11+e1*a1.11*pi.1,  a2.12+e1*a1.11*pi.2,          e1*a1.12*pi.3,        e1*a1.12*pi.4,
                                   a2.21+e2*a1.11*pi.1,  a2.22+e2*a1.11*pi.2,          e2*a1.12*pi.3,        e2*a1.12*pi.4,
                                   e3*a1.21*pi.1,        e3*a1.21*pi.2,          a2.33+e3*a1.22*pi.3,  a2.34+e3*a1.22*pi.4,
                                   e4*a1.21*pi.1,        e4*a1.21*pi.2,          a2.43+e4*a1.22*pi.3,  a2.44+e4*a1.22*pi.4),
                                 nrow=4,byrow=T)

        } else if (type=="FHMM"){
            #          p1          (1-p1)          (1-p2)        p2
            Gamma.1 = matrix(c(prob.i(1,1), 1-prob.i(1,1), 1-prob.i(1,2), prob.i(1,2)),nrow=2,byrow=T)
            #          q1          (1-q1)          (1-q2)        q2
            Gamma.2 = matrix(c(prob.i(1,3), 1-prob.i(1,3), 1-prob.i(1,4), prob.i(1,4)),nrow=2,byrow=T)
            # Produit de Kronecker entre les 2 matrices
            Gamma.Build=kronecker(Gamma.1,Gamma.2);
        }

    } else {
        if (Transition.Type=="Diebold"){

            # if (nbStepsBack!=(dim(prob.i)[2]-1)){stop(paste("Not enough parameters for the number of steps back.\n nbStepsBack =",nbStepsBack,
            #                                                 "while you have given",dim(prob.i)[2]-1,"parameters, excluding beta_0"))}

            if (is.null(data)){stop("Please provide the past observation for the non-homogeneous transition probabilities model.")}

            Gamma.Build <- matrix(0,nrow=2,ncol=2)

            vec.1 <- c(1, data)
            vec.2 <- c(1, data)

            Gamma.Build[1,1] <- exp(prob.i[1,] %*% vec.1)/(1+exp(prob.i[1,] %*% vec.1))
            Gamma.Build[2,2] <- exp(prob.i[2,] %*% vec.2)/(1+exp(prob.i[2,] %*% vec.2))

            Gamma.Build <- Gamma.Build +  (diag(1,2) - Gamma.Build) %*% matrix(c(0, 1, 1, 0), byrow=T, ncol=2)

        } else if (Transition.Type=="Diebold.w.filter"){

            # if (2*nbStepsBack!=(dim(prob.i)[2]-1)){stop(paste("Not enough parameters for the number of steps back.\n nbStepsBack =",2*nbStepsBack,
            #                                                 "while you have given",dim(prob.i)[2]-1,"parameters, excluding beta_0"))}

            if (is.null(data)){stop("Please provide the past observation for the non-homogeneous transition probabilities model.")}

            Gamma.Build <- matrix(0,nrow=2,ncol=2)

            vec.1 <- c(1, data[1], data[2])
            vec.2 <- c(1, data[1], data[3])

            Gamma.Build[1,1] <- exp(prob.i[1,] %*% vec.1)/(1+exp(prob.i[1,] %*% vec.1))
            Gamma.Build[2,2] <- exp(prob.i[2,] %*% vec.2)/(1+exp(prob.i[2,] %*% vec.2))

            Gamma.Build <- Gamma.Build +  (diag(1,2) - Gamma.Build) %*% matrix(c(0, 1, 1, 0), byrow=T, ncol=2)
        }
    }
    return(list(Gamma=Gamma.Build,
                pers.Gamma=diag(Gamma.Build)))
}


# --------------------------------------------------------
# OBTENTION PARAMÈTRES INITIAUX construction des vecteurs de valeurs
# initiales pour mu et sigma
# --------------------------------------------------------
#       Utilisation du vecteur de données, mis en ordre puis segmenté en
#       nbRegime vecteurs. Ce sont les moyennes et écarts types de ces vecteurs
#       qui initialise les paramètres du modèle
# --------------------------------------------------------
init.Params = function(data,nbRegime){

    # Version ordonnée du vecteur de données
    sortData <- sort(data)

    # Initialisation des vecteurs résultats
    mu0 <- rep(0,nbRegime)
    sigma0 <- rep(0,nbRegime)

    # Nombre d'observations par tranche de données
    Tranche <- floor(length(sortData)/nbRegime)

    # Construction vecteur
    cat("Initial distributions\n")
    for (i in 1:nbRegime){
        dataBlock <- sortData[((i-1)*Tranche+1):(i*Tranche)]
        mu0[i] <- mean(dataBlock)
        sigma0[i] <- sd(dataBlock)
        cat(i,": N(",s(mu0[i],2),", ",s(sigma0[i],2),")\n",sep='')
    }
    return(list(mu0=mu0,
                sigma0=sigma0))
}



# ------------------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////
# FILTERING ALGORITHMS
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# ------------------------------------------------------------------------------------

# --------------------------------------------------------
# HAMILTON FILTER Implementation du filtre d'Hamilton
#    Permet le calcul des probabilités filtrées : omega_t = P[C_t | X_1:t]
# --------------------------------------------------------
normal.HMM.HamiltonFilter = function(mu,
                                     sigma,
                                     Gamma,
                                     initial.Distribution=NULL,
                                     data,
                                     distribution="Normal",
                                     Transition.Type,
                                     nbStepsBack,
                                     nu=4){

    nbRegime <- length(mu)
    n <- length(data)

    # Verification that in the event of non-homogeneous transitions, the initial distribution
    #  is correctly specified
    if (Transition.Type != "Homogeneous"){
        if (is.null(initial.Distribution)){
            stop(paste("The initial distribution must be provided explicitly if the Transitions are non-Homogeneous.\nPlease specify a vector of length ",
                       nbRegime," for the 'initial.Distribution' parameter."))
        }
    }

    if (is.null(initial.Distribution)){
        # Utilisation de la distribution stationnaire
        initial.Distribution <- solve(t(diag(nbRegime)-Gamma+1),rep(1,nbRegime))
    }

    # -------------------------------------------------
    # HAMILTON FORWARD FILTERING
    # -------------------------------------------------
    p.ct.x1t <- p.ct.x1tm1 <- matrix(nrow=nbRegime,ncol=n)
    pers.Gamma <- matrix(nrow=nbRegime,ncol=(n-1))

    # FIRST STEP OF THE PROCESS
    p.ct.x1tm1[,1] <- initial.Distribution

    if (distribution=="Normal"){
        a.j <- p.ct.x1tm1[,1]*(1/sqrt(2*pi*sigma^2))*exp((-(data[1]-mu)^2)/(2*sigma^2))

    } else if  (distribution=="Student"){
        x <- (data[1]-mu)/sigma
        a.j <- p.ct.x1tm1[,1]*((gamma((nu+1)/2)/gamma(nu/2))/sqrt(nu*pi)/((1+(x^2)/nu)^((nu+1)/2)))/sigma

    }else {
        stop("The distribution provided is not accepted. \n  Please provide one of the two accepted choices provided here : \n     Normal, Student")
    }

    a <- sum(a.j)
    llk <- log(a)
    omega.t <- a.j/a
    p.ct.x1t[,1] <- omega.t


    # STEPS (2:n) OF THE PROCESS
    for (i in 2:n){
        if (Transition.Type != "Homogeneous"){
            if (Transition.Type=="Diebold"){
                vec <- data[i-1]
            } else if (Transition.Type=="Diebold.w.filter"){
                vec <- c(data[i-1],p.ct.x1t[,(i-1)])
            }
            Gamma.Build = Gamma.Build(prob.i = Gamma,
                                      Transition.Type=Transition.Type,
                                      nbStepsBack=nbStepsBack,
                                      data=vec,
                                      type=type)

            pers.Gamma[,i-1] <- Gamma.Build$pers.Gamma
            Gamma.Build <- Gamma.Build$Gamma

        } else {
            Gamma.Build <- Gamma
        }


        p.ct.x1tm1[,i] <- omega.t%*%Gamma.Build

        if (distribution=="Normal"){
            a.j <- p.ct.x1tm1[,i]*(1/sqrt(2*pi*sigma^2))*exp((-(data[i]-mu)^2)/(2*sigma^2))

        } else if  (distribution=="Student"){
            x <- (data[i]-mu)/sigma
            a.j <- p.ct.x1tm1[,i]*((gamma((nu+1)/2)/gamma(nu/2))/sqrt(nu*pi)/((1+(x^2)/nu)^((nu+1)/2)))/sigma

        }

        a <- sum(a.j)
        llk <- llk + log(a)
        omega.t <- a.j/a
        p.ct.x1t[,i] <- omega.t
    }


    return(list(llk=llk,
                p.ct.x1tm1=p.ct.x1tm1,
                p.ct.x1t=p.ct.x1t,
                pers.Gamma=pers.Gamma))
}



# --------------------------------------------------------
# KIM FILTER - ALGORITHME DE LISSAGE Implementation du filtre de Kim
#    Permet le calcul des probabilités lisées : phi_t = P[C_t | X_1:T]
#
# Either the filtered probs are given in the function call or the
# parameters mu and sigma required to calculate these probabilities
# --------------------------------------------------------
normal.HMM.KimFilter = function(data,
                                type,
                                distribution,
                                mu=NULL,
                                sigma=NULL,
                                Gamma,
                                p.ct.x1t=NULL,
                                initial.Distribution=NULL,
                                Transition.Type,
                                nbStepsBack,
                                nu){

    # Testing if the parameters are given or if
    # the chain of filtered probabilities is given
    # Gamma must be specified either way
    if ((is.null(mu) || is.null(sigma)) && is.null(p.ct.x1t)){
        stop(cat("mu =",mu,'\nsigma =',sigma,'\np.ct.x1t =',p.ct.x1t,'\n'))
    }

    if (type=="HHMM"){
        Gamma <- Gamma.Build(prob.i=Gamma,
                             type=type,
                             Transition.Type=Transition.Type,
                             nbStepsBack=nbStepsBack)

        Gamma <- Gamma$Gamma
    }

    n <- length(data)
    nbRegime <- length(mu)

    if (Transition.Type != "Homogeneous"){
        if (is.null(initial.Distribution)){
            stop(paste("The initial distribution must be provided explicitly if the Transitions are non-Homogeneous.\nPlease specify a vector of length ",
                       nbRegime," for the 'initial.Distribution' parameter."))
        }
    }

    # If the initial distribution P[C_1] isn't specified
    # the stationary distribution of Gamma is used
    if (is.null(initial.Distribution)){
        initial.Distribution <- solve(t(diag(nbRegime)-Gamma+1),rep(1,nbRegime))
    }

    # Calcul des probabilités filtrées par l'algorithme du filtre d'Hamilton
    # et retour de la diagonale de la matrice de transition si les transitions
    # ne sont pas homogènes
    Hamilton.Filter <- normal.HMM.HamiltonFilter(mu=mu,
                                                 sigma=sigma,
                                                 Gamma=Gamma,
                                                 initial.Distribution=initial.Distribution,
                                                 data=data,
                                                 distribution=distribution,
                                                 Transition.Type=Transition.Type,
                                                 nbStepsBack=nbStepsBack,
                                                 nu=nu)
    # We keep the filtered probabilities from the Hamilton filter
    p.ct.x1t <- Hamilton.Filter$p.ct.x1t

    # Initialization of the matrix to keep the diagonal of the transition matrix
    # step-by-step if the transitions aren't homogeneous
    pers.Gamma <- matrix(nrow=nbRegime,ncol=(n-1))
    if (Transition.Type!="Homogeneous"){
        pers.Gamma <- Hamilton.Filter$pers.Gamma
        if (sum(is.na(pers.Gamma))>0){stop(paste("NAs in pers.Gamma are at",which(is.na(pers.Gamma))," in Kim filter\n"))}
    }

    # Alert message and information on localization of mistake if NAs are returned
    # by the Hamilton filter
    if (sum(is.na(p.ct.x1t))>0){stop(paste("NAs in p.ct.x1t are at",paste(which(is.na(p.ct.x1t)), collapse=""),"in Kim filter\n"))}

    # Initialisation de la matrice de probabilités lissées
    p.ct.x1T <- matrix(NA,nbRegime,n)

    # Probabilite de la derniere observation de la variable latente
    p.ct.x1T[,n] <- p.ct.x1t[,n]

    # Algorithme de lissage (Filtre de Kim)
    for (t in (n-1):1){
        phi.t.ij <- matrix(NA,nbRegime,nbRegime)

        if (Transition.Type != "Homogeneous"){
            Gamma.Built <- diag(pers.Gamma[,t]) +  (diag(1,2) - diag(pers.Gamma[,t])) %*% matrix(c(0, 1, 1, 0), byrow=T, ncol=2)
            if (sum(is.na(Gamma.Built))){stop(paste("Problem with Gamma at step",t,"in the Kim filter function."))}
        } else {
            Gamma.Built <- Gamma
        }

        for (i in 1:nbRegime){
            for (j in 1:nbRegime){
                # Calcul des probabilités conjointes
                phi.t.ij[i,j] <- p.ct.x1T[j,t+1]*(p.ct.x1t[i,t]*Gamma.Built[i,j]/sum(p.ct.x1t[,t]*Gamma.Built[,j]))
            }
        }

        # Throw error message if NAs in the intermediate step of kim Smoother
        if (sum(is.na(phi.t.ij))){
            print("pers.Gamma[,t]")
            print(pers.Gamma[,t])
            print("p.ct.x1T[,t+1]")
            print(p.ct.x1T[,t+1])
            print("Gamma.Built")
            print(Gamma.Built)
            print("phi.t.ij")
            print(phi.t.ij)
            stop(paste("Problem with phi.t.ij at step t =",t,"in the Kim filter function."))
        }

        p.ct.x1T[,t] <- rowSums(phi.t.ij)
    }

    return(list(p.ct.x1T=p.ct.x1T,
                p.ct.x1t=p.ct.x1t,
                pers.Gamma=pers.Gamma))
}

# ------------------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////
# OPTIMIZATION FUNCTIONS
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# ------------------------------------------------------------------------------------

# --------------------------------------------------------
# FONCTION OBJECTIF afin de calculer la vraisemblance et la log-vraisemblance
#     Le calcul est effectué à l'aide du filtre d'hamilton (Probabilités
#     filtrées P[C_t | X_1:t])
# --------------------------------------------------------
normal.HMM.mllk = function(parvect, data, type, distribution,Transition.Type,
                           nbStepsBack, nu,
                           initial.Distribution=initial.Distribution){
    n <- length(data)

    natural.par <- normal.HMM.W2N(parvect=parvect,
                                  type=type,
                                  Transition.Type=Transition.Type)
    mu <- natural.par$mu
    sigma <- natural.par$sigma
    Gamma <- natural.par$Gamma

    if (type=="HHMM"){
        Gamma.vec <- Gamma.Build(prob.i=Gamma,
                             type=type,
                             Transition.Type=Transition.Type,
                             nbStepsBack=nbStepsBack)

        Gamma <- Gamma.vec$Gamma
    }

    optim.results <- normal.HMM.HamiltonFilter(mu=mu,
                                               sigma=sigma,
                                               Gamma=Gamma,
                                               data=data,
                                               distribution=distribution,
                                               nu=nu,
                                               Transition.Type=Transition.Type,
                                               nbStepsBack=nbStepsBack,
                                               initial.Distribution=initial.Distribution)
    # print(optim.results$llk)
    llk <- -optim.results$llk
    return(llk)
}


# --------------------------------------------------------
#FONCTION OPTIMISATION Fonction afin d'optimiser numériquement la log-vraisemblance
#   Optimisation à l'aide du filtre d'Hamilton (Probabilités filtrées
#   P[C_t | X_1:t])
# --------------------------------------------------------
normal.HMM.mle <- function(data, mu, sigma, Gamma, type, distribution,Transition.Type,
                           nbStepsBack, nu,
                           initial.Distribution=initial.Distribution){
    nbRegime <- length(mu)
    parvect0       <- normal.HMM.N2W(mu, sigma, Gamma, type,
                                     Transition.Type=Transition.Type)

    #Calcul du EVM
    start_time <- Sys.time()
    mod            <- nlm(normal.HMM.mllk,
                          parvect0,
                          data=data,
                          type=type,
                          distribution=distribution,
                          nu=nu,
                          Transition.Type=Transition.Type,
                          nbStepsBack=nbStepsBack,
                          initial.Distribution=initial.Distribution,
                          iterlim=1000,
                          gradtol = 1e-20)
    print(mod)

    natural.par.optim <- normal.HMM.W2N(mod$estimate,type,
                                        Transition.Type=Transition.Type)

    if (distribution=='Student'){
        optim.sigma=natural.par.optim$sigma*sqrt(nu/(nu-2))
    } else {
        optim.sigma=natural.par.optim$sigma
    }

    end_time <- Sys.time()
    cat("Time for optimization :",end_time - start_time,"\n ----------------------------- \n")

    #Log-vraisemblance
    mllk           <- -mod$minimum
    np             <- length(parvect0)
    AIC            <- 2*(-mllk+np)
    n              <- sum(!is.na(data))
    BIC            <- np*log(n)-2*mllk
    list(mu=natural.par.optim$mu, sigma=optim.sigma, Gamma=natural.par.optim$Gamma, mllk=mllk,
         AIC=AIC, BIC=BIC)
}



# ------------------------------------------------------------------------------------
# ////////////////////////////////////////////////////////////////////
# GRAPHICAL OUTPUT FUNCITONS
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# ------------------------------------------------------------------------------------


# --------------------------------------------------------
# GRAPHICAL PRINT OF TIME-SERIES
# --------------------------------------------------------
ts.plot.dev = function(y,nbTicks=30, xlabels = NULL,custom.Colors=NULL, axis.label.size=1.5,max.Prob.States=NULL,lwd=4,nbStepsBack){
    par(mar=c(20,16,12,12),cex.axis=axis.label.size,cex.lab=6,las=2)
    # Getting length of the time series
    n <- length(y)
    # print(n)
    # print(c(max.Prob.States))

    # Spacing between the ticks on the x axis
    xTicksSpace <- n/nbTicks

    # plot(y, type='n',xaxt="n",xlab="", ylab=deparse(substitute(y)),cex.main=axis.label.size, xlim=c(1,n),ylim=c(min(y)-1,max(y)+1))
    if (is.null(max.Prob.States)){
        plot(y, type='h',lwd=lwd, xaxt="n",xlab="", ylab="",cex.main=axis.label.size, xlim=c(1,n),ylim=c(min(y)-1,max(y)+1), xaxs="i")
    } else {
        plot(y, type='h', lwd=lwd, xaxt="n",xlab="", ylab="",cex.main=axis.label.size, xlim=c(1,n),ylim=c(min(y)-1,max(y)+1), xaxs="i",col=custom.Colors[max.Prob.States])
    }

    # lines(y, col="black")
    title(ylab="Timeserie", line=7)

    # Adding the correct ticks on the x axis
    axis(side=1,at=c(seq(from=1, to=n, by=xTicksSpace),n),labels = xlabels[c(seq(from=1, to=n, by=xTicksSpace),n)],lwd=4,lwd.ticks=2,col.ticks="red",cex.axis=axis.label.size)

    # Adding grid to the plot
    grid(lty=1, col=gray(.9))
    # Adding an horizontal line at y=0
    abline(h=0, col="blue")
}


# --------------------------------------------------------
# PROBABILITIES FILL COLOR
# --------------------------------------------------------
HMM.Stack.Plot = function(series,xlabels,nbTicks=30,custom.Colors, axis.label.size=1.5, nbStepsBack,custom.y.range=c(0,1)){

    par(mar=c(20,16,12,12),cex.axis=axis.label.size,cex.lab=6,las=2)

    nbRegime <- dim(series)[1] # Nb of regimes
    n <- dim(series)[2] # Nb of observations
    cat(nbRegime,"regimes and",n,"observations.\n")

    # Spacing between the ticks on the x axis
    xTicksSpace <- n/nbTicks

    seq.labels <- c(seq(from=1, to=n, by=xTicksSpace),n)


    stacked.series <- matrix(0,nbRegime+1,n)
    # Stacking
    for (i in 1:nbRegime){
        stacked.series[i+1,] <- stacked.series[i,] + series[i,]
    }

    plot(1:n,stacked.series[1,],type = 'n',xaxt="n", ylim = custom.y.range,
         ylab =  "", xlab = '', col = "red", xlim=c(1,n), xaxs="i", yaxs="i")
    title(ylab = 'probabilities', line=9)

    for (i in 2:nbRegime+1){
        lines(stacked.series[i,], col = "black")
    }

    # Adding the correct ticks on the x axis
    axis(side=1,at=seq.labels,labels = xlabels[seq.labels],lwd=4,lwd.ticks=2,col.ticks="red",cex.axis=axis.label.size)

    for (i in 1:nbRegime){
        polygon(c(1:n, rev(1:n)), c(stacked.series[i+1,], rev(stacked.series[i,])),
                col = custom.Colors[i], border = NA)
    }
    abline(h=0.5, col="white")
    # print(t(stacked.series))
}






# --------------------------------------------------------
# MAIN PLOT FUNCTION
# --------------------------------------------------------
Grid.Plot.HMM.Full = function(timeseries, dates, state.Probs, type, mu,
                              sigma,nu,distribution, Gamma, name,path,
                              resolution,nbTicks,Transition.Type,
                              nbStepsBack,optimResults,
                              diag.pers.Gamma,
                              Time.String=NULL){

    print("In Grid.Plot.HMM.Full")
    # Colors
    custom.Colors.vec <- matrix(c("firebrick3","blue4",NA,NA,
                              "firebrick3","deepskyblue1","blue4",NA,
                              "firebrick3","darksalmon","deepskyblue1","blue4"),ncol=4,byrow=T)

    nbRegime <- length(mu)
    custom.Colors <- custom.Colors.vec[nbRegime-1,]

    n <- length(dates)

    # Title size
    title.size <- 6
    axis.label.size <- 5
    line.width <- 8
    cex.lab <- 7.5

    # If the model is Hierarchical, the Gamma matrix is built here
    if (type=="HHMM"){
        Gamma.Built <- Gamma.Build(prob.i=Gamma,
                                   type=type,
                                   Transition.Type)

    }

    # Getting dates in correct format while keeping only year and month
    dates.split <- matrix(unlist(strsplit(dates,split="\\-")),ncol=3,byrow=T)
    dates.Plot <- paste(month.abb[as.numeric(dates.split[,2])],dates.split[,1],sep='-')

    # Testing for the validity of the probabilities received
    if (sum(is.na(state.Probs))>0){stop(paste("These steps of the smoothed probabilities are NA : ",
                                              paste(which(is.na(state.Probs)),collapse=", "),"\n",sep=""))}

    # Getting the state with the maximum probability
    max.Prob.States <- colMax(state.Probs)

    # JPEG DIMENSIONS FOR OUTPUTED FILE
    image.width <- 7000
    image.heigth <- 3000

    # Defines the layout of the figure output
    if (Transition.Type=="Homogeneous"){
        # Creating the empty image to be filled with plots
        jpeg(paste(path,"/",name,".jpg",sep=""), width = image.width, height = image.heigth)
        layout(mat = matrix(c(1,1,
                              2,2,
                              3,4,
                              3,4,
                              5,6,
                              7,6),ncol=2,nrow = 6,byrow=TRUE))
    } else {
        # Creating the empty image to be filled with plots
        jpeg(paste(path,"/",name,".jpg",sep=""), width = image.width, height = image.heigth*7/6)
        layout(mat = matrix(c(1,1,
                              2,2,
                              8,8,
                              3,4,
                              3,4,
                              5,6,
                              7,6),ncol=2,nrow = 7,byrow=TRUE))
    }

    # -----------------------------------------------------------------------------------------
    print("Time series printing")
    # -----------------------------------------------------------------------------------------

    ts.plot.dev(y=timeseries[nbStepsBack:length(timeseries)],
                xlabels=dates.Plot[nbStepsBack:length(timeseries)],
                nbTicks=nbTicks,
                custom.Colors=custom.Colors,
                axis.label.size=axis.label.size,
                max.Prob.States=max.Prob.States,
                lwd = 8,
                nbStepsBack=nbStepsBack)
    title(main = paste(Sys.time(), "\n",name, " series - ",Time.String, sep=""),cex.main=title.size)

    # -----------------------------------------------------------------------------------------
    print("Stacking")
    # -----------------------------------------------------------------------------------------

    xlabels <- dates.Plot[nbStepsBack:length(timeseries)]

    HMM.Stack.Plot(series=state.Probs,
                   xlabels=xlabels,
                   nbTicks=nbTicks,
                   custom.Colors=custom.Colors,
                   axis.label.size=axis.label.size)

    title(main = "Smooth probabilities - P[C_t | X_1:T]",cex.main=title.size)

    # -----------------------------------------------------------------------------------------
    print("Distribution plotting")
    # -----------------------------------------------------------------------------------------

    par(mar=c(20,24,24,24),cex.axis=axis.label.size,cex.lab=cex.lab,las=0)

    # Calculation of the frequency of each regime
    state.Realisation <- matrix(ncol=nbRegime)
    for (i in 1:nbRegime){
        state.Realisation[i] <- sum(max.Prob.States==i)
    }
    if (nbRegime<4){
        state.Realisation.pad <- c(state.Realisation,rep(0,4-nbRegime))
    } else {
        state.Realisation.pad <- state.Realisation
    }

    # Keeping only the parameters of the realized states for the distribution plotting
    nbRegime.vector <- 1:nbRegime
    mu.Realized <- mu[state.Realisation>0]
    sigma.Realized <- sigma[state.Realisation>0]
    custom.Colors.Realized <- custom.Colors[state.Realisation.pad>0]

    regimes.Realized <- nbRegime.vector[state.Realisation>0]

    nbRegime.Realized <- length(mu.Realized)

    step.distr = 2*(max(mu.Realized)+3.5*max(sigma.Realized))/100000
    x <- seq(from= -(max(mu.Realized)+3.5*max(sigma.Realized)),
             to=  (max(mu.Realized)+3.5*max(sigma.Realized)),
             by= step.distr)
    y <- matrix(nrow=nbRegime.Realized,ncol=length(x))

    if (distribution=="Normal"){
        for (i in 1:nbRegime.Realized){
            y[i,] <- dnorm(x, mean=mu.Realized[i], sd=sigma.Realized[i])
        }
    } else if (distribution=="Student"){
        for (i in 1:nbRegime.Realized){
            y[i,] <- dt((x-mu.Realized[i])/sigma.Realized[i],df=nu)/sigma.Realized[i]
        }
    }

    plot(x,y[1,],col=custom.Colors.Realized[1],
         type="l",lwd=line.width,
         xlab="",
         ylab="",xlim=c(min(x)-1,max(x)+1),ylim=c(0,max(y)+0.1),xaxt="n")

    if (nbRegime.Realized>1){
        for (i in 2:nbRegime.Realized){
            lines(x=x,y=y[i,],col=custom.Colors.Realized[i],lwd=line.width,xaxt="n")
        }
    }

    abline(v=mu.Realized, col=custom.Colors.Realized,lwd=line.width, lty=2)

    title(main = "Statistical Distributions\nof realized states",
          cex.main=title.size)
    axis(side=1,line=2, tick=F)

    # -----------------------------------------------------------------------------------------
    print("BoxPlotting")
    # -----------------------------------------------------------------------------------------

    par(mar=c(20,24,24,24),cex.axis=axis.label.size,cex.lab=cex.lab,las=0)
    obs.plus.state <- data.frame(rbind(t(timeseries[nbStepsBack:length(timeseries)]),t(max.Prob.States)),
                                 row.names = c("timeseries","max.Prob.States"))

    boxplot(timeseries[nbStepsBack:length(timeseries)]~max.Prob.States,
            data=obs.plus.state,
            col=custom.Colors.Realized,
            border="brown",
            lwd=4,
            varwidth = TRUE, xaxt="n")

    title(main = "Boxplots of observations by\nmaximum probability",
          cex.main=title.size)
    title(xlab="Regimes",
          ylab="Observations",
          line=10)

    boxlabels <- paste("state ",regimes.Realized,
                       " - [ ",
                       as.character(round(100*state.Realisation[state.Realisation>0]/sum(state.Realisation[state.Realisation>0]),2)),"% ]",sep="")

    axis(side=1,at=1:length(mu.Realized),
         labels=boxlabels,
         line=3, tick=F)

    # -----------------------------------------------------------------------------------------
    print("mu and sigma")
    # -----------------------------------------------------------------------------------------

    rowNames <- c("mu","sigma")
    if (type=="HMM"){
        colNames <- paste("state",as.character(1:nbRegime),sep="")
    } else if (type=="HHMM"){
        colNames <- as.character(paste(paste("state",c(1,1,2,2),sep=""),1:2,sep=""))
    }

    #Set up the plot area and plot the matrix
    par(mar=c(16, 60, 16, 60))

    color2D.matplot(rbind(mu,sigma),cellcolors="#ffffff",
                    main="Parameters",
                    cex.main=title.size,
                    show.values=3,
                    vcol=rgb(0,0,0),
                    axes=FALSE,
                    vcex=cex.lab, xlab="", ylab="")

    # ROW LABELS
    axis(2, at=seq(2, 1, -1)-0.5, labels=rowNames, padj=0.5,line = 2, las=2, tick=F,cex.axis=cex.lab)

    # COLUMN LABELS
    for (i in 1:nbRegime){
        axis(1, at=i-0.5, labels=colNames[i], padj=0.5,line = 2, las=0, tick=F,cex.axis=cex.lab,col.axis=custom.Colors[i])
    }

    # -----------------------------------------------------------------------------------------
    print("Gamma")
    # -----------------------------------------------------------------------------------------
    if (Transition.Type=="Homogeneous"){
        if (type=="HMM"){
            rowNames <- paste("state",as.character(1:nbRegime),sep="")
            colNames <- paste("state",as.character(1:nbRegime),sep="")
        } else if (type=="HHMM"){
            rowNames <- c("a^2.11 - a^2.12 - e1",
                          "a^2.22 - a^2.21 - e2",
                          "a^2.33 - a^2.34 - e3",
                          "a^2.44 - a^2.43 - e4",
                          "a^1.11 - a^1.12",
                          "a^1.22 - a^1.21",
                          "pi.1 - pi.2",
                          "pi.3 - pi.4")
            colNames <- rep("",3)
        }
        xmin <- 0
        xmax <- 1
        #Generate the palette for the matrix and the legend.  Generate labels for the legend
        palmat <- color.scale(Gamma, c(1, 0.4), c(1, 0.4), c(0.96, 1))

    } else if (Transition.Type=="Diebold"){
        rowNames <- paste("Beta",as.character(1:nbRegime),sep="")

        if (nbStepsBack>1){
            colNames <- colNames <- c("1",paste("x_(t-",nbStepsBack:1,")",sep=""))
        } else {
            colNames <- c("1","x_(t-1)")
        }

        #Generate the palette for the matrix and the legend.  Generate labels for the legend
        palmat <- "#ffffff"

    } else if (Transition.Type=="Diebold.w.filter"){
        rowNames <- paste("Beta",as.character(1:nbRegime),sep="")

        if (nbStepsBack>1){
            colNames <- c("1",
                          paste("x_(t-",nbStepsBack:1,")",sep=""),
                          paste("P[C_(t-",nbStepsBack:1,")|x_1:(t-",nbStepsBack:1,")]",sep=""))
        } else {
            colNames <- c("1",
                          "x_(t-1)",
                          "P[C_(t-1)|x_1:(t-1)]")
        }

        #Generate the palette for the matrix and the legend.  Generate labels for the legend
        palmat <- "#ffffff"
    }

    #Set up the plot area and plot the matrix
    par(mar=c(16, 60, 16, 60))

    color2D.matplot(Gamma, cellcolors=palmat,
                    main="Gamma",
                    cex.main=title.size,
                    show.values=3,
                    vcol=rgb(0,0,0),
                    axes=FALSE,
                    vcex=cex.lab, xlab="", ylab="")

    # ROW LABELS
    if (type=="HMM"){
        if (Transition.Type=="Diebold" || Transition.Type=="Diebold.w.filter"){
            for (i in 1:nbRegime){
                axis(2, at=(nbRegime-i+1)-0.5, labels=rowNames[i], padj=0.5,line = 2, las=0, tick=F,cex.axis=cex.lab,col.axis=custom.Colors[i])
            }
        } else if (Transition.Type=="Homogeneous"){
            axis(2, at=seq(nbRegime,1, -1)-0.5, labels=rowNames, padj=0.5,line = 2, las=0, tick=F,cex.axis=cex.lab)
        }
    } else if (type=="HHMM"){
        axis(2, at=seq(from=8, to=1, by=-1), labels=rowNames, padj=1.5,line = 2, las=2, tick=F,cex.axis=cex.lab)
    }

    # COLUMN LABELS
    if (type=="HMM"){
        if (Transition.Type=="Diebold"){
            axis(1, at=c(1:(2+(nbStepsBack-1)))-0.5, labels=colNames, padj=0.5,line = 2, las=0, tick=F,cex.axis=cex.lab)
        } else if (Transition.Type=="Diebold.w.filter"){
            axis(1, at=c(1:(3+2*(nbStepsBack-1)))-0.5, labels=colNames, padj=0.5,line = 2, las=0, tick=F,cex.axis=cex.lab)
        } else {
            axis(1, at=seq(1, nbRegime, 1)-0.5, labels=colNames, padj=0.5,line = 2, las=0, tick=F,cex.axis=cex.lab)
        }
    } else if (type=="HHMM"){

    }

    # -----------------------------------------------------------------------------------------
    print("mllk, AIC, BIC")
    # -----------------------------------------------------------------------------------------

    rowNames <- c("Results")
    colNames <- c("mllk", "AIC", "BIC")

    #Set up the plot area and plot the matrix
    par(mar=c(16, 60, 16, 60))

    likelihood.values <- matrix(c(optimResults$mllk, optimResults$AIC, optimResults$BIC),nrow=1)

    color2D.matplot(likelihood.values,cellcolors="#ffffff",
                    main="Parameters",
                    cex.main=title.size,
                    show.values=3,
                    vcol=rgb(0,0,0),
                    axes=FALSE,
                    vcex=cex.lab, xlab="", ylab="")

    # ROW LABELS
    axis(2, at=0.5, labels=rowNames, padj=0.5,line = 4, las=0, tick=F,cex.axis=cex.lab)

    # COLUMN LABELS
    for (i in 1:3){
        axis(1, at=i-0.5, labels=colNames[i], padj=0.5,line = 2, las=0, tick=F,cex.axis=cex.lab)
    }

    # -----------------------------------------------------------------------------------------
    if (Transition.Type!="Homogeneous"){
        print("Transition probabilities on diagonal of Gamma")
    # -----------------------------------------------------------------------------------------

        HMM.Stack.Plot(series=diag.pers.Gamma,
                       xlabels=xlabels,
                       nbTicks=nbTicks,
                       custom.Colors=custom.Colors,
                       axis.label.size=axis.label.size,
                       custom.y.range=c(0,2))

        # par(mar=c(20,16,12,12),cex.axis=axis.label.size,cex.lab=6,las=2)
        # plot(1:(n-1),diag.pers.Gamma[1,],col=custom.Colors[1,1],type="l",lwd=line.width, xaxt="n",xlab="", ylab="",cex.main=axis.label.size, xlim=c(1,n),ylim=c(0,1), xaxs="i")
        # lines(diag.pers.Gamma[2,],col=custom.Colors[1,2],lwd=line.width)
        #
        # # Spacing between the ticks on the x axis
        # xTicksSpace <- n/nbTicks

        title(main="Time-dependent transition probabilities matrix diagonal", line=2,cex.main=title.size)

        # # Adding the correct ticks on the x axis
        # axis(side=1,
        #      at=c(seq(from=1, to=n, by=xTicksSpace),n),
        #      labels = xlabels[c(seq(from=1, to=n, by=xTicksSpace),n)],
        #      lwd=line.width,
        #      lwd.ticks=2,
        #      col.ticks="red",
        #      cex.axis=axis.label.size)
        #
        # # Adding grid to the plot
        # grid(lty=1, col=gray(.9))
        # # Adding an horizontal line at y=0
        # abline(h=0, col="blue")
    }

    # --------------------------------------------------------
    dev.off()
    print(name)
}



