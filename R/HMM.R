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
# Last version : april 16th, 2019
# ----------------------------------------------------------------------------------------------------
#set working directory
path <- '~/Documents/GitHub/HHMM'
setwd(path.expand(path)) # Setting path

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

# ——————————————————————————————————————————————————————————————————————————
# FONCTION PERMETTANT D'IMPRIMER LE RÉSULTATS EN FORMAT LATEX
# ——————————————————————————————————————————————————————————————————————————
Make.LaTeX.Table = function(R.Matrix.Object,
                            caption=TRUE,
                            title = "",
                            Col.Titles = NULL,
                            Row.Titles = NULL,
                            n.dec = 3,
                            type = "Table",
                            Row.Pos="c",
                            Cross.Lines=FALSE,
                            print.Cons=TRUE,
                            copy.CB=TRUE,...){ # Ecrit le code LaTeX dun environnement Tabular a partir dune matrice R
    
    n.col = max(ncol(R.Matrix.Object),1)
    n.Row = max(nrow(R.Matrix.Object),length(R.Matrix.Object))
    
    Row.Titles.Ind = FALSE
    Col.Titles.Ind = FALSE
    if(!is.null(Row.Titles)){Row.Titles.Ind = TRUE}
    if(!is.null(Col.Titles)){Col.Titles.Ind = TRUE}
    if((Row.Titles.Ind | Cross.Lines)  & !(n.Row==1)){
        pos.col = "r |"
        nb.col.pos = n.col
        if(Cross.Lines & !Row.Titles.Ind){
            nb.col.pos = nb.col.pos - 1
        }
    }else{
        pos.col=""
        nb.col.pos = n.col
    }
    
    pos.col = paste(pos.col," ",chr(42),'{',nb.col.pos,"}","{",Row.Pos,"} ",sep="")
    
    str = paste("% ------------------------ \n", 
                chr(92),"begin{minipage}{",
                chr(92),"linewidth} \n",
                chr(92),"centering",
                if(caption){paste(chr(92),"captionof{table}{",title,"} \n",sep="")},
                chr(92),"label{",title,"} \n",
                chr(92),"begin{tabular}[t]{",pos.col,"} \n",sep="")
    
    Print.Nth.Element = function(Element,String,Col,last=FALSE){
        if (is.numeric(Element)){
            String=paste(String," $",round(Element,n.dec),"$ ",if(!last){chr(38)}else{""},sep="")
        }else{
            String=paste(String," ",Element," ",if(!last){chr(38)}else{""},sep="")
        }
        Print.Nth.Element = String # Returns the input string modified
    }
    
    if (Col.Titles.Ind){
        lenoff = length(Col.Titles)
        if (Row.Titles.Ind){
            if (length(Col.Titles) == n.col){
                Col.Titles = append(Col.Titles,"",after=0)
            }
        }
        for (i in 1:(lenoff-1)){
            str = paste(str,Col.Titles[i],chr(38))
        }
        str = paste(str," ",Col.Titles[lenoff]," ", chr(92),chr(92),"\n ", chr(92),"hline \n", sep="")
    }
    
    for (i in 1:n.Row){
        if (!is.null(Row.Titles[i])){ # Adds the Row titles if there are any
            str = paste(str,Row.Titles[i],chr(38)) # Adds the ith one
        }
        
        if (n.col!=1){ # If there is more than one Column
            for (k in 1:n.col){ # Does all the Columns and doesnt add & on the last one
                str = Print.Nth.Element(R.Matrix.Object[i,k],str,k,last=isTRUE(k==n.col))
            }
        } else {
            str = Print.Nth.Element(R.Matrix.Object[i],str,n.col,last=TRUE) # If one Column only, does that one only
        }
        str=paste(str," ",chr(92),chr(92),sep="") # Adds the double backslash at the end of the line
        if(!(i==n.Row)){str=paste(str,"\n",sep="")} # If it is the last line, it wont add line jump
        if (i==1){ # If the first line was already printed
            if (Cross.Lines & !Col.Titles.Ind & !(n.Row==1)){ # If the line is to be printed but no col.titles giving expressively
                str = paste(str, chr(92),"hline \n", sep="")
            }
        }
    }
    
    
    
    # Fin de lexecution et retour du resultat
    if(copy.CB){
        Copie.Presse.Papier(paste(str,"\n",
                                  chr(92),"end{tabular} \n",
                                  chr(92),"end{minipage}","\n ~",chr(92),chr(92),
                                  "\n % ------------------------ \n",sep="")) # Copie le string concatene au presse-papier
    }
    
    if(print.Cons){
        str=cat(paste(str,"\n",
                      chr(92),"end{tabular} \n",
                      chr(92),"end{minipage}","\n ~",chr(92),chr(92),
                      "\n % ------------------------ \n",sep="")) # Limprime aussi dans la console R
    }
}

chr <- function(n) { rawToChar(as.raw(n)) }

Copie.Presse.Papier <- function(string) {
    os <- Sys.info()[['sysname']]
    if (os == "Windows") { # Si systeme dexploitation windows
        return(utils::writeClipboard(string))
    } else if (os == "Darwin") { # Si systeme dexploitation iOS
        Mac.Copie.Presse.Papier <- function(string){
            presse.papier <- pipe("pbcopy", "w")
            cat(string, file = presse.papier, sep = "\n")
            close(presse.papier)	# Fermer lobjet presse-papier
        }
        return(Mac.Copie.Presse.Papier(string))
    }
}

# ——————————————————————————————————————————————————————————————————————————
# Rounding with n zeros function <- Returns a string and not a numeric type
# ——————————————————————————————————————————————————————————————————————————
s <- function(x,n){
    sprintf(paste("%.",n,"f",sep=""), round(x,n))
}

# ——————————————————————————————————————————————————————————————————————————
# GETTING COLUMNWISE MAXIMUM
# ——————————————————————————————————————————————————————————————————————————
colMax <- function(data){
	apply(data,2, which.max)
}

# ————————————————————————————————————————————————————————————————————————————————————
# ////////////////////////////////////////////////////////////////////
# TRANSFORMATION DES PARAMÈTRES
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# ————————————————————————————————————————————————————————————————————————————————————

# ————————————————————————————————————————————————————————
# PARAMÈTRES CONTRAINTS -> PARAMÈTRES NON-CONTRAINTS
# ————————————————————————————————————————————————————————
# Function to convert from natural parameters to their working equivalents
#   The idea is that by changing the constrained parameters to
#   unconstrained versions, the optimization can be done without
#   constraints
# ————————————————————————————————————————————————————————
normal.HMM.N2W = function(mu, sigma, Gamma, type){
    
    tsigma <- log(sigma)
    tmu <- mu
    tGamma <- NULL
    
    if (type=='HMM'){
        nbRegime <- length(mu)
        trans.Gamma <- log(Gamma/Gamma[,nbRegime])
        # Retrait de la dernière colonne de notre matrice de transition
        tGamma      <- as.vector(t(trans.Gamma[,-nbRegime]))
        
    } else if (type=='FHMM'){
        p1 <- Gamma[1,1]+Gamma[1,2]
        q1 <- Gamma[3,1]+Gamma[3,3]
        p2 <- Gamma[3,3]+Gamma[3,4]
        p2 <- Gamma[2,2]+Gamma[2,4]
        tGamma <- log(c(p1, p2, q1, q2)/(1-c(p1, p2, q1, q2)))
        
    } else if (type=='HHMM'){
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

# ————————————————————————————————————————————————————————
# PARAMÈTRES NON-CONTRAINTS -> PARAMÈTRES CONTRAINTS
# ————————————————————————————————————————————————————————
# Function to convert from Working parameters to their natural equivalents
# ————————————————————————————————————————————————————————
normal.HMM.W2N = function(parvect, type){
    if (type=="HMM"){
        
        nparvect=length(parvect)
        nbRegime=(-1+sqrt(1+4*nparvect))/2
        
        nat.param  <- exp(parvect)
        mu         <- log(nat.param[1:nbRegime])
        sigma      <- nat.param[(nbRegime+1):(2*nbRegime)]
        
        Gamma   <- t(matrix(nat.param[(2*nbRegime+1):(2*nbRegime+nbRegime*(nbRegime-1))], 
                            nrow = nbRegime-1, ncol = nbRegime))
        Gamma   <- cbind(Gamma, rep(1, nbRegime))
        Gamma   <- Gamma/apply(Gamma,1,sum)
        
        # Just to fit the format of the HHMM model below
        Gamma.Built = Gamma;
        
    } else if (type=="FHMM"){ # Factorial Hidden Markov Model
        nbRegime <- 4;
        mu <- parvect[1:4];
        
        nat.param  <- exp(parvect)
        sigma <- nat.param[5:8];
        
        Gamma <- matrix(nat.param[9:12],ncol=2,byrow=T)
        Gamma <- cbind(Gamma, rep(1, nbRegime))
        Gamma.temp.1 <- Gamma/apply(Gamma,1,sum)
        
        # Building the transition matrix Gamma only for the purpose of delta
        # calculation
        Gamma.Built = Gamma.Build.HHMM(prob.i,type);
        
    } else if (type=="HHMM"){ # Hierchical Hidden Markov Model
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
        
        # Building the transition matrix Gamma only for the purpose of delta
        # calculation
        Gamma.Built = Gamma.Build.HHMM(prob.i=Gamma,type);
        
    }
    
    #distribution stationnaire
    delta    <- solve(t(diag(nbRegime)-Gamma.Built+1),rep(1,nbRegime))
    
    return(list(mu=mu, sigma=sigma, Gamma=Gamma))
}

# ————————————————————————————————————————————————————————
#   Construction de la matrice de transition (Gamma)
# ————————————————————————————————————————————————————————
#   Cette fonction construit la matrice de transition d'un modèle HHMM à 2
#   états parents avec 2 états enfants chaques par le biais des
#   probabilités de transitions horizontales (p1-p2) et transitions
#   verticales (q1-q2)
# ———————————————————————————————
# STRUCTURE POUR LE HIERARCHICAL HMM
# ———————————————————————————————
#     Gamma.i = [a^2.11  +  a^2.12  +  e1     ;  = 1
#                a^2.22  +  a^2.21  +  e2     ;  = 1
#                a^2.33  +  a^2.34  +  e3     ;  = 1
#                a^2.44  +  a^2.43  +  e4     ;  = 1
#                a^1.11  +  a^1.12  +  0      ;  = 1
#                a^1.22  +  a^1.21  +  0      ;  = 1
#                pi.1    +  pi.2    +  0      ;  = 1
#                pi.3    +  pi.4    +  0    ] ;  = 1
# ———————————————————————————————
# STRUCTURE POUR LE FACTORIAL HMM
# ———————————————————————————————
# #          p1          (1-p1)          (1-p2)        p2
# Gamma.1 = [prob.i(1,1) 1-prob.i(1,1) ; 1-prob.i(1,2) prob.i(1,2)];
# #          q1          (1-q1)          (1-q2)        q2
# Gamma.2 = [prob.i(1,3) 1-prob.i(1,3) ; 1-prob.i(1,4) prob.i(1,4)];
# # Produit de Kronecker entre les 2 matrices
# Gamma=kron(Gamma.1,Gamma.2);
# —————————————————————————
# ————————————————————————————————————————————————————————
Gamma.Build.HHMM = function(prob.i, type){
    
    if (type=="HHMM"){
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
        
        Gamma.Build.HHMM = matrix(c(a2.11+e1*a1.11*pi.1,  a2.12+e1*a1.11*pi.2,          e1*a1.12*pi.3,        e1*a1.12*pi.4,
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
        Gamma.Build.HHMM=kronecker(Gamma.1,Gamma.2);
    }
    
    
}


# ————————————————————————————————————————————————————————
# OBTENTION PARAMÈTRES INITIAUX construction des vecteurs de valeurs
# initiales pour mu et sigma
# ————————————————————————————————————————————————————————
#       Utilisation du vecteur de données, mis en ordre puis segmenté en
#       nbRegime vecteurs. Ce sont les moyennes et écarts types de ces vecteurs
#       qui initialise les paramètres du modèle
# ————————————————————————————————————————————————————————
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



# ————————————————————————————————————————————————————————————————————————————————————
# ////////////////////////////////////////////////////////////////////
# FILTERING ALGORITHMS
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# ————————————————————————————————————————————————————————————————————————————————————

# ————————————————————————————————————————————————————————
# HAMILTON FILTER Implementation du filtre d'Hamilton
#    Permet le calcul des probabilités filtrées : omega_t = P[C_t | X_1:t]
# ————————————————————————————————————————————————————————
normal.HMM.HamiltonFilter = function(mu,
                                     sigma,
                                     Gamma,
                                     initial.Distribution=NULL,
                                     data,
                                     distribution="Normal",
                                     nu=4){
    
    nbRegime <- length(mu)
    n <- length(data)
    
    if (is.null(initial.Distribution)){
        initial.Distribution <- solve(t(diag(nbRegime)-Gamma+1),rep(1,nbRegime))
    }
    
    # —————————————————————————————————————————————————
    # HAMILTON FORWARD FILTERING
    # —————————————————————————————————————————————————
    p.ct.x1t <- p.ct.x1tm1 <- matrix(NA,nbRegime,n)
    
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
        # if (sum(is.na(omega.t))>0){
        #     cat("i = ",i,"\n ---- \n",omega.t,'\n')
        # }
        
        p.ct.x1tm1[,i] <- omega.t%*%Gamma
        
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
                p.ct.x1t=p.ct.x1t))
}



# ————————————————————————————————————————————————————————
# KIM FILTER - ALGORITHME DE LISSAGE Implementation du filtre de Kim
#    Permet le calcul des probabilités lisées : phi_t = P[C_t | X_1:T]
# 
# Either the filtered probs are given in the function call or the 
# parameters mu and sigma required to calculate these probabilities
# ————————————————————————————————————————————————————————
normal.HMM.KimFilter = function(data,
                                type,
                                distribution,
                                mu=NULL,
                                sigma=NULL,
                                Gamma,
                                p.ct.x1t=NULL,
                                initial.Distribution=NULL,
                                nu=nu){
    
    # Testing if the parameters are given or if 
    # the chain of filtered probabilities is given
    # Gamma must be specified either way
    if ((is.null(mu) || is.null(sigma)) && is.null(p.ct.x1t)){
        stop(cat("mu =",mu,'\nsigma =',sigma,'\np.ct.x1t =',p.ct.x1t,'\n'))
    }
    
    if (type=="HHMM"){
        Gamma <- Gamma.Build.HHMM(prob.i=Gamma,
                                  type=type)
    }
    
    n <- length(data)
    nbRegime <- length(mu)
    
    # If the initial distribution P[C_1] isn't specified
    # the stationary distribution of Gamma is used
    if (is.null(initial.Distribution)){
        initial.Distribution <- solve(t(diag(nbRegime)-Gamma+1),rep(1,nbRegime))
    }
    
    # Calcul des probabilités filtrées par l'algorithme du filtre d'Hamilton
    if (is.null(p.ct.x1t)){
        Hamilton.Filter <- normal.HMM.HamiltonFilter(mu=mu,
                                  sigma=sigma,
                                  Gamma=Gamma,
                                  initial.Distribution=initial.Distribution,
                                  data=data,
                                  distribution=distribution,
                                  nu=nu)
        p.ct.x1t <- Hamilton.Filter$p.ct.x1t
    }
    
    # Initialisation de la matrice de probabilités lissées
    p.ct.x1T <- matrix(NA,nbRegime,n)
    
    # Probabilite de la derniere observation de la variable latente
    p.ct.x1T[,n] <- p.ct.x1t[,n]
    
    # Algorithme de lissage (Filtre de Kim)
    for (t in (n-1):1){
        phi.t.ij <- matrix(NA,nbRegime,nbRegime)
        
        for (i in 1:nbRegime){
            for (j in 1:nbRegime){
                # Calcul des probabilités conjointes
                phi.t.ij[i,j] <- p.ct.x1T[j,t+1]*(p.ct.x1t[i,t]*Gamma[i,j]/sum(p.ct.x1t[,t]*Gamma[,j]))
            }
        }
        p.ct.x1T[,t] <- rowSums(phi.t.ij)
    }
    
    return(list(p.ct.x1T=p.ct.x1T,
                p.ct.x1t=p.ct.x1t))
}

# ————————————————————————————————————————————————————————————————————————————————————
# ////////////////////////////////////////////////////////////////////
# OPTIMIZATION FUNCTIONS
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# ————————————————————————————————————————————————————————————————————————————————————

# ————————————————————————————————————————————————————————
# FONCTION OBJECTIF afin de calculer la vraisemblance et la log-vraisemblance 
#     Le calcul est effectué à l'aide du filtre d'hamilton (Probabilités
#     filtrées P[C_t | X_1:t])
# ————————————————————————————————————————————————————————
normal.HMM.mllk = function(parvect, data, type, distribution, nu){
    n <- length(data)
    
    natural.par <- normal.HMM.W2N(parvect=parvect,
                                  type=type)
    mu <- natural.par$mu
    sigma <- natural.par$sigma
    Gamma <- natural.par$Gamma
    
    if (type=="HHMM"){
        Gamma <- Gamma.Build.HHMM(prob.i=Gamma,
                                  type=type)
    }
    
    optim.results <- normal.HMM.HamiltonFilter(mu=mu,
                                               sigma=sigma,
                                               Gamma=Gamma,
                                               data=data,
                                               distribution=distribution,
                                               nu=nu)
    # print(optim.results$llk)
    llk <- -optim.results$llk
    return(llk)
}


# ————————————————————————————————————————————————————————
#FONCTION OPTIMISATION Fonction afin d'optimiser numériquement la log-vraisemblance
#   Optimisation à l'aide du filtre d'Hamilton (Probabilités filtrées 
#   P[C_t | X_1:t])
# ————————————————————————————————————————————————————————
normal.HMM.mle <- function(data, mu, sigma, Gamma, type, distribution, nu){
    nbRegime <- length(mu)
    parvect0       <- normal.HMM.N2W(mu, sigma, Gamma, type)
    
    #Calcul du EVM
    start_time <- Sys.time()
    mod            <- nlm(normal.HMM.mllk, 
                          parvect0, 
                          data=data, 
                          type=type, 
                          distribution=distribution, 
                          nu=nu,
                          iterlim=1000,
                          gradtol = 1e-20,
                          hessian=TRUE)
    print(mod)
    natural.par.optim <- normal.HMM.W2N(mod$estimate,type)
    end_time <- Sys.time()
    cat("Time for optimzation :",end_time - start_time,"\n ----------------------------- \n")
    
    #Log-vraisemblance
    mllk           <- -mod$minimum
    np             <- length(parvect0)
    AIC            <- 2*(-mllk+np)
    n              <- sum(!is.na(data))
    BIC            <- np*log(n)-2*mllk
    list(mu=natural.par.optim$mu, sigma=natural.par.optim$sigma, Gamma=natural.par.optim$Gamma, mllk=mllk,
         AIC=AIC, BIC=BIC)
}



# ————————————————————————————————————————————————————————————————————————————————————
# ////////////////////////////////////////////////////////////////////
# TRAINING THE MODEL
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# ————————————————————————————————————————————————————————————————————————————————————

# ————————————————————————————————————————————————————————
# FULL HMM EXPERIENCE : TRAINING FOR THE PARAMETERS
# ————————————————————————————————————————————————————————
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
                     nu=4,
                     nbRegime=4,
                     type="HMM",
                     auto.assign=FALSE,
                     data.Origin="MATLAB",
                     PlotName = "test",
					 path="",
					 nbTicks=20){
    
    if (type=="HHMM"){
        nbRegime <- 4
    }
    
    n <- length(data[,1])
    # --------------------------------------------
    # Data manipulation
    # --------------------------------------------
    # Date bounds for training
    (if (is.null(start.date)){start.date <- data[1,1]})
    (if (is.null(start.date)){end.date <- data[n,1]})
    
    if (data.Origin=="MATLAB"){
        # converting the dates from MATLAB's format to R's
        dates <- as.Date(as.POSIXct((data[,1] - 719529)*86400, origin = "1970-01-01", format="%Y-%m-%d"))
        print("MATLAB")
    }else{
        dates <- as.Date(as.character(Data[,1]))
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
    
    # dim()
    ts.plot.dev(y=logR, xlabels=dates)
    
    
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
    # TESTING CORRECT PARAMETER TRANSFORMATION
    # --------------------------------------------
    parvect <- normal.HMM.N2W(mu0,sigma0,Gamma0,type)
    natpar <- normal.HMM.W2N(parvect,type)
    
    
    
    # --------------------------------------------
    # TRAINING FOR OPTIMAL PARAMETERS (HAMILTON)
    # --------------------------------------------
    
    start_time <- Sys.time()
    HMM.Train <- normal.HMM.mle(data=logR, 
                   mu=mu0,
                   sigma=sigma0,
                   Gamma=Gamma0, 
                   type=type, 
                   distribution=distribution, 
                   nu=nu)
    end_time <- Sys.time()
    cat("Time for full estimation :",end_time - start_time,"\n ----------------------------- \n")
    
    
    print(HMM.Train)
    
    # If the model is the basic HMM, the parameters are ordered in decreasing
    # order of the volatility values
    if (type=="HMM"){
        # Getting sigma order
        augmented.sigma <- cbind(HMM.Train$sigma,c(1:nbRegime))
        order.matrix <- augmented.sigma[order(augmented.sigma[1,], decreasing = TRUE),]
        order.vector <- order.matrix[,2]
        
        # Ordering the parameters in decreasing order of volatility
        mu.F <- HMM.Train$mu[order.vector]
        sigma.F <- HMM.Train$sigma[order.vector]
        Gamma.F <- HMM.Train$Gamma[order.vector,order.vector]
        
        # Reassignment of the parameters for the returned parameters in 
        # HMM.Train to match the ordered values passed in the other functions
        HMM.Train$mu <- mu.F
        HMM.Train$sigma <- sigma.F
        HMM.Train$Gamma <- Gamma.F
    }
    
    # --------------------------------------------
    # KIM FILTER -> SMOOTHED PROBABILITIES P[C_t|x_1:T]
    # --------------------------------------------
    
    smooth.Prob <- normal.HMM.KimFilter(data=logR,
                                        type=type,
                                        distribution=distribution,
                                        mu=HMM.Train$mu,
                                        sigma=HMM.Train$sigma,
                                        Gamma=HMM.Train$Gamma,
                                        nu=nu)
    
    # Uniform color:
    # barplot(my_vector, col=rgb(0.2,0.4,0.6,0.6) )
    
    # HMM.Stack.Plot(smooth.Prob,t_[-1])
    # --------------------------------------------
    
    Grid.Plot.HMM.Full(timeseries=logR,
    				   dates=dates,
    				   state.Probs=smooth.Prob$p.ct.x1T,
    				   mu=HMM.Train$mu,
    				   sigma=HMM.Train$sigma,
    				   Gamma=HMM.Train$Gamma,
    				   name=PlotName,
    				   path=path,
    				   nbTicks=nbTicks)
    	
    return(list(HMM.Train=HMM.Train,
                smooth.Prob=smooth.Prob$p.ct.x1T,
                filtered.Prob=smooth.Prob$p.ct.x1t,
                dates.char=dates.Char[-1],
                dates=dates))

}



# ————————————————————————————————————————————————————————————————————————————————————
# ////////////////////////////////////////////////////////////////////
# GRAPHICAL OUTPUT FUNCITONS
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# ————————————————————————————————————————————————————————————————————————————————————


# ————————————————————————————————————————————————————————
# GRAPHICAL PRINT OF TIME-SERIES
# ————————————————————————————————————————————————————————
ts.plot.dev = function(y,nbTicks=30, xlabels = NULL,Colors, axis.label.size=1.5){
	par(mar=c(10,8,6,6),cex.axis=axis.label.size,cex.lab=2,las=2)
	# Getting length of the time series
	n <- length(y)
	
	# Spacing between the ticks on the x axis
	xTicksSpace <- n/nbTicks
	
	plot(y, type='n',xaxt="n",xlab="", ylab=deparse(substitute(y)),cex.main=axis.label.size, xlim=c(1,n),ylim=c(min(y)-1,max(y)+1))
	lines(y, col="black")
	
	# Adding the correct ticks on the x axis
	axis(side=1,at=c(seq(from=1, to=n, by=xTicksSpace),n),labels = xlabels[c(seq(from=1, to=n, by=xTicksSpace),n)],lwd=4,lwd.ticks=2,col.ticks="red",cex.axis=axis.label.size)
	
	# Adding grid to the plot
	grid(lty=1, col=gray(.9))
	# Adding an horizontal line at y=0
	abline(h=0, col="blue")
}


# ————————————————————————————————————————————————————————
# PROBABILITIES FILL COLOR
# ————————————————————————————————————————————————————————
HMM.Stack.Plot = function(series,xlabels,nbTicks=30,Colors, axis.label.size=1.5){
    
    par(mar=c(10,8,6,6),cex.axis=axis.label.size,cex.lab=2,las=2)
    
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
    
    plot(1:n,stacked.series[1,],type = 'n',xaxt="n", ylim = c(0,1),
         ylab = 'probabilities', xlab = '', col = "red", xlim=c(1,n))
    
    for (i in 2:nbRegime+1){
        lines(stacked.series[i,], col = "black")
    }
    
    # Adding the correct ticks on the x axis
    axis(side=1,at=seq.labels,labels = xlabels[seq.labels],lwd=4,lwd.ticks=2,col.ticks="red",cex.axis=axis.label.size)
    
    for (i in 1:nbRegime){
        polygon(c(1:n, rev(1:n)), c(stacked.series[i+1,], rev(stacked.series[i,])),
                col = Colors[(nbRegime-1),i], border = NA)
    }
    # print(t(stacked.series))
}




Grid.Plot.HMM.Full = function(timeseries, dates, state.Probs, type="HMM", mu, sigma, Gamma, name,path,resolution=300,nbTicks=20){
    # Colors
    custom.Colors <- matrix(c("firebrick3","blue4",NA,NA,
                              "firebrick3","deepskyblue1","blue4",NA,
                              "firebrick3","darksalmon","deepskyblue1","blue4"),ncol=4,byrow=T)
    
    # Title size
    title.size <- 3
    axis.label.size <- 2.5
    
    # If the model is Hierarchical, the Gamma matrix is built here
    if (type=="HHMM"){
        Gamma <- Gamma.Build.HHMM(prob.i=Gamma,
                                  type=type)
    }
    
    # Getting dates in correct format while keeping only year and month
    dates.split <- matrix(unlist(strsplit(dates,split="\\-")),ncol=3,byrow=T)
    dates.Plot <- paste(month.abb[as.numeric(dates.split[,2])],dates.split[,1],sep='-')
    
	# Getting the state with the maximum probability
	max.Prob.States <- colMax(state.Probs)
	
	# Creating the empty image to be filled with plots
	# bitmap(paste(path,'/2 - Graphiques/',name,".png",sep=""), res=resolution)
	jpeg(paste(path,'/2 - Graphiques/',name,".jpg",sep=""), width = 2500, height = 1500) 
	
	# Defines the layout of the figure output
	layout(mat = matrix(c(1,1,2,2,3,4,5,6),ncol=2,nrow = 4,byrow=TRUE))
	
	print("Time series printing")
	ts.plot.dev(y=timeseries, 
	            xlabels=dates.Plot,
				nbTicks=nbTicks,
				Colors=custom.Colors,
				axis.label.size=axis.label.size)
	title(main = paste(name, "series - ",Sys.time()),cex.main=title.size)
	
	print("Stacking")
	HMM.Stack.Plot(series=state.Probs,
	               xlabels=dates.Plot,
				   nbTicks=nbTicks,
				   Colors=custom.Colors,
				   axis.label.size=axis.label.size)
	title(main = "Smooth probabilities - P[C_t | X_1:T]",cex.main=title.size)
	
	# Gamma.df <- data.frame(Gamma,row.names = paste("state.",as.character(1:nbRegime),sep=""))
	grid.table(Gamma)
	
	dev.off()
}


# test <- HMM.Train(index="SPXIndex",
#                   Data,
#                   start.date,
#                   end.date,
#                   frequency=5,
#                   mult=52,
#                   Gamma0=Gamma.HMM.3,
#                   nbRegime=3,
#                   type="HMM",
#                   distribution='Normal',
#                   PlotName = "test1",
#                   data.Origin="MATLAB",
#                   auto.assign=T,
#                   path=path,
#                   nbTicks=50)



