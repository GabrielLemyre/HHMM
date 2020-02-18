# ----------------------------------------------------------------------------------------------------
# Hidden Markov Models
# Unconstraining parameters
# ----------------------------------------------------------------------------------------------------
# written
# Gabriel LEMYRE
# ----------------------------------------------------------------------------------------------------
# Under the supervision of :
# Maciej AUGUSTYNIAK
# ----------------------------------------------------------------------------------------------------
# Last version : january 8th, 2020
# Last version : january 9th, 2020
# ----------------------------------------------------------------------------------------------------

# --------------------------------------------------------
# PARAMÈTRES CONTRAINTS -> PARAMÈTRES NON-CONTRAINTS
# --------------------------------------------------------
# Function to convert from natural parameters to their working counterparts
#     The idea is that by changing the constrained parameters to
#     unconstrained versions, the optimization can be done without
#     constraints
# --------------------------------------------------------
normal.HMM.N2W = function(mu, sigma, Gamma,
                          type, Transition.Type=NULL){
  
  if (length(mu)!=length(sigma)){
    stop("Les dimensions des vecteurs mu et sigma ne concordent pas.")
  }
  
  if ((length(mu)!=dim(Gamma)[1])&&(Transition.Type!="HHMM")){
    stop("Les dimensions de la matrice Gamma et des vecteurs mu et sigma ne concordent pas.")
  }
  
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
    # ----------------
    # HIERARCHICAL HMM
    # ----------------
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
