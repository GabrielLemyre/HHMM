# ----------------------------------------------------------------------------------------------------
# Hidden Markov Models
# Moments of the Markov Chain
# ----------------------------------------------------------------------------------------------------
# written
# Gabriel LEMYRE
# ----------------------------------------------------------------------------------------------------
# Under the supervision of :
# Maciej AUGUSTYNIAK
# ----------------------------------------------------------------------------------------------------
# Last version : january 10th, 2020
# Last version : january 13th, 2020
# ----------------------------------------------------------------------------------------------------

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

Gamma.test=matrix(c(9.640156e-01, 1.956318e-30, 0.03598444,
                     3.039625e-11, 9.777013e-01, 0.02229874,
                     4.921709e-03, 2.455015e-02, 0.97052814),
                   byrow=T,nrow=3)

test <- MomentsMC(Gamma.test,initial.Distribution=c(0,0,1),nbStepsEsp=4000,nbStepsACF=100)
# test