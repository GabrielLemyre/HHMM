# ----------------------------------------------------------------------------------------------------
# Hidden Markov Models
# Kim Filtering
# ----------------------------------------------------------------------------------------------------
# written
# Gabriel LEMYRE
# ----------------------------------------------------------------------------------------------------
# Under the supervision of :
# Maciej AUGUSTYNIAK
# ----------------------------------------------------------------------------------------------------
# Last version : january 9th, 2020
# Last version : january 9th, 2020
# ----------------------------------------------------------------------------------------------------

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