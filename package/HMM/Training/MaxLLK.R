# ----------------------------------------------------------------------------------------------------
# Hidden Markov Models
# Maximum Likelihood computation
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