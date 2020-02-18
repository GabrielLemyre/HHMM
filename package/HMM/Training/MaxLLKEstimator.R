# ----------------------------------------------------------------------------------------------------
# Hidden Markov Models
# Maximum Likelihood Estimator
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
# OPTIMISATION Fonction afin d'optimiser numériquement la log-vraisemblance
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