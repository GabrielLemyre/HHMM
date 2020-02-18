# ----------------------------------------------------------------------------------------------------
# Hidden Markov Models
# Hamilton Filtering
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