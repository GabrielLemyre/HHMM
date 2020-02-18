# ----------------------------------------------------------------------------------------------------
# Hidden Markov Models
# Constraining parameters
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