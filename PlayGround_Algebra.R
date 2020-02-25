# K <- 3
# vecun <- rep(1,K)
# (M <- matrix(seq(1,K,1),ncol=K))
# temp <- matrix(NA, K,K)
# (b <- lower.tri(temp, diag = TRUE))
# 
# (triangulaire.inf <- vecun%*%M * b)
# 
# M%*%triangulaire.inf%*%vecun


# mu.t <- matrix(c(-0.14491789,  0.01172504,  0.05939242),ncol=1)
mu.t <- c(-0.14491789,  0.01172504,  0.05939242)
sigma.g.test <- c(2.5575702, 1.0183486, 0.5295502)

gamma.g.test <- matrix(c(9.557654e-01, 0.04423460, 2.404844e-30,
                         6.046790e-03, 0.97682473, 1.712848e-02,
                         2.103557e-11, 0.01543173, 9.845683e-01),nrow=3,byrow=T)

test <- Diagnosis(gamma.g.test, 
                  additional.info = rbind(mu.t,sigma.g.test), 
                  special.Weights=sigma.g.test,
                  network.print = FALSE,
                  docPath=path)

dist.stat.g <- test$dist.asympt

variance <- function(u.t,mu.g,sigma.g){
  K <- length(mu.g) # Nombre de régimes
  
  vecun <- matrix(rep(1,K),ncol=1) # Vecteur de 1
  
  Omega <- matrix(u.t*t(mu.t),nrow=1,byrow=T) # Vecteur ligne des moyennes conditionnelles pondérées
  
  # ————————————————————————————————
  part.1 <- u.t %*%diag(sigma.g)^2 %*%vecun
  cat("part.1 =",part.1,"\n")
  # ————————————————————————————————
  part.2 <- u.t %*%diag(mu.g)^2%*%(vecun-t(u.t))
  cat("part.2 =",part.2,"\n")
  # ————————————————————————————————
  temp <- matrix(NA, K,K)
  triSupNoDiag <- upper.tri(temp, diag = FALSE)
  
  part.3 <- Omega%*%((vecun%*%Omega)*triSupNoDiag)%*%vecun
  cat("part.3 =",part.3,"\n")
  # ————————————————————————————————
    
  var <- part.1 + part.2 - 2 * part.3
}

var.Memoire.Gabriel <- variance(u.t=dist.stat.g,
                mu.g=mu.t,
                sigma.g.test)
var.Memoire.Gabriel

var.Notes.Maciej <- dist.stat.g%*%diag(mu.t)%*%mu.t - (dist.stat.g%*%mu.t)^2
var.Notes.Maciej

M <- u.t*t(mu.t)
K <- length(u.t)
vecun <- rep(1,K)

temp <- matrix(NA, K,K)
(b <- lower.tri(temp, diag = TRUE))


First.part <- (M%*%vecun)^2

(triangulaire.inf <- vecun%*%M * b)
Second.part <- M%*%triangulaire.inf%*%vecun

(Result.A <- First.part-Second.part)

# Au long

Result.Sum <- function(u,mu){
  K <- length(u)
  res <- 0
  for (i in 2:K){
    for (j in 1:(i-1))
      # print(paste("i =",i," j =",j))
      res <- res+mu[i]*u[i]*mu[j]*u[j]
  }
  return(res)
}

(Result.S <- Result.Sum(u.t,mu.t))

cat(paste("Algebra :",Result.A,"\nSummation :",Result.S))

