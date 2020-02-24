# K <- 3
# vecun <- rep(1,K)
# (M <- matrix(seq(1,K,1),ncol=K))
# temp <- matrix(NA, K,K)
# (b <- lower.tri(temp, diag = TRUE))
# 
# (triangulaire.inf <- vecun%*%M * b)
# 
# M%*%triangulaire.inf%*%vecun


u.t <- matrix(c(0.1,0.4,0.5),nrow=1)
mu.t <- matrix(c(3,7,2),ncol=1)

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