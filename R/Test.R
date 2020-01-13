# ----------------------------------------------------------------------------------------------------
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
# Last version : january 13th, 2020
# ----------------------------------------------------------------------------------------------------
options(max.print=1000000)

# -------------------------------------------------------
# 
# 
# -------------------------------------------------------
# Sourcing private libraries
# set Sourcing directory
# -------------------------------------------------------
path <- '~/Documents/GitHub/HHMM'
setwd(path.expand(path)) # Setting Sourcing path
source("R/HMM.R")


# -------------------------------------------------------
# 
# 
# -------------------------------------------------------
# Initialisation des variables et tableaux de résultats
# -------------------------------------------------------
Analyse <- NULL
GraphPath <- paste(path,"R/1 - Test Graphs",sep="/")

# -------------------------------------------------------
# 
# 
# -------------------------------------------------------
# Préparation des données
# -------------------------------------------------------
data.Origin="R"
nbTicks=30

#load data
# Data <- read.csv("1 - Data/Table_Daily.csv") # MATLAB
Data <- read.csv("1 - Data/DATA_INDEX.csv") # R
head(Data)
data.freq <- "daily"
mult <- 100          #multiply all returns by -mult-

index <- "SP500"
start.date <- as.Date("1999-12-31")
end.date   <- as.Date("2016-12-30")
#weekdays(start.date); weekdays(end.date);

dates <- as.Date(as.character(Data[,1]))
st <- which(dates==start.date)
en <- which(dates==end.date)
dates <- dates[st:en]
Pt <- as.numeric(as.character(Data[st:en,index]))
names(Pt) <- dates
Pt <- Pt[!is.na(Pt)]
dates <- names(Pt)

logR <- log( Pt[-1] / Pt[-length(Pt)] )
logR <- logR * mult
#the length of logR is one less than Pt

#convert dates into a numeric vector (useful for figures)
t_ <- as.POSIXlt(c(as.character(start.date),names(logR)), format="%Y-%m-%d")
n.days <- rep(365, length(t_))
n.days[((1900+t_$year) %% 4)==0] <- 366 #leap years
t_  <- 1900 + t_$year + (t_$yday+1)/n.days
#when yday=0, we are jan. 1st
#as.POSIXlt("2012-01-01",format="%Y-%m-%d")$yday

# no factors in a data.frame : stringsAsFactors=F

#descriptive statistics
des_stats <- c("Mean" = mean(logR),
               "StDev" = sd(logR),
               "Skewness" = timeDate::skewness(logR, method="moment"),
               "Kurtosis" = timeDate::kurtosis(logR, method="moment"), "Minimum" = min(logR),
               "Maximum" = max(logR), "n" = length(logR)
)
des_stats

# ————————————————————————————————————————————————————————————————————————————————————
# ////////////////////////////////////////////////////////////////////
# TESTING THE CODE
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# -———————————————————————————————————————————————————————————————————————————————————


mu.Test.2 <- c(0.08,-0.1)
sigma.Test.2 <- c(2,0.6)
Gamma.HMM.2=matrix(c(0.98, 0.02,
                     0.02, 0.98),
                   byrow=T,nrow=2)


mu.Test.3 <- c(-0.12771163,  0.08938697, -0.03075640)
sigma.Test.3 <- c(2.6139977, 0.5848203, 1.1616007)
Gamma.HMM.3=matrix(c(9.640156e-01, 1.956318e-30, 0.03598444,
                     3.039625e-11, 9.777013e-01, 0.02229874,
                     4.921709e-03, 2.455015e-02, 0.97052814),
                   byrow=T,nrow=3)


mu.Test.4 <- c(-0.06,  0.09, -0.002, -0.29)
sigma.Test.4 <- c(1.59, 0.54, 0.99, 3.54)
Gamma.HMM.4=matrix(c(9.791353e-01, 2.695902e-07, 1.527788e-02, 5.586534e-03,
                     6.668266e-11, 9.682687e-01, 3.173053e-02, 7.292526e-07,
                     9.760572e-03, 3.632286e-02, 9.532491e-01, 6.674701e-04,
                     3.702408e-02, 1.075279e-09, 1.519283e-09, 9.629759e-01),
                   byrow=T,nrow=4)

# Gamma.HMM.4 <- matrix(1/300,ncol=4,nrow=4)
# diag(Gamma.HMM.4) <- 0.99

Gamma.HHMM <- matrix(c(0.95  ,  0.04  , 0.01 ,
                       0.95  ,  0.04 ,  0.01 ,
                       0.95   , 0.04 ,  0.01 ,
                       0.95  ,  0.04  , 0.01 ,
                       0.85 ,   0.15  , 0    ,
                       0.85  ,  0.15  , 0  ,
                       0.5  ,   0.5  ,  0  ,
                       0.5  ,   0.5 ,   0 ),ncol=3,byrow=T)

Gamma.HMM.Diebold <- matrix(c( 2,-4.928700,  # beta_0.0, beta_0.1
                               1, 7.773485), # beta_1.0, beta_1.1
                            nrow=2, byrow=T)

Gamma.HMM.Diebold.w.filter <- matrix(c(2, -5, 2, # beta_0.0, beta_0.1, beta_0.2
                                       1,  1, 6), # beta_1.0, beta_1.1, beta_1.2
                                     nrow=2, byrow=T)

Gamma.Test=matrix(c(0.90229999, 0.0977, 0, 0.00000001,
                    0.0096, 0.9778, 0, 0.0126,
                    0.0109, 0.0093, 0.8739, 0.1059,
                    0, 0.0005, 0.0947, 0.9048),
                  byrow=T,nrow=4)




# ————————————————————————————————————————————————————————————————————————————————————
# ////////////////////////////////////////////////////////////////////
# TRANSFORMATION DES PARAMÈTRES
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# ————————————————————————————————————————————————————————————————————————————————————

parvect <- normal.HMM.N2W(mu.Test[c(1,2)],sigma.Test[c(1,2)],Gamma.HMM.Diebold,type="HMM",Transition.Type="Diebold")
natpar <- normal.HMM.W2N(parvect,type="HMM",Transition.Type="Diebold")







# ————————————————————————————————————————————————————————————————————————————————————
# ////////////////////////////////////////////////////////////////////
# STRAIGHT FORWARD HMM WITH K=2,3,4
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# ————————————————————————————————————————————————————————————————————————————————————

test.2 <- HMM.Train(index <- index,
                    logR,
                    mu0=mu.Test.2,
                    sigma0=sigma.Test.2,
                    Gamma0=Gamma.HMM.2,
                    type="HMM",
                    distribution='Normal',
                    Transition.Type="Homogeneous",
                    PlotName = NULL,
                    data.Origin=data.Origin,
                    auto.assign=F,
                    path=GraphPath,
                    nbTicks=nbTicks)
real.mllk.2 <- -6190.254


test.3 <- HMM.Train(index <- index,
                    logR,
                    mu0=mu.Test.3,
                    sigma0=sigma.Test.3,
                    Gamma0=Gamma.HMM.3,
                    type="HMM",
                    distribution='Normal',
                    Transition.Type="Homogeneous",
                    PlotName = NULL,
                    data.Origin=data.Origin,
                    auto.assign=F,
                    path=GraphPath,
                    nbTicks=nbTicks)
real.mllk.3 <- -6020.947


test.4 <- HMM.Train(index <- index,
                    logR,
                    mu0=mu.Test.4,
                    sigma0=sigma.Test.4,
                    Gamma0=Gamma.HMM.4,
                    type="HMM",
                    distribution='Normal',
                    Transition.Type="Homogeneous",
                    PlotName = NULL,
                    data.Origin=data.Origin,
                    auto.assign=F,
                    path=GraphPath,
                    nbTicks=nbTicks)
real.mllk.4 <- -5970.916




# ————————————————————————————————————————————————————————————————————————————————————
# ////////////////////////////////////////////////////////////////////
# HHMM models
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# ————————————————————————————————————————————————————————————————————————————————————
test.HHMM <- HMM.Train(index <- index,
                       logR,
                       mu0=mu.Test.4,
                       sigma0=sigma.Test.4,
                       Gamma0=Gamma.HHMM,
                       type="HHMM",
                       distribution='Normal',
                       PlotName = NULL,
                       data.Origin=data.Origin,
                       auto.assign=F,
                       path=GraphPath,
                       nbTicks=nbTicks)




# ————————————————————————————————————————————————————————————————————————————————————
# ////////////////////////////////////////////////////////////////////
# HMM - DIEBOLD
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# ————————————————————————————————————————————————————————————————————————————————————
test.Diebold.1.0 <- HMM.Train(index <- index,
                              logR,
                              mu0=mu.Test.2,
                              sigma0=sigma.Test.2,
                              Gamma0=Gamma.HMM.Diebold,
                              type="HMM",
                              Transition.Type="Diebold",
                              distribution='Normal',
                              PlotName = NULL,
                              data.Origin=data.Origin,
                              auto.assign=F,
                              path=GraphPath,
                              nbTicks=nbTicks,
                              initial.Distribution = c(1,0))

test.Diebold.0.1 <- HMM.Train(index <- index,
                              logR,
                              mu0=mu.Test.2,
                              sigma0=sigma.Test.2,
                              Gamma0=Gamma.HMM.Diebold,
                              type="HMM",
                              Transition.Type="Diebold",
                              distribution='Normal',
                              PlotName = NULL,
                              data.Origin=data.Origin,
                              auto.assign=F,
                              path=GraphPath,
                              nbTicks=nbTicks,
                              initial.Distribution = c(0,1))

test.Diebold.w.Filter.1.0 <- HMM.Train(index <- index,
                                       logR,
                                       mu0=mu.Test.2,
                                       sigma0=sigma.Test.2,
                                       Gamma0=Gamma.HMM.Diebold.w.filter,
                                       type="HMM",
                                       Transition.Type="Diebold.w.filter",
                                       distribution='Normal',
                                       PlotName = NULL,
                                       data.Origin=data.Origin,
                                       auto.assign=F,
                                       path=GraphPath,
                                       nbTicks=nbTicks,
                                       initial.Distribution = c(1,0))

test.Diebold.w.Filter.0.1 <- HMM.Train(index <- index,
                                       logR,
                                       mu0=mu.Test.2,
                                       sigma0=sigma.Test.2,
                                       Gamma0=Gamma.HMM.Diebold.w.filter,
                                       type="HMM",
                                       Transition.Type="Diebold.w.filter",
                                       distribution='Normal',
                                       PlotName = NULL,
                                       data.Origin=data.Origin,
                                       auto.assign=F,
                                       path=GraphPath,
                                       nbTicks=nbTicks,
                                       initial.Distribution = c(0,1))

# ————————————————————————————————————————————————————————————————————————————————————
# ////////////////////////////////////////////////////////////////////
# TEST RESULTS
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# ————————————————————————————————————————————————————————————————————————————————————
Analyse <- NULL

Analyse <- rbind(Analyse,
                 c("HMM-2",              
                   s(test.2$HMM.Train$mllk,3),                        
                   s(test.2$HMM.Train$AIC,3),          
                   s(test.2$HMM.Train$BIC,3),                 
                   paste(s(abs((test.2$HMM.Train$mllk-real.mllk.2)/real.mllk.2)*100,6),"%",sep=""),
                   test.2$timeSpent.String))         # Time to run Full model

Analyse <- rbind(Analyse,
                 c("HMM-3",    
                   s(test.3$HMM.Train$mllk,3),               
                   s(test.3$HMM.Train$AIC,3),          
                   s(test.3$HMM.Train$BIC,3),   
                   paste(s(abs((test.3$HMM.Train$mllk-real.mllk.3)/real.mllk.3)*100,6),"%",sep=""),
                   test.3$timeSpent.String))         # Time to run Full model

Analyse <- rbind(Analyse,
                 c("HMM-4",               # Methode
                   s(test.4$HMM.Train$mllk,3), 
                   s(test.4$HMM.Train$AIC,3),
                   s(test.4$HMM.Train$BIC,3), 
                   paste(s(abs((test.4$HMM.Train$mllk-real.mllk.4)/real.mllk.4)*100,6),"%",sep=""),
                   test.4$timeSpent.String))         # Time to run Full model

Analyse <- rbind(Analyse,
                 c("HMM-Diebold (1,0)",               # Methode
                   s(test.Diebold.1.0$HMM.Train$mllk,3), 
                   s(test.Diebold.1.0$HMM.Train$AIC,3),
                   s(test.Diebold.1.0$HMM.Train$BIC,3), 
                   "-",
                   test.Diebold.1.0$timeSpent.String))         # Time to run Full model

Analyse <- rbind(Analyse,
                 c("HMM-Diebold (0,1)",               # Methode
                   s(test.Diebold.0.1$HMM.Train$mllk,3), 
                   s(test.Diebold.0.1$HMM.Train$AIC,3),
                   s(test.Diebold.0.1$HMM.Train$BIC,3), 
                   "-",
                   test.Diebold.0.1$timeSpent.String))         # Time to run Full model

Analyse <- rbind(Analyse,
                 c("HMM-Diebold.w.Filter (1,0)",               # Methode
                   s(test.Diebold.w.Filter.1.0$HMM.Train$mllk,3), 
                   s(test.Diebold.w.Filter.1.0$HMM.Train$AIC,3),
                   s(test.Diebold.w.Filter.1.0$HMM.Train$BIC,3), 
                   "-",
                   test.Diebold.w.Filter.1.0$timeSpent.String))         # Time to run Full model

Analyse <- rbind(Analyse,
                 c("HMM-Diebold.w.Filter (0,1)",               # Methode
                   s(test.Diebold.w.Filter.0.1$HMM.Train$mllk,3), 
                   s(test.Diebold.w.Filter.0.1$HMM.Train$AIC,3),
                   s(test.Diebold.w.Filter.0.1$HMM.Train$BIC,3), 
                   "-",
                   test.Diebold.w.Filter.0.1$timeSpent.String))         # Time to run Full model

Analyse <- rbind(Analyse,
                 c("HHMM",              
                   s(test.HHMM$HMM.Train$mllk,3),                        
                   s(test.HHMM$HMM.Train$AIC,3),          
                   s(test.HHMM$HMM.Train$BIC,3),                 
                   "-",
                   test.HHMM$timeSpent.String))         # Time to run Full model

colnames(Analyse) <- c('Modèle', "mllk", "AIC", "BIC", "Erreur", "Temps")
as.data.frame(Analyse)


