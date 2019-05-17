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
# Last version : april 16th, 2019
# ----------------------------------------------------------------------------------------------------
options(max.print=1000000)

#set working directory
path <- '~/Documents/GitHub/HHMM'
setwd(path.expand(path)) # Setting path

source("R/HMM.R")

index <- "SPXIndex"
start.date <- as.Date("2018-01-01")
end.date   <- as.Date("2019-12-30")

frequency <- 5
mult <- 52

data.Origin="MATLAB"

nbTicks=30

Data <- read.csv("1 - Data/Table_Daily.csv") # MATLAB
# Data <- read.csv("1 - Data/DATA_INDEX.csv") # R

head(Data)
data.freq <- "daily"



# ————————————————————————————————————————————————————————————————————————————————————
# ////////////////////////////////////////////////////////////////////
# TESTING THE CODE
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# -———————————————————————————————————————————————————————————————————————————————————

Gamma.HMM.2=matrix(c(0.98, 0.02,
                     0.02, 0.98),
                   byrow=T,nrow=2)

Gamma.HMM.3=matrix(c(0.9048, 0.0952, 0.0000000001,
                     0.0122, 0.9485, 0.0393,
                     0.0000000001, 0.0429, 0.9571),
                   byrow=T,nrow=3)

Gamma.HMM.4=matrix(c(9.791353e-01, 2.695902e-07, 1.527788e-02, 5.586534e-03,
                     6.668266e-11, 9.682687e-01, 3.173053e-02, 7.292526e-07,
                     9.760572e-03, 3.632286e-02, 9.532491e-01, 6.674701e-04,
                     3.702408e-02, 1.075279e-09, 1.519283e-09, 9.629759e-01),
                   byrow=T,nrow=4)

Gamma.HMM.4 <- matrix(1/300,ncol=4,nrow=4)
diag(Gamma.HMM.4) <- 0.99

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

mu.Test <- c(-0.4882, 0.0862, 0.0170, 0.1967)
sigma.Test <- c(2.9746, 1.3410, 0.9686, 0.5169)

Gamma.Test=matrix(c(0.90229999, 0.0977, 0, 0.00000001,
                   0.0096, 0.9778, 0, 0.0126,
                   0.0109, 0.0093, 0.8739, 0.1059,
                   0, 0.0005, 0.0947, 0.9048),
                 byrow=T,nrow=4)

# parvect <- normal.HMM.N2W(mu.Test,sigma.Test,Gamma.Test,type="HMM")
# natpar <- normal.HMM.W2N(parvect,type="HMM")
#
# normal.HMM.mllk(normal.HMM.N2W(mu.Test, sigma.Test, Gamma.Test, type="HMM"),
#                 logR,
#                 type="HMM",
#                 distribution="Normal",
#                 nu=4)



# (Gamma <- Gamma.Build.HHMM(prob.i=test$HMM.Train$Gamma,
#                           type="HHMM"))


# HMM.Stack.Plot(test$smooth.Prob,
#                test$dates)
#
# HMM.Stack.Plot(test$filtered.Prob,
#                test$dates)


# test <- HMM.Train(index <- index,
#                   Data,
#                   start.date <- start.date,
#                   end.date   <- end.date,
#                   frequency=frequency,
#                   mult=mult,
#                   # mu0=mu.Test,
#                   # sigma0=sigma.Test,
#                   Gamma0=Gamma.HMM.Diebold.w.filter,
#                   type="HMM",
#                   distribution='Normal',
#                   Transition.Type="Diebold.w.filter",
#                   PlotName = NULL,
#                   data.Origin=data.Origin,
#                   auto.assign=T,
#                   path=path,
#                   nbTicks=nbTicks,
#                   initial.Distribution = c(1,0))
#
#
# test <- HMM.Train(index <- index,
#                   Data,
#                   start.date <- start.date,
#                   end.date   <- end.date,
#                   frequency=frequency,
#                   mult=mult,
#                   # mu0=mu.Test,
#                   # sigma0=sigma.Test,
#                   Gamma0=Gamma.HMM.Diebold.w.filter,
#                   type="HMM",
#                   distribution='Student',
#                   Transition.Type="Diebold.w.filter",
#                   PlotName = NULL,
#                   data.Origin=data.Origin,
#                   auto.assign=T,
#                   path=path,
#                   nbTicks=nbTicks,
#                   initial.Distribution = c(1,0))
#
#
# test <- HMM.Train(index <- index,
#                   Data,
#                   start.date <- start.date,
#                   end.date   <- end.date,
#                   frequency=frequency,
#                   mult=mult,
#                   # mu0=mu.Test,
#                   # sigma0=sigma.Test,
#                   Gamma0=Gamma.HMM.Diebold.w.filter,
#                   type="HMM",
#                   distribution='Normal',
#                   Transition.Type="Diebold.w.filter",
#                   PlotName = NULL,
#                   data.Origin=data.Origin,
#                   auto.assign=T,
#                   path=path,
#                   nbTicks=nbTicks,
#                   initial.Distribution = c(0,1))
#
#
# test <- HMM.Train(index <- index,
#                   Data,
#                   start.date <- start.date,
#                   end.date   <- end.date,
#                   frequency=frequency,
#                   mult=mult,
#                   # mu0=mu.Test,
#                   # sigma0=sigma.Test,
#                   Gamma0=Gamma.HMM.Diebold.w.filter,
#                   type="HMM",
#                   distribution='Student',
#                   Transition.Type="Diebold.w.filter",
#                   PlotName = NULL,
#                   data.Origin=data.Origin,
#                   auto.assign=T,
#                   path=path,
#                   nbTicks=nbTicks,
#                   initial.Distribution = c(0,1))
#
#
# test <- HMM.Train(index <- index,
#                   Data,
#                   start.date <- start.date,
#                   end.date   <- end.date,
#                   frequency=frequency,
#                   mult=mult,
#                   # mu0=mu.Test,
#                   # sigma0=sigma.Test,
#                   Gamma0=Gamma.HMM.Diebold,
#                   type="HMM",
#                   distribution='Normal',
#                   Transition.Type="Diebold",
#                   PlotName = NULL,
#                   data.Origin=data.Origin,
#                   auto.assign=T,
#                   path=path,
#                   nbTicks=nbTicks,
#                   initial.Distribution = c(1,0))
#
#
# test <- HMM.Train(index <- index,
#                   Data,
#                   start.date <- start.date,
#                   end.date   <- end.date,
#                   frequency=frequency,
#                   mult=mult,
#                   # mu0=mu.Test,
#                   # sigma0=sigma.Test,
#                   Gamma0=Gamma.HMM.Diebold,
#                   type="HMM",
#                   distribution='Student',
#                   Transition.Type="Diebold",
#                   PlotName = NULL,
#                   data.Origin=data.Origin,
#                   auto.assign=T,
#                   path=path,
#                   nbTicks=nbTicks,
#                   initial.Distribution = c(1,0))
#
#
# test <- HMM.Train(index <- index,
#                   Data,
#                   start.date <- start.date,
#                   end.date   <- end.date,
#                   frequency=frequency,
#                   mult=mult,
#                   # mu0=mu.Test,
#                   # sigma0=sigma.Test,
#                   Gamma0=Gamma.HMM.Diebold,
#                   type="HMM",
#                   distribution='Normal',
#                   Transition.Type="Diebold",
#                   PlotName = NULL,
#                   data.Origin=data.Origin,
#                   auto.assign=T,
#                   path=path,
#                   nbTicks=nbTicks,
#                   initial.Distribution = c(0,1))
#
#
# test <- HMM.Train(index <- index,
#                   Data,
#                   start.date <- start.date,
#                   end.date   <- end.date,
#                   frequency=frequency,
#                   mult=mult,
#                   # mu0=mu.Test,
#                   # sigma0=sigma.Test,
#                   Gamma0=Gamma.HMM.Diebold,
#                   type="HMM",
#                   distribution='Student',
#                   Transition.Type="Diebold",
#                   PlotName = NULL,
#                   data.Origin=data.Origin,
#                   auto.assign=T,
#                   path=path,
#                   nbTicks=nbTicks,
#                   initial.Distribution = c(0,1))
#
#
# test <- HMM.Train(index <- index,
#                   Data,
#                   start.date <- start.date,
#                   end.date   <- end.date,
#                   frequency=frequency,
#                   mult=mult,
#                   # mu0=mu.Test,
#                   # sigma0=sigma.Test,
#                   Gamma0=Gamma.HMM.2,
#                   type="HMM",
#                   distribution='Normal',
#                   Transition.Type="Homogeneous",
#                   PlotName = NULL,
#                   data.Origin=data.Origin,
#                   auto.assign=T,
#                   path=path,
#                   nbTicks=nbTicks)
#
#
# test <- HMM.Train(index <- index,
#                   Data,
#                   start.date <- start.date,
#                   end.date   <- end.date,
#                   frequency=frequency,
#                   mult=mult,
#                   # mu0=mu.Test,
#                   # sigma0=sigma.Test,
#                   Gamma0=Gamma.HMM.2,
#                   type="HMM",
#                   distribution='Student',
#                   Transition.Type="Homogeneous",
#                   PlotName = NULL,
#                   data.Origin=data.Origin,
#                   auto.assign=T,
#                   path=path,
#                   nbTicks=nbTicks)
#
#
# test <- HMM.Train(index <- index,
#                   Data,
#                   start.date <- start.date,
#                   end.date   <- end.date,
#                   frequency=frequency,
#                   mult=mult,
#                   # mu0=mu.Test,
#                   # sigma0=sigma.Test,
#                   Gamma0=Gamma.HMM.3,
#                   type="HMM",
#                   distribution='Normal',
#                   Transition.Type="Homogeneous",
#                   PlotName = NULL,
#                   data.Origin=data.Origin,
#                   auto.assign=T,
#                   path=path,
#                   nbTicks=nbTicks)
#
#
# test <- HMM.Train(index <- index,
#                   Data,
#                   start.date <- start.date,
#                   end.date   <- end.date,
#                   frequency=frequency,
#                   mult=mult,
#                   # mu0=mu.Test,
#                   # sigma0=sigma.Test,
#                   Gamma0=Gamma.HMM.3,
#                   type="HMM",
#                   distribution='Student',
#                   Transition.Type="Homogeneous",
#                   PlotName = NULL,
#                   data.Origin=data.Origin,
#                   auto.assign=T,
#                   path=path,
#                   nbTicks=nbTicks)
#
#
# test <- HMM.Train(index <- index,
#                   Data,
#                   start.date <- start.date,
#                   end.date   <- end.date,
#                   frequency=frequency,
#                   mult=mult,
#                   # mu0=mu.Test,
#                   # sigma0=sigma.Test,
#                   Gamma0=Gamma.HMM.4,
#                   type="HMM",
#                   distribution='Normal',
#                   Transition.Type="Homogeneous",
#                   PlotName = NULL,
#                   data.Origin=data.Origin,
#                   auto.assign=T,
#                   path=path,
#                   nbTicks=nbTicks)
#
#
test <- HMM.Train(index <- index,
                  Data,
                  start.date <- start.date,
                  end.date   <- end.date,
                  frequency=frequency,
                  mult=mult,
                  # mu0=mu.Test,
                  # sigma0=sigma.Test,
                  Gamma0=Gamma.HMM.4,
                  type="HMM",
                  distribution='Student',
                  Transition.Type="Homogeneous",
                  PlotName = NULL,
                  data.Origin=data.Origin,
                  auto.assign=T,
                  path=path,
                  nbTicks=nbTicks)


test <- HMM.Train(index <- index,
                  Data,
                  start.date <- start.date,
                  end.date   <- end.date,
                  frequency=frequency,
                  mult=mult,
                  # mu0=mu.Test,
                  # sigma0=sigma.Test,
                  Gamma0=Gamma.HHMM,
                  type="HHMM",
                  distribution='Normal',
                  Transition.Type="Homogeneous",
                  PlotName = NULL,
                  data.Origin=data.Origin,
                  auto.assign=T,
                  path=path,
                  nbTicks=nbTicks)


# test <- HMM.Train(index <- index,
#                   Data,
#                   start.date <- start.date,
#                   end.date   <- end.date,
#                   frequency=frequency,
#                   mult=mult,
#                   # mu0=mu.Test,
#                   # sigma0=sigma.Test,
#                   Gamma0=Gamma.HHMM,
#                   type="HHMM",
#                   distribution='Student',
#                   Transition.Type="Homogeneous",
#                   PlotName = NULL,
#                   data.Origin=data.Origin,
#                   auto.assign=T,
#                   path=path,
#                   nbTicks=nbTicks)

