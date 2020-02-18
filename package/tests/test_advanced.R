Analyse <- NULL
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
#                   distribution='Normal',
#                   Transition.Type="Homogeneous",
#                   PlotName = NULL,
#                   data.Origin=data.Origin,
#                   auto.assign=T,
#                   path=path,
#                   nbTicks=nbTicks)


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

