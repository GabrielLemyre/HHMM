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
GraphPath <- paste(path,"R/1 - Graphs",sep="/")

# -------------------------------------------------------
# 
# 
resdev=300
# -------------------------------------------------------
# Préparation des données
# -------------------------------------------------------
data.Origin="R"
nbTicks=10

#load data
Data <- read.csv("1 - Data/Table_Daily.csv") # MATLAB
# Data <- read.csv("1 - Data/DATA_INDEX.csv") # R
head(Data)
data.freq <- "daily"
mult <- 100          #multiply all returns by -mult-


index <- "SPXIndex"
start.date <- as.Date("1999-12-31")
end.date   <- as.Date("2016-12-28")
#weekdays(start.date); weekdays(end.date);

dates <- as.Date(as.character(Matlab2Rdate(Data[,1])))
st <- which(dates==start.date)
en <- min(which(dates==end.date),length(dates))
dates <- dates[st:en]
Pt <- as.numeric(as.character(Data[st:en,index]))

names(Pt) <- dates
Pt <- Pt[!is.na(Pt)]
dates <- names(Pt)

#  Correction valeur éronnée dans la table du S&P500
if (index=="SPXIndex"){
  Pt["2015-05-14"] <- 2121.10
}


# AR.Order.Original <- ar(Pt)


logR <- log( Pt[-1] / Pt[-length(Pt)] )
logR <- logR * mult
#the length of logR is one less than Pt


#descriptive statistics
des_stats <- c("Moyenne" = mean(logR),
               "Écart-type" = sd(logR),
               "Asymétrie" = timeDate::skewness(logR, method="moment"),
               "Aplatissement" = timeDate::kurtosis(logR, method="moment"), "Minimum" = min(logR),
               "Maximum" = max(logR), "n" = length(logR)
)
des_stats

Make.LaTeX.Table(matrix(des_stats,nrow=1),
                 Col.Titles = names(des_stats),
                 Cross.Lines = T, 
                 Row.Pos = 'c',
                 title=paste("Statistiques descriptives sur ",index,sep=""))

# JPEG DIMENSIONS FOR OUTPUTED FILE
image.width <- 1250
image.heigth <- 666

# Plot the dataset and export resulting plot
jpeg(paste(GraphPath,"/1 - DataStats/Dataset_",index,".png",sep=""), width = image.width, height = image.heigth)
layout(matrix(c(1,2),2,1,byrow=TRUE))
par(mar = rep(6, 4)) # Set the margin on all sides to 2
xtick<-seq(1, length(dates), by=floor(length(dates)/nbTicks))


plot(Pt,type="l",xlab="Dates",ylab="Value", xaxt="n", cex.axis=1.8, cex.lab=2)
axis(side=1, at=xtick, labels=dates[xtick], cex.axis=1.8)
title(main=paste("Valeur de l'indice ",index,sep=""), cex.main=2)

plot(logR,type="l",xlab="Dates", xaxt="n", cex.axis=1.8, cex.lab=2)
axis(side=1, at=xtick, labels=dates[xtick], cex.axis=1.8)
title(main=paste("Log-rendement sur l'indice ",index,sep=""), cex.main=2)
abline(h=0, col="blue")


dev.off()



# Plot the normality tests and export resulting plot
jpeg(paste(GraphPath,"/1 - DataStats/Dataset_Normality_",index,".png",sep=""), width = image.width, height = image.heigth)
layout(matrix(c(1,2),2,1,byrow=TRUE))
par(mar = rep(6, 4)) # Set the margin on all sides to 2

z.norm<-(logR-mean(logR))/sd(logR) ## standardized data 

# Estimation de la densité des données
plot(density(logR),main=paste("Estimation de la densité non-paramétrique des log-rendements de",index), cex.main=2)

# Comparaison des quantiles de l'échantillon par rapport à ceux d'une normale
qqnorm(z.norm, cex.main=2) ## drawing the QQplot 
abline(0,1) ## drawing a 45-degree reference line

dev.off()

# # Test de Normalité
# qqnorm(logR)
# qqline(logR, col=2)





# ————————————————————————————————————————————————————————————————————————————————————
# ////////////////////////////////////////////////////////////////////
# FAITS STYLISÉS DES RENDEMENTS FINANCIERS
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# -———————————————————————————————————————————————————————————————————————————————————

epsilon <- logR-mean(logR)

# Plot de l'autocorrélation des rendements et du carré des erreurs au carrée
jpeg(paste(GraphPath,"/1 - DataStats/Dataset_ACF_",index,".png",sep=""), width = image.width, height = image.heigth)
layout(matrix(c(1,2,3),3,1,byrow=TRUE))
par(mar = rep(6, 4)) # Set the margin on all sides to 2

# Fonction d'autocorrélation des rendements
AutoCorrelation <- acf(logR, plot = FALSE)
plot(AutoCorrelation, main="")
title(main = paste("ACF des log-rendements de ",index,sep=""), cex.main=2)

# Fonction d'autocorrélation des rendements centrés au carré
AutoCorrelationRES <- acf(epsilon^2, plot = FALSE)
plot(AutoCorrelationRES, main="")
title(main = paste("ACF des log-rendements centrés au carré de ",index,sep=""), cex.main=2)

# Fonction d'autocorrélation des rendements centrés au carré
AutoCorrelationAbsRES <- acf(abs(epsilon), plot = FALSE)
plot(AutoCorrelationAbsRES, main="")
title(main = paste("ACF de la valeur absolue des log-rendements centrés de ",index,sep=""), cex.main=2)

dev.off()


maxf <- function(x){max(x,0)}
negImpact <- sapply(-epsilon,FUN=maxf)
posImpact <- sapply(epsilon,FUN=maxf)

# Plot afin d'illustrer l'effet de levier
jpeg(paste(GraphPath,"/1 - DataStats/Dataset_LEVIER_",index,".png",sep=""), width = image.width, height = image.heigth)
layout(matrix(c(1,2,3),3,1,byrow=TRUE))
par(mar = rep(6, 4)) # Set the margin on all sides to 2

# Corrélation entre la valeur absolue des erreurs et un impact négatif
AutoCorrelationLEVIER.NEG <- ccf(abs(epsilon), negImpact, ylab = "Effet négatif", plot=F)
# On ne conserve que les lag positifs
AutoCorrelationLEVIER.NEG.lag <- AutoCorrelationLEVIER.NEG$acf[AutoCorrelationLEVIER.NEG$lag>=0]
names(AutoCorrelationLEVIER.NEG.lag) <- AutoCorrelationLEVIER.NEG$lag[AutoCorrelationLEVIER.NEG$lag>=0]
# Diagramme à bande des corrélations croisées
barplot(height=c(AutoCorrelationLEVIER.NEG.lag), main="", xlab="lags",names.arg=names(AutoCorrelationLEVIER.NEG.lag))
title(main = "Corrélation avec un impact négatif", cex.main=2)

# Corrélation entre la valeur absolue des erreurs et un impact positif
AutoCorrelationLEVIER.POS <- ccf(abs(epsilon), posImpact, ylab = "Effet positif", plot=F)
# On ne conserve que les lag positifs
AutoCorrelationLEVIER.POS.lag <- AutoCorrelationLEVIER.POS$acf[AutoCorrelationLEVIER.POS$lag>=0]
names(AutoCorrelationLEVIER.POS.lag) <- AutoCorrelationLEVIER.POS$lag[AutoCorrelationLEVIER.POS$lag>=0]
# Diagramme à bande des corrélations croisées
barplot(height=c(AutoCorrelationLEVIER.POS.lag), main="", xlab="lags",names.arg=names(AutoCorrelationLEVIER.POS.lag))
title(main = "Corrélation avec un impact positif", cex.main=2)

# Diagramme à bande de la différence entre les corrélations d'un impact négatif et d'un impact positif
barplot(height=c(AutoCorrelationLEVIER.NEG.lag-AutoCorrelationLEVIER.POS.lag),names.arg=names(AutoCorrelationLEVIER.NEG.lag), col="red")
title(main = "correlation impact négatif - correlation impact positif", cex.main=2)

dev.off()

boxplot(logR)

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


