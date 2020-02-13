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
# First version : april 15th, 2019
# Last version : february 13th, 2020
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
layout(matrix(c(1,1,1,1,2,3,2,3),4,2,byrow=TRUE))
par(mar = rep(6, 4)) # Set the margin on all sides to 2

# Quantité pour taille et options affichage graphique
cex.axis   <- 0.85
cex.mtext  <- 2
cex.legend <- 2.5
cex.pch    <- 1

z.norm<-(logR-mean(logR))/sd(logR) ## standardized data 

# Estimation de la densité des données et histogramme plus courbe normale
hist_plot <- hist(logR, freq=FALSE, breaks=100,
                  plot=TRUE, col="light gray", border="white",
                  xlim=c(-6,6), xaxs="i", xlab="", ylab="",main="",axes=FALSE)
# Ajout de la densité non-paramétrique
lines(density(logR), col="black", lwd=1)
# Ajout de la courbe d'une loi normale avec paramètre empirique
curve(dnorm(x, mean=mu, sd=sigma),
      col="black", lwd=1, lty=2, add=TRUE)
# Ajout des axes, orientation du texte et ajustement de l'écart en
#   l'axe et les étiquettes
axis(side=1, cex.axis=cex.axis, padj=-0.6)
axis(side=2, cex.axis=cex.axis, las=1)
# Ajout d'un titre sur le graphique, modification de son emplacement,
#   de sa taille et de sa police
mtext(paste(index,"\nDensité des rendements",sep=""), side=3, line=0.5, font=4, las=0, cex=cex.mtext)
# Ajout d'une légende sur le graphique
legend("topleft", inset=0.05, legend=c("Non-paramétrique", "Normale"),
       lty=c(1,2), lwd=1, pch=NA, pt.cex=cex.pch, col=c("black","black"),
       bty="n", cex=cex.legend, ncol=1, x.intersp=1, y.intersp=1.5)

# Comparaison des quantiles de l'échantillon par rapport à ceux d'une normale
qqnorm(z.norm, main="") ## drawing the QQplot 
mtext("Q-Q Plot", side=3, line=2, font=4, las=0, cex=cex.mtext)
abline(0,1) ## drawing a 45-degree reference line

# box plot
boxplot(logR,
        border="black",
        lwd=1,
        varwidth = FALSE, xaxt="n", yaxt="n",xlab="",ylab="")
title(main = "Box Plot",
      cex.main=cex.mtext*1.5)
title(xlab="log-rendements",
      line=1)
axis(side=2, las=1)
title(ylab="Observations",
      line=2)

dev.off()

# # Test de Normalité
# qqnorm(logR)
# qqline(logR, col=2)





# ————————————————————————————————————————————————————————————————————————————————————
# ////////////////////////////////////////////////////////////////////
# FAITS STYLISÉS DES RENDEMENTS FINANCIERS
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# -———————————————————————————————————————————————————————————————————————————————————
title.size.var <- 12
axis.lab.size.var <- 8

epsilon <- logR-mean(logR)
conf.level <- 0.95
ciline <- qnorm((1 - conf.level)/2)/sqrt(length(epsilon))

# Fonction d'autocorrélation des rendements
AutoCorrelation <- Acf(logR, plot = FALSE,lag.max = 200)
# # Fonction d'autocorrélation des rendements centrés au carré
AutoCorrelationCarreRES <- Acf(epsilon^2, plot = FALSE,lag.max = 200)
# # Fonction d'autocorrélation des rendements centrés en valeur absolue
AutoCorrelationAbsRES <- Acf(abs(epsilon), plot = FALSE,lag.max = 200)

A.AutoCorrelation <- with(AutoCorrelation, data.frame(lag, acf))
B.AutoCorrelation <- with(AutoCorrelationCarreRES, data.frame(lag, acf))
C.AutoCorrelation <- with(AutoCorrelationAbsRES, data.frame(lag, acf))

# Différence entre carré et val absolue
D.AutoCorrelation <- data.frame(B.AutoCorrelation$lag, (C.AutoCorrelation-B.AutoCorrelation)$acf)
names(D.AutoCorrelation) <- c("lag","acf")

A <- ggplot(data = A.AutoCorrelation, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + 
  geom_hline(aes(yintercept = ciline), linetype = 3, color = 'darkblue') + 
  geom_hline(aes(yintercept = -ciline), linetype = 3, color = 'darkblue') +
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  labs(title=paste("ACF des log-rendements de ",index,sep=""),
       x ="Lags", 
       y = "Auto-Corrélation") + 
  theme(
    plot.title = element_text(color="Black", size=title.size.var, face="bold.italic"),
    axis.title.x = element_text(color="black", size=axis.lab.size.var, face="bold"),
    axis.title.y = element_text(color="black", size=axis.lab.size.var, face="bold")
  )

B <- ggplot(data = B.AutoCorrelation, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + 
  geom_hline(aes(yintercept = ciline), linetype = 3, color = 'darkblue') + 
  geom_hline(aes(yintercept = -ciline), linetype = 3, color = 'darkblue') +
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  labs(title=paste("ACF des log-rendements centrés au carré de ",index,sep=""),
       x ="Lags", 
       y = "Auto-Corrélation") + 
  theme(
    plot.title = element_text(color="Black", size=title.size.var, face="bold.italic"),
    axis.title.x = element_text(color="black", size=axis.lab.size.var, face="bold"),
    axis.title.y = element_text(color="black", size=axis.lab.size.var, face="bold")
  )

C <- ggplot(data = C.AutoCorrelation, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + 
  geom_hline(aes(yintercept = ciline), linetype = 3, color = 'darkblue') + 
  geom_hline(aes(yintercept = -ciline), linetype = 3, color = 'darkblue') +
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  labs(title=paste("ACF de la valeur absolue des log-rendements centrés de ",index,sep=""),
       x ="Lags", 
       y = "Auto-Corrélation") + 
  theme(
    plot.title = element_text(color="Black", size=title.size.var, face="bold.italic"),
    axis.title.x = element_text(color="black", size=axis.lab.size.var, face="bold"),
    axis.title.y = element_text(color="black", size=axis.lab.size.var, face="bold")
  )

D <- ggplot(data = D.AutoCorrelation, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + 
  geom_hline(aes(yintercept = ciline), linetype = 3, color = 'darkblue') + 
  geom_hline(aes(yintercept = -ciline), linetype = 3, color = 'darkblue') +
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  labs(title="Différence entre Corrélation valeur absolue et au carré",
       x ="Lags", 
       y = "Différence") + 
  theme(
    plot.title = element_text(color="Black", size=title.size.var, face="bold.italic"),
    axis.title.x = element_text(color="black", size=axis.lab.size.var, face="bold"),
    axis.title.y = element_text(color="black", size=axis.lab.size.var, face="bold")
  )

p <- plot_grid(A,B,C,D, labels = "AUTO", ncol = 1)
ggsave(paste(GraphPath,"/1 - DataStats/Dataset_ACF_",index,".png",sep=""), p)






# --------------------------------------------------------
# 
# --------------------------------------------------------

maxf <- function(x){max(x,0)}
negImpact <- sapply(-epsilon,FUN=maxf)
posImpact <- sapply(epsilon,FUN=maxf)


# Corrélation entre la valeur absolue des erreurs et un choc négatif
AutoCorrelationLEVIER.NEG <- Ccf(abs(epsilon), negImpact, ylab = "Effet négatif", plot=F,lag.max = 200,level=0.95)
# Corrélation entre la valeur absolue des erreurs et un choc positif
AutoCorrelationLEVIER.POS <- Ccf(abs(epsilon), posImpact, ylab = "Effet positif", plot=F,lag.max = 200)

NEG.AutoCorrelation <- with(AutoCorrelationLEVIER.NEG, data.frame(lag, acf))
NEG.AutoCorrelation.df <- subset(NEG.AutoCorrelation, lag >= 0)
POS.AutoCorrelation <- with(AutoCorrelationLEVIER.POS, data.frame(lag, acf))
POS.AutoCorrelation.df <- subset(POS.AutoCorrelation, lag >= 0)

# Corrélation entre la valeur absolue des erreurs et un choc positif
AutoCorrelationLEVIER.DIFF <- data.frame(NEG.AutoCorrelation.df$lag, (NEG.AutoCorrelation.df-POS.AutoCorrelation.df)$acf)
names(AutoCorrelationLEVIER.DIFF) <- c("lag","acf")

NEG <- ggplot(data = NEG.AutoCorrelation.df, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + 
  geom_hline(aes(yintercept = ciline), linetype = 3, color = 'darkblue') + 
  geom_hline(aes(yintercept = -ciline), linetype = 3, color = 'darkblue') +
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  labs(title="Corrélation avec un choc négatif",
       x ="Lag", 
       y = "Auto-Corrélation") + 
  theme(
    plot.title = element_text(color="Black", size=title.size.var, face="bold.italic"),
    axis.title.x = element_text(color="black", size=axis.lab.size.var, face="bold"),
    axis.title.y = element_text(color="black", size=axis.lab.size.var, face="bold")
  )

POS <- ggplot(data = POS.AutoCorrelation.df, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + 
  geom_hline(aes(yintercept = ciline), linetype = 3, color = 'darkblue') + 
  geom_hline(aes(yintercept = -ciline), linetype = 3, color = 'darkblue') +
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  labs(title="Corrélation avec un choc positif",
       x ="Lag", 
       y = "Auto-Corrélation") + 
  theme(
    plot.title = element_text(color="Black", size=title.size.var, face="bold.italic"),
    axis.title.x = element_text(color="black", size=axis.lab.size.var, face="bold"),
    axis.title.y = element_text(color="black", size=axis.lab.size.var, face="bold")
  )

DIFF <- ggplot(data = AutoCorrelationLEVIER.DIFF, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + 
  geom_hline(aes(yintercept = ciline), linetype = 3, color = 'darkblue') + 
  geom_hline(aes(yintercept = -ciline), linetype = 3, color = 'darkblue') +
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  labs(title="Corrélation choc négatif - Corrélation choc positif",
       x ="Lag", 
       y = "Différence") + 
  theme(
    plot.title = element_text(color="Black", size=title.size.var, face="bold.italic"),
    axis.title.x = element_text(color="black", size=axis.lab.size.var, face="bold"),
    axis.title.y = element_text(color="black", size=axis.lab.size.var, face="bold")
  )

p <- plot_grid(NEG,POS,DIFF, labels = "AUTO", ncol = 1)
ggsave(paste(GraphPath,"/1 - DataStats/Dataset_LEVIER_",index,".png",sep=""), p)

# maxf <- function(x){max(x,0)}
# negImpact <- sapply(-epsilon,FUN=maxf)
# posImpact <- sapply(epsilon,FUN=maxf)
# 
# # Plot afin d'illustrer l'effet de levier
# jpeg(paste(GraphPath,"/1 - DataStats/Dataset_LEVIER_",index,".png",sep=""), width = image.width, height = image.heigth)
# layout(matrix(c(1,2,3),3,1,byrow=TRUE))
# par(mar = rep(6, 4)) # Set the margin on all sides to 2
# 
# # Corrélation entre la valeur absolue des erreurs et un impact négatif
# AutoCorrelationLEVIER.NEG <- ccf(abs(epsilon), negImpact, ylab = "Effet négatif", plot=F)
# # On ne conserve que les lag positifs
# AutoCorrelationLEVIER.NEG.lag <- AutoCorrelationLEVIER.NEG$acf[AutoCorrelationLEVIER.NEG$lag>=0]
# names(AutoCorrelationLEVIER.NEG.lag) <- AutoCorrelationLEVIER.NEG$lag[AutoCorrelationLEVIER.NEG$lag>=0]
# # Diagramme à bande des corrélations croisées
# barplot(height=c(AutoCorrelationLEVIER.NEG.lag), main="", xlab="lags",names.arg=names(AutoCorrelationLEVIER.NEG.lag))
# title(main = "Corrélation avec un impact négatif", cex.main=2)
# 
# # Corrélation entre la valeur absolue des erreurs et un impact positif
# AutoCorrelationLEVIER.POS <- ccf(abs(epsilon), posImpact, ylab = "Effet positif", plot=F)
# # On ne conserve que les lag positifs
# AutoCorrelationLEVIER.POS.lag <- AutoCorrelationLEVIER.POS$acf[AutoCorrelationLEVIER.POS$lag>=0]
# names(AutoCorrelationLEVIER.POS.lag) <- AutoCorrelationLEVIER.POS$lag[AutoCorrelationLEVIER.POS$lag>=0]
# # Diagramme à bande des corrélations croisées
# barplot(height=c(AutoCorrelationLEVIER.POS.lag), main="", xlab="lags",names.arg=names(AutoCorrelationLEVIER.POS.lag))
# title(main = "Corrélation avec un impact positif", cex.main=2)
# 
# # Diagramme à bande de la différence entre les corrélations d'un impact négatif et d'un impact positif
# barplot(height=c(AutoCorrelationLEVIER.NEG.lag-AutoCorrelationLEVIER.POS.lag),names.arg=names(AutoCorrelationLEVIER.NEG.lag), col="red")
# title(main = "correlation impact négatif - correlation impact positif", cex.main=2)
# 
# dev.off()

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


