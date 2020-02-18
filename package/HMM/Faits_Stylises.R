# ----------------------------------------------------------------------------------------------------
# Hidden Markov Models
# Faits stylisés des rendements financiers
# ----------------------------------------------------------------------------------------------------
# written
# Gabriel LEMYRE
# ----------------------------------------------------------------------------------------------------
# Under the supervision of :
# Maciej AUGUSTYNIAK
# ----------------------------------------------------------------------------------------------------
# First version : january 8th, 2020
# Last version : february 13th, 2020
# ----------------------------------------------------------------------------------------------------

# --------------------------------------------------------
# ANALYSE D'UNE SÉRIE ET RÉSULTATS SUR LES FAITS STYLISÉS
# --------------------------------------------------------
# Fonction produisant des tables et graphiques
#     L'idée est de montrer qu'une série présente les
#     caractéristiques des faits stylisés des rendements
#     financiers
# --------------------------------------------------------

# Instalation du package forecast permettant d'utiliser les fonctions ACF, PAcf et CCf
# install.packages('forecast', dependencies = TRUE) # Installation du package
library(forecast)

# Permet de créer des graphiques très complexes
# install.packages('ggplot2', dependencies = TRUE) # Installation du package
library(ggplot2)

# Instalation du package cowplot opermettant de combiner des graphiques
# install.packages('cowplot', dependencies = TRUE) # Installation du package
library(cowplot)

# --------------------------------------------------------
# Paramètres graphiques
# --------------------------------------------------------
title.size.var <- 12
axis.lab.size.var <- 8
lty.ciline <- 2
lwd.ciline <- 0.3
lwd.predict <- 0.2

# Nombre de points à afficher sur l'axe horizontal
nbTicks=10

# JPEG DIMENSIONS FOR OUTPUTED FILE
image.width <- 1250
image.heigth <- 666

# Faits stylisés
title.size.var <- 12
axis.lab.size.var <- 8
resdev=300
# --------------------

# --------------------------------------------------------
# Point (1) et (2) : Absence d'acf pour rendement mais 
# présence pour erreurs^2 et abs(erreurs)
# --------------------------------------------------------
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
x <- A.AutoCorrelation$lag
y <- A.AutoCorrelation$acf
lo.A <- loess(y~x)

B.AutoCorrelation <- with(AutoCorrelationCarreRES, data.frame(lag, acf))
x <- B.AutoCorrelation$lag
y <- B.AutoCorrelation$acf
lo.B <- loess(y~x)

C.AutoCorrelation <- with(AutoCorrelationAbsRES, data.frame(lag, acf))
x <- C.AutoCorrelation$lag
y <- C.AutoCorrelation$acf
lo.C <- loess(y~x)

# Différence entre carré et val absolue
D.AutoCorrelation <- data.frame(B.AutoCorrelation$lag, (C.AutoCorrelation-B.AutoCorrelation)$acf)
names(D.AutoCorrelation) <- c("lag","acf")

A <- ggplot(data = A.AutoCorrelation, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + 
  geom_hline(aes(yintercept = ciline), linetype = lty.ciline, color = 'darkblue', size=lwd.ciline) + 
  geom_hline(aes(yintercept = -ciline), linetype = lty.ciline, color = 'darkblue', size=lwd.ciline) +
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  labs(title=paste("ACF des log-rendements de ",index,sep=""),
       x ="Lags", 
       y = "Auto-Corrélation") + 
  theme(
    plot.title = element_text(color="Black", size=title.size.var, face="bold.italic"),
    axis.title.x = element_text(color="black", size=axis.lab.size.var, face="bold"),
    axis.title.y = element_text(color="black", size=axis.lab.size.var, face="bold")
  ) +
  geom_line(aes(x = lag, y = predict(lo.A)), color = "red", linetype = 1,lwd=lwd.predict)

B <- ggplot(data = B.AutoCorrelation, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + 
  geom_hline(aes(yintercept = ciline), linetype = lty.ciline, color = 'darkblue', size=lwd.ciline) + 
  geom_hline(aes(yintercept = -ciline), linetype = lty.ciline, color = 'darkblue', size=lwd.ciline) +
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  labs(title=paste("ACF des log-rendements centrés au carré de ",index,sep=""),
       x ="Lags", 
       y = "Auto-Corrélation") + 
  theme(
    plot.title = element_text(color="Black", size=title.size.var, face="bold.italic"),
    axis.title.x = element_text(color="black", size=axis.lab.size.var, face="bold"),
    axis.title.y = element_text(color="black", size=axis.lab.size.var, face="bold")
  ) +
  geom_line(aes(x = lag, y = predict(lo.B)), color = "red", linetype = 1,lwd=lwd.predict)

C <- ggplot(data = C.AutoCorrelation, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + 
  geom_hline(aes(yintercept = ciline), linetype = lty.ciline, color = 'darkblue', size=lwd.ciline) + 
  geom_hline(aes(yintercept = -ciline), linetype = lty.ciline, color = 'darkblue', size=lwd.ciline) +
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  labs(title=paste("ACF de la valeur absolue des log-rendements centrés de ",index,sep=""),
       x ="Lags", 
       y = "Auto-Corrélation") + 
  theme(
    plot.title = element_text(color="Black", size=title.size.var, face="bold.italic"),
    axis.title.x = element_text(color="black", size=axis.lab.size.var, face="bold"),
    axis.title.y = element_text(color="black", size=axis.lab.size.var, face="bold")
  ) +
  geom_line(aes(x = lag, y = predict(lo.C)), color = "red", linetype = 1,lwd=lwd.predict)

D <- ggplot(data = D.AutoCorrelation, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + 
  geom_hline(aes(yintercept = ciline), linetype = lty.ciline, color = 'darkblue', size=lwd.ciline) + 
  geom_hline(aes(yintercept = -ciline), linetype = lty.ciline, color = 'darkblue', size=lwd.ciline) +
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
# Point (3) : Clustering des volatilités
# --------------------------------------------------------
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



# --------------------------------------------------------
# Point (4) : Asymétrie négative et grand aplatissement
# --------------------------------------------------------
#descriptive statistics
des_stats <- c("Moyenne" = mean(logR),
               "Écart-type" = sd(logR),
               "Asymétrie" = timeDate::skewness(logR, method="moment"),
               "Aplatissement" = timeDate::kurtosis(logR, method="moment"), 
               "Minimum" = min(logR),
               "Maximum" = max(logR), "n" = length(logR)
)
des_stats



# --------------------------------------------------------
# Point (5) : Effet de levier 
# --------------------------------------------------------
maxf <- function(x){max(x,0)}
negChoc <- sapply(-epsilon,FUN=maxf)
posChoc <- sapply(epsilon,FUN=maxf)


# Corrélation entre la valeur absolue des erreurs et un choc négatif
AutoCorrelationLEVIER.NEG <- Ccf(abs(epsilon), negChoc, plot=F,lag.max = 200,level=0.95)
# Corrélation entre la valeur absolue des erreurs et un choc positif
AutoCorrelationLEVIER.POS <- Ccf(abs(epsilon), posChoc, plot=F,lag.max = 200)

NEG.AutoCorrelation <- with(AutoCorrelationLEVIER.NEG, data.frame(lag, acf))
NEG.AutoCorrelation.df <- subset(NEG.AutoCorrelation, lag >= 0)
x <- NEG.AutoCorrelation.df$lag
y <- NEG.AutoCorrelation.df$acf
lo.LEVIER.NEG <- loess(y~x)

POS.AutoCorrelation <- with(AutoCorrelationLEVIER.POS, data.frame(lag, acf))
POS.AutoCorrelation.df <- subset(POS.AutoCorrelation, lag >= 0)
x <- POS.AutoCorrelation.df$lag
y <- POS.AutoCorrelation.df$acf
lo.LEVIER.POS <- loess(y~x)

# Corrélation entre la valeur absolue des erreurs et un choc positif
AutoCorrelationLEVIER.DIFF <- data.frame(NEG.AutoCorrelation.df$lag, (NEG.AutoCorrelation.df-POS.AutoCorrelation.df)$acf)
names(AutoCorrelationLEVIER.DIFF) <- c("lag","acf")

NEG <- ggplot(data = NEG.AutoCorrelation.df, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + 
  geom_hline(aes(yintercept = ciline), linetype = lty.ciline, color = 'darkblue') + 
  geom_hline(aes(yintercept = -ciline), linetype = lty.ciline, color = 'darkblue') +
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  labs(title="Corrélation avec un choc négatif",
       x ="Lag", 
       y = "Auto-Corrélation") + 
  theme(
    plot.title = element_text(color="Black", size=title.size.var, face="bold.italic"),
    axis.title.x = element_text(color="black", size=axis.lab.size.var, face="bold"),
    axis.title.y = element_text(color="black", size=axis.lab.size.var, face="bold")
  ) +
  geom_line(aes(x = lag, y = predict(lo.LEVIER.NEG)), color = "red", linetype = 1,lwd=lwd.predict)

POS <- ggplot(data = POS.AutoCorrelation.df, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + 
  geom_hline(aes(yintercept = ciline), linetype = lty.ciline, color = 'darkblue') + 
  geom_hline(aes(yintercept = -ciline), linetype = lty.ciline, color = 'darkblue') +
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  labs(title="Corrélation avec un choc positif",
       x ="Lag", 
       y = "Auto-Corrélation") + 
  theme(
    plot.title = element_text(color="Black", size=title.size.var, face="bold.italic"),
    axis.title.x = element_text(color="black", size=axis.lab.size.var, face="bold"),
    axis.title.y = element_text(color="black", size=axis.lab.size.var, face="bold")
  ) +
  geom_line(aes(x = lag, y = predict(lo.LEVIER.POS)), color = "red", linetype = 1,lwd=lwd.predict)

DIFF <- ggplot(data = AutoCorrelationLEVIER.DIFF, mapping = aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) + 
  geom_hline(aes(yintercept = ciline), linetype = lty.ciline, color = 'darkblue') + 
  geom_hline(aes(yintercept = -ciline), linetype = lty.ciline, color = 'darkblue') +
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
