# ----------------------------------------------------------------------------------------------------
# Hidden Markov Models
# Unconstraining parameters
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
# TEST DE NORMALITÉ
# --------------------------------------------------------
# Fonctions permettant de tester la normalité de l'échantillon
#     En premier lieu on trace la densité des observations,
#     puis on trace le graphique des quantiles de la loi normale
#     en comparaison avec ceux de notre échantillon.
# --------------------------------------------------------

# Quantité pour taille et options affichage graphique
cex.axis   <- 0.85
cex.mtext  <- 0.85
cex.legend <- 0.6
cex.pch    <- 0.5

# Moyenne et Écart-type échantillonale pour densité normale univariée
mu <- mean(logR)
sigma <- sd(logR)
# Données standardisées
z.norm<-(logR-mu)/sigma


# --------------------------------------------------------
# GRAPHIQUE : Densité paramétrique, non-paramétrique et histogramme
# --------------------------------------------------------
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
mtext("Densité des rendements", side=3, line=0.5, font=4, las=0, cex=cex.mtext)
# Ajout d'une légende sur le graphique
legend("topleft", inset=0.05, legend=c("Non-paramétrique", "Normale"),
       lty=c(1,2), lwd=1, pch=NA, pt.cex=cex.pch, col=c("black","black"),
       bty="n", cex=cex.legend, ncol=1, x.intersp=1, y.intersp=1.5)


# --------------------------------------------------------
# GRAPHIQUE : Comparaison des quantilles de la loi normale et de l'échantillon
# --------------------------------------------------------
# Comparaison des quantiles de l'échantillon par rapport à ceux d'une normale
qqnorm(z.norm, main="") ## drawing the QQplot 
mtext("Q-Q Plot", side=3, line=2, font=4, las=0, cex=cex.mtext)
abline(0,1) ## drawing a 45-degree reference line


# --------------------------------------------------------
# GRAPHIQUE : Box Plot de l'échantillon
# --------------------------------------------------------
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