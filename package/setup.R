# ----------------------------------------------------------------------------------------------------
# Hidden Markov Models
# Setup environment for testing and running code
# ----------------------------------------------------------------------------------------------------
# written
# Gabriel LEMYRE
# ----------------------------------------------------------------------------------------------------
# Under the supervision of :
# Maciej AUGUSTYNIAK
# ----------------------------------------------------------------------------------------------------
# First version : 17 février 2020
# Last version :  17 février 2020
# ----------------------------------------------------------------------------------------------------


# -------------------------------------------------------
# Sourcing private libraries
# set Sourcing directory
# -------------------------------------------------------
path <- '~/Documents/GitHub/R-Tools'
setwd(path.expand(path)) # Setting Sourcing path
source("LaTeXTable.R")
source("RunTime.R")
source("BasicFunctions.R")
# -------------------------------------------------------
# 
# 
# -------------------------------------------------------
# set working directory
# -------------------------------------------------------
path <- '~/Documents/GitHub/Memoire_Maitrise/1 - Babebones Code/HMM'
setwd(path.expand(path)) # Setting path
# -------------------------------------------------------
# 
# 
# -------------------------------------------------------
# Sourcing all files in path folder
# -------------------------------------------------------
source.all( path, grepstring="\\.R",  print.source=TRUE)
# -------------------------------------------------------
# 
# 
# -------------------------------------------------------
# Chemin vers les figures
# -------------------------------------------------------
GraphPath <- paste(path,"R/1 - Graphs",sep="/")
# -------------------------------------------------------
# 
# 
# --------------------------------------------------------
# Préparation de l'échantillon
# --------------------------------------------------------
# Fonction manipulant les données
#     L'idée est de prendre l'échantillon et de préparer
#     les données afin d'être utilisées directement dans
#     les modèles
# --------------------------------------------------------

# Provenance des données (Langage de programmation [R ou Matlab])
data.Origin="R"

#load data
Data <- read.csv("1 - Data/Table_Daily.csv") # MATLAB
# Data <- read.csv("1 - Data/DATA_INDEX.csv") # R

head(Data)

# Paramètre de transformation des données
data.freq <- "daily"
mult <- 100          #multiply all returns by -mult-

# Nom de la variable à conserver
index <- "SPXIndex"

# Transformation de la variable date
dates <- as.Date(as.character(Matlab2Rdate(Data[,1])))

# Borne de dates pour l'échantillon final
start.date <- as.Date("1960-01-04")
end.date   <- as.Date("2019-02-15")
weekdays(start.date); weekdays(end.date);

st <- max(which(dates==start.date),1)
en <- min(which(dates==end.date),length(dates))

dates <- dates[st:en]

# Construction de la variable de prix et ajout des dates comme noms
Pt <- as.numeric(as.character(Data[st:en,index]))
names(Pt) <- dates
Pt <- Pt[!is.na(Pt)]
dates <- names(Pt)

#  Correction valeur éronnée dans la table du S&P500
if (index=="SPXIndex"){
  Pt["2015-05-14"] <- 2121.10
}

# Calcul des log-rendements
logR <- log( Pt[-1] / Pt[-length(Pt)] )
# Multipliction des log-rendements par un facteur 'mult'
logR <- logR * mult


