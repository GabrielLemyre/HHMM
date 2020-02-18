# ----------------------------------------------------------------------------------------------------
# Hidden Markov Models
# Initialisation
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

# install.packages('numDeriv')
# install.packages('Rsolnp')
# install.packages('rootSolve')
# install.packages('data.table')
# install.packages('matrixcalc') # To test weither matrix is positive definite
# install.packages('magic') #combiner matrice par blocs
# install.packages('gridExtra') # add tables in plots
# install.packages('grid') # add tables in plots
# install.packages('plotrix') # add tables in plots
# install.packages('timeDate') # add tables in plots
# install.packages('rootSolve') # Outil d'analyse de séries chronologiques
# install.packages('MASS') # permet d'entrainer une normale sans spécifier la llk explicitement
# install.packages('forecast', dependencies = TRUE) # Installation du package
# install.packages('ggplot2', dependencies = TRUE) # Installation du package
# install.packages('cowplot', dependencies = TRUE) # Installation du package
# install.packages("miceadds")

#load numDeriv package for computing numerical derivatives
#enables the use of functions grad and hessian
library(numDeriv)
#load Rsolnp package for numerical optimization with solnp function
library(Rsolnp)
library(rootSolve)
library(data.table) # Permet d'utiliser rbindlist
library(matrixcalc) # Permet d'utiliser rbindlist
library(magic)
library(gridExtra)
library(grid)
library(plotrix)
library(timeDate)
library("rootSolve", lib.loc="~/Library/R/3.1/library")
library(MASS) ## loading package MASS

# Instalation du package forecast permettant d'utiliser les fonctions ACF, PAcf et CCf
library(forecast)

# Permet de créer des graphiques très complexes
library(ggplot2)

# Instalation du package cowplot opermettant de combiner des graphiques
library(cowplot)

# Allows the use of source.all to source all files in a folder
library(miceadds)
