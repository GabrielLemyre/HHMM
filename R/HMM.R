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
# Last version : april 15th, 2019
# ----------------------------------------------------------------------------------------------------
#set working directory
path <- '~/Documents/GitHub/HHMM'
setwd(path.expand(path)) # Setting path



# ———————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# TRANSFORMATION DES PARAMÈTRES
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# ———————————————————————————————————————————————————————————————————————————————————————————————————————————————————————

# ————————————————————————————————————————————————————————
# PARAMÈTRES CONTRAINTS -> PARAMÈTRES NON-CONTRAINTS
# ————————————————————————————————————————————————————————
normal_HMM_N2W = function(mu, sigma, gamma, type){
# Function to convert from natural parameters to their working equivalents
#   The idea is that by changing the constrained parameters to
#   unconstrained versions, the optimization can be done without
#   constraints
# ———————————————————————————————————————————
# Creation              : 29 novembre 2018
# ———————————————————————————————————————————
# Dernière modification : 25 janvier 2019
# ———————————————————————————————————————————
    
    tsigma ->
    
}

# ————————————————————————————————————————————————————————
# MATRICE GAMMA
# ————————————————————————————————————————————————————————
gamma_Build_HHMM = function(prob_i, type="HMM"){
#   Construction de la matrice de transition (gamma)
#   Cette fonction construit la matrice de transition d'un modèle HHMM à 2
#   états parents avec 2 états enfants chaques par le biais des
#   probabilités de transitions horizontales (p1-p2) et transitions
#   verticales (q1-q2)
#    
# ———————————————————————————————
# STRUCTURE POUR LE HIERARCHICAL HMM
# ———————————————————————————————
#     gamma_i = [a^2_11  +  a^2_12  +  e1     ;  = 1
#                a^2_22  +  a^2_21  +  e2     ;  = 1
#                a^2_33  +  a^2_34  +  e3     ;  = 1
#                a^2_44  +  a^2_43  +  e4     ;  = 1
#                a^1_11  +  a^1_12  +  0      ;  = 1
#                a^1_22  +  a^1_21  +  0      ;  = 1
#                pi_1    +  pi_2    +  0      ;  = 1
#                pi_3    +  pi_4    +  0    ] ;  = 1
#
# ———————————————————————————————
# STRUCTURE POUR LE FACTORIAL HMM
# ———————————————————————————————
# #          p1          (1-p1)          (1-p2)        p2
# gamma_1 = [prob_i(1,1) 1-prob_i(1,1) ; 1-prob_i(1,2) prob_i(1,2)];
# #          q1          (1-q1)          (1-q2)        q2
# gamma_2 = [prob_i(1,3) 1-prob_i(1,3) ; 1-prob_i(1,4) prob_i(1,4)];
# # Produit de Kronecker entre les 2 matrices
# gamma=kron(gamma_1,gamma_2);
# —————————————————————————


t = num2cell(prob_i(1,:));
[a2_11, a2_12, e1] = deal(t{:});

t = num2cell(prob_i(2,:));
[a2_22, a2_21, e2] = deal(t{:});

t = num2cell(prob_i(3,:));
[a2_33, a2_34, e3] = deal(t{:});

t = num2cell(prob_i(4,:));
[a2_44, a2_43, e4] = deal(t{:});

t = num2cell(prob_i(5,:));
[a1_11, a1_12, ~] = deal(t{:});

t = num2cell(prob_i(6,:));
[a1_22, a1_21, ~] = deal(t{:});

t = num2cell(prob_i(7,:));
[pi1,   pi2, ~] = deal(t{:});

t = num2cell(prob_i(8,:));
[pi3,   pi4, ~] = deal(t{:});

gamma = [a2_11+e1*a1_11*pi1  a2_12+e1*a1_11*pi2          e1*a1_12*pi3        e1*a1_12*pi4 ;
a2_21+e2*a1_11*pi1  a2_22+e2*a1_11*pi2          e2*a1_12*pi3        e2*a1_12*pi4 ;
e3*a1_21*pi1        e3*a1_21*pi2    a2_33+e3*a1_22*pi3  a2_34+e3*a1_22*pi4 ;
e4*a1_21*pi1        e4*a1_21*pi2    a2_43+e4*a1_22*pi3  a2_44+e4*a1_22*pi4];
}