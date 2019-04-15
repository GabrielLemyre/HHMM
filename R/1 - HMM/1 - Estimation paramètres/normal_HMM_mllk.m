function [llk] = normal_HMM_mllk(distribution,parvect, data,type)
%FONCTION OBJECTIF afin de calculer la vraisemblance et la log-vraisemblance 
%   Le calcul est effectué à l'aide du filtre d'hamilton (Probabilités
%   filtrées P[C_t | X_1:t])
% ———————————————————————————————————————————
% Creation              : 30 novembre 2018
% ———————————————————————————————————————————
% Dernière modification : 28 janvier 2019
% ———————————————————————————————————————————

n=length(data);

[mu,sigma,gamma,delta]=normal_HMM_W2N(parvect,type);

if strcmp(type,'HHMM')
    gamma = gamma_Build_HHMM(gamma);
end

[llk,~,~] = normal_HMM_HamiltonFilter(distribution,mu,sigma,gamma,delta,data);

llk=-llk;

end