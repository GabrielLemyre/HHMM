function [errMu,errSigma,errGamma,tabMU,tabSIGMA,tabGAMMA,tabGAMMAperc] = EstimationErreurType(funMLLK,funWN,wp)
%CALCUL ERREUR TYPE Cette fonction calcul les erreurs types des paramètres
%   Cette fonction reçoit le 'function handle' de la fonction de calcul de la
%   LOG-VRAISEMBLANCE (funMLLK) ainsi que le 'function handle' de la
%   fonction de transformation des paramètres NC en paramètres C (funWN)
%
%   Le calcul passe par l'estimation numérique du gradient selon les
%   paramètres de travail Non-contraints (NC) puis par la multiplication
%   avec le Jacobien numérique de la transformation de ces paramètres en
%   paramètres naturels Contraints (C)


% Transformation des paramètres NC en paramètres C
[muFinal,sigmaFinal,gammaFinal,~]=normal_HMM_W2N(wp,'HMM');

nb_wp=length(wp);
nbRegime=(-1+sqrt(1+4*nb_wp))/2; % Récupération du nombre de régimes du modèle

% Initialisation du vecteur d'erreurs
errGamma=zeros(nbRegime,nbRegime);

% Approximation numérique de la matrice Hessienne du modèle selon les
% paramètres NC
hess=hessian(funMLLK,wp);

% Calcul du Jacobien de la transformation des paramètres
gj=jacobianest(funWN,wp);

% Transformation de la matrice hessienne NC en hessienne C
mat_var_cov=gj*inv(hess)*gj.';

% Extraction des erreurs types de la diagonale de la matrice hessienne
err=sqrt(diag(mat_var_cov));

% SÉPARATION DES RÉSULTATS PAR PARAMÈTRE
errMu=err(1:nbRegime).';
errSigma=err(nbRegime+1:2*nbRegime).';

for i=1:nbRegime
    errGamma(i,:)=err(2*nbRegime+(i-1)*nbRegime+1:2*nbRegime+(i)*nbRegime);
end


% PRÉSENTATION DES RÉSULTATS
tabMU=table(round(muFinal,4).',round(errMu,4).',abs(100*round((errMu./muFinal),4)).');
tabMU.Properties.RowNames = cellstr(strcat('mu_',string(1:nbRegime)));
tabMU.Properties.VariableNames = {'mu','err','err_perc'}

tabSIGMA=table(round(sigmaFinal,4).',round(errSigma,4).',abs(100*round((errSigma./sigmaFinal),4)).');
tabSIGMA.Properties.RowNames = cellstr(strcat('sigma_',string(1:nbRegime)));
tabSIGMA.Properties.VariableNames = {'sigma','err','err_perc'}

tabGAMMA=table(strcat(string(round(gammaFinal,3)),'±(',string(round(errGamma,3)),')'));
tabGAMMA.Properties.RowNames = cellstr(strcat('etat_',string(1:nbRegime)));
tabGAMMA.Properties.VariableNames = {'Erreur_Gamma'}

tabGAMMAperc=table(abs(100*round((round(errGamma,5)./round(gammaFinal,5)),4)));
tabGAMMAperc.Properties.RowNames = cellstr(strcat('etat_',string(1:nbRegime)));
tabGAMMAperc.Properties.VariableNames = {'Erreur_Gamma_percentage'}

end

