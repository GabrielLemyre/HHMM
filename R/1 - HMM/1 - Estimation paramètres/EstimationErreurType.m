function [errMu,errSigma,errGamma,tabMU,tabSIGMA,tabGAMMA,tabGAMMAperc] = EstimationErreurType(funMLLK,funWN,wp)
%CALCUL ERREUR TYPE Cette fonction calcul les erreurs types des param�tres
%   Cette fonction re�oit le 'function handle' de la fonction de calcul de la
%   LOG-VRAISEMBLANCE (funMLLK) ainsi que le 'function handle' de la
%   fonction de transformation des param�tres NC en param�tres C (funWN)
%
%   Le calcul passe par l'estimation num�rique du gradient selon les
%   param�tres de travail Non-contraints (NC) puis par la multiplication
%   avec le Jacobien num�rique de la transformation de ces param�tres en
%   param�tres naturels Contraints (C)


% Transformation des param�tres NC en param�tres C
[muFinal,sigmaFinal,gammaFinal,~]=normal_HMM_W2N(wp,'HMM');

nb_wp=length(wp);
nbRegime=(-1+sqrt(1+4*nb_wp))/2; % R�cup�ration du nombre de r�gimes du mod�le

% Initialisation du vecteur d'erreurs
errGamma=zeros(nbRegime,nbRegime);

% Approximation num�rique de la matrice Hessienne du mod�le selon les
% param�tres NC
hess=hessian(funMLLK,wp);

% Calcul du Jacobien de la transformation des param�tres
gj=jacobianest(funWN,wp);

% Transformation de la matrice hessienne NC en hessienne C
mat_var_cov=gj*inv(hess)*gj.';

% Extraction des erreurs types de la diagonale de la matrice hessienne
err=sqrt(diag(mat_var_cov));

% S�PARATION DES R�SULTATS PAR PARAM�TRE
errMu=err(1:nbRegime).';
errSigma=err(nbRegime+1:2*nbRegime).';

for i=1:nbRegime
    errGamma(i,:)=err(2*nbRegime+(i-1)*nbRegime+1:2*nbRegime+(i)*nbRegime);
end


% PR�SENTATION DES R�SULTATS
tabMU=table(round(muFinal,4).',round(errMu,4).',abs(100*round((errMu./muFinal),4)).');
tabMU.Properties.RowNames = cellstr(strcat('mu_',string(1:nbRegime)));
tabMU.Properties.VariableNames = {'mu','err','err_perc'}

tabSIGMA=table(round(sigmaFinal,4).',round(errSigma,4).',abs(100*round((errSigma./sigmaFinal),4)).');
tabSIGMA.Properties.RowNames = cellstr(strcat('sigma_',string(1:nbRegime)));
tabSIGMA.Properties.VariableNames = {'sigma','err','err_perc'}

tabGAMMA=table(strcat(string(round(gammaFinal,3)),'�(',string(round(errGamma,3)),')'));
tabGAMMA.Properties.RowNames = cellstr(strcat('etat_',string(1:nbRegime)));
tabGAMMA.Properties.VariableNames = {'Erreur_Gamma'}

tabGAMMAperc=table(abs(100*round((round(errGamma,5)./round(gammaFinal,5)),4)));
tabGAMMAperc.Properties.RowNames = cellstr(strcat('etat_',string(1:nbRegime)));
tabGAMMAperc.Properties.VariableNames = {'Erreur_Gamma_percentage'}

end

