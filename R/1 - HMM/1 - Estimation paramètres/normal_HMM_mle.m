function [muSort,sigmaSort,gammaSort,deltaSort,tempsMoyenSejour,prob,mllk,AIC,BIC] = normal_HMM_mle(data,distribution,mu,sigma,gamma,maxiter,type)
% FONCTION OPTIMISATION Fonction afin d'optimiser numériquement la log-vraisemblance
%   Optimisation à l'aide du filtre d'Hamilton (Probabilités filtrées 
%   P[C_t | X_1:t])
% ———————————————————————————————————————————
% Creation              : 30 novembre 2018
% ———————————————————————————————————————————
% Dernière modification : 25 février 2019
% ———————————————————————————————————————————

% Choix du nombre de degré de liberté pour la distribution de Student
% Aussi défini dans : normal_HMM_HamiltonFilter -- normal_HMM_Viterbi -- plotHMM_sequence
nu=4;

% Modification du nombre maximal d'itérations lors de l'optimisation
% puisque certains modèles prennent plus de temps à converger
options =optimoptions('fminunc','MaxFunctionEvaluations', maxiter);

% Transformation des paramètres naturels en paramètres de travail
parvect0=normal_HMM_N2W(mu,sigma,gamma,type);

% Définition de la fonction objective pour la minimisation
f=@(parvect)normal_HMM_mllk(distribution, parvect, data,type);

% Minimisation de la fonction de -logvraisemblance
[x, mllk]=fminunc(f,parvect0,options);
mllk=-mllk;

% Transformation des paramètres de travail en paramètres naturels
[mu,sigma,gamma,delta]=normal_HMM_W2N(x,type);

% Si c'est une distribution de Student, alors on doit multiplier par la
% variance du modèle pour obtenir la vraie variance des rendements
% R_t = mu + sigma * eta_t, où eta_t ~ student(ddl = nu)
if strcmp(distribution,'Student')
    sigma=sigma*sqrt(nu/(nu-2));
end

prob=zeros(1,4);

switch type
    case 'HMM'
        % Si modèle HMM, réorganisation des paramètres en ordre décroissant
        % de volatilité
        [sigmaSort,sigmaSortORDER]=sort(sigma,'descend'); % Get the order of sigma and sort the values
        muSort=mu(sigmaSortORDER); % sort sigma according to sigma order
        gammaSort=gamma(sigmaSortORDER,sigmaSortORDER); % sort gamma according to sigma order
        deltaSort=delta(sigmaSortORDER); % sort delta according to sigma order
        
    case 'HHMM'
        muSort = mu;
        sigmaSort = sigma;
        
        % —————————————————————————
        % STRUCTURE POUR LE FACTORIAL HMM
        % —————————————————————————
        % gamma = [ p1*q1        p1(1-q1)      (1-p1)q1      (1-p1)(1-q1);
        %           p1(1-q2)     p1*q2         (1-p1)(1-q2)  (1-p1)q2    ;
        %          (1-p2)q1     (1-p2)(1-q1)    p2*q1         p2(1-q1)   ;
        %          (1-p2)(1-q2) (1-p2)q2        p2(1-q2)      p2*q2      ]
        %         prob(1)=gamma(1,1)+gamma(1,2); % p1
        %         prob(2)=gamma(3,3)+gamma(3,4); % p2
        %         prob(3)=gamma(3,1)+gamma(3,3); % q1
        %         prob(4)=gamma(2,2)+gamma(2,4); % q2
        % —————————————————————————

        % Si modèle HHMM, c'est les probabilités telles que définies 
        % ci-haut qui sont enregistrées dans gammaSort au lieu des 
        % probabilités asymptotiques puisque gamma peut être reconstruite à
        % partir de ces probabilités grâce à la fonction 'gamma_Build_HHMM'
        % qui se trouve dans 'PVP Direction des risques - udem\Travaux\A -
        % Outils'
        prob = gamma;
        gammaSort = gamma_Build_HHMM(gamma);
        deltaSort = delta;
end

% Calcul du temps moyen dans un état tm_i = 1/(1-gamma_ii) puisque
% distribution géométrique du temps de séjour dans un état
tempsMoyenSejour = 1./(1-diag(gammaSort));

% Nombre de paramètres du modèle
np=length(parvect0);
AIC=2*(-mllk+np);
BIC=-2*mllk+log(length(data))*np;

end

