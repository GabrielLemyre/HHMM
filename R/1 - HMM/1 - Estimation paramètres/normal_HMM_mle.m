function [muSort,sigmaSort,gammaSort,deltaSort,tempsMoyenSejour,prob,mllk,AIC,BIC] = normal_HMM_mle(data,distribution,mu,sigma,gamma,maxiter,type)
% FONCTION OPTIMISATION Fonction afin d'optimiser num�riquement la log-vraisemblance
%   Optimisation � l'aide du filtre d'Hamilton (Probabilit�s filtr�es 
%   P[C_t | X_1:t])
% �������������������������������������������
% Creation              : 30 novembre 2018
% �������������������������������������������
% Derni�re modification : 25 f�vrier 2019
% �������������������������������������������

% Choix du nombre de degr� de libert� pour la distribution de Student
% Aussi d�fini dans : normal_HMM_HamiltonFilter -- normal_HMM_Viterbi -- plotHMM_sequence
nu=4;

% Modification du nombre maximal d'it�rations lors de l'optimisation
% puisque certains mod�les prennent plus de temps � converger
options =optimoptions('fminunc','MaxFunctionEvaluations', maxiter);

% Transformation des param�tres naturels en param�tres de travail
parvect0=normal_HMM_N2W(mu,sigma,gamma,type);

% D�finition de la fonction objective pour la minimisation
f=@(parvect)normal_HMM_mllk(distribution, parvect, data,type);

% Minimisation de la fonction de -logvraisemblance
[x, mllk]=fminunc(f,parvect0,options);
mllk=-mllk;

% Transformation des param�tres de travail en param�tres naturels
[mu,sigma,gamma,delta]=normal_HMM_W2N(x,type);

% Si c'est une distribution de Student, alors on doit multiplier par la
% variance du mod�le pour obtenir la vraie variance des rendements
% R_t = mu + sigma * eta_t, o� eta_t ~ student(ddl = nu)
if strcmp(distribution,'Student')
    sigma=sigma*sqrt(nu/(nu-2));
end

prob=zeros(1,4);

switch type
    case 'HMM'
        % Si mod�le HMM, r�organisation des param�tres en ordre d�croissant
        % de volatilit�
        [sigmaSort,sigmaSortORDER]=sort(sigma,'descend'); % Get the order of sigma and sort the values
        muSort=mu(sigmaSortORDER); % sort sigma according to sigma order
        gammaSort=gamma(sigmaSortORDER,sigmaSortORDER); % sort gamma according to sigma order
        deltaSort=delta(sigmaSortORDER); % sort delta according to sigma order
        
    case 'HHMM'
        muSort = mu;
        sigmaSort = sigma;
        
        % �������������������������
        % STRUCTURE POUR LE FACTORIAL HMM
        % �������������������������
        % gamma = [ p1*q1        p1(1-q1)      (1-p1)q1      (1-p1)(1-q1);
        %           p1(1-q2)     p1*q2         (1-p1)(1-q2)  (1-p1)q2    ;
        %          (1-p2)q1     (1-p2)(1-q1)    p2*q1         p2(1-q1)   ;
        %          (1-p2)(1-q2) (1-p2)q2        p2(1-q2)      p2*q2      ]
        %         prob(1)=gamma(1,1)+gamma(1,2); % p1
        %         prob(2)=gamma(3,3)+gamma(3,4); % p2
        %         prob(3)=gamma(3,1)+gamma(3,3); % q1
        %         prob(4)=gamma(2,2)+gamma(2,4); % q2
        % �������������������������

        % Si mod�le HHMM, c'est les probabilit�s telles que d�finies 
        % ci-haut qui sont enregistr�es dans gammaSort au lieu des 
        % probabilit�s asymptotiques puisque gamma peut �tre reconstruite �
        % partir de ces probabilit�s gr�ce � la fonction 'gamma_Build_HHMM'
        % qui se trouve dans 'PVP Direction des risques - udem\Travaux\A -
        % Outils'
        prob = gamma;
        gammaSort = gamma_Build_HHMM(gamma);
        deltaSort = delta;
end

% Calcul du temps moyen dans un �tat tm_i = 1/(1-gamma_ii) puisque
% distribution g�om�trique du temps de s�jour dans un �tat
tempsMoyenSejour = 1./(1-diag(gammaSort));

% Nombre de param�tres du mod�le
np=length(parvect0);
AIC=2*(-mllk+np);
BIC=-2*mllk+log(length(data))*np;

end

