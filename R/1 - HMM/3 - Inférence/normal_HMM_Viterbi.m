function [sequence,epsilon] = normal_HMM_Viterbi(distribution,parvect,data,type)
%VITERBI ALGORITHM Implementation de l'algorithme de Viterbi
%   Permet d'obtenir la suite d'騁ats la plus probable, 騁ant donn馥
%   l'馗hantillon complet des donn馥s.
%   ARGMAX(c(1:T)|x(1:T))
%   c(1:T)
% 覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧�
% Creation              : 7 d馗embre 2018
% 覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧�
% Derni鑽e modification : 28 janvier 2019
% 覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧�

% Choix du nombre de degr� de libert� pour la distribution de Student
% Aussi d馭ini dans : normal_HMM_HamiltonFilter -- normal_HMM_mle -- plotHMM_sequence
nu=4;

% Transformation des param鑼res non-contraints en param鑼res contraints
[mu,sigma,Gamma,delta]=normal_HMM_W2N(parvect,type);

if strcmp(type,'HHMM')
    Gamma = gamma_Build_HHMM(Gamma);
end

% Obtention longeur de la cha�ne et nombre de r馮imes du mod鑞e
n=length(data);
nbRegime=length(Gamma);

% Initialisation du vecteur interm馘iaire d'espilon
epsilon=zeros(n,nbRegime);

% Calcul valeur de E_1 non normalis馥 puis normalisation (ﾉviter probl鑪e
% de underflow et overflow)
switch distribution
    case 'Normale'
        temp_epsilon=delta.*(1./sqrt(2*pi*sigma.^2)).*exp((-(data(1)-mu).^2)./(2*sigma.^2));
    case 'Student'
        % "Frequency distribution of standard deviations of samples drawn
        % from a normal population." -  William Sealy Gosset's 1908 paper in Biometrika
        x = (data(1)-mu)./sigma;
        temp_epsilon=delta.*((gamma((nu+1)/2)/gamma(nu/2))/sqrt(nu*pi)./((1+x.^2/nu).^((nu+1)/2)))./sigma;
end
epsilon(1,:)=temp_epsilon/sum(temp_epsilon);

for i=2:n
    switch distribution
    case 'Normale'
        temp_epsilon=max(epsilon(i-1,:)*Gamma,[],1).*(1./sqrt(2*pi*sigma.^2)).*exp((-(data(i)-mu).^2)./(2*sigma.^2));
    case 'Student'
        % "Frequency distribution of standard deviations of samples drawn
        % from a normal population." -  William Sealy Gosset's 1908 paper in Biometrika
        x = (data(i)-mu)./sigma;
        temp_epsilon=max(epsilon(i-1,:)*Gamma,[],1).*((gamma((nu+1)/2)/gamma(nu/2))/sqrt(nu*pi)./((1+x.^2/nu).^((nu+1)/2)))./sigma;
    end
    
    epsilon(i,:)=temp_epsilon/sum(temp_epsilon);
end

% Initialisation de la s駲uence d'騁ats
sequence=zeros(1,n);

% Construction de la s駲uence d'騁ats la plus probable
[~,whichMax1]=max(epsilon(n,:));
sequence(n)=whichMax1;

for i=(n-1):-1:1
    [~,whichMax]=max(Gamma(:,sequence(i+1)).*epsilon(i,:).');
    sequence(i)=whichMax;
end

end

