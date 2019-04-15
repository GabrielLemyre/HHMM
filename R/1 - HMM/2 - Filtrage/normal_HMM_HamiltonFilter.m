function [llk,p_ct_x1t,p_ct_x1tm1] = normal_HMM_HamiltonFilter(distribution,mu,sigma,Gamma,delta,data)
%HAMILTON FILTER Implementation du filtre d'Hamilton
%   Permet le calcul des probabilités filtrées : omega_t = P[C_t | X_1:t]
% ———————————————————————————————————————————
% Creation              : 7 décembre 2018
% ———————————————————————————————————————————
% Dernière modification : 31 janvier 2019
% ———————————————————————————————————————————

% Choix du nombre de degré de liberté pour la distribution de Student
% Aussi défini dans : normal_HMM_Viterbi -- normal_HMM_mle -- plotHMM_sequence
nu=4;

n=length(data);
p_ct_x1t=zeros(length(mu),n);
%  —————————————————————————————————————————
%  Probabilités prédictive à un pas de temps
%  —————————————————————————————————————————
p_ct_x1tm1=zeros(length(mu),n);

% —————————————————————————————————————————————
% HAMILTON FORWARD FILTERING
% —————————————————————————————————————————————
% Calcul récursif des termes de la log vraisemblance
p_ct_x1tm1(:,1)=delta;
switch distribution
    case 'Normale'
        a_j=p_ct_x1tm1(:,1).'.*(1./sqrt(2*pi*sigma.^2)).*exp((-(data(1)-mu).^2)./(2*sigma.^2));
%         a_j=delta.*normpdf((data(1)-mu)./sigma,0,1)./sigma;
    case 'Student'
        % "Frequency distribution of standard deviations of samples drawn
        % from a normal population." -  William Sealy Gosset's 1908 paper in Biometrika
        % Utilisation de sigma car lors de la transformation en paramètres
        % de travail, c'est cette valeur qui est tenue positive ce qui
        % correspond aux conditions sur nu dans la distribution de student
%         a_j=delta.*tpdf((data(1)-mu)./sigma,nu)./sigma;
        x = (data(1)-mu)./sigma;
        a_j=p_ct_x1tm1(:,1).'.*((gamma((nu+1)/2)/gamma(nu/2))/sqrt(nu*pi)./((1+x.^2/nu).^((nu+1)/2)))./sigma;
end
a=sum(a_j);
llk=log(a);
omega_t=a_j./a;
p_ct_x1t(:,1)=omega_t;

for i=2:n
    p_ct_x1tm1(:,i)=omega_t*Gamma;
    switch distribution
        case 'Normale'
            a_j=p_ct_x1tm1(:,i).'.*(1./sqrt(2*pi*sigma.^2)).*exp((-(data(i)-mu).^2)./(2*sigma.^2));
    %         a_j=omega_t*gamma.*normpdf((data(i)-mu)./sigma)./sigma;
        case 'Student'
            % "Frequency distribution of standard deviations of samples drawn
            % from a normal population." -  William Sealy Gosset's 1908 paper in Biometrika
            % Utilisation de sigma car lors de la transformation en paramètres
            % de travail, c'est cette valeur qui est tenue positive ce qui
            % correspond aux conditions sur nu dans la distribution de student
    %         a_j=omega_t*Gamma.*tpdf((data(i)-mu)./sigma,nu)./sigma;
            x = (data(i)-mu)./sigma;
            a_j=p_ct_x1tm1(:,i).'.*((gamma((nu+1)/2)/gamma(nu/2))/sqrt(nu*pi)./((1+x.^2/nu).^((nu+1)/2)))./sigma;
    end
    a=sum(a_j);
    llk=llk+log(a);
    omega_t=a_j./a;
    p_ct_x1t(:,i)=omega_t;
end

end

