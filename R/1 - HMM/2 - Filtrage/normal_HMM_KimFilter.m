function [p_ct_x1t,p_ct_x1T] = normal_HMM_KimFilter(distribution,parvect,data,type)
%KIM FILTER - ALGORITHME DE LISSAGE Implementation du filtre de Kim
%   Permet le calcul des probabilit�s lis�es : phi_t = P[C_t | X_1:T]
% �������������������������������������������
% Creation              : 7 d�cembre 2018
% �������������������������������������������
% Derni�re modification : 28 janvier 2019
% �������������������������������������������

% Transformation des param�tres non-contraints en param�tres contraints
[mu,sigma,gamma,delta]=normal_HMM_W2N(parvect,type);

if strcmp(type,'HHMM')
    gamma = gamma_Build_HHMM(gamma);
end

n=length(data);
nbRegime=length(gamma);

% Calcul des probabilit�s filtr�es par l'algorithme du filtre d'Hamilton
[~,p_ct_x1t,~] = normal_HMM_HamiltonFilter(distribution,mu,sigma,gamma,delta,data);

% Initialisation de la matrice de probabilit�s liss�es
p_ct_x1T=zeros(nbRegime,n);

% Probabilite de la derniere observation de la variable latente
p_ct_x1T(:,n)=p_ct_x1t(:,n);

% Algorithme de lissage (Filtre de Kim)
for t=(n-1):-1:1
    phi_t_ij=zeros(nbRegime,nbRegime);
    
    for i=1:nbRegime
        for j=1:nbRegime
            % Calcul des probabilit�s conjointes
            phi_t_ij(i,j)=p_ct_x1T(j,t+1)*(p_ct_x1t(i,t)*gamma(i,j)/sum(p_ct_x1t(:,t).'*gamma(:,j)));
        end
    end
    p_ct_x1T(:,t)=sum(phi_t_ij,2);
end

end

