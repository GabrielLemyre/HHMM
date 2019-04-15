function [] = Model_Test_Allocation(distribution,data,type,nbRegime)
%FONCTION TEST Calcul des rendements d'une stratégie dynamique
%   La fonction entraîne le modèle tous les 'n' pas et calcul ensuite 'n'
%   valeurs de l'état prédit un pas plus loin. Ainsi le modèle est entrainé
%   puis on fait 'n' pas de prédiction puis on réentraine le modèle et on
%   refait 'n' pas et ainsi de suite.
% ———————————————————————————————————————————
% Creation              : 7 décembre 2018
% ———————————————————————————————————————————
% Dernière modification : 28 janvier 2019
% ———————————————————————————————————————————
n=length(data);
n_before_Train=52;
n_step=floor(n/n_before_Train);

etatsPredits=zeros(nbRegime,n);
% Choix du nombre de degré de liberté pour la distribution de Student
% Aussi défini dans : normal_HMM_HamiltonFilter -- normal_HMM_mle -- plotHMM_sequence
nu=4;

% On initialise d'abord les paramètres du modèle
switch type
    case 'HMM'
        switch nbRegime
            case 2
                %mu_i=[-mWD mWD];
                %sigma_i=[stdWD*2 stdWD/2];
                [mu_i, sigma_i]=ParametresInitiaux(WorkData,2);
                gamma_i=[ 0.97 0.03; 0.03 0.97];
            case 3
        %         mu_i=[-mWD mWD 2*mWD];
        %         sigma_i=[stdWD/2 stdWD stdWD*2];
                [mu_i, sigma_i]=ParametresInitiaux(WorkData,3);
                gamma_i=[0.99 0.05 0.05; 0.005 0.99 0.005; 0.05 0.05 0.99];
            case 4
        %         mu_i=[-4*mWD -mWD mWD 2*mWD];
        %         sigma_i=[stdWD*4 stdWD stdWD/2 stdWD/4];
                [mu_i, sigma_i]=ParametresInitiaux(WorkData,4);
                gamma_i=[ 0.97 0.01 0.01 0.01; 0.01 0.97 0.01 0.01; 0.01 0.01 0.97 0.01; 0.01 0.01 0.01 0.97];
        end
    case 'HHMM'
        [mu_i, sigma_i]=ParametresInitiaux(WorkData,4);
        gamma_i=gamma_Build_HHMM([0.9 0.9 0.9 0.9]);
end

for i=1:n_step
    if i>1
        mu_i=mu_f;
        sigma_i=sigma_f;
        gamma_i=gamma_f;
    end
    [mu_f,sigma_f,gamma_f,delta_f,~,~,~,~]=normal_HMM_mle(data(1:(i*n_before_Train)),distribution,mu_i,sigma_i,gamma_i,2500,type{index})
   
    % Calcul des probabilités filtrées par l'algorithme du filtre d'Hamilton
    [~,p_ct_x1t,~] = normal_HMM_HamiltonFilter(distribution,mu_f,sigma_f,gamma_f,delta_f,data(1:(i-1)*n_before_Train));
    etatsPredits(:,(i-1)*n_before_Train+1:(i*n_before_Train))=p_ct_x1t
end

% —————————————————————————————————————————
%  Filtre d'Hamilton et algorithme EM m=3
%  ————————————————————————————————————————————
[mu_f,sigma_f,gamma_f,delta_f,prob_f,~,~,~]=normal_HMM_mle(WorkData,distribution{index},mu_i,sigma_i,gamma_i,2500,type{index})

% Transformation en paramètres de travail des valeurs obtenues
workparam_f=normal_HMM_N2W(mu_f,sigma_f,gamma_f,type{index});

%  —————————————————————————————————————————
%  Filtre de Kim
%  —————————————————————————————————————————
[~,smoothProb_f] = normal_HMM_KimFilter(distribution{index},workparam_f,WorkData,type{index});

%  —————————————————————————————————————————
%  Inférence séquence plus probable avec algorithme de Viterbi
%  —————————————————————————————————————————
[sequence,~] = normal_HMM_Viterbi(distribution{index},workparam_f,WorkData,type{index});


% Calcul des probabilités filtrées par l'algorithme du filtre d'Hamilton
[~,p_ct_x1t,~] = normal_HMM_HamiltonFilter(distribution,mu,sigma,gamma,delta,data);

end

