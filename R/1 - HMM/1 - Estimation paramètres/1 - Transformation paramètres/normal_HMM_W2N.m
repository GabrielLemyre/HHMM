function [mu,sigma,gamma,delta] = normal_HMM_W2N(parvect,type)
% Function to convert from WORKING parameters to their NATURAL equivalents
%   Retransforms the WORKING paramaters (unconstrained) to the parameters
%   of the model 
% ———————————————————————————————————————————
% Creation              : 30 novembre 2018
% ———————————————————————————————————————————
% Dernière modification : 25 février  2019
% ———————————————————————————————————————————

switch type 
    case 'HMM' % Hidden Markov Model
        nparvect=length(parvect);
        nbRegime=(-1+sqrt(1+4*nparvect))/2; % Récupération du nombre de régimes du modèle

        mu=parvect(1:nbRegime);
        sigma=exp(parvect((nbRegime+1):(2*nbRegime)));

        gamma=exp(parvect((2*nbRegime+1):length(parvect))).';
        gamma=reshape(gamma,[(nbRegime-1),nbRegime]).';
        gamma=[gamma ones(nbRegime,1)];
        gammaRows = num2cell(gamma, 2);       % Collect the rows into cells
        rowSums = cellfun(@sum, gammaRows);   % Calculate the sum by row
        gamma=gamma./rowSums;
        
        % Just to fit the format of the HHMM model below
        gamma_Built = gamma;
        
    case 'HHMM' % Hierchical Hidden Markov Model
        nbRegime = 4;
        mu=parvect(1:4);
        sigma=exp(parvect(5:8));

        gamma=exp(parvect(9:16)).';
        gamma=reshape(gamma,[2,4]).';
        gamma=[gamma ones(4,1)];
        gammaRows = num2cell(gamma, 2);       % Collect the rows into cells
        rowSums = cellfun(@sum, gammaRows);   % Calculate the sum by row
        gamma1=gamma./rowSums;

        gamma=exp(parvect(17:20)).';
        gamma=reshape(gamma,[1,4]).';
        gamma=[gamma ones(4,1)];
        gammaRows = num2cell(gamma, 2);       % Collect the rows into cells
        rowSums = cellfun(@sum, gammaRows);   % Calculate the sum by row
        gamma2=gamma./rowSums;
        
        % Padding the second part to fit dimensions of the fisrt one and
        % then stacking them
        gamma = [gamma1 ; gamma2 zeros(4,1)];
        
        % Building the transition matrix gamma only for the purpose of delta
        % calculation
        gamma_Built = gamma_Build_HHMM(gamma);
        
        
% —————————————————————————
% STRUCTURE POUR LE FACTORIAL HMM
% —————————————————————————
%         gamma_1 = [gamma_vec(1,1) 1-gamma_vec(1,1) ; 1-gamma_vec(1,2) gamma_vec(1,2)];
%         gamma_2 = [gamma_vec(1,3) 1-gamma_vec(1,3) ; 1-gamma_vec(1,4) gamma_vec(1,4)];
%         gamma=kron(gamma_1,gamma_2);
% —————————————————————————
end

delta=ones(nbRegime).'/(eye(nbRegime)-gamma_Built+1); % Formule p.68 notes
delta=delta(1,:);
end