function [natpar] = normal_HMM_W2N_JAC(parvect)
% Function to convert from WORKING parameters to their NATURAL equivalents
%   Retransforms the WORKING paramaters (unconstrained) to the parameters
%   of the model 
% 覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧�
% Creation              : 6 d馗embre 2018
% 覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧�
% Derni鑽e modification : 6 d馗embre 2018
% 覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧覧�

nparvect=length(parvect);
nbRegime=(-1+sqrt(1+4*nparvect))/2; % R馗up駻ation du nombre de r馮imes du mod鑞e

mu=parvect(1:nbRegime);
sigma=exp(parvect((nbRegime+1):(2*nbRegime)));
gamma_vec=[];

if nbRegime>1
    gamma=exp(parvect((2*nbRegime+1):length(parvect))).';
    gamma=reshape(gamma,[(nbRegime-1),nbRegime]).';
    gamma=[gamma ones(nbRegime,1)];
    
    gammaRows = num2cell(gamma, 2);       % Collect the rows into cells
    rowSums = cellfun(@sum, gammaRows);   % Calculate the sum by row
    
    gamma=gamma./rowSums;
    
    for i=1:nbRegime
        gamma_vec=[gamma_vec gamma(i,:)];
    end

natpar=[mu sigma gamma_vec];
end