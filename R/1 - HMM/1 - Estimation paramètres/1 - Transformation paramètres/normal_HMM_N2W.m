function [parvect] = normal_HMM_N2W(mu,sigma,gamma,type)
% Function to convert from natural parameters to their working equivalents
%   The idea is that by changing the constrained parameters to
%   unconstrained versions, the optimization can be done without
%   constraints
% ———————————————————————————————————————————
% Creation              : 29 novembre 2018
% ———————————————————————————————————————————
% Dernière modification : 25 janvier 2019
% ———————————————————————————————————————————

tsigma=log(sigma);
tmu=mu;
tgamma_vec=[];

switch type
    case 'HMM' % Hidden Markov Model
        nbRegime=length(mu);
        trans_gamma=log(gamma./gamma(:,nbRegime));
        tgamma=trans_gamma(:,1:(nbRegime-1));
        for i=1:nbRegime
            tgamma_vec=[tgamma_vec tgamma(i,:)];
        end
        
    case 'HHMM' % Hierchical Hidden Markov Model
        
% —————————————————————————
% STRUCTURE POUR LE FACTORIAL HMM
% —————————————————————————
%         p1=gamma(1,1)+gamma(1,2);
%         q1=gamma(3,1)+gamma(3,3);
%         p2=gamma(3,3)+gamma(3,4);
%         q2=gamma(2,2)+gamma(2,4);
% 
%         tgamma_vec=log([p1 p2 q1 q2]./(1-[p1 p2 q1 q2]));
% —————————————————————————
%         
        gamma1 = gamma(1:4,:);
        trans_gamma=log(gamma1./gamma1(:,3));
        tgamma=trans_gamma(:,1:2);
        for i=1:4
            tgamma_vec=[tgamma_vec tgamma(i,:)];
        end
        
        gamma2 = gamma(5:8,1:2);
        trans_gamma=log(gamma2./gamma2(:,2));
        tgamma=trans_gamma(:,1);
        for i=1:4
            tgamma_vec=[tgamma_vec tgamma(i,:)];
        end
        
end

parvect = [tmu	tsigma	tgamma_vec];
end

