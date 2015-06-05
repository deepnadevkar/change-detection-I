function [perfmodel1_fitted, perfmodel_fitted, K_estimate, alpha_estimate, J1_estimate,temp,BIC_ll,AIC_ll,BMC] = SRmodelfitting_J(humandata)

% Simulate an experiment in which a monkey decides at which of two
% locations an orientation change occurred
% set size (N) and magnitude of change (Delta) were variable
% This code fits the SR model

load SRmodelpred_J

Ndata = humandata(:,1);
deltadata = humandata(:,2);
corrdata = humandata(:,3);
    
Nvec = 2:5;
deltavec = 10:10:90;
perfdata1 = zeros(length(Nvec),1);
perfdata2 = zeros(length(Nvec),length(deltavec));

%% Calculating monkey performance
for Nind = 1:length(Nvec)
    N = Nvec(Nind);
    perfdata1(Nind) = mean(corrdata(Ndata == N));
    for deltaind = 1:length(deltavec)
            delta = deltavec(deltaind);
            perfdata2(Nind,deltaind) = mean(corrdata(find(Ndata == N & deltadata == delta)));
    end
end

%% Fitting model to data (by maximizing likelihood)
% Look at the contributions of each of the stimulus condition (defined by set size N and change magnitude delta)
  
for alphaind = 1:length(alphavec)
    alpha = alphavec(alphaind);
    for Jind = 1:length(Jvec)
        J1 = Jvec(Jind);
        for Kind = 1:length(Kvec)
            K = Kvec(Kind);
    
        loglike = zeros(length(Nvec), length(deltavec));
        for Nind = 1:length(Nvec)
            N = Nvec(Nind);
            for deltaind = 1:length(deltavec)
                delta = deltavec(deltaind);
                idx = find(Ndata == N & deltadata == delta);
                loglike(Nind, deltaind) = sum(corrdata(idx)==1) * log(perfmodel(Kind, Nind, alphaind, Jind, deltaind)) + sum(corrdata(idx)==0) * log(1-perfmodel(Kind, Nind, alphaind, Jind, deltaind));
            end
        end
        totalloglike (Kind,alphaind,Jind) = sum(sum(loglike));
        end
    end
end 

[temp, idx] = max(totalloglike(:));
[Kmax_idx,alphamax_idx, J1max_idx] = ind2sub(size(totalloglike),idx); %lets you calculate totalloglike for a matrix 
K_estimate = Kvec(Kmax_idx);
alpha_estimate = alphavec(alphamax_idx);
J1_estimate = Jvec(J1max_idx);

%calculate BIC and AIC
p = 2; % number of parameters
sample = 11520; %number of trials
AIC = -2*temp + 2*p;
AIC_ll = AIC/-2;
BIC = -2*temp + p*log(sample);
BIC_ll = BIC/-2;

%calculate BMC
Kmax = 5;
Ksteps = 1;
alphamax = 3;
alphasteps = 0.5;
Jmax = 25;
Jsteps = 2.5;

BMC = log((Ksteps/Kmax)*(alphasteps/alphamax)*(Jsteps/Jmax))+temp+log(sum(exp(totalloglike(:)-temp)));

perfmodel1 = squeeze(mean(perfmodel,5));
perfmodel1_fitted = squeeze(perfmodel1(Kmax_idx,:,alphamax_idx,J1max_idx));
perfmodel_fitted = squeeze(perfmodel(Kmax_idx,:,alphamax_idx,J1max_idx,:));

%save SRmodelfitting_J
