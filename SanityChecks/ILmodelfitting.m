function [perfmodel1_fitted, perfmodel_fitted, eps_estimate ,K_estimate, temp, BIC_ll,AIC_ll,BMC]  = ILmodelfitting(humandata)


load ILmodelpred

% Simulate an experiment in which a monkey decides at which of two
% locations an orientation change occurred
% set size (N) and magnitude of change (Delta) were variable
% This code fits the IL model

Nvec = 2:5;
deltavec = 10:10:90;
Ndata = humandata(:,1);
deltadata = humandata(:,2);
corrdata = humandata(:,3); 

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

%% IL model fitting

for epsind = 1:length(epsvec)
    eps = epsvec(epsind);
    
    for Kind = 1:length(Kvec)
        K = Kvec(Kind);
        
        loglike = zeros(length(Nvec),1);
        for Nind = 1:length(Nvec)
            N = Nvec(Nind);
            idx = find(Ndata == N);
            loglike(Nind) = sum(corrdata(idx)==1) * log(perfmodel(Nind,epsind,Kind)) + sum(corrdata(idx)==0) * log(1-perfmodel(Nind,epsind,Kind));
        end
        totalloglike (epsind,Kind) = sum(loglike);
    end
end

[temp, idx] = max(totalloglike(:));
[epsmax_idx,Kmax_idx] = ind2sub(size(totalloglike),idx); %lets you calculate totalloglike for a matrix
eps_estimate = epsvec(epsmax_idx);
K_estimate = Kvec(Kmax_idx);

%calculate BIC and AIC
p = 2; % number of parameters
sample = 11520; % number of trials
AIC = -2*temp + 2*p;
AIC_ll = AIC/-2;
BIC = -2*temp + p*log(sample);
BIC_ll = BIC/-2;

%calculate BMC
epsmax = 0.3;
epssteps = 0.01;
Kmax = 5;
Ksteps = 1;

BMC = log((epssteps/epsmax)*(Ksteps/Kmax))+temp+log(sum(exp(totalloglike(:)-temp)));


perfmodel1 = squeeze(mean(perfmodel,1));

perfmodel_fitted = perfmodel1(epsmax_idx,Kmax_idx);
perfmodel1_fitted = squeeze(perfmodel(:,epsmax_idx,Kmax_idx));
perfmodel_fitted = repmat(perfmodel1_fitted,1,9);

