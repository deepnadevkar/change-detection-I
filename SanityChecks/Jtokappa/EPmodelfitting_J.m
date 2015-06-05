function [perfmodel1_fitted, perfmodel_fitted, alpha_estimate, J1_estimate, temp] = EPmodelfitting_J(m1data)

% Simulate an experiment in which a monkey decides at which of two
% locations an orientation change occurred
% set size (N) and magnitude of change (Delta) were variable
% This code fits the EP model

load EPmodelpred_J

mdata = ' ';
switch mdata
    case {'m1data'}
        Ndata = m1data(:,1);
        deltadata = m1data(:,2);
        corrdata = m1data(:,3);      
    case {'m2data'}
        Ndata = m2data(:,1);
        deltadata = m2data(:,2);
        corrdata = m2data(:,3);
    case {'fakedataEP'}
        Ndata = fakedataEP(:,1);
        deltadata = fakedataEP(:,2);
        corrdata = fakedataEP(:,3);
end

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
    
        loglike = zeros(length(Nvec), length(deltavec));
        for Nind = 1:length(Nvec)
            N = Nvec(Nind);
            for deltaind = 1:length(deltavec)
                delta = deltavec(deltaind);
                idx = find(Ndata == N & deltadata == delta);
                loglike(Nind, deltaind) = sum(corrdata(idx)==1) * log(perfmodel(Nind,deltaind, alphaind,Jind)) + sum(corrdata(idx)==0) * log(1-perfmodel(Nind,deltaind, alphaind,Jind));
            end
        end
        totalloglike (alphaind,Jind) = sum(sum(loglike));

    end
end 

[temp, idx] = max(totalloglike(:));
[alphamax_idx, J1max_idx] = ind2sub(size(totalloglike),idx); %lets you calculate totalloglike for a matrix 
alpha_estimate = alphavec(alphamax_idx);
J1_estimate = Jvec(J1max_idx);

% figure;
% imagesc(kappa1vec,alphavec,totalloglike);
% xlabel('\kappa_1'); ylabel('alpha'); axis xy; colorbar

perfmodel1 = squeeze(mean(perfmodel,2));

perfmodel1_fitted = squeeze(perfmodel1(:,alphamax_idx,J1max_idx));
perfmodel_fitted = squeeze(perfmodel(:,:,alphamax_idx,J1max_idx));

save EPmodelfitting_J
  