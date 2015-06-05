function [perfmodel1_fitted, perfmodel_fitted,alpha_estimate,J_bar_estimate,tau_estimate,temp, BIC_ll,AIC_ll,BMC] = VPmodelfitting_J(humandata)

% Simulate an experiment in which a monkey decides at which of two
% locations an orientation change occurred
% set size (N) and magnitude of change (Delta) were variable
% This code fits the VP model


load VPmodelpred_J

Ndata = humandata(:,1);
deltadata = humandata(:,2);
corrdata = humandata(:,3);
sample = 11520;

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
% Look at the contributions of each of the stimulus condition (defined
% by set size N and change magnitude delta)


totalloglike = zeros(length(alphavec), length(Jbarvec), length(tauvec));

for Jbarind = 1:length(Jbarvec)
    Jbarind
    for tauind = 1:length(tauvec)
        
        for alphaind = 1:length(alphavec)
            loglike = zeros(length(Nvec), length(deltavec));
            for Nind = 1:length(Nvec)
                N = Nvec(Nind);
                
                for deltaind = 1:length(deltavec)
                    delta = deltavec(deltaind);
                    idx = find(Ndata == N & deltadata == delta);
                    loglike(Nind, deltaind) = sum(corrdata(idx)==1) * log(perfmodel(Nind,deltaind, alphaind,Jbarind,tauind)) + sum(corrdata(idx)==0) * log(1-perfmodel(Nind,deltaind, alphaind,Jbarind,tauind));
                end
            end
            totalloglike (alphaind,Jbarind,tauind) = sum(sum(loglike));
        end
        
    end
    
end

[temp, idx] = max(totalloglike(:));
[alphamax_idx, J_barmax_idx,taumax_idx] = ind2sub(size(totalloglike),idx); %lets you calculate totalloglike for a matrix
alpha_estimate = alphavec(alphamax_idx);
J_bar_estimate = Jbarvec(J_barmax_idx);
tau_estimate = tauvec(taumax_idx);

perfmodel1 = squeeze(mean(perfmodel,2));
perfmodel1_fitted = squeeze(perfmodel1(:,alphamax_idx,J_barmax_idx,taumax_idx));
perfmodel_fitted = squeeze(perfmodel(:,:,alphamax_idx,J_barmax_idx,taumax_idx));

%calculate BIC and AIC
p = 2; % p = number of parameters
AIC = -2*temp + 2*p;
AIC_ll = AIC/-2;
BIC = -2*temp + p*log(sample);
BIC_ll = BIC/-2;

%calculate BMC
alphamax = 3;
alphasteps = 0.33;
Jmax = 60;
Jsteps = 6.67;
taumax = 60;
tausteps = 7.56; 
BMC = log((alphasteps/alphamax)*(Jsteps/Jmax)*(tausteps/taumax))+temp+log(sum(exp(totalloglike(:)-temp)));


% figure;
% imagesc(alphavec,Jbarvec, totalloglike(:,:,taumax_idx));
% xlabel('alpha'); ylabel('Jbar'); axis xy; colorbar
% 
% figure;
% imagesc(Jbarvec,tauvec,squeeze(totalloglike(alphamax_idx,:,:)));
% xlabel('Jbar'); ylabel('tau'); axis xy; colorbar
% 
% figure;
% imagesc(alphavec,tauvec,squeeze(totalloglike(:,J_barmax_idx,:)));
% xlabel('alpha'); ylabel('tau'); axis xy; colorbar

% figure;
% plot(Nvec, perfdata1); hold on;
% plot(Nvec, perfmodel1(:,alphamax_idx,J_barmax_idx,taumax_idx),'--');
% xlabel('Set size'); ylabel('Proportion correct');
% ylim([0.5 1])

