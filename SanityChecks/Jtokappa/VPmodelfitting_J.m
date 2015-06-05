function [perfmodel1_fitted, perfmodel_fitted,alpha_estimate,J_bar_estimate,tau_estimate,temp] = VPmodelfitting_J(m1data)

% Simulate an experiment in which a monkey decides at which of two
% locations an orientation change occurred
% set size (N) and magnitude of change (Delta) were variable
% This code fits the VP model


load VPmodelpred_J

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
    case {'fakedataVP'}
        Ndata = fakedataVP(:,1);
        deltadata = fakedataVP(:,2);
        corrdata = fakedataVP(:,3);
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
% Look at the contributions of each of the stimulus condition (defined
% by set size N and change magnitude delta)


totalloglike = zeros(length(alphavec), length(Jvec), length(tauvec));

for Jind = 1:length(Jvec)
    Jind
    for tauind = 1:length(tauvec)
        
        for alphaind = 1:length(alphavec)
            loglike = zeros(length(Nvec), length(deltavec));
            for Nind = 1:length(Nvec)
                N = Nvec(Nind);
                
                for deltaind = 1:length(deltavec)
                    delta = deltavec(deltaind);
                    idx = find(Ndata == N & deltadata == delta);
                    loglike(Nind, deltaind) = sum(corrdata(idx)==1) * log(perfmodel(Nind,deltaind, alphaind,Jind,tauind)) + sum(corrdata(idx)==0) * log(1-perfmodel(Nind,deltaind, alphaind,Jind,tauind));
                end
            end
            totalloglike (alphaind,Jind,tauind) = sum(sum(loglike));
        end
        
    end
    
end

[temp, idx] = max(totalloglike(:));
[alphamax_idx, J_barmax_idx,taumax_idx] = ind2sub(size(totalloglike),idx); %lets you calculate totalloglike for a matrix
alpha_estimate = alphavec(alphamax_idx);
J_bar_estimate = Jvec(J_barmax_idx);
tau_estimate = tauvec(taumax_idx);

perfmodel1 = squeeze(mean(perfmodel,2));
perfmodel1_fitted = squeeze(perfmodel1(:,alphamax_idx,J_barmax_idx,taumax_idx));
perfmodel_fitted = squeeze(perfmodel(:,:,alphamax_idx,J_barmax_idx,taumax_idx));

%figure;
%imagesc(alphavec,kappabarvec, totalloglike(:,:,taumax_idx));
%xlabel('alpha'); ylabel('kappabar'); axis xy; colorbar

%figure;
%imagesc(kappabarvec,tauvec,squeeze(totalloglike(alphamax_idx,:,:)));
%xlabel('\kappabar'); ylabel('tau'); axis xy; colorbar

%figure;
%imagesc(alphavec,tauvec,squeeze(totalloglike(:,kappabarmax_idx,:)));
%xlabel('alpha'); ylabel('tau'); axis xy; colorbar

%figure;
%plot(Nvec, perfdata1); hold on;
%plot(Nvec, perfmodel1(:,alphamax_idx,kappabarmax_idx,taumax_idx),'--');
%xlabel('Set size'); ylabel('Proportion correct');
%ylim([0.5 1])

save VPmodelfitting_J.mat
