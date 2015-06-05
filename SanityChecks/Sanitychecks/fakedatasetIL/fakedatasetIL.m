Nvec = 2:5;
deltavec = 10:10:90; % change magnitude in degrees
Kvec = 1:5;
epsvec = 0:0.01:0.3;
ndatasets = 10;


epsest = zeros(1,ndatasets);
Kest = zeros(1,ndatasets);
totalloglike = zeros(ndatasets);
perfdata1 = zeros(ndatasets,length(Nvec));
perfdata2  = zeros(ndatasets, length(Nvec),length(deltavec));


for i = 1:ndatasets
    
    eps(i) = epsvec(randi(length(epsvec),1));   % randomly generating eps from epsvec
    K(i)= Kvec(randi(length(Kvec),1)); % randomly generating K from Kvec
    
    fakedataIL = ILfakedata(eps(i),K(i));
    
    modelfit = 'SR';
    switch modelfit
        case {'IL'}
            [perfmodel1_fitted, perfmodel_fitted, eps_estimate ,K_estimate, temp]  = ILmodelfitting(fakedataIL)
            epsest(i) = eps_estimate;
            Kest(i) = K_estimate;
        case {'EP'}
            [perfmodel1_fitted, perfmodel_fitted, alpha_estimate, kappa1_estimate, temp] = EPmodelfitting(fakedataIL)
        case {'SR'}
            [perfmodel1_fitted, perfmodel_fitted, K_estimate, alpha_estimate, kappa1_estimate,temp] = SRmodelfitting(fakedataIL)
        case {'VP'}
            [perfmodel1_fitted, perfmodel_fitted,alpha_estimate,kappabar_estimate,tau_estimate,temp] = VPmodelfitting(fakedataIL)
    end
    
    
    totalloglike (i) = temp;
    Ndata = fakedataIL(:,1);
    deltadata = fakedataIL(:,2);
    corrdata = fakedataIL(:,3);
    
    %% calculating and plotting fakedata performance across set sizes and change magnitudes
    
    
    for Nind = 1:length(Nvec)
        N = Nvec(Nind);
        perfdata1(i, Nind) = mean(corrdata(Ndata == N)); % perfdata1 calculates performance for different set sizes
        
        for deltaind = 1:length(deltavec)
            delta = deltavec(deltaind);
            perfdata2(i,Nind,deltaind) = mean(corrdata(Ndata == N & deltadata == delta)); % perfdata2 calculates performance for different set sizes and change in magnitude
        end
        
        perfmodel1(i,:) = perfmodel1_fitted;
        perfmodel2(i,:,:) = perfmodel_fitted;
        
    end
    
    
end

%% Calculates and plots monkey data and IL model fits
dataN_mean = nanmean(perfdata1,1);
dataN_ste = nanstd(perfdata1,0,1)/sqrt(length(ndatasets));

datadelta_mean = squeeze(nanmean(perfdata2,1));
datadelta_ste = squeeze(nanstd(perfdata2,0,1)/sqrt(length(ndatasets)));

modelN_mean = mean(perfmodel1,1);
modelN_ste = std(perfmodel1,0,1)/sqrt(length(ndatasets));

modeldelta_mean = squeeze(mean(perfmodel2,1));
modeldelta_ste = squeeze(std(perfmodel2,0,1))/sqrt(length(ndatasets));

RMSE = sqrt(mean((datadelta_mean(:)- modeldelta_mean(:)).^2));
%Rsquared = 1-((datadelta_mean(:)- modeldelta_mean(:)).^2)/(((datadelta_mean(:)- mean(datadelta_mean(:))).^2));

SE = ((datadelta_mean(:)- modeldelta_mean(:)).^2);
SSE = sum(SE);
ST = ((datadelta_mean(:)- mean(datadelta_mean(:))).^2);
SST = sum(ST);
SSEoverSST = SSE/SST;
Rsquared = 1 - SSEoverSST;


%% Set size plot -- temporary
figure; %graphs monkey performance across set sizes with error bars
errorbar(Nvec, dataN_mean, dataN_ste); hold on;
errorbar(Nvec, modelN_mean, modelN_ste ,'--');
xlabel('Set size'); ylabel('Proportion correct');axis([1.8 5.2 0.4 1.0])
%ylim([0.5 1])
%xlim([1.8 5.5])
set(gca,'YTick',0.4:.1:1.0)
set(gca,'XTick', 2:1:5)

%% Change magnitude plot
% define colors (optional)
colorvec = get(gca, 'ColorOrder');
colorvec = min(colorvec+.65,1);

figure;

% draw shade first
for pp = 1:length(Nvec)
    patch([deltavec, wrev(deltavec)],[modeldelta_mean(pp, :) - modeldelta_ste(pp, :), wrev(modeldelta_mean(pp, :) + modeldelta_ste(pp, :))], colorvec(pp, :),'Linestyle','None');
end
hold on;

% draw error bars on top of it
set(gca,'YTick',0.4:.1:1)
set(gca,'XTick',10:10:90)
errorbar(repmat(deltavec,4,1)', datadelta_mean', datadelta_ste','o'); hold on;
% errorbar(repmat(deltavec, 4, 1)', modeldelta_mean', modeldelta_ste', 'o');
xlabel('Change magnitude'); ylabel('Proportion correct');axis([7 93 0.4 1]);
legend(strcat('N= ',int2str(Nvec')), 4);
%text(40, .2, ['RMSE' '=' num2str("insert value here", 3)]);


% % Parameter recovery plots
% figure;
% scatter(eps, epsest)
% xlabel('eps'); ylabel('eps estimate');
% 
% hold on;
% plot([0 0.3], [0 0.3], 'k--');% plots a diagonal with 0 to 3 being the range of alpha
% 
% figure;
% scatter(K, Kest)
% xlabel('K'); ylabel('K estimate');
% 
% hold on;
% plot([1 5], [1 5], 'k--'); % plots a diagonal with 0 to 10 being the range of kappa1
% 

