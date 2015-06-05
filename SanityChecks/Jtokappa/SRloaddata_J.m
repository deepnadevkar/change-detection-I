clear all; close all;
% Simulate an experiment in which a monkey decides at which of two
% locations an orientation change occurred
% set size (N) and magnitude of change (Delta) were variable
% This code loads and plots monkey data by sessions and calculates mean and
% standard deviation across sessions 

mdata = ' ';
switch mdata
    case {'m1data'}
        load m1databysession
        m1data.data.x1051(147,4)=1; % correct a manual NaN problem
        sessions = struct2cell(m1data.data);
    case {'m2data'}
        load m2databysession
        sessions = struct2cell(m2data.data);
        
end

Nvec = 2:5;
deltavec = 10:10:90; % change magnitude in degrees

s_perfdata1 = zeros(length(sessions),length(Nvec));
s_perfdata2  = zeros(length(sessions),length(Nvec),length(deltavec));

%% Loading data and plotting monkey's performance across sessions

for sind = 1:1:length(sessions)
    sind
    
    s_Ndata = sessions{sind}(:,2);
    s_deltadata = sessions{sind}(:,3);
    s_corrdata = sessions{sind}(:,4);

    for Nind = 1:length(Nvec)
        N = Nvec(Nind);
        s_perfdata1(sind,Nind) = mean(s_corrdata(s_Ndata == N)); % s_perfdata1 calculates performance for different set sizes across sessions
        
        for deltaind = 1:length(deltavec)
            delta = deltavec(deltaind);
            s_perfdata2(sind,Nind,deltaind) = mean(s_corrdata(s_Ndata == N & s_deltadata == delta)); % s_perfdata2 calculates performance for different set sizes and change in magnitude across sessions
            
        end
    end
    [temp1, temp2, temp3, temp4, temp5] = SRmodelfitting_J(sessions{sind}(:,2:4));
    s_perfmodel1(sind,:) = temp1;
    s_perfmodel2(sind,:,:) = temp2;
    s_K(sind,:,:) = temp3;
    s_alpha(sind,:,:) = temp4;
    s_J1(sind,:,:) = temp5;
end

%% Calculates and plots monkey data and EP model fits 
dataN_mean = nanmean(s_perfdata1,1);
dataN_ste = nanstd(s_perfdata1,0,1)/sqrt(length(sessions));

datadelta_mean = squeeze(nanmean(s_perfdata2,1));
datadelta_ste = squeeze(nanstd(s_perfdata2,0,1)/sqrt(length(sessions)));

modelN_mean = mean(s_perfmodel1,1);
modelN_ste = std(s_perfmodel1,0,1)/sqrt(length(sessions));

modeldelta_mean = squeeze(mean(s_perfmodel2,1));
modeldelta_ste = squeeze(std(s_perfmodel2,0,1))/sqrt(length(sessions));

K_mean = mean(s_K,1);
K_ste = std(s_K,0,1)/sqrt(length(sessions));

alpha_mean = mean(s_alpha,1);
alpha_ste = std(s_alpha,0,1)/sqrt(length(sessions));

J1_mean = mean(s_J1,1);
J1_ste = std(s_J1,0,1)/sqrt(length(sessions));

RMSE = sqrt(mean((datadelta_mean(:)- modeldelta_mean(:)).^2));
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
xlabel('Set size'); ylabel('Proportion correct');axis([1.8 5.2 0.5 0.8])
%ylim([0.5 1])
%xlim([1.8 5.5])
set(gca,'YTick',0.4:.1:0.8)
set(gca,'XTick', 2:1:5)

%figure; %graphs monkey performance across set sizes and change magnitude with error bars
%errorbar(repmat(deltavec,4,1)', datadelta_mean', datadelta_ste'); hold on; 
%errorbar(repmat(deltavec,4,1)', modeldelta_mean', modeldelta_ste' ,'--');
%xlabel('Change magnitude'); ylabel('Proportion correct');

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

% draw error bars on top of it (if only run this part, it graphs behavioral
% data only)
set(gca,'YTick',0.4:.1:1)
set(gca,'XTick',10:10:90)
errorbar(repmat(deltavec,4,1)', datadelta_mean', datadelta_ste','o'); hold on; 
% errorbar(repmat(deltavec, 4, 1)', modeldelta_mean', modeldelta_ste', 'o');
xlabel('Change magnitude'); ylabel('Proportion correct');axis([7 93 0.4 1]);
legend(strcat('N= ',int2str(Nvec')), 4);
%text(40, .2, ['RMSE' '=' num2str("insert value here", 3)]);

%% Set size plot 
% define colors (optional)
colorvec = get(gca, 'ColorOrder');
colorvec = min(colorvec+.65,1);

figure;

% draw shade first
for pp = 1:length(Nvec)
    patch([Nvec, wrev(Nvec)],[modelN_mean(pp, :) - modelN_ste(pp, :), wrev(modelN_mean(pp, :) + modelN_ste(pp, :))], colorvec(pp, :),'Linestyle','None');
end
hold on;

% draw error bars on top of it
set(gca,'YTick',0:.2:1)
set(gca,'XTick',2:1:5)
errorbar(Nvec', dataN_mean', dataN_ste', 'o'); 
%errorbar(Nvec, modelN_mean, modelN_ste);
xlabel('Set size'); ylabel('Proportion correct');axis([2 5 0 1]);
legend(strcat('N= ',int2str(Nvec')), 4);
%text(40, .2, ['RMSE' '=' num2str(RMSE(1), 3)]);

save SRloaddata_J
