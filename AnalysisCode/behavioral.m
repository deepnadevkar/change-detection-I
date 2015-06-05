clear all; close all;
% Simulate an experiment in which a monkey decides at which of two
% locations an orientation change occurred
% set size (N) and magnitude of change (Delta) were variable
% This code loads and plots the data by sessions and calculates mean and
% standard deviation across sessions

% load dataset
load datasets_m1
sessions = datasets;
   

Nvec = 2:5; %set size
deltavec = 10:10:90; % change magnitude in degrees

s_perfdata1 = zeros(length(sessions),length(Nvec));
s_perfdata2  = zeros(length(sessions),length(Nvec),length(deltavec));
s_perfdata3 = zeros(length(sessions),length(deltavec));

% Loading data and plotting monkey's performance across sessions

for sind = 1:1:length(sessions)
    sind
    
    s_Ndata = sessions{sind}(:,1);
    s_deltadata = sessions{sind}(:,2);
    s_corrdata = sessions{sind}(:,3);
    
    for Nind = 1:length(Nvec)
        N = Nvec(Nind);
        s_perfdata1(sind,Nind) = mean(s_corrdata(s_Ndata == N)); % s_perfdata1 calculates performance for different set sizes across sessions
        
        for deltaind = 1:length(deltavec)
            delta = deltavec(deltaind);
            s_perfdata2(sind,Nind,deltaind) = mean(s_corrdata(s_Ndata == N & s_deltadata == delta)); % s_perfdata2 calculates performance for different set sizes and change in magnitude across sessions
            
        end
    end
    for deltaind = 1:length(deltavec)
        delta = deltavec(deltaind);
        s_perfdata3(sind,deltaind) = mean(s_corrdata(s_deltadata == delta)); % s_perfdata2 calculates performance for different set sizes and change in magnitude across sessions
        
    end
    
end

SS210= s_perfdata2(:,1,1);
SS310= s_perfdata2(:,2,1);
SS410= s_perfdata2(:,3,1);
SS510= s_perfdata2(:,4,1);
SS220= s_perfdata2(:,1,2);
SS320= s_perfdata2(:,2,2);
SS420= s_perfdata2(:,3,2);
SS520= s_perfdata2(:,4,2);
SS230= s_perfdata2(:,1,3);
SS330= s_perfdata2(:,2,3);
SS430= s_perfdata2(:,3,3);
SS530= s_perfdata2(:,4,3);
SS240= s_perfdata2(:,1,4);
SS340= s_perfdata2(:,2,4);
SS440= s_perfdata2(:,3,4);
SS540= s_perfdata2(:,4,4);
SS250= s_perfdata2(:,1,5);
SS350= s_perfdata2(:,2,5);
SS450= s_perfdata2(:,3,5);
SS550= s_perfdata2(:,4,5);
SS260= s_perfdata2(:,1,6);
SS360= s_perfdata2(:,2,6);
SS460= s_perfdata2(:,3,6);
SS560= s_perfdata2(:,4,6);
SS270= s_perfdata2(:,1,7);
SS370= s_perfdata2(:,2,7);
SS470= s_perfdata2(:,3,7);
SS570= s_perfdata2(:,4,7);
SS280= s_perfdata2(:,1,8);
SS380= s_perfdata2(:,2,8);
SS480= s_perfdata2(:,3,8);
SS580= s_perfdata2(:,4,8);
SS290= s_perfdata2(:,1,9);
SS390= s_perfdata2(:,2,9);
SS490= s_perfdata2(:,3,9);
SS590= s_perfdata2(:,4,9);

setsizebymag = [SS210 SS220 SS230 SS240 SS250 SS260 SS270 SS280 SS290 SS310 SS320 SS330 SS340 SS350 SS360 SS370 SS380 SS390 SS410 SS420 SS430 SS440 SS450 SS460 SS470 SS480 SS490 SS510 SS520 SS530 SS540 SS550 SS560 SS570 SS580 SS590]; 

%% Calculates monkey data & runs tests of significance

% dataN_mean = mean(s_perfdata1,1);
% dataN_ste = std(s_perfdata1,0,1)/sqrt(length(sessions));
% [p_N,tbl_N,stats_N] = anova1(s_perfdata1);
% % 
% datadelta_mean = mean(s_perfdata2,1);
% datadelta_ste = std(s_perfdata2,0,1)/sqrt(length(sessions));
% [p_delta,tbl_delta,stats_delta] = anovan(s_perfdata2);
% 

dataN_mean = nanmean(s_perfdata1,1);
dataN_ste = nanstd(s_perfdata1,0,1)/sqrt(length(sessions));

datadelta_mean = squeeze(nanmean(s_perfdata2,1));
datadelta_ste = squeeze(nanstd(s_perfdata2,0,1)/sqrt(length(sessions)));

%% plots the data
figure;
set(gca,'FontSize',11,'FontName','Arial');
xfac = 1;
yfac = 1;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac yfac]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac yfac]);
hold on;
set(gca,'YTick',0.4:0.1:1.0)
set(gca,'XTick',10:10:90)
errorbar(repmat(deltavec,4,1)', datadelta_mean', datadelta_ste','-o'); hold on; 
xlabel('Change magnitude'); ylabel('Proportion correct');axis([7 93 0.4 1.0]);
legend(strcat('N= ',int2str(Nvec')), 4);


