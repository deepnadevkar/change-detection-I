clear all close all;

load fakedataILfits_1
load fakedataILfits_2
load fakedataILfits_3
load fakedataILfits_4
load fakedataILfits_5

Nvec = 2:1:5;
deltavec = 10:10:90;

perfdata1 = [data1_1'; data1_2'; data1_3'; data1_4'; data1_5'];
perfdata2(1,:,:) = [data2_1];
perfdata2(2,:,:) = [data2_2];
perfdata2(3,:,:) = [data2_3];
perfdata2(4,:,:) = [data2_4];
perfdata2(5,:,:) = [data2_5];

perfIL1 = [IL1_1'; IL1_2'; IL1_3'; IL1_4'; IL1_5'];
perfIL2(1,:,:) = [IL2_1];
perfIL2(2,:,:) = [IL2_2];
perfIL2(3,:,:) = [IL2_3];
perfIL2(4,:,:) = [IL2_4];
perfIL2(5,:,:) = [IL2_5];

perfEP1 = [EP1_1'; EP1_2'; EP1_3'; EP1_4'; EP1_5'];
perfEP2(1,:,:) = [EP2_1];
perfEP2(2,:,:) = [EP2_2];
perfEP2(3,:,:) = [EP2_3];
perfEP2(4,:,:) = [EP2_4];
perfEP2(5,:,:) = [EP2_5];

perfSR1 = [SR1_1; SR1_2; SR1_3; SR1_4; SR1_5];
perfSR2(1,:,:) = [SR2_1];
perfSR2(2,:,:) = [SR2_2];
perfSR2(3,:,:) = [SR2_3];
perfSR2(4,:,:) = [SR2_4];
perfSR2(5,:,:) = [SR2_5];

perfVP1 = [VP1_1'; VP1_2'; VP1_3'; VP1_4'; VP1_5'];
perfVP2(1,:,:) = [VP2_1];
perfVP2(2,:,:) = [VP2_2];
perfVP2(3,:,:) = [VP2_3];
perfVP2(4,:,:) = [VP2_4];
perfVP2(5,:,:) = [VP2_5];


%% Calculates and plots fake data and model fits
dataN_mean = nanmean(perfdata1,1);
dataN_ste = nanstd(perfdata1,0,1)/sqrt(5);

datadelta_mean = squeeze(nanmean(perfdata2,1));
datadelta_ste = squeeze(nanstd(perfdata2,0,1)/sqrt(5));

ILN_mean = mean(perfIL1,1);
ILN_ste = nanstd(perfIL1,0,1)/sqrt(5);
ILdelta_mean = squeeze(mean(perfIL2,1));
ILdelta_ste = squeeze(nanstd(perfIL2,0,1))/sqrt(5);

%% SS plot
figure; %graphs monkey performance across set sizes with error bars
set(gca,'FontSize',11,'FontName','Arial');
xfac = 1;
yfac = 1;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac yfac]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac yfac]);

patch([Nvec fliplr(Nvec)],[ILN_mean+ILN_ste fliplr(ILN_mean-ILN_ste)],[0.7 0.7 0.7], 'Linestyle','None'); hold on;
errorbar(Nvec, dataN_mean, dataN_ste, 'ok'); 
xlabel('Set size'); ylabel('Proportion correct');axis([1.8 5.2 0.6 1.0])
%ylim([0.5 1])
%xlim([1.8 5.5])
set(gca,'YTick',0.6:.1:1.0)
set(gca,'XTick', 2:1:5)

%Change magnitude plot
% define colors (optional)
colorvec = get(gca, 'ColorOrder');
colorvec = min(colorvec+.65,1);

figure;

% draw shade first
for pp = 1:length(Nvec)
    patch([deltavec, wrev(deltavec)],[ILdelta_mean(pp, :) - ILdelta_ste(pp, :), wrev(ILdelta_mean(pp, :) + ILdelta_ste(pp, :))], colorvec(pp, :),'Linestyle','None');
end
hold on;

% draw error bars on top of it
set(gca,'YTick',0.3:.1:1)
set(gca,'XTick',10:10:90)
errorbar(repmat(deltavec,4,1)', datadelta_mean', datadelta_ste','o'); hold on;
% errorbar(repmat(deltavec, 4, 1)', modeldelta_mean', modeldelta_ste', 'o');
xlabel('Change magnitude'); ylabel('Proportion correct');axis([7 93 0.3 1]);
legend(strcat('N= ',int2str(Nvec')), 4);
%text(40, .2, ['RMSE' '=' num2str("insert value here", 3)]);

%%EP model fits
EPN_mean = mean(perfEP1,1);
EPN_ste = std(perfEP1,0,1)/sqrt(5);
EPdelta_mean = squeeze(mean(perfEP2,1));
EPdelta_ste = squeeze(std(perfEP2,0,1))/sqrt(5);
%% SS plot
figure; %graphs monkey performance across set sizes with error bars
set(gca,'FontSize',11,'FontName','Arial');
xfac = 1;
yfac = 1;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac yfac]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac yfac]);

patch([Nvec fliplr(Nvec)],[EPN_mean+EPN_ste fliplr(EPN_mean-EPN_ste)],[0.7 0.7 0.7], 'Linestyle','None'); hold on;
errorbar(Nvec, dataN_mean, dataN_ste, 'ok'); 
xlabel('Set size'); ylabel('Proportion correct');axis([1.8 5.2 0.6 1.0])
%ylim([0.5 1])
%xlim([1.8 5.5])
set(gca,'YTick',0.6:.1:1.0)
set(gca,'XTick', 2:1:5)
%Change magnitude plot
% define colors (optional)
colorvec = get(gca, 'ColorOrder');
colorvec = min(colorvec+.65,1);

figure;

% draw shade first
for pp = 1:length(Nvec)
    patch([deltavec, wrev(deltavec)],[EPdelta_mean(pp, :) - EPdelta_ste(pp, :), wrev(EPdelta_mean(pp, :) + EPdelta_ste(pp, :))], colorvec(pp, :),'Linestyle','None');
end
hold on;

% draw error bars on top of it
set(gca,'YTick',0.3:.1:1)
set(gca,'XTick',10:10:90)
errorbar(repmat(deltavec,4,1)', datadelta_mean', datadelta_ste','o'); hold on;
% errorbar(repmat(deltavec, 4, 1)', modeldelta_mean', modeldelta_ste', 'o');
xlabel('Change magnitude'); ylabel('Proportion correct');axis([7 93 0.3 1]);
legend(strcat('N= ',int2str(Nvec')), 4);
%text(40, .2, ['RMSE' '=' num2str("insert value here", 3)]);

SRN_mean = mean(perfSR1,1);
SRN_ste = std(perfSR1,0,1)/sqrt(5);
SRdelta_mean = squeeze(mean(perfSR2,1));
SRdelta_ste = squeeze(std(perfSR2,0,1))/sqrt(5);

%% SS plot
figure; %graphs monkey performance across set sizes with error bars
set(gca,'FontSize',11,'FontName','Arial');
xfac = 1;
yfac = 1;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac yfac]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac yfac]);

patch([Nvec fliplr(Nvec)],[SRN_mean+SRN_ste fliplr(SRN_mean-SRN_ste)],[0.7 0.7 0.7], 'Linestyle','None'); hold on;
errorbar(Nvec, dataN_mean, dataN_ste, 'ok'); 
xlabel('Set size'); ylabel('Proportion correct');axis([1.8 5.2 0.6 1.0])
%ylim([0.5 1])
%xlim([1.8 5.5])
set(gca,'YTick',0.6:.1:1.0)
set(gca,'XTick', 2:1:5)

%Change magnitude plot
% define colors (optional)
colorvec = get(gca, 'ColorOrder');
colorvec = min(colorvec+.65,1);

figure;

% draw shade first
for pp = 1:length(Nvec)
    patch([deltavec, wrev(deltavec)],[SRdelta_mean(pp, :) - SRdelta_ste(pp, :), wrev(SRdelta_mean(pp, :) + SRdelta_ste(pp, :))], colorvec(pp, :),'Linestyle','None');
end
hold on;

% draw error bars on top of it
set(gca,'YTick',0.3:.1:1)
set(gca,'XTick',10:10:90)
errorbar(repmat(deltavec,4,1)', datadelta_mean', datadelta_ste','o'); hold on;
% errorbar(repmat(deltavec, 4, 1)', modeldelta_mean', modeldelta_ste', 'o');
xlabel('Change magnitude'); ylabel('Proportion correct');axis([7 93 0.3 1]);
legend(strcat('N= ',int2str(Nvec')), 4);
%text(40, .2, ['RMSE' '=' num2str("insert value here", 3)]);

VPN_mean = mean(perfVP1,1);
VPN_ste = std(perfVP1,0,1)/sqrt(5);
VPdelta_mean = squeeze(mean(perfVP2,1));
VPdelta_ste = squeeze(std(perfVP2,0,1))/sqrt(5);
%% SS plot
figure; %graphs monkey performance across set sizes with error bars
set(gca,'FontSize',11,'FontName','Arial');
xfac = 1;
yfac = 1;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac yfac]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac yfac]);

patch([Nvec fliplr(Nvec)],[VPN_mean+VPN_ste fliplr(VPN_mean-VPN_ste)],[0.7 0.7 0.7], 'Linestyle','None'); hold on;
errorbar(Nvec, dataN_mean, dataN_ste, 'ok'); 
xlabel('Set size'); ylabel('Proportion correct');axis([1.8 5.2 0.6 1.0])
%ylim([0.5 1])
%xlim([1.8 5.5])
set(gca,'YTick',0.6:.1:1.0)
set(gca,'XTick', 2:1:5)%Change magnitude plot
% define colors (optional)
colorvec = get(gca, 'ColorOrder');
colorvec = min(colorvec+.65,1);

figure;

% draw shade first
for pp = 1:length(Nvec)
    patch([deltavec, wrev(deltavec)],[VPdelta_mean(pp, :) - VPdelta_ste(pp, :), wrev(VPdelta_mean(pp, :) + VPdelta_ste(pp, :))], colorvec(pp, :),'Linestyle','None');
end
hold on;

% draw error bars on top of it
set(gca,'YTick',0.3:.1:1)
set(gca,'XTick',10:10:90)
errorbar(repmat(deltavec,4,1)', datadelta_mean', datadelta_ste','o'); hold on;
% errorbar(repmat(deltavec, 4, 1)', modeldelta_mean', modeldelta_ste', 'o');
xlabel('Change magnitude'); ylabel('Proportion correct');axis([7 93 0.3 1]);
legend(strcat('N= ',int2str(Nvec')), 4);
%text(40, .2, ['RMSE' '=' num2str("insert value here", 3)]);