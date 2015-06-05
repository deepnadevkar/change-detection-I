%% This script loads and plots the human behavioral data and model fits of each of the 5 models. 
 % It also calculates R-squared values and model comparison metrics such as
 % AIC, AICc, BIC, and BMC (log bayes factor)

clear all; close all;

load datasets_h
sessions = struct2cell(humandata.data);

load modelfitting_H
    
Nvec = 2:5;
deltavec = 10:10:90; % change magnitude in degrees

perfdata1 = zeros(length(sessions),length(Nvec));
perfdata2  = zeros(length(sessions),length(Nvec),length(deltavec));

%% Loading data and plotting monkey's performance across sessions

for sind = 1:1:length(sessions)
    sind
    
    Ndata = sessions{sind}(:,1);
    deltadata = sessions{sind}(:,2);
    corrdata = sessions{sind}(:,3);
    
    for Nind = 1:length(Nvec)
        N = Nvec(Nind);
        perfdata1(sind,Nind) = mean(corrdata(Ndata == N)); 
        for deltaind = 1:length(deltavec)
            delta = deltavec(deltaind);
            perfdata2(sind,Nind,deltaind) = mean(corrdata(Ndata == N & deltadata == delta)); 
        end
    end
    
end  
%criterion diffs from VP model
AIC_diff = [AIC_all(:,1)-AIC_all(:,4) AIC_all(:,2)-AIC_all(:,4) AIC_all(:,3)-AIC_all(:,4) AIC_all(:,5)-AIC_all(:,4)];
AICc_diff = [AICc_all(:,1)-AICc_all(:,4) AICc_all(:,2)-AICc_all(:,4) AICc_all(:,3)-AICc_all(:,4)  AICc_all(:,5)-AICc_all(:,4)];
BIC_diff = [BIC_all(:,1)-BIC_all(:,4) BIC_all(:,2)-BIC_all(:,4) BIC_all(:,3)-BIC_all(:,4) BIC_all(:,5)-BIC_all(:,4)];
BMC_diff = [BMC_all(:,1)-BMC_all(:,4) BMC_all(:,2)-BMC_all(:,4) BMC_all(:,3)-BMC_all(:,4) BMC_all(:,5)-BMC_all(:,4)];

%% means 
dataN_mean = nanmean(perfdata1,1);
datadelta_mean = squeeze(nanmean(perfdata2,1));
IL1_mean = mean(IL1_all,1);
IL2_mean = squeeze(mean(IL2_all,1));
EP1_mean = mean(EP1_all,1);
EP2_mean = squeeze(mean(EP2_all,1));
EPF1_mean = mean(EPF1_all,1);
EPF2_mean = squeeze(mean(EPF2_all,1));
VP1_mean = mean(VP1_all,1);
VP2_mean = squeeze(mean(VP2_all,1));
VPF1_mean = mean(VPF1_all,1);
VPF2_mean = squeeze(mean(VPF2_all,1));
pars_IL_mean = mean(pars_IL_all,1);
pars_EP_mean = mean(pars_EP_all,1);
pars_EPF_mean = mean(pars_EPF_all,1);
pars_VP_mean = mean(pars_VP_all,1);
pars_VPF_mean = mean(pars_VPF_all,1);
BMC_mean = mean(BMC_all,1);
BIC_mean = mean(BIC_all,1);
AIC_mean = mean(AIC_all,1);
AICc_mean = mean(AICc_all,1); 
Rs_IL_mean = mean(Rs_IL_all,1);
Rs_EP_mean = mean(Rs_EP_all,1);
Rs_EPF_mean = mean(Rs_EPF_all,1);
Rs_VP_mean = mean(Rs_VP_all,1);
Rs_VPF_mean = mean(Rs_VPF_all,1);
AIC_meandiff = mean(AIC_diff,1);
AICc_meandiff = mean(AICc_diff,1);
BIC_meandiff = mean(BIC_diff,1);
BMC_meandiff = mean(BMC_diff,1);


%% standard errors
   
    dataN_ste = nanstd(perfdata1,0,1)/sqrt(length(sessions));
    datadelta_ste = squeeze(nanstd(perfdata2,0,1))/sqrt(length(sessions));
    IL1_ste = std(IL1_all,0,1)/sqrt(length(sessions));
    IL2_ste = squeeze(std(IL2_all,0,1))/sqrt(length(sessions));
    EP1_ste = std(EP1_all,0,1)/sqrt(length(sessions));
    EP2_ste = squeeze(std(EP2_all,0,1))/sqrt(length(sessions));
    EPF1_ste = std(EPF1_all,0,1)/sqrt(length(sessions));
    EPF2_ste = squeeze(std(EPF2_all,0,1))/sqrt(length(sessions));
    VP1_ste = std(VP1_all,0,1)/sqrt(length(sessions));
    VP2_ste = squeeze(std(VP2_all,0,1))/sqrt(length(sessions));
    VPF1_ste = std(VPF1_all,0,1)/sqrt(length(sessions));
    VPF2_ste = squeeze(std(VPF2_all,0,1))/sqrt(length(sessions));
    pars_IL_ste = std(pars_IL_all,0,1)/sqrt(length(sessions));
    pars_EP_ste = std(pars_EP_all,0,1)/sqrt(length(sessions));
    pars_EPF_ste = std(pars_EPF_all,0,1)/sqrt(length(sessions));
    pars_VP_ste = std(pars_VP_all,0,1)/sqrt(length(sessions));
    pars_VPF_ste = std(pars_VPF_all,0,1)/sqrt(length(sessions));
    BMC_ste = std(BMC_all,0,1)/sqrt(length(sessions));
    BIC_ste = std(BIC_all,0,1)/sqrt(length(sessions));
    AIC_ste = std(AIC_all,0,1)/sqrt(length(sessions));
    AICc_ste = std(AICc_all,0,1)/sqrt(length(sessions));
    BMC_diffste = std(BMC_diff,0,1)/sqrt(length(sessions));
    BIC_diffste = std(BIC_diff,0,1)/sqrt(length(sessions));
    AIC_diffste = std(AIC_diff,0,1)/sqrt(length(sessions));
    AICc_diffste = std(AICc_diff,0,1)/sqrt(length(sessions));
    Rs_IL_ste = std(Rs_IL_all,0,1)/sqrt(length(sessions));
    Rs_EP_ste = std(Rs_EP_all,0,1)/sqrt(length(sessions));
    Rs_EPF_ste = std(Rs_EPF_all,0,1)/sqrt(length(sessions));
    Rs_VP_ste = std(Rs_VP_all,0,1)/sqrt(length(sessions));
    Rs_VPF_ste = std(Rs_VPF_all,0,1)/sqrt(length(sessions));
    

% behavioral plots

figure;
set(gca,'FontSize',11,'FontName','Arial');
xfac = 1;
yfac = 1;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac yfac]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac yfac]);
hold on;
set(gca,'YTick',0.4:0.2:1.0)
set(gca,'XTick',10:10:90)
errorbar(repmat(deltavec,4,1)', datadelta_mean', datadelta_ste','-o'); hold on; 
xlabel('Change magnitude'); ylabel('Proportion correct');axis([7 93 0.4 1.0]);
legend(strcat('N= ',int2str(Nvec')), 4);

%% SS plots

%IL 
figure; %graphs monkey performance across set sizes with error bars
set(gca,'FontSize',11,'FontName','Arial');
xfac = 1;
yfac = 1;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac yfac]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac yfac]);
patch([Nvec fliplr(Nvec)],[IL1_mean+IL1_ste fliplr(IL1_mean-IL1_ste)],[0.7 0.7 0.7], 'Linestyle','None'); hold on;
errorbar(Nvec, dataN_mean, dataN_ste, 'ok');
xlabel('Set size'); ylabel('Proportion correct');axis([1.8 5.2 0.4 1.0])
set(gca,'YTick',0.4:.2:1.0)
set(gca,'XTick', 2:1:5)

%EP 
figure; %graphs monkey performance across set sizes with error bars
set(gca,'FontSize',11,'FontName','Arial');
xfac = 1;
yfac = 1;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac yfac]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac yfac]);
patch([Nvec fliplr(Nvec)],[EP1_mean+EP1_ste fliplr(EP1_mean-EP1_ste)],[0.7 0.7 0.7], 'Linestyle','None'); hold on;
errorbar(Nvec, dataN_mean, dataN_ste, 'ok');
xlabel('Set size'); ylabel('Proportion correct');axis([1.8 5.2 0.4 1.0])
set(gca,'YTick',0.4:.2:1.0)
set(gca,'XTick', 2:1:5)

%EPF 
figure; %graphs monkey performance across set sizes with error bars
set(gca,'FontSize',11,'FontName','Arial');
xfac = 1;
yfac = 1;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac yfac]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac yfac]);
patch([Nvec fliplr(Nvec)],[EPF1_mean+EPF1_ste fliplr(EPF1_mean-EPF1_ste)],[0.7 0.7 0.7], 'Linestyle','None'); hold on;
errorbar(Nvec, dataN_mean, dataN_ste, 'ok');
xlabel('Set size'); ylabel('Proportion correct');axis([1.8 5.2 0.4 1.0])
set(gca,'YTick',0.4:.2:1.0)
set(gca,'XTick', 2:1:5)

%VP 
figure; %graphs monkey performance across set sizes with error bars
set(gca,'FontSize',11,'FontName','Arial');
xfac = 1;
yfac = 1;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac yfac]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac yfac]);
patch([Nvec fliplr(Nvec)],[VP1_mean+VP1_ste fliplr(VP1_mean-VP1_ste)],[0.7 0.7 0.7], 'Linestyle','None'); hold on;
errorbar(Nvec, dataN_mean, dataN_ste, 'ok');
xlabel('Set size'); ylabel('Proportion correct');axis([1.8 5.2 0.4 1.0])
set(gca,'YTick',0.4:.2:1.0)
set(gca,'XTick', 2:1:5)

%VPF 
figure; %graphs monkey performance across set sizes with error bars
set(gca,'FontSize',11,'FontName','Arial');
xfac = 1;
yfac = 1;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac yfac]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac yfac]);
patch([Nvec fliplr(Nvec)],[VPF1_mean+VPF1_ste fliplr(VPF1_mean-VPF1_ste)],[0.7 0.7 0.7], 'Linestyle','None'); hold on;
errorbar(Nvec, dataN_mean, dataN_ste, 'ok');
xlabel('Set size'); ylabel('Proportion correct');axis([1.8 5.2 0.4 1.0])
set(gca,'YTick',0.4:.2:1.0)
set(gca,'XTick', 2:1:5)



%% Change magnitude plot
colorvec = get(gca, 'ColorOrder');
colorvec = min(colorvec+.65,1);

%IL
figure;
set(gca,'FontSize',11,'FontName','Arial');
xfac = 1;
yfac = 1;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac yfac]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac yfac]);
for pp = 1:length(Nvec)
    patch([deltavec, wrev(deltavec)],[IL2_mean(pp, :) - IL2_ste(pp, :), wrev(IL2_mean(pp, :) + IL2_ste(pp, :))], colorvec(pp, :),'Linestyle','None');
end
hold on;
set(gca,'YTick',0.4:.2:1.0)
set(gca,'XTick',0:30:90)
errorbar(repmat(deltavec,4,1)', datadelta_mean', datadelta_ste','o'); hold on;
xlabel('Change magnitude'); ylabel('Proportion correct');axis([7 93 0.4 1.0]);
legend(strcat('N= ',int2str(Nvec')), 4);

%EP
figure;
set(gca,'FontSize',11,'FontName','Arial');
xfac = 1;
yfac = 1;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac yfac]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac yfac]);
for pp = 1:length(Nvec)
    patch([deltavec, wrev(deltavec)],[EP2_mean(pp, :) - EP2_ste(pp, :), wrev(EP2_mean(pp, :) + EP2_ste(pp, :))], colorvec(pp, :),'Linestyle','None');
end
hold on;
set(gca,'YTick',0.4:.2:1.0)
set(gca,'XTick',0:30:90)
errorbar(repmat(deltavec,4,1)', datadelta_mean', datadelta_ste','o'); hold on;
xlabel('Change magnitude'); ylabel('Proportion correct');axis([7 93 0.4 1.0]);
legend(strcat('N= ',int2str(Nvec')), 4);

%EPF
figure;
set(gca,'FontSize',11,'FontName','Arial');
xfac = 1;
yfac = 1;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac yfac]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac yfac]);
for pp = 1:length(Nvec)
    patch([deltavec, wrev(deltavec)],[EPF2_mean(pp, :) - EPF2_ste(pp, :), wrev(EPF2_mean(pp, :) + EPF2_ste(pp, :))], colorvec(pp, :),'Linestyle','None');
end
hold on;
set(gca,'YTick',0.4:.2:1.0)
set(gca,'XTick',0:30:90)
errorbar(repmat(deltavec,4,1)', datadelta_mean', datadelta_ste','o'); hold on;
xlabel('Change magnitude'); ylabel('Proportion correct');axis([7 93 0.4 1.0]);
legend(strcat('N= ',int2str(Nvec')), 4);

%VP
figure;
set(gca,'FontSize',11,'FontName','Arial');
xfac = 1;
yfac = 1;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac yfac]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac yfac]);
for pp = 1:length(Nvec)
    patch([deltavec, wrev(deltavec)],[VP2_mean(pp, :) - VP2_ste(pp, :), wrev(VP2_mean(pp, :) + VP2_ste(pp, :))], colorvec(pp, :),'Linestyle','None');
end
hold on;
set(gca,'YTick',0.4:.2:1.0)
set(gca,'XTick',0:30:90)
errorbar(repmat(deltavec,4,1)', datadelta_mean', datadelta_ste','o'); hold on;
xlabel('Change magnitude'); ylabel('Proportion correct');axis([7 93 0.4 1.0]);
legend(strcat('N= ',int2str(Nvec')), 4);

%VPF
figure;
set(gca,'FontSize',11,'FontName','Arial');
xfac = 1;
yfac = 1;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac yfac]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac yfac]);
for pp = 1:length(Nvec)
    patch([deltavec, wrev(deltavec)],[VPF2_mean(pp, :) - VPF2_ste(pp, :), wrev(VPF2_mean(pp, :) + VPF2_ste(pp, :))], colorvec(pp, :),'Linestyle','None');
end
hold on;
set(gca,'YTick',0.4:.2:1.0)
set(gca,'XTick',0:30:90)
errorbar(repmat(deltavec,4,1)', datadelta_mean', datadelta_ste','o'); hold on;
xlabel('Change magnitude'); ylabel('Proportion correct');axis([7 93 0.4 1.0]);
legend(strcat('N= ',int2str(Nvec')), 4);


% BMC plot
figure
set(gca,'FontSize',11,'FontName','Arial');
xfac = 1;
yfac = 1;
set(gcf,'Position',get(gcf,'Position').*[.1 .1 xfac yfac]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[.1 .1 xfac yfac]);

mnames = {'IL','EP','EPF','VPF'};
bar(1:4,BMC_meandiff,'k');
hold on
errorbar(1:4,BMC_meandiff,BMC_diffste,'k','LineStyle','none');
set(gca,'XTick',1:4,'XTickLabel',mnames);
set(gca,'YTick',-200:25:0);
ylabel('Model log likelihood relative to VP');


%save loadbysession_H