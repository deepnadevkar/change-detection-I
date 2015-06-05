clear all close all;

load ILloaddata_hd
%%change AICc depending on parameters!!
sample = 1152; %number of trials
p = 2; % number of parameters

subtotalll = s_total_ll;
%calculate AIC, BIC and BMC
subAIC = (-2*subtotalll + 2*p);
subAIC_ll = (subAIC/-2);
%subAICc = subAIC +(24/(sample)-2);%for SR and VP
subAICc = subAIC +(12/(sample)-1);%for IL and EP

subAICc_ll = (subAICc/-2);
subBIC = (-2*s_total_ll + p*log(sample));
subBIC_ll = (subBIC/-2);
subBMC = s_BMC;

%means and ste of criterions
subtotalll_mean = mean(subtotalll,1);
subtotalll_ste = std(subtotalll,0,1)/sqrt(length(sessions));

subBICll_mean = mean(subBIC_ll,1);
subBICll_ste = std(subBIC_ll,0,1)/sqrt(length(sessions));

subAICll_mean = mean(subAIC_ll,1);
subAICll_ste = std(subAIC_ll,0,1)/sqrt(length(sessions));

subAICcll_mean = mean(subAICc_ll,1);
subAICcll_ste = std(subAICc_ll,0,1)/sqrt(length(sessions));

subBMC_mean = mean(s_BMC,1);
subBMC_ste = std(s_BMC,0,1)/sqrt(length(sessions));

%give
criterions_means = [subAICll_mean subAICcll_mean subBICll_mean subBMC_mean]
criterions_ste = [subAICll_ste subAICcll_ste subBICll_ste subBMC_ste]

% %Calculate overall R-squared and RMSE for set size and CM (uses entire
% %dataset)
% RMSE_SS_all = sqrt(mean((dataN_mean(:)- modelN_mean(:)).^2));
% RMSE_CM_all = sqrt(mean((datadelta_mean(:)- modeldelta_mean(:)).^2));
% 
% Rsquared_SS_all = 1-(sum((dataN_mean(:)- modelN_mean(:)).^2))/sum((dataN_mean(:)- mean(dataN_mean(:))).^2);
% Rsquared_CM_all = 1-(sum((datadelta_mean(:)- modeldelta_mean(:)).^2))/sum((datadelta_mean(:)- mean(datadelta_mean(:))).^2);
% 
% %R-squared and RMSE per subject and then averaged across

for sind = 1:1:length(sessions)
    sind
    
    RMSE_SS(sind) = sqrt(mean((s_perfdata1(sind,:)- s_perfmodel1(sind,:)).^2));
    RMSE_CM(sind) = sqrt(mean((s_perfdata2(sind,:)- s_perfmodel2(sind,:)).^2));
    
    Rsquared_SS (sind)= 1-(sum((s_perfdata1(sind,:)- s_perfmodel1(sind,:)).^2))/sum((s_perfdata1(sind,:)- mean(s_perfdata1(sind,:))).^2);
    Rsquared_CM (sind)= 1-(sum((s_perfdata2(sind,:)- s_perfmodel2(sind,:)).^2))/sum((s_perfdata2(sind,:)- mean(s_perfdata2(sind,:))).^2);
    
end

Rsquared_SS_mean = mean(Rsquared_SS);
Rsquared_SS_ste = std(Rsquared_SS)/sqrt(length(sessions));

Rsquared_CM_mean = mean(Rsquared_CM);
Rsquared_CM_ste = std(Rsquared_CM)/sqrt(length(sessions));

RMSE_SS_mean = mean(RMSE_SS);
RMSE_SS_ste = std(RMSE_SS)/sqrt(length(sessions));

RMSE_CM_mean = mean(RMSE_CM);
RMSE_CM_ste = std(RMSE_CM)/sqrt(length(sessions));

