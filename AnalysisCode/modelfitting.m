%% This function runs the model fits for each of the 5 models using the corresponding model's predictions. 
% The model predictions are fit to the behavioral data of humans and monkeys 
%(Note: the model parameters are different for humans and monkeys.You can comment and uncomment accordingly)

function [IL1, IL2, EP1, EP2, SR1, SR2, VP1, VP2, VPF1, VPF2, pars_IL, pars_EP, pars_SR, pars_VP, pars_VPF, LL, BMC, BIC, AIC, AICc, Rs_IL, Rs_EP, Rs_SR, Rs_VP Rs_VPF] = modelfitting(humandata)

%for humandata
%load humandata
load ILmodelpred_H
load EPmodelpred_H
load SRmodelpred_H
load VPmodelpred_H
load VPFmodelpred_H
Ndata = humandata(:,1);
deltadata = humandata(:,2);
corrdata = humandata(:,3);

% %for monkey data
% load EPmodelpred_M
% load ILmodelpred_M
% load SRmodelpred_M
% load VPmodelpred_M
% load VPFmodelpred_M
% Ndata = m3data(:,1);
% deltadata = m3data(:,2);
% corrdata = m3data(:,3);

Nvec = 2:5;
deltavec = 10:10:90;

for Nind = 1:length(Nvec)
    N = Nvec(Nind);
    idx = find(Ndata == N);
    corrdataidx_1_N(Nind)= sum(corrdata(idx) ==1);
    corrdataidx_0_N(Nind)= sum(corrdata(idx) ==0);
    perfdata1(Nind) = mean(corrdata(Ndata == N)); 
    for deltaind = 1:length(deltavec)
        delta = deltavec(deltaind);
        idx2 = find(Ndata == N & deltadata == delta);
        corrdataidx_1(Nind, deltaind) = sum(corrdata(idx2)==1);
        corrdataidx_0(Nind, deltaind)= sum(corrdata(idx2)==0);
        perfdata2(Nind,deltaind) = mean(corrdata(Ndata == N & deltadata == delta)); 
    end
end

%% Fitting models (by maximizing likelihood)   

%IL model fit
loglike_IL = bsxfun(@times, corrdataidx_1_N', log(perfmodel_IL))+ bsxfun(@times, corrdataidx_0_N', log(1 - perfmodel_IL));
totalloglike_IL = squeeze(sum(loglike_IL)); 
[LL_IL, idx] = max(totalloglike_IL(:));
[epsmax_idx_IL,Kmax_idx_IL] = ind2sub(size(totalloglike_IL),idx); 
pars_IL = [epsvec(epsmax_idx_IL) Kvec(Kmax_idx_IL)]; 
IL1 = squeeze(perfmodel_IL(:,epsmax_idx_IL,Kmax_idx_IL));
IL2 = repmat(IL1,1,9);
epsmax = 0.3; epssteps = 0.01; Kmax = 5; Ksteps = 1;
BMC_IL = log((epssteps/epsmax)*(Ksteps/Kmax))+LL_IL+log(sum(exp(totalloglike_IL(:)-LL_IL)));
BIC_IL = ((-2*LL_IL + 2*log(length(Ndata))))/-2; 
AIC_IL = ((-2*LL_IL + 2*2))/-2;
AICc_IL = AIC_IL +(12/(length(Ndata))-1);

RMSE_SS_IL = sqrt(mean((perfdata1(:)- IL1(:)).^2));
RMSE_CM_IL = sqrt(mean((perfdata2(:)- IL2(:)).^2));
Rsquared_SS_IL = 1-(sum((perfdata1(:)- IL1(:)).^2))/sum((perfdata1(:)- mean(perfdata1(:))).^2);
Rsquared_CM_IL = 1-(sum((perfdata2(:)- IL2(:)).^2))/sum((perfdata2(:)- mean(perfdata2(:))).^2);
Rs_IL = [RMSE_SS_IL RMSE_CM_IL Rsquared_SS_IL Rsquared_CM_IL];

%EP model fit
loglike_EP = bsxfun(@times, corrdataidx_1, log(perfmodel_EP))+ bsxfun(@times, corrdataidx_0, log(1 - perfmodel_EP));
totalloglike_EP = squeeze(sum(sum(loglike_EP,2),1)); 
[LL_EP, idx2] = max(totalloglike_EP(:));
[alphamax_idx_EP, J1max_idx_EP] = ind2sub(size(totalloglike_EP),idx2);  
pars_EP = [alphavec(alphamax_idx_EP) Jvec(J1max_idx_EP)];
perfmodel1_EP = squeeze(mean(perfmodel_EP,2));
EP1 = squeeze(perfmodel1_EP(:,alphamax_idx_EP,J1max_idx_EP));
EP2 = squeeze(perfmodel_EP(:,:,alphamax_idx_EP,J1max_idx_EP));
alphamax = 3; alphasteps = 0.1; Jmax = 25; Jsteps = 0.2;%humans
% alphamax = 3; alphasteps = 0.1; Jmax = 15; Jsteps = 0.2;%monkeys
BMC_EP = log((alphasteps/alphamax)*(Jsteps/Jmax))+LL_EP+log(sum(exp(totalloglike_EP(:)-LL_EP)));
BIC_EP = ((-2*LL_EP + 2*log(length(Ndata))))/-2; 
AIC_EP = ((-2*LL_EP + 2*2))/-2; 
AICc_EP = AIC_EP +(12/(length(Ndata))-1);

RMSE_SS_EP = sqrt(mean((perfdata1(:)- EP1(:)).^2));
RMSE_CM_EP = sqrt(mean((perfdata2(:)- EP2(:)).^2));
Rsquared_SS_EP = 1-(sum((perfdata1(:)- EP1(:)).^2))/sum((perfdata1(:)- mean(perfdata1(:))).^2);
Rsquared_CM_EP = 1-(sum((perfdata2(:)- EP2(:)).^2))/sum((perfdata2(:)- mean(perfdata2(:))).^2);
Rs_EP = [RMSE_SS_EP RMSE_CM_EP Rsquared_SS_EP Rsquared_CM_EP];

%SR model fit
loglike_SR = bsxfun(@times, corrdataidx_1, log(perfmodel_SR))+ bsxfun(@times, corrdataidx_0, log(1 - perfmodel_SR));
totalloglike_SR = squeeze(sum(sum(loglike_SR,2),1)); 
[LL_SR, idx2] = max(totalloglike_SR(:));
[Kmax_idx_SR,alphamax_idx_SR, J1max_idx_SR] = ind2sub(size(totalloglike_SR),idx2); 
pars_SR = [Kvec(Kmax_idx_SR) alphavec(alphamax_idx_SR) Jvec(J1max_idx_SR)]; 
perfmodel1_SR = squeeze(mean(perfmodel_SR,2));
SR1 = squeeze(perfmodel1_SR(:,Kmax_idx_SR,alphamax_idx_SR,J1max_idx_SR));
SR2 = squeeze(perfmodel_SR(:,:,Kmax_idx_SR,alphamax_idx_SR,J1max_idx_SR));
Kmax = 5; Ksteps = 1;alphamax = 3;alphasteps = 0.1;Jmax = 25;Jsteps = 0.2;%humans
% Kmax = 5; Ksteps = 1;alphamax = 3;alphasteps = 0.1;Jmax = 15;Jsteps = 0.2;%monkeys
BMC_SR = log((Ksteps/Kmax)*(alphasteps/alphamax)*(Jsteps/Jmax))+LL_SR+log(sum(exp(totalloglike_SR(:)-LL_SR)));
BIC_SR = ((-2*LL_SR + 3*log(length(Ndata))))/-2; 
AIC_SR = ((-2*LL_SR + 2*3))/-2; 
AICc_SR = AIC_SR +(24/(length(Ndata))-2);

RMSE_SS_SR = sqrt(mean((perfdata1(:)- SR1(:)).^2));
RMSE_CM_SR = sqrt(mean((perfdata2(:)- SR2(:)).^2));
Rsquared_SS_SR = 1-(sum((perfdata1(:)- SR1(:)).^2))/sum((perfdata1(:)- mean(perfdata1(:))).^2);
Rsquared_CM_SR = 1-(sum((perfdata2(:)- SR2(:)).^2))/sum((perfdata2(:)- mean(perfdata2(:))).^2);
Rs_SR = [RMSE_SS_SR RMSE_CM_SR Rsquared_SS_SR Rsquared_CM_SR];


%VP model fit 
loglike_VP = bsxfun(@times, corrdataidx_1, log(perfmodel_VP))+ bsxfun(@times, corrdataidx_0, log(1 - perfmodel_VP));
totalloglike_VP = squeeze(sum(sum(loglike_VP,2),1)); 
[LL_VP, idx2] = max(totalloglike_VP(:));
[alphamax_idx_VP, J1barmax_idx_VP,taumax_idx_VP] = ind2sub(size(totalloglike_VP),idx2); 
pars_VP = [alphavec(alphamax_idx_VP) J1barvec(J1barmax_idx_VP) tauvec(taumax_idx_VP)]; 
perfmodel1_VP = squeeze(mean(perfmodel_VP,2));
VP1 = squeeze(perfmodel1_VP(:,alphamax_idx_VP,J1barmax_idx_VP,taumax_idx_VP));
VP2 = squeeze(perfmodel_VP(:,:,alphamax_idx_VP,J1barmax_idx_VP,taumax_idx_VP));
% J1barmax = 30; J1barsteps = 0.5;alphamax = 3;alphasteps = 0.1;taumax = 30;tausteps = 0.5; %monkeys
J1barmax = 100; J1barsteps = 1.67;alphamax = 3;alphasteps = 0.1;taumax = 100;tausteps = 1.67; %humans
BMC_VP = log((J1barsteps/J1barmax)*(alphasteps/alphamax)*(tausteps/taumax))+LL_VP+log(sum(exp(totalloglike_VP(:)-LL_VP)));

BIC_VP = ((-2*LL_VP + 3*log(length(Ndata))))/-2; 
AIC_VP = ((-2*LL_VP + 2*3))/-2; 
AICc_VP = AIC_VP +(24/(length(Ndata))-2);

RMSE_SS_VP = sqrt(mean((perfdata1(:)- VP1(:)).^2));
RMSE_CM_VP = sqrt(mean((perfdata2(:)- VP2(:)).^2));
Rsquared_SS_VP = 1-(sum((perfdata1(:)- VP1(:)).^2))/sum((perfdata1(:)- mean(perfdata1(:))).^2);
Rsquared_CM_VP = 1-(sum((perfdata2(:)- VP2(:)).^2))/sum((perfdata2(:)- mean(perfdata2(:))).^2);
Rs_VP = [RMSE_SS_VP RMSE_CM_VP Rsquared_SS_VP Rsquared_CM_VP];

%VPF model fit 
loglike_VPF = bsxfun(@times, corrdataidx_1, log(perfmodel_VPF))+ bsxfun(@times, corrdataidx_0, log(1 - perfmodel_VPF));
totalloglike_VPF = squeeze(sum(sum(loglike_VPF,2),1)); 
[LL_VPF, idx2] = max(totalloglike_VPF(:));
[alphamax_idx_VPF, J1barmax_idx_VPF,taumax_idx_VPF,Kmax_idx_VPF] = ind2sub(size(totalloglike_VPF),idx2); 
pars_VPF = [alphavec(alphamax_idx_VPF) J1barvec(J1barmax_idx_VPF) tauvec(taumax_idx_VPF) Kvec(Kmax_idx_VPF)]; 
perfmodel1_VPF = squeeze(mean(perfmodel_VPF,2));
VPF1 = squeeze(perfmodel1_VPF(:,alphamax_idx_VPF,J1barmax_idx_VPF,taumax_idx_VPF, Kmax_idx_VPF));
VPF2 = squeeze(perfmodel_VPF(:,:,alphamax_idx_VPF,J1barmax_idx_VPF,taumax_idx_VPF,Kmax_idx_VPF));
% Kmax = 5; Ksteps = 1; J1barmax = 30; J1barsteps = 0.5;alphamax = 3;alphasteps = 0.1;taumax = 30;tausteps = 0.5; %monkeys
Kmax = 5; Ksteps = 1; J1barmax = 100; J1barsteps = 1.67;alphamax = 3;alphasteps = 0.1;taumax = 100;tausteps = 1.67; %humans
BMC_VPF = log((Ksteps/Kmax)*(J1barsteps/J1barmax)*(alphasteps/alphamax)*(tausteps/taumax))+LL_VPF+log(sum(exp(totalloglike_VPF(:)-LL_VPF)));

BIC_VPF = ((-2*LL_VPF + 4*log(length(Ndata))))/-2; 
AIC_VPF = ((-2*LL_VPF + 2*4))/-2; 
AICc_VPF = AIC_VPF +(72/(length(Ndata))-3);

RMSE_SS_VPF = sqrt(mean((perfdata1(:)- VPF1(:)).^2));
RMSE_CM_VPF = sqrt(mean((perfdata2(:)- VPF2(:)).^2));
Rsquared_SS_VPF = 1-(sum((perfdata1(:)- VPF1(:)).^2))/sum((perfdata1(:)- mean(perfdata1(:))).^2);
Rsquared_CM_VPF = 1-(sum((perfdata2(:)- VPF2(:)).^2))/sum((perfdata2(:)- mean(perfdata2(:))).^2);
Rs_VPF = [RMSE_SS_VPF RMSE_CM_VPF Rsquared_SS_VPF Rsquared_CM_VPF];

%Report all
LL = [LL_IL LL_EP LL_SR LL_VP LL_VPF];
BMC = [BMC_IL BMC_EP BMC_SR BMC_VP BMC_VPF];
BIC = [BIC_IL BIC_EP BIC_SR BIC_VP BIC_VPF];
AIC = [AIC_IL AIC_EP AIC_SR AIC_VP AIC_VPF];
AICc = [AICc_IL AICc_EP AICc_SR AICc_VP AICc_VPF];

