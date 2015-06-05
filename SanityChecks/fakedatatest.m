function [data1, data2, IL1, IL2, EP1, EP2, SR1, SR2, VP1, VP2, pars, parsest_IL, parsest_EP, parsest_SR, parsest_VP, LL, BIC_LL, AIC_LL, AICc_LL, BIC, AIC, AICc, BMC] = fakedatatest(generatingmodel)

%% GENERATING DATA
%clear all close all

Nvec = 2:5;
deltavec = 10:10:90; % change magnitude in degrees
fakedatatrials = 10000; % number of trials in each fake data set

%generatingmodel = 2;
switch generatingmodel
    case 1
        Kvec = 1:5; % capacity
        epsvec = 0:0.01:0.3;
        eps = epsvec(randi(length(epsvec),1));   % randomly generating eps from epsvec
        K = Kvec(randi(length(Kvec),1)); % randomly generating K from Kvec
        pars = [eps K];
        fakedata = ILfakedata(eps,K);
    case 2
        Jvec = 0:1:25;   %precision when set size = 1
        alphavec = 0:0.2:3;  %power for relationship between set size and precision
        J1 = Jvec(randi(length(Jvec),1)); % randomly generating kappa1 from kappa1vec
        alpha =  alphavec(randi(length(alphavec),1));   % randomly generating alpha from alphavec
        pars = [alpha J1];
        fakedata = EPfakedata_J(alpha,J1);
    case 3
        Jvec = 0:2.5:25;   %precision when set size = 1
        Kvec = 1:5; % capacity
        alphavec = 0:0.5:3;  %power for relationship between set size and precision
        K = Kvec(randi(length(Kvec),1)); % randomly generating K from Kvec
        alpha =  alphavec(randi(length(alphavec),1));   % randomly generating alpha from alphavec
        J1  = Jvec(randi(length(Jvec),1)); % randomly generating kappa1 from kappa1vec
        pars = [alpha K J1];
        fakedata = SRfakedata_J(K,alpha,J1);
    case 4
        nsteps = 10;
        Jbarvec = linspace(0,60,nsteps);
        alphavec = linspace(0,3,nsteps);
        tauvec = linspace(1,80,nsteps);
        Jbar = Jbarvec(randi(length(Jbarvec),1)); % randomly generating kappa_bar from kappabarvec
        tau = tauvec(randi(length(tauvec),1)); % randomly generating tau from tauvec
        alpha = alphavec(randi(length(alphavec),1));   % randomly generating alpha from alphavec
        pars = [alpha Jbar tau];
        fakedata = VPfakedata_J(alpha,Jbar,tau);
        
end

%% calculating fakedata performance across set sizes and change magnitudes
Ndata = fakedata(:,1);
deltadata = fakedata(:,2);
corrdata = fakedata(:,3);

data1 = zeros(length(Nvec),1);
data2 = zeros(length(Nvec),length(deltavec));
data1_ste = zeros(length(Nvec),1);
data2_ste = zeros(length(Nvec),length(deltavec));

for Nind = 1:length(Nvec)
    N = Nvec(Nind);
    data1(Nind) = mean(corrdata(Ndata == N)); % data1 calculates ormance for different set sizes
    data1_ste (Nind) = std(corrdata(Ndata == N));
    for deltaind = 1:length(deltavec)
        delta = deltavec(deltaind);
        data2(Nind,deltaind) = mean(corrdata(Ndata == N & deltadata == delta)); % data2 calculates ormance for different set sizes and change in magnitude
        data2_ste (Nind, deltaind) = std(corrdata(Ndata == N & deltadata == delta));
    end
end

%% TESTING MODELS ON THE GENERATED DATA

%IL model
[IL1, IL2, epsest ,Kest, LL_IL, BMC_IL] = ILmodelfitting(fakedata);
parsest_IL = [epsest Kest];
BIC_LL_IL = (-2*LL_IL + 2*log(fakedatatrials))/-2; %calculate BIC-corrected LL
AIC_LL_IL = (-2*LL_IL + 2*2)/-2; %calculate AIC-corrected LL

BIC_IL = (-2*LL_IL + 2*log(fakedatatrials)); %calculate BIC-corrected LL
AIC_IL = (-2*LL_IL + 2*2); %calculate AIC-corrected LL
AICc_IL = AIC_IL +(12/(fakedatatrials)-1);
AICc_LL_IL = AICc_IL/-2;

% epsmax = 0.3; epssteps = 0.01; Kmax = 5; Ksteps = 1;
% BMC_IL = log((epssteps/epsmax)*(Ksteps/Kmax))+LL_IL+log(sum(exp(totalloglike(:)-LL_IL))); % calculate BMC

%EP model
[EP1, EP2, alphaest, J1est, LL_EP, BMC_EP] = EPmodelfitting_J(fakedata);
parsest_EP = [alphaest J1est];
BIC_LL_EP = (-2*LL_EP + 2*log(fakedatatrials))/-2; %calculate BIC-corrected LL
AIC_LL_EP = (-2*LL_EP + 2*2)/-2; %calculate AIC-corrected LL

BIC_EP = (-2*LL_EP + 2*log(fakedatatrials)); %calculate BIC-corrected LL
AIC_EP = (-2*LL_EP + 2*2); %calculate AIC-corrected LL
AICc_EP = AIC_EP +(12/(fakedatatrials)-1);
AICc_LL_EP = AICc_EP/-2;


% alphamax = 3; alphasteps = 0.1; Jmax = 10; Jsteps = 0.1;
% BMC_EP = log((alphasteps/alphamax)*(Jsteps/Jmax))+LL_EP+log(sum(exp(totalloglike(:)-LL_EP))); % calculate BMC

%SR model
[SR1, SR2, Kest, alphaest, J1est,LL_SR, BMC_SR] = SRmodelfitting_J(fakedata);
parsest_SR = [Kest alphaest J1est];
BIC_LL_SR = (-2*LL_SR + 3*log(fakedatatrials))/-2; %calculate BIC-corrected LL
AIC_LL_SR = (-2*LL_SR + 2*3)/-2; %calculate AIC-corrected LL

BIC_SR = (-2*LL_SR + 3*log(fakedatatrials)); %calculate BIC-corrected LL
AIC_SR = (-2*LL_SR + 2*3); %calculate AIC-corrected LL
AICc_SR = AIC_SR +(24/(fakedatatrials)-2);
AICc_LL_SR = AICc_SR/-2;


% Kmax = 5; Ksteps = 1; alphamax = 3; alphasteps = 0.1; Jmax = 10; Jsteps = 0.1;
% BMC_SR = log((Ksteps/Kmax)*(alphasteps/alphamax)*(Jsteps/Jmax))+LL_SR+log(sum(exp(totalloglike(:)-LL_SR))); %calculate BMC

%VP model
[VP1, VP2,alphaest,Jbarest,tauest, LL_VP, BMC_VP] = VPmodelfitting_J(fakedata);
parsest_VP = [alphaest Jbarest tauest];
BIC_LL_VP = (-2*LL_VP + 3*log(fakedatatrials))/-2; %calculate BIC-corrected LL
AIC_LL_VP = (-2*LL_VP + 2*3)/-2; %calculate AIC-corrected LL

BIC_VP = (-2*LL_VP + 3*log(fakedatatrials)); %calculate BIC-corrected LL
AIC_VP = (-2*LL_VP + 2*3); %calculate AIC-corrected LL
AICc_VP = AIC_VP +(24/(fakedatatrials)-2);
AICc_LL_VP = AICc_VP/-2;


%alphamax = 3; alphasteps = 0.75; Jmax = 10; Jsteps = 2.5; taumax = 20; tausteps = 5.75;
% BMC_VP = log((alphasteps/alphamax)*(Jsteps/Jmax)*(tausteps/taumax))+LL_VP+log(sum(exp(totalloglike(:)-LL_VP)));

%Report all
LL = [LL_IL LL_EP LL_SR LL_VP];
BIC_LL = [BIC_LL_IL BIC_LL_EP BIC_LL_SR BIC_LL_VP];
AIC_LL = [AIC_LL_IL AIC_LL_EP AIC_LL_SR AIC_LL_VP];
AICc_LL = [AICc_LL_IL AICc_LL_EP AICc_LL_SR AICc_LL_VP];

BIC = [BIC_IL BIC_EP BIC_SR BIC_VP];
AIC = [AIC_IL AIC_EP AIC_SR AIC_VP];
AICc = [AICc_IL AICc_EP AICc_SR AICc_VP];
BMC = [BMC_IL BMC_EP BMC_SR BMC_VP];