function [data1, data2, IL1, IL2, EP1, EP2, SR1, SR2, VP1, VP2, pars, parsest_IL, parsest_EP, parsest_SR, parsest_VP, LL, BIC_LL] = fakedatatest(generatingmodel, testingmodel)

%% GENERATING DATA

Nvec = 2:5;
deltavec = 10:10:90; % change magnitude in degrees
fakedatatrials = 10000; % number of trials in each fake data set

generatingmodel = '1';
switch generatingmodel
    case 1
        Kvec = 1:5; % capacity
        epsvec = 0:0.01:0.3;
        eps = epsvec(randi(length(epsvec),1));   % randomly generating eps from epsvec
        K = Kvec(randi(length(Kvec),1)); % randomly generating K from Kvec
        pars = [eps K];
        fakedata = ILfakedata(eps,K);
    case 2
        Jvec = 0:0.1:10;   %precision when set size = 1
        alphavec = 0:0.1:3;  %power for relationship between set size and precision
        J1 = Jvec(randi(length(Jvec),1)); % randomly generating kappa1 from kappa1vec
        alpha =  alphavec(randi(length(alphavec),1));   % randomly generating alpha from alphavec
        pars = [alpha J1];
        fakedata = EPfakedata_J(alpha,J1);
    case 3
        Jvec = 0:0.1:10;   %precision when set size = 1
        Kvec = 1:5; % capacity
        alphavec = 0:0.1:3;  %power for relationship between set size and precision
        K = Kvec(randi(length(Kvec),1)); % randomly generating K from Kvec
        alpha =  alphavec(randi(length(alphavec),1));   % randomly generating alpha from alphavec
        J1  = Jvec(randi(length(Jvec),1)); % randomly generating kappa1 from kappa1vec
        pars = [alpha K J1];
        fakedata = SRfakedata_J(K,alpha,J1);
    case 4
        nsteps = 100;
        Jvec = linspace(0,10,nsteps);
        alphavec = linspace(0,3,nsteps);
        tauvec = linspace(1,20,nsteps);
        Jbar1 = Jvec(randi(length(Jvec),1)); % randomly generating kappa_bar from kappabarvec
        tau = tauvec(randi(length(tauvec),1)); % randomly generating tau from tauvec
        alpha = alphavec(randi(length(alphavec),1));   % randomly generating alpha from alphavec
        pars = [alpha Jbar1 tau];
        fakedata = VPfakedata_J(alpha,Jbar1,tau);
        
end

%% calculating fakedata performance across set sizes and change magnitudes
Ndata = fakedata(:,1);
deltadata = fakedata(:,2);
corrdata = fakedata(:,3);

for Nind = 1:length(Nvec)
    N = Nvec(Nind);
    data1(Nind) = mean(corrdata(Ndata == N)); % data1 calculates ormance for different set sizes
    
    for deltaind = 1:length(deltavec)
        delta = deltavec(deltaind);
        data2(Nind,deltaind) = mean(corrdata(Ndata == N & deltadata == delta)); % data2 calculates ormance for different set sizes and change in magnitude
    end
end

dataN_mean = mean(data1,1);
dataN_ste = std(data2,0,1);
datadelta_mean = squeeze(mean(data2,1));
datadelta_ste = squeeze(std(data2,0,1));

%% TESTING MODELS ON THE GENERATED DATA

testingmodel = '1';
switch testingmodel
    case 1
        [IL1, IL2, epsest ,Kest, LL_IL] = ILmodelfitting(fakedata);
        parsest_IL = [epsest Kest];
        modelN_mean = mean(IL1,1); % calculates IL model fitted mean for all set sizes
        modelN_ste = std(IL1,0,1);  % calculates IL model fitted standard deviation for all set sizes
        modeldelta_mean = squeeze(mean(IL2,1)); % calculates IL model fitted mean for all change magnitudes
        modeldelta_ste = squeeze(std(IL2,0,1)); % calculates IL model fitted standard deviation for all change magnitudes
        BIC_LL_IL = (-2*LL_IL + 2*log(fakedatatrials))/-2; %calculate BIC-corrected LL
        
    case 2
        [EP1, EP2, alphaest, J1est, LL_EP] = EPmodelfitting_J(fakedata);
        parsest_EP = [alphaest J1est];
        modelN_mean = mean(EP1,1); % calculates EP model fitted mean for all set sizes
        modelN_ste = std(EP1,0,1);  % calculates EP model fitted standard deviation for all set sizes
        modeldelta_mean = squeeze(mean(EP2,1)); % calculates EP model fitted mean for all change magnitudes
        modeldelta_ste = squeeze(std(EP2,0,1)); % calculates EP model fitted standard deviation for all change magnitudes
        BIC_LL_EP = (-2*LL_EP + 2*log(fakedatatrials))/-2; %calculate BIC-corrected LL
        
    case 3
        [SR1, SR2, Kest, alphaest, J1est,LL_SR] = SRmodelfitting_J(fakedata);
        parsest_SR = [alphaest Kest J1est];
        modelN_mean = mean(SR1,1); % calculates SR model fitted mean for all set sizes
        modelN_ste = std(SR1,0,1);  % calculates SR model fitted standard deviation for all set sizes
        modeldelta_mean = squeeze(mean(SR2,1)); % calculates SR model fitted mean for all change magnitudes
        modeldelta_ste = squeeze(std(SR2,0,1)); % calculates SR model fitted standard deviation for all change magnitudes
        BIC_LL_SR = (-2*LL_SR + 3*log(fakedatatrials))/-2; %calculate BIC-corrected LL
        
    case 4
        [VP1, VP2,alphaest,Jbar1est,tauest, LL_VP] = VPmodelfitting_J(fakedata);
        parsest_VP = [alphaest Jbar1est tauest];
        modelN_mean = mean(VP1,1); % calculates VP model fitted mean for all set sizes
        modelN_ste = std(VP1,0,1);  % calculates VP model fitted standard deviation for all set sizes
        modeldelta_mean = squeeze(mean(VP2,1)); % calculates VP model fitted mean for all change magnitudes
        modeldelta_ste = squeeze(std(VP2,0,1)); % calculates VP model fitted standard deviation for all change magnitudes
        BIC_LL_VP = (-2*LL_VP + 3*log(fakedatatrials))/-2; %calculate BIC-corrected LL
        
end

LL = [LL_IL LL_EP LL_SR LL_VP];
BIC_LL = [BIC_LL_IL BIC_LL_EP BIC_LL_SR BIC_LL_VP];