function VPmodelpred_forcluster(J1barind)
% Simulate an experiment in which a monkey decides at which of two
% locations an orientation change occurred
% set size (N) and magnitude of change (Delta) were variable
% This code create model predictions for the VP model
tic
Nvec = 2:5;
deltavec = 10:10:90; % change magnitude in degrees
Ntrials = 10000; % number of trials per set size
nsteps = 100;

alphavec = linspace(0,3,nsteps);

%Parameter ranges for monkey data
J1barvec = linspace(0,30,nsteps);
tauvec = linspace(0.1,30,nsteps);

%Parameter ranges for human data
%J1barvec = linspace(0,100,nsteps);
%tauvec = linspace(0.1,100,nsteps);

% interpolation of kappa over J
upper_bound = 3000;
k_vec = linspace(0,upper_bound,1e4);
J_vec = k_vec .* besseli(1,k_vec,1)./besseli(0,k_vec,1);

%% Calculate model predictions and fitting model

perfmodel = NaN(length(Nvec), length(deltavec), length(alphavec), length(tauvec));

J1bar = J1barvec(J1barind);

for tauind = 1:length(tauvec)
    tau = tauvec(tauind)
    
    for alphaind = 1:length(alphavec)
        alpha = alphavec(alphaind);
        
        for Nind = 1:length(Nvec)
            N = Nvec(Nind);
            
            for deltaind = 1: length(deltavec)
                delta = deltavec(deltaind) * 2* pi/180;
                
                Jbar = J1bar/N^alpha; %Jbar at set size 1 that changes as set size changes
                Jvec = gamrnd(Jbar/tau, tau, Ntrials, 2);
                
                % J to k transformation
                kappavec = interp1(J_vec,k_vec,Jvec);
                
                % Generative model: simulating observer's measurements
                x = circ_vmrnd(zeros(Ntrials,2),kappavec); % noisy measurements from first display

                % Decision variable
                d = -log (besseli(0, kappavec(:,2)))+ kappavec(:,2) .* cos(x(:,2)) + log (besseli(0, kappavec(:,1)))- kappavec(:,1) .* cos(x(:,1)-delta);
                
                % Proportion correct
                perfmodel(Nind, deltaind, alphaind,tauind) = (sum(d>0) + 0.5 * sum(d==0))/Ntrials;
                
            end
        end
    end
end

perfmodel(perfmodel==0) = 1/Ntrials;
perfmodel(perfmodel==1) = 1 - 1/Ntrials;
perfmodel_VP = perfmodel;
toc

filename = strcat('VPmodelpred',num2str(J1barind))
save(filename, 'perfmodel_VP','J1bar','alphavec','tauvec')