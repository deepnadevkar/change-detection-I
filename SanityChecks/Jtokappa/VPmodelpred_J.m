clear all; close all;
% Simulate an experiment in which a monkey decides at which of two
% locations an orientation change occurred
% set size (N) and magnitude of change (Delta) were variable
%This code creates
tic
Nvec = 2:5;
deltavec = 10:10:90; % change magnitude in degrees

Ntrials = 10000; % number of trials per set size
nsteps = 60;
Jvec = linspace(0,10,nsteps);
alphavec = linspace(0,3,nsteps);
tauvec = linspace(1,20,nsteps);

% interpolation of kappa over J
upper_bound = 3000;
k_vec = linspace(0,upper_bound,1e4);
J_vec = k_vec .* besseli(1,k_vec,1)./besseli(0,k_vec,1);

%% Calculate model predictions and fitting model

totalloglike = zeros(length(alphavec), length(Jvec));

perfmodel = zeros(length(Nvec), length(deltavec), length(alphavec), length(Jvec),length(tauvec));

for Jind = 1:length(Jvec)
    J_bar = Jvec(Jind);
    
    for tauind = 1:length(tauvec)
        tau = tauvec(tauind);
        
        for alphaind = 1:length(alphavec)
            alpha = alphavec(alphaind);
            
            for Nind = 1:length(Nvec)
                N = Nvec(Nind);
                
                Jvec = gamrnd(J_bar/tau, tau, Ntrials, 2)/N^alpha;
                % J to k transformation
                kappavec = interp1(J_vec,k_vec,J);
                
                % Generative model: simulating observer's measurements
                x = circ_vmrnd(zeros(Ntrials,2),kappavec); % noisy measurements from first display
                phi1 = deltavec( randi(length(deltavec),[1 Ntrials]))' * 2* pi/180;    % orientation at first location (where change occurs) in second display
                phi2 = zeros(Ntrials,1);                                               % orientation at second location (remains the same as first display) in second display
                
                for deltaind = 1: length(deltavec)
                    delta = deltavec(deltaind) * 2* pi/180;
                    trials = find(phi1== delta);
                    
                    d = -log (besseli(0, kappavec(:,2)))+ kappavec(:,2) .* cos (x(:,2)-phi2) + log (besseli(0, kappavec(:,1)))- kappavec(:,1) .* cos(x(:,1)-phi1); % decision variable
                    
                    perfmodel(Nind, deltaind, alphaind,Jind,tauind) = (sum(d(trials)>0) + 0.5 * sum(d(trials)==0))/length(trials);
                    
                end
            end
        end
    end
end
toc
save VPmodelpred_J perfmodel Jvec alphavec tauvec