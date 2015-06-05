clear all; close all;

% Simulate an experiment in which a monkey decides at which of two
% locations an orientation change occurred
% set size (N) and magnitude of change (Delta) were variable
% This code creates EP model predictions

Nvec = 2:5;
deltavec = 10:10:90; % change magnitude in degrees
Jvec = 0:0.2:25;   %precision when set size = 1
alphavec = 0:0.1:3;  %power for relationship between set size and precision
Ntrials = 10000; % number of trials per set size

% interpolation of kappa over J
upper_bound = 3000;
k_vec = linspace(0,upper_bound,1e4);
J_vec = k_vec .* besseli(1,k_vec,1)./besseli(0,k_vec,1);

perfmodel = zeros(length(Nvec), length(deltavec), length(alphavec), length(Jvec));

for alphaind = 1:length(alphavec)
    alpha = alphavec(alphaind);
    
    for Jind = 1:length(Jvec)
        J1 = Jvec(Jind);
        
        for Nind = 1:length(Nvec)
            N = Nvec(Nind);
                      
            J = J1/N^alpha* ones(Ntrials,2);
            % J to k transformation
            kappa = interp1(J_vec,k_vec,J);
            x = circ_vmrnd(zeros(Ntrials,2),kappa); % noisy measurements from first display                    
                        
            for deltaind = 1: length(deltavec)
                delta = deltavec(deltaind);
                % Generative model: simulating observer's measurements
                
                phi1 = delta * ones(Ntrials,1)* 2* pi/180; % orientation at first location (where change occurs) in second display
                
                % Decision rule
                d = -log (besseli(0, kappa(:,2)))+ kappa(:,2) .* cos (x(:,2)) + log (besseli(0, kappa(:,1)))- kappa(:,1) .* cos(x(:,1)-phi1); % decision variable
                
                perfmodel(Nind, deltaind, alphaind, Jind) = (sum(d>0) + 0.5 * sum(d==0))/Ntrials;
               
            end
        end
    end
end

ones = find(perfmodel== 1);
perfmodel(ones)= (1/100000);
zeros = find(perfmodel == 0);
perfmodel(zeros) = (0/100000);

save EPmodelpred_J perfmodel Jvec alphavec
