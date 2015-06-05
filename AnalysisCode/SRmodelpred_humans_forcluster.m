function SRmodelpred_humans_forcluster(Jind)

% Simulate an experiment in which a monkey decides at which of two
% locations an orientation change occurred
% set size (N) and magnitude of change (Delta) were variable
% This code creates SR model predictions
tic
Nvec = 2:5;
deltavec = 10:10:90; % change magnitude in degrees
nsteps = 100;

%Parameter ranges for human data
Kvec = 1:5; % capacity
Jvec = linspace(0,25,nsteps);   %precision when set size = 1
alphavec = linspace(0,3,nsteps);  %power for relationship between set size and precision
Ntrials = 10000; % number of trials per set size

% interpolation of kappa over J
upper_bound = 3000;
k_vec = linspace(0,upper_bound,1e4);
J_vec = k_vec .* besseli(1,k_vec,1)./besseli(0,k_vec,1);

perfmodel = NaN(length(Nvec), length(deltavec), length(Kvec), length(alphavec));

J1 = Jvec(Jind);

for alphaind = 1:length(alphavec)
    alpha = alphavec(alphaind)
    
    for Kind = 1:length(Kvec)
        K = Kvec(Kind);
        
        for deltaind = 1:length(deltavec)
            delta = deltavec(deltaind) * 2* pi/180;
            
            for Nind = 1:length(Nvec)
                N = Nvec(Nind);
                if N <= K
                    J = J1/N^alpha* ones(Ntrials,2);
                    % J to k transformation
                    kappa = interp1(J_vec,k_vec,J);
                    
                    x = circ_vmrnd(zeros(Ntrials,2),kappa); % noisy measurements from first display
                    % Decision rule
                    d = -log (besseli(0, kappa(:,2)))+ kappa(:,2) .* cos (x(:,2)) + log (besseli(0, kappa(:,1)))- kappa(:,1) .* cos(x(:,1)-delta); % decision variable
                    perfmodel(Nind, deltaind, Kind, alphaind) = (sum(d>0) + 0.5 * sum(d==0))/Ntrials;
                    
                elseif N > K
                    J = J1/K^alpha* ones(Ntrials,2);
                    % J to k transformation
                    kappa = interp1(J_vec,k_vec,J);
                    
                    x = circ_vmrnd(zeros(Ntrials,2),kappa); % noisy measurements from first display
                    
                    % case 1: phi1 and phi2 are both encoded
                    d = -log (besseli(0, kappa(:,2)))+ kappa(:,2) .* cos (x(:,2)) + log (besseli(0, kappa(:,1)))- kappa(:,1) .* cos(x(:,1)-delta); % decision variable
                    temp1 = (sum(d>0) + 0.5 * sum(d==0))/Ntrials;
                    
                    % case 2: phi1 is encoded and phi2 is not encoded
                    d = log (besseli(0, kappa(:,1))) - kappa(:,1) .* cos(x(:,1)-delta); % decision variable
                    temp2 = (sum(d>0) + 0.5 * sum(d==0))/Ntrials;
                    
                    % case 3: phi2 is encoded and phi1 is not encoded
                    d = -log (besseli(0, kappa(:,2))) + kappa (:,2) .* cos(x(:,2)); % decision variable
                    temp3 = (sum(d>0) + 0.5 * sum(d==0))/Ntrials;
                    
                    % case 4: phi1 and phi2 are both not encoded -->
                    % then proportion correct is simply 0.5
                    
                    perfmodel(Nind, deltaind, Kind, alphaind) = K/N *(K-1)/(N-1) * temp1  + K/N*(N-K)/(N-1) * temp2 + K/N*(N-K)/(N-1) * temp3 + (N-K)/N*(N-K-1)/(N-1)* 0.5;
                end
            end
        end
    end
end


perfmodel(perfmodel==0) = 1/Ntrials;
perfmodel(perfmodel==1) = 1 - 1/Ntrials;
perfmodel_SR = perfmodel;
toc

filename = strcat('SRmodelpred_H',num2str(Jind))
save(filename, 'perfmodel_SR','J1','Kvec','alphavec')
