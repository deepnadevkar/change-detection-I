function EPmodelpred_humans_forcluster(Jind)

% Simulate an experiment in which a monkey decides at which of two
% locations an orientation change occurred
% set size (N) and magnitude of change (Delta) were variable
% This code creates EP model predictions

Nvec = 2:5;
deltavec = 10:10:90; % change magnitude in degrees
nsteps = 100;

%Paramter ranges for human data
Jvec = linspace(0,25,nsteps);   %precision when set size = 1
alphavec = linspace(0,3,nsteps);  %power for relationship between set size and precision
Ntrials = 10000; % number of trials per set size

% interpolation of kappa over J
upper_bound = 3000;
k_vec = linspace(0,upper_bound,1e4);
J_vec = k_vec .* besseli(1,k_vec,1)./besseli(0,k_vec,1);

perfmodel = NaN(length(Nvec), length(deltavec), length(alphavec));
J1 = Jvec(Jind);

for alphaind = 1:length(alphavec)
    alpha = alphavec(alphaind)
    
    for Nind = 1:length(Nvec)
        N = Nvec(Nind);
        
        J = J1/N^alpha* ones(Ntrials,2);
        % J to k transformation
        kappa = interp1(J_vec,k_vec,J);
        x = circ_vmrnd(zeros(Ntrials,2),kappa); % noisy measurements from first display
        
        for deltaind = 1: length(deltavec)
            delta = deltavec(deltaind) * 2* pi/180;
                                             
            % Decision rule
            d = -log (besseli(0, kappa(:,2)))+ kappa(:,2) .* cos (x(:,2)) + log (besseli(0, kappa(:,1)))- kappa(:,1) .* cos(x(:,1)-delta); % decision variable
            
            perfmodel(Nind, deltaind, alphaind) = (sum(d>0) + 0.5 * sum(d==0))/Ntrials;
            
        end
    end
end


perfmodel(perfmodel==1) = 1-1/Ntrials;
perfmodel(perfmodel==0) = 1/Ntrials;
perfmodel_EP = perfmodel;

filename = strcat('EPmodelpred_H',num2str(Jind))
save(filename, 'perfmodel_EP','J1','Kvec','alphavec')