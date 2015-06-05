function VPFmodelpredH_forcluster(J1barind)
% This code creates VPF (VP + IL) model predictions
tic
Nvec = 2:5; %set size
deltavec = 10:10:90; % change magnitude in degrees
nsteps = 100;

Kvec = 1:5;
alphavec = linspace(0,3,50);

%Parameter ranges for monkey data
% J1barvec = linspace(0,30,nsteps);
% tauvec = linspace(0.1,30,nsteps);

%Parameter ranges for human data
J1barvec = linspace(0,200,nsteps);
tauvec = linspace(0.1,200,nsteps);

Ntrials = 10000; % number of trials per set size

% interpolation of kappa over J
upper_bound = 3000;
k_vec = linspace(0,upper_bound,1e4);
J_vec = k_vec .* besseli(1,k_vec,1)./besseli(0,k_vec,1);

perfmodel = NaN(length(Nvec), length(deltavec), length(alphavec), length(tauvec), length(Kvec));

J1bar = J1barvec(J1barind);

for tauind = 1:length(tauvec)
    tau = tauvec(tauind);
    tauind
    
    for alphaind = 1:length(alphavec)
        alpha = alphavec(alphaind);
        
        for Kind = 1:length(Kvec)
            K = Kvec(Kind);
            
            for deltaind = 1:length(deltavec)
                delta = deltavec(deltaind) * 2* pi/180;
                
                for Nind = 1:length(Nvec)
                    N = Nvec(Nind);
                    
                    if N <= K
                        Jbar = J1bar/N^alpha;
                        Jvec = gamrnd(Jbar/tau, tau, Ntrials, 2);
                        % J to k transformation
                        kappavec = interp1(J_vec,k_vec,Jvec);
                        
                        % Generative model: simulating observer's measurements
                        x = circ_vmrnd(zeros(Ntrials,2),kappavec); % noisy measurements from first display
                        phi1 = delta;               % orientation at first location (where change occurs) in second display
                        phi2 = zeros(Ntrials,1);    % orientation at second location (remains the same as first display) in second display
                        
                        d = -log (besseli(0, kappavec(:,2)))+ kappavec(:,2) .* cos (x(:,2)-phi2) + log (besseli(0, kappavec(:,1)))- kappavec(:,1) .* cos(x(:,1)-phi1); % decision variable
                        
                        perfmodel(Nind, deltaind, alphaind,tauind,Kind) = (sum(d>0) + 0.5 * sum(d==0))/Ntrials;
                        
                    elseif N > K
                        Jbar = J1bar/K^alpha;
                        Jvec = gamrnd(Jbar/tau, tau, Ntrials, 2);
                        % J to k transformation
                        kappavec = interp1(J_vec,k_vec,Jvec);
                        
                        x = circ_vmrnd(zeros(Ntrials,2),kappavec); % noisy measurements from first display
                        phi1 = delta;               % orientation at first location (where change occurs) in second display
                        phi2 = zeros(Ntrials,1);    % orientation at second location (remains the same as first display) in second display
                        
                        % case 1: phi1 and phi2 are both encoded
                        d = -log (besseli(0, kappavec(:,2)))+ kappavec(:,2) .* cos (x(:,2)-phi2) + log (besseli(0, kappavec(:,1)))- kappavec(:,1) .* cos(x(:,1)-phi1); % decision variable
                        temp1 = (sum(d>0) + 0.5 * sum(d==0))/Ntrials;
                        
                        % case 2: phi1 is encoded and phi2 is not encoded
                        d = log (besseli(0, kappavec(:,1))) - kappavec(:,1) .* cos(x(:,1)-phi1); % decision variable
                        temp2 = (sum(d>0) + 0.5 * sum(d==0))/Ntrials;
                        
                        % case 3: phi2 is encoded and phi1 is not encoded
                        d = -log (besseli(0, kappavec(:,2))) + kappavec (:,2) .* cos(x(:,2)-phi2); % decision variable
                        temp3 = (sum(d>0) + 0.5 * sum(d==0))/Ntrials;
                        
                        % case 4: phi1 and phi2 are both not encoded -->
                        % then proportion correct is simply 0.5
                        
                        perfmodel(Nind, deltaind, alphaind,tauind,Kind) = K/N *(K-1)/(N-1) * temp1  + K/N*(N-K)/(N-1) * temp2 + K/N*(N-K)/(N-1) * temp3 + (N-K)/N*(N-K-1)/(N-1)* 0.5;
                    end
                end
            end
        end
    end
end

perfmodel(perfmodel==1) = 1-1/Ntrials;
perfmodel(perfmodel==0) = 1/Ntrials;
perfmodel_VPF = perfmodel;
filename = strcat('VPFmodelpred_bigrange_',num2str(J1barind))
save(filename, 'perfmodel_VPF' ,'J1bar','alphavec', 'tauvec','kappavec','Kvec')
toc