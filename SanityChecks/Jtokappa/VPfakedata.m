function fakedataVP = fakedatasetVP(alpha, kappabar, tau)

% This code creates a n fake data sets for the VP model

Ntrials = 2500;
Nvec = 2:5;
deltavec = 10:10:90; % change magnitude in degrees


%% Calculate model predictions and fitting model

fakedataVP = [];
for Nind = 1:length(Nvec)
    N = Nvec(Nind);
    
    fakedataN = zeros(Ntrials,3);
    
    
    kappavec = gamrnd(kappabar/tau, tau, Ntrials, 2)/N^alpha;
        
    % Generative model: simulating observer's measurements
    x = circ_vmrnd(zeros(Ntrials,2),kappavec); % noisy measurements from first display
    phi1 = deltavec( randi(length(deltavec),[1 Ntrials]))'* 2* pi/180;    % orientation at first location (where change occurs) in second display
    phi2 = zeros(Ntrials,1);                                               % orientation at second location (remains the same as first display) in second display
    
    d = -log (besseli(0, kappavec(:,2)))+ kappavec(:,2) .* cos (x(:,2)-phi2) + log (besseli(0, kappavec(:,1)))- kappavec(:,1) .* cos(x(:,1)-phi1); % decision variable
    
    fakedataN(:,1) = N;
    fakedataN(:,2) = round(phi1/(2* pi/180));
    fakedataN(:,3) = d>0 + (d==0).* round(rand(Ntrials,1));
    
    fakedataVP = [fakedataVP ; fakedataN];
end

save fakedataVP.mat
