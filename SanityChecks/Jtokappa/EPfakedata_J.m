function fakedataEP = fakedatasetEP(alpha, J1)

% This code generates a fake EP data set

% We create fake data with these parameters (and then hope that the model
% fitting code can reproduce those values)

Ntrials = 2500; %number of trials per set size
Nvec = 2:5;
deltavec = 10:10:90; % change magnitude in degrees

% interpolation of kappa over J
upper_bound = 3000;
k_vec = linspace(0,upper_bound,1e4);
J_vec = k_vec .* besseli(1,k_vec,1)./besseli(0,k_vec,1);

fakedataEP = [];

for Nind = 1:length(Nvec)
    N = Nvec(Nind);
    
    fakedataN = zeros(Ntrials,3);
    
    J = J1/N^alpha* ones(Ntrials,2);
    % J to k transformation
    kappa = interp1(J_vec,k_vec,J);
          
    % Generative model: simulating observer's measurements
    x = circ_vmrnd(zeros(Ntrials,2),kappa); % noisy measurements from first display
    phi1 = deltavec(randi(length(deltavec),[1 Ntrials]))'* 2* pi/180;    % orientation at first location (where change occurs) in second display
    
    % Decision rule
    d = -log (besseli(0, kappa(:,2)))+ kappa(:,2) .* cos (x(:,2)) + log (besseli(0, kappa(:,1)))- kappa(:,1) .* cos(x(:,1)-phi1); % decision variable
    
    fakedataN(:,1) = N;
    fakedataN(:,2) = round(phi1/(2* pi/180));
    fakedataN(:,3) = d>0 + (d==0).* round(rand(Ntrials,1));
    
    fakedataEP = [fakedataEP ; fakedataN];
    
end

save fakedataEP_J.mat