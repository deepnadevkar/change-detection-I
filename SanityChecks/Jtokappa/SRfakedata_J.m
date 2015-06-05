function fakedataSR = fakedatasetSR_J(K, alpha, J1)

% This code creates a n fake data sets for the SR model


Ntrials = 2500; % number of trials per set size
Nvec = 2:5;
deltavec = 10:10:90; % change magnitude in degrees
Kvec = 1:5; % capacity
Jvec = 0:0.1:10;   %precision when set size = 1
alphavec = 0:0.1:3;  %power for relationship between set size and precision

% interpolation of kappa over J
upper_bound = 3000;
k_vec = linspace(0,upper_bound,1e4);
J_vec = k_vec .* besseli(1,k_vec,1)./besseli(0,k_vec,1);

fakedataSR = [];
for Nind = 1:length(Nvec)
    N = Nvec(Nind);
    
    fakedataN = zeros(Ntrials,3);
    % Generative model: simulating observer's measurements
    phi1 = deltavec(randi(length(deltavec),[1 Ntrials]))'* 2*pi/180; % randomly generating change magnitude from deltavec
    
    fakedataN(:,1) = N;
    
    if N <= K
        J = J1/N^alpha* ones(Ntrials,2);
        % J to k transformation
        kappa = interp1(J_vec,k_vec,J);
        
        x = circ_vmrnd(zeros(Ntrials,2),kappa); % noisy measurements from first display
        % Decision rule
        d = -log (besseli(0, kappa(:,2)))+ kappa(:,2) .* cos (x(:,2)) + log (besseli(0, kappa(:,1)))- kappa(:,1) .* cos(x(:,1)-phi1); % decision variable
        fakedataN(:,2) = round(phi1/(2*pi/180));
        fakedataN(:,3) = d>0;
        
    elseif N > K
        
        caseprobs = [K*(K-1)  K*(N-K)   K*(N-K)  (N-K)*(N-K-1)] /N/(N-1);
        casecount = mnrnd(Ntrials, caseprobs); % multinomial draw from the four trial types; gives four trial counts, one for each type
        cumcasecount = cumsum(casecount);   % cumulative sum of these counts
        
        for casenum=1:4         % loop over the four trial types
            if casecount(casenum)>0 % check if there are any trials of this type
                J = J1/K^alpha* ones(casecount(casenum),2);
                % J to k transformation
                kappa = interp1(J_vec,k_vec,J);
                x = circ_vmrnd(zeros(casecount(casenum),2),kappa); % noisy measurements from first display
                if casenum == 1
                    % case 1: phi1 and phi2 are both encoded
                    d = -log (besseli(0, kappa(:,2)))+ kappa(:,2) .* cos (x(:,2)) + log (besseli(0, kappa(:,1)))- kappa(:,1) .* cos(x(:,1)-phi1(1:cumcasecount(1))); % decision variable
                    fakedataN(1:cumcasecount(1),3) = d>0;
                elseif casenum == 2
                    % case 2: phi1 is encoded and phi2 is not encoded
                    d = log (2 * pi * besseli(0, kappa(:,1))) - kappa(:,1) .* cos(x(:,1)-phi1(cumcasecount(1)+1:cumcasecount(2))); % decision variable
                    fakedataN(cumcasecount(1)+1:cumcasecount(2),3) = d>0 + (d==0).* round(rand(casecount(casenum),1));
                elseif casenum == 3
                    % case 3: phi2 is encoded and phi1 is not encoded
                    d = -log (2 * pi * besseli(0, kappa(:,2))) + kappa (:,2) .* cos(x(:,2)); % decision variable
                    fakedataN(cumcasecount(2)+1:cumcasecount(3),3) = d>0 + (d==0).* round(rand(casecount(casenum),1));
                elseif casenum == 4
                    % case 4: phi1 and phi2 are both not encoded -->
                    % then proportion correct is simply 0.5
                    fakedataN(cumcasecount(3)+1:Ntrials,3) = round(rand(casecount(casenum),1));
                end
            end
        end
        fakedataN(:,2) = round(phi1/(2*pi/180));
    end
    
    fakedataSR = [fakedataSR ; fakedataN];
end

% %% calculating and plotting fakedata performance across set sizes and change magnitudes
%
% perfdata1 = zeros(1,length(Nvec));
% perfdata2  = zeros(length(Nvec),length(deltavec));
%
% Ndata =fakedataSR(:,1);
% deltadata = fakedataSR(:,2);
% corrdata = fakedataSR(:,3);
%
% for Nind = 1:length(Nvec)
%     N = Nvec(Nind);
%     perfdata1(Nind) = mean(corrdata(Ndata == N)); % perfdata1 calculates performance for different set sizes
%
%     for deltaind = 1:length(deltavec)
%         delta = deltavec(deltaind);
%         perfdata2(Nind,deltaind) = mean(corrdata(Ndata == N & deltadata == delta)); % perfdata2 calculates performance for different set sizes and change in magnitude
%     end
% end
%
%
% dataN_mean = perfdata1;
% datadelta_mean = perfdata2;
%
% figure;
% plot(Nvec, dataN_mean);
% xlabel('Set size'); ylabel('Proportion correct');axis([1.8 5.2 0.5 0.8])
% set(gca,'YTick',0.4:.1:0.8)
% set(gca,'XTick', 2:1:5)
%
% colorvec = get(gca, 'ColorOrder');
% colorvec = min(colorvec+.65,1);
%
% figure;
% set(gca,'YTick',0.4:.1:1)
% set(gca,'XTick',10:10:90)
% plot(repmat(deltavec,4,1)', datadelta_mean','-o'); hold on;
% xlabel('Change magnitude'); ylabel('Proportion correct');axis([10 90 0.4 1]);
% legend(strcat('N= ',int2str(Nvec')), 4);
%
save fakedataSR.mat
