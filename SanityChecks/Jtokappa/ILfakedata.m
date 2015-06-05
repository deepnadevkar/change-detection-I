function fakedataIL = ILfakedata(eps, K)

% Simulate an experiment in which a monkey decides at which of two
% locations an orientation change occurred
% set size (N) and magnitude of change (Delta) were variable
% This code creates EP model predictions

Nvec = 2:5;
deltavec = 10:10:90; % change magnitude in degrees

Ntrials = 2500; %Number of trials per set size

fakedataIL = [];
for Nind = 1:length(Nvec)
    N = Nvec(Nind);
    
    fakedataN = zeros(Ntrials,3);
    % Generative model: simulating observer's measurements
    phi1 = deltavec(randi(length(deltavec),[1 Ntrials]))'; % randomly generating change magnitude from deltavec
    
    fakedataN(:,1) = N;
    fakedataN(:,2) = phi1;
    fakedataN(:,3) = ones(Ntrials,1);
    
    if N <= K
        idx = rand(Ntrials,1)<eps;
        fakedataN(idx,3) = 0;
        
    elseif N > K
        caseprobs = [K*(K-1)+K*(N-K)+K*(N-K)  (N-K)*(N-K-1)] /N/(N-1);
        casecount = mnrnd(Ntrials, caseprobs); % multinomial draw from the four trial types; gives four trial counts, one for each type
        
        for casenum=1:2         % loop over the four trial types
            if casecount(casenum)>0 % check if there are any trials of this type
                
                if casenum == 1
                    % case 1: either phi1 or phi2 or both are encoded -->
                    % then proportion correct is 1-eps
                    idx = find(rand(casecount(1),1)<eps);
                    fakedataN(idx,3) = 0;
                    
                elseif casenum == 2
                    % case 2: phi1 and phi2 are both not encoded -->
                    % then proportion correct is simply 0.5
                    idx = find(rand(casecount(2),1) <0.5);
                    fakedataN(casecount(1)+idx,3) = 0;
                end
            end
        end
    end
    fakedataIL = [fakedataIL ; fakedataN];
end


% calculating and plotting fakedata performance across set sizes and change magnitudes
% 
% perfdata1 = zeros(1,length(Nvec));
% perfdata2  = zeros(length(Nvec),length(deltavec));
% 
% Ndata = fakedataIL(:,1);
% deltadata = fakedataIL(:,2);
% corrdata = fakedataIL(:,3);
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
% xlabel('Set size'); ylabel('Proportion correct');axis([1.8 5.2 0.5 1.0])
% set(gca,'YTick',0.4:.1:1.0)
% set(gca,'XTick', Nvec)
% 
% colorvec = get(gca, 'ColorOrder');
% colorvec = min(colorvec+.65,1);
% 
% figure;
% set(gca,'YTick',0.4:.1:1)
% set(gca,'XTick',10:10:90)
% plot(repmat(deltavec,length(Nvec),1)', datadelta_mean','-o'); hold on;
% xlabel('Change magnitude'); ylabel('Proportion correct');axis([10 90 0.4 1]);
% legend(strcat('N= ',int2str(Nvec')), 4);
% 
% 
