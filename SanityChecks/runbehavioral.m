clear all close all;

load fakedataSR

Nvec = 2:5;
deltavec = 10:10:90; % change magnitude in degrees

Ndata = fakedata(:,1);
deltadata = fakedata(:,2);
corrdata = fakedata(:,3); 

perfdata1 = zeros(length(Nvec),1);
perfdata2 = zeros(length(Nvec),length(deltavec));

%% Calculating monkey performance
for Nind = 1:length(Nvec)
    N = Nvec(Nind);
    perfdata1(Nind) = mean(corrdata(Ndata == N));
    for deltaind = 1:length(deltavec)
            delta = deltavec(deltaind);
            perfdata2(Nind,deltaind) = mean(corrdata(find(Ndata == N & deltadata == delta)));
    end
end

%set size plot
figure;
bar(Nvec, perfdata1); hold on;
xlabel('Set size'); ylabel('Proportion correct');axis([1.5 5.5 0.5 0.8])
set(gca,'YTick',0.5:.1:0.8)
set(gca,'XTick', 2:1:5)

%change mgnitude plot
figure;
set(gca,'YTick',0.2:.1:1)
set(gca,'XTick', 10:10:90)
plot(repmat(deltavec,4,1)', perfdata2','-o'); hold on; 
xlabel('Change magnitude'); ylabel('Proportion correct');axis([7 93 0.2 1]);
legend(strcat('N= ',int2str(Nvec')), 4);