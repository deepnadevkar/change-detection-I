clear all; close all;

% Simulate an experiment in which a monkey decides at which of two
% locations an orientation change occurred
% set size (N) and magnitude of change (Delta) were variable
% This code creates EP model predictions

Nvec = 2:5;
deltavec = 10:10:90; % change magnitude in degrees
Kvec = 1:5;
epsvec = 0:0.01:0.3;

perfmodel = zeros(length(Nvec), length(epsvec), length(Kvec));

for epsind = 1:length(epsvec)
    eps = epsvec(epsind);
    
    for Kind = 1:length(Kvec)
        K = Kvec(Kind);
        for Nind = 1:length(Nvec)
            N = Nvec(Nind);
            if N <= K
                perfmodel(Nind,epsind,Kind) = (1 - eps);
            elseif N > K
                perfmodel(Nind,epsind,Kind) = (1 - eps)-(N-K)*(N-K-1)/(N*(N-1))*(0.5-eps);
            end
        end
    end
end

ones = find(perfmodel== 1);
perfmodel(ones)= (1/100000);
zeros = find(perfmodel == 0);
perfmodel(zeros) = (0/100000);

save ILmodelpred perfmodel Kvec epsvec
