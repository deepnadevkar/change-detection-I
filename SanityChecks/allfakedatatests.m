clear all; close all;
Nvec = 2:1:5;
deltavec = 10:10:90;

generatingmodel = 4;
tic
[data1, data2, IL1, IL2, EP1, EP2, SR1, SR2, VP1, VP2, pars, parsest_IL, parsest_EP, parsest_SR, parsest_VP, LL, BIC_LL, AIC_LL, AICc_LL, BIC, AIC, AICc, BMC] = fakedatatest(generatingmodel);
% Parameter recovery & difference in criterions from winning model 
switch generatingmodel
    case 1
        parsest = parsest_IL;
        
    case 2
        parsest = parsest_EP;
        
    case 3
        parsest = parsest_SR;
              
    case 4
        parsest = parsest_VP;
              
        
end
[pars; parsest]
% BIC_LL
BIC_LL
AIC_LL
BIC
AIC
AICc
BMC

% IL Model set size and change magnitude fits
figure; plot(Nvec,data1,'-o')
hold on; plot(Nvec,IL1, '--');
xlabel('Set size'); ylabel('Proportion correct');axis([1.8 5.2 0.3 1.1])
set(gca,'YTick',0.3:.1:1.0)
set(gca,'XTick', 2:1:5)
title_IL_ = title('IL model fit');
set(title_IL_,'String','IL model fit')

figure;
plot(deltavec,data2'); hold on;
plot(deltavec, IL2','--');
xlabel('Change magnitude'); ylabel('Proportion correct');
ylim([0.3 1.3])
legend(strcat('N =',num2str(Nvec')),'Location','Best')
title_IL_ = title('IL model fit');
set(title_IL_,'String','ILmodel fit')

% EP Model set size and change magnitude fits
figure; plot(Nvec,data1,'-o')
hold on; plot(Nvec,EP1, '--');
xlabel('Set size'); ylabel('Proportion correct');axis([1.8 5.2 0.3 1.1])
set(gca,'YTick',0.3:.1:1.0)
set(gca,'XTick', 2:1:5)
title_IL_ = title('EPmodel fit');
set(title_IL_,'String','EP model fit')

figure;
plot(deltavec,data2'); hold on;
plot(deltavec, EP2','--');
xlabel('Change magnitude'); ylabel('Proportion correct');
ylim([0.3 1.3])
legend(strcat('N =',num2str(Nvec')),'Location','Best')
title_IL_ = title('EP model fit');
set(title_IL_,'String','EP model fit')

% SR Model set size and change magnitude fits
figure; plot(Nvec,data1,'-o')
hold on; plot(Nvec,SR1, '--');
xlabel('Set size'); ylabel('Proportion correct');axis([1.8 5.2 0.3 1.1])
set(gca,'YTick',0.3:.1:1.0)
set(gca,'XTick', 2:1:5)
title_IL_ = title('SR model fit');
set(title_IL_,'String','SR model fit')

figure;
plot(deltavec,data2'); hold on;
plot(deltavec, SR2','--');
xlabel('Change magnitude'); ylabel('Proportion correct');
ylim([0.3 1.3])
legend(strcat('N =',num2str(Nvec')),'Location','Best')
title_IL_ = title('SR model fit');
set(title_IL_,'String','SR model fit')

% VP Model set size and change magnitude fits
figure; plot(Nvec,data1,'-o')
hold on; plot(Nvec,VP1, '--');
xlabel('Set size'); ylabel('Proportion correct');axis([1.8 5.2 0.3 1.1])
set(gca,'YTick',0.3:.1:1.0)
set(gca,'XTick', 2:1:5)
title_IL_ = title('VP model fit');
set(title_IL_,'String','VP model fit')

figure;
plot(deltavec,data2'); hold on;
plot(deltavec, VP2','--');
xlabel('Change magnitude'); ylabel('Proportion correct');
ylim([0.3 1.3])
legend(strcat('N =',num2str(Nvec')),'Location','Best')
title_IL_ = title('VP model fit');
set(title_IL_,'String','VP model fit')

toc






