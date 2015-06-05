clear all close all;

load fakedataEPfits_1
load fakedataEPfits_2
load fakedataEPfits_3
load fakedataEPfits_4
load fakedataEPfits_5

diff_AIC_LL = [diff_AIC_LL_1; diff_AIC_LL_2; diff_AIC_LL_3; diff_AIC_LL_4; diff_AIC_LL_5];
diff_AICc_LL = [diff_AICc_LL_1; diff_AICc_LL_2; diff_AICc_LL_3; diff_AICc_LL_4; diff_AICc_LL_5];
diff_BIC_LL = [diff_BIC_LL_1; diff_BIC_LL_2; diff_BIC_LL_3; diff_BIC_LL_4; diff_BIC_LL_5];
diff_BMC = [diff_BMC_1; diff_BMC_2; diff_BMC_3; diff_BMC_4; diff_BMC_5];

mean_diff_AIC_LL = mean(diff_AIC_LL);
ste_diff_AIC_LL = (std(diff_AIC_LL)/sqrt(5));

mean_diff_AICc_LL = mean(diff_AICc_LL);
ste_diff_AICc_LL = (std(diff_AICc_LL)/sqrt(5));

mean_diff_BIC_LL = mean(diff_BIC_LL);
ste_diff_BIC_LL = (std(diff_BIC_LL)/sqrt(5));

mean_diff_BMC = mean(diff_BMC);
ste_diff_BMC = (std(diff_BMC)/sqrt(5));
