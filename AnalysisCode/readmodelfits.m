clear all close all;

for ii = 1:100
    fname = sprintf('Z:/dthakkar/CDTCode - 9.1.14_nolapse/modelfits_M3/modelfitsM3_%i.mat', ii);
    load(fname);
    AIC_all(ii,:) = AIC;
    AICc_all(ii,:) = AICc;
    BIC_all(ii,:) = BIC;
    BMC_all(ii,:) = BMC;
    IL1_all(ii,:) = IL1;
    IL2_all(ii,:,:) = IL2;
    EP1_all(ii,:) = EP1;
    EP2_all(ii,:,:) = EP2;
    EPF1_all(ii,:) = EPF1;
    EPF2_all(ii,:,:) = EPF2;
    VP1_all(ii,:) = VP1;
    VP2_all(ii,:,:) = VP2;
    VPF1_all(ii,:) = VPF1;
    VPF2_all(ii,:,:) = VPF2;
    LL_all(ii,:) = LL;
    Rs_IL_all(ii,:) = Rs_IL;
    Rs_EP_all(ii,:) = Rs_EP;
    Rs_EPF_all(ii,:) = Rs_EPF;
    Rs_VP_all(ii,:) = Rs_VP;
    Rs_VPF_all(ii,:) = Rs_VPF;
    pars_IL_all(ii,:) = pars_IL;
    pars_EP_all(ii,:) = pars_EP;
    pars_EPF_all(ii,:) = pars_EPF;
    pars_VP_all(ii,:) = pars_VP;
    pars_VPF_all(ii,:) = pars_VPF;
end

save modelfitting_m3 