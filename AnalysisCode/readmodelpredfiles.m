clear all close all; 
perfmodel = zeros(4,9,50,100,5,100);

for ii = 1:100
   fname = sprintf('C:/Users/dthakkar/Google Drive/Deepna cluster results/Exp 1/modelpreds_from_cluster_20140830/modelpredictions_VPF_range200/VPFmodelpred_bigrange_%i.mat', ii);
   % you will have to change the first part of the folder path and make
   % sure that they are all forward slashes and not back slashes
   load(fname);
   perfmodel(:,:,:,:,:,ii) = perfmodel_VPF;
end

clear perfmodel_VPF
perfmodel_VPF = perfmodel; 

save modelpred_VPF_range200 perfmodel_VPF 

