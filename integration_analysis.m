%Integration calculated, including part coeff, zmodule, modularity
    %loader_id, with subject ID information and run number
    %integration calculated details https://arxiv.org/abs/1511.02976
for ii = 1:60
    if ii < 23
        abc = sprintf('%s%d%s%d%s','load(''sub-0',loader(ii,1),'_run-',loader(ii,2),'_mtd_raw_4.mat'');');
        filename_partcoeff = sprintf('%s%d%s%d%s','sub-0',loader(ii,1),'_run-',loader(ii,2),'_partcoeff_mtd.mat');
        eval(abc) %loads raw mtd file
    else
        abc = sprintf('%s%d%s%d%s','load(''sub-',loader(ii,1),'_run-',loader(ii,2),'_mtd_raw_4.mat'');');
        filename_partcoeff = sprintf('%s%d%s%d%s','sub-',loader(ii,1),'_run-',loader(ii,2),'_partcoeff_mtd.mat');
        eval(abc) %loads raw mtd file
    end
    mtd2 = mtd(:,:,6:end -5);%removes first 5 and bottom 5 end timepoints
    [~,~,part,~,~] = integration(mtd2,1); % integration function from https://sites.google.com/site/bctnet/
    save(filename_partcoeff,'part');
    sprintf('%f%s',ii,'completed')
end
%% Fit to GLM 
for ii = 1:60
    if ii < 23
        dsmtx = sprintf('%s%d%s%d%s','load(''sub-0',loader(ii,1),'_run-',loader(ii,2),'_dsmtx_normalvsplankvsfreeze.mat'');');
        part = sprintf('%s%d%s%d%s','load(''sub-0',loader(ii,1),'_run-',loader(ii,2),'_partcoeff_mtd.mat'');');
        filename_part =sprintf('%s%d%s%d%s','sub-0',loader(ii,1),'_run-',loader(ii,2),'_pc_beta.mat');
        eval(dsmtx) %loads design matrix
        eval(part)%loads partcoeff file
    else
        dsmtx = sprintf('%s%d%s%d%s','load(''sub-',loader(ii,1),'_run-',loader(ii,2),'_dsmtx_normalvsplankvsfreeze.mat'');');
        part = sprintf('%s%d%s%d%s','load(''sub-',loader(ii,1),'_run-',loader(ii,2),'_partcoeff_mtd.mat'');');
        filename_paty =sprintf('%s%d%s%d%s','sub-',loader(ii,1),'_run-',loader(ii,2),'_pc_beta.mat');
        eval(dsmtx) %loads design matrix
        eval(part)%loads partcoeff file
    end
    nROI = 150; %number of ROIs
    for nn = 1:nROI
        pc_beta(:,nn) = glmfit(dsmtx_full(6:end -5,:),z(nn,:)');%run glmfit on partcoeff with timing of conditions
        sprintf('%d',nn)
    end
save(filename_part,'pc_beta');
sprintf('%f%s',ii,'completed')
end



%% Analysis of integration participation coefficients
% beta values are the values from GLM, defining relationship
for ii = 1:60
if ii < 23
abc = sprintf('%s%d%s%d%s','load(''sub-0',loader(ii,1),'_run-',loader(ii,2),'_pc_beta.mat'');');
eval(abc)
normal_beta_pc(ii,:) = pc_beta(2,:);
plank_beta_pc(ii,:) = pc_beta(3,:);
freeze_beta_pc(ii,:) = pc_beta(4,:);
cope_beta_pc_all(ii,:) = pc_beta(3,:) - pc_beta(2,:); %calculate difference in beta pc values for plank vs normal
else
abc = sprintf('%s%d%s%d%s','load(''sub-',loader(ii,1),'_run-',loader(ii,2),'_pc_beta.mat'');');
eval(abc)
normal_beta_pc(ii,:) = pc_beta(2,:);
plank_beta_pc(ii,:) = pc_beta(3,:);
freeze_beta_pc(ii,:) = pc_beta(4,:);
cope_beta_pc_all(ii,:) = pc_beta(3,:) - pc_beta(2,:); %calculate difference in beta pc values for plank vs normal
end
end

% run cope_beta_pc_all through permutation.m function - get pval for relationship
% between plank vs normal for 5000 iterations
pval_betapc_bin = double(pval<0.05);

%figure of difference between mean values of beta pc for plank vs normal,
%significance determined from permutation
figure
imagesc((cope_beta_pc_mean(:,sort_mat_id)).*pval_betapc_bin(:,sort_mat_id)) %sort_mat_id is id of networks sorted into motor, affective, cognitive and subcort regions
hold on
xline(54,'-','LineWidth',3);
xline(78,'--','LineWidth',3);
xline(140, '-.','LineWidth',3);
colorbar
title('significant mean beta pc plank vs normal')

%[h_pc_plankvsnormal,p_pc_plankvsnormal] = ttest(cope_beta_pc_all,zeros(60,150),'dim',1);

mean_plank_beta_pc = mean(plank_beta_pc);
mean_normal_beta_pc = mean(normal_beta_pc);
mean_cope_beta_pc = mean(cope_beta_pc_all);

%figure created using python code - thanks to https://github.com/RainCloudPlots/RainCloudPlots 
% mean(beta values for plank + intercept), mean(beta values for normal +
% intecept)




