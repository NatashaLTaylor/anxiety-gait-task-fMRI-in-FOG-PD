%% DYNAMIC CONNECTIVITY ANALYSIS
%this analysis involves the addition of time-series from PAG
function d = dFC_analysis(filename_ts,subnum,run,window)
load(filename_ts);
% files create from ts_extraction function
sub_ts_2_name = sprintf('%s_%s-ts_parc2.mat',subnum, run); %additional region of interest
sub_ts_2 = load(sub_ts_2_name);
% this part of code, creates the specific ROIs time series
load('special_ROIs_parc3.mat');
ts_roi_1 = sub_ts(special_ROIs,:);
ts_roi_2 = sub_ts_2.sub_ts;
ts_roi_all = vertcat(ts_roi_1,ts_roi_2); %concatenate all ROIs, PAG was added as a region
save([subnum '_' run '_ts_roi_all_5.mat'],'ts_roi_all');
%% Design Matrix
%load design matrix of normal, plank,freeze conditions for subject
dsmtx_name = sprintf('%s_%s_dsmtx_normalvsplankvsfreeze.mat',subnum,run);
load(dsmtx_name);
%figure
%plot(dsmtx_full)
%% BOLD Signal
nROI = 150; %number of ROIs
for nn = 1:nROI
bold_beta(:,nn) = glmfit(dsmtx_full(6:end -5,:),ts_roi_all(nn,6:end -5)');%run glmfit on bold signal with epochs of all conditions
sprintf('%d',nn)
end
save([subnum '_' run '_beta_bold_5.mat'],'bold_beta');
%cope = contrast of parameter estimates, cope for difference between plank
%and normal (plank - normal)
sub_cope = bold_beta(3,:)-bold_beta(2,:); %row 2 is normal and row 3 is plank, row 4 is freeze 
save([subnum '_' run '_cope_bold_5.mat'],'sub_cope'); %save based on subject and run
%% MTD
% calculate multiplication of temporal derivatives, coupling function https://github.com/macshine/coupling/
mtd = coupling(ts_roi_all',window); %remember transpose matrix to time x nodes, remove first and last 5 timepoints
template = find(triu(ones(nROI))-eye(nROI)); %finding unique combination of pairs
save([subnum '_' run '_mtd_raw_4.mat'],'mtd');

mtd = mtd(:,:,6:end -5);%removes first 5 and bottom 5 end timepoints
for tt = 1:140
temp = mtd(:,:,tt);
mtd_flat(:,tt) = temp(template); %flattens mtd across timepoints, 2D matrix
end
save([subnum '_' run '_mtd_flat_4.mat'],'mtd_flat');
% calculates general linear model for relationship between MTD and design
% matrix conditions
nPairs = 11175; %number of unique pairs from template
nTime = 140;
for qq = 1:nPairs
beta_mtd(:,qq) = glmfit(dsmtx_full(6:end -5,:),mtd_flat(qq,:)'); %removes first 5 and bottom 5 end timepoints from design matrix
sprintf('%d',qq)
end
save([subnum '_' run '_beta_mtd_4.mat'],'beta_mtd');
% cope for mtd between plank and walk
sub_cope_mtd = beta_mtd(3,:)-beta_mtd(2,:);
save([subnum '_' run '_cope_mtd_4.mat'],'sub_cope_mtd');
end
