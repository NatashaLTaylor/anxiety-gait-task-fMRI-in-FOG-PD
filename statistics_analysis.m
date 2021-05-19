%% Stats + Figs section
%% load data
loader = [1;2;2;2;3;3;3;4;4;4;5;5;5;6;6;7;8;8;9;9;9;10;10;10;11;11;11;13;13;13;15;15;16;16;16;17;18;19;19;20;20;21;21;22;22;23;23;24;24;25;25;26;26;28;28;29;29;30;30];
loader(:,2) = [1;1;2;3;1;2;3;1;2;3;1;2;3;1;2;2;1;2;3;1;2;3;1;2;3;1;2;3;1;2;3;1;1;2;3;1;1;1;2;1;2;1;2;1;2;1;2;1;2;1;2;1;2;1;2;1;2;1;2];
find(loader(:,1)==10,1)
%load data
for ii = 1:60
if ii < 23
abc = sprintf('%s%d%s%d%s','load(''sub-0',loader(ii,1),'_run-',loader(ii,2),'_cope_mtd_5.mat'');');
eval(abc)
cope_mtd_plankvsnormal_all(ii,:) = sub_cope_mtd;
else
abc = sprintf('%s%d%s%d%s','load(''sub-',loader(ii,1),'_run-',loader(ii,2),'_cope_mtd_5.mat'');');
eval(abc)
cope_mtd_plankvsnormal_all(ii,:) = sub_cope_mtd;
end
end
%load individual cope (difference between plank vs normal) for bold
for ii = 1:60
if ii < 23
abc = sprintf('%s%d%s%d%s','load(''sub-0',loader(ii,1),'_run-',loader(ii,2),'_cope_bold_5.mat'');');
eval(abc)
cope_bold_plankvsnormal_all(ii,:) = sub_cope;
else
abc = sprintf('%s%d%s%d%s','load(''sub-',loader(ii,1),'_run-',loader(ii,2),'_cope_bold_5.mat'');');
eval(abc)
cope_bold_plankvsnormal_all(ii,:) = sub_cope;
end
end
% load bold and mtd beta values
for ii = 1:60
if ii < 23
abc = sprintf('%s%d%s%d%s','load(''sub-0',loader(ii,1),'_run-',loader(ii,2),'_beta_bold_5.mat'');');
eval(abc)
normal_bold_all(ii,:) = bold_beta(2,:);
plank_bold_all(ii,:) = bold_beta(3,:);
freeze_bold_all(ii,:) = bold_beta(4,:);
abc = sprintf('%s%d%s%d%s','load(''sub-0',loader(ii,1),'_run-',loader(ii,2),'_beta_mtd_5.mat'');');
eval(abc)
normal_mtd_all(ii,:) = beta_mtd(2,:);
plank_mtd_all(ii,:) = beta_mtd(3,:);
freeze_mtd_all(ii,:) = beta_mtd(4,:);
else
abc = sprintf('%s%d%s%d%s','load(''sub-',loader(ii,1),'_run-',loader(ii,2),'_beta_bold_5.mat'');');
eval(abc)
normal_bold_all(ii,:) = bold_beta(2,:);
plank_bold_all(ii,:) = bold_beta(3,:);
freeze_bold_all(ii,:) = bold_beta(4,:);
abc = sprintf('%s%d%s%d%s','load(''sub-',loader(ii,1),'_run-',loader(ii,2),'_beta_mtd_5.mat'');');
eval(abc)
normal_mtd_all(ii,:) = beta_mtd(2,:);
plank_mtd_all(ii,:) = beta_mtd(3,:);
freeze_mtd_all(ii,:) = beta_mtd(4,:);
end
end

%% 1. MTD beta analysis
[h_mtd,p_mtd] = ttest(cope_mtd_plankvsnormal_all,zeros(60,11175),'dim',1);
sum(h_mtd)
sum(p_mtd<.0001)
sig_mtd = matify(mean(cope_mtd_plankvsnormal_all)').*matify(h_mtd'); %see matify function

pval_mat_bin = double(matify(pval'<0.05));

%% 2. BOLD beta analysis
% initial stats test analysis, plank vs normal 150 ROIs
[h_bold,p_bold] = ttest(cope_bold_all,zeros(60,150),'dim',1);
%permuted pvalue use permutation.m function
pval_bold_bin = double(pval_bold_permute<0.05);


%% 3. finding unique roi pairs that were hyperconnected in the cope mtd (plank - normal, without freeze), using sig_mtd (which is mean of individual cope)
roi_hyperconnected = find(sig_mtd(sort_mat_id,sort_mat_id)>0);
%remove row 507 because it does not exist 
roi_hyperconnected(507) = [];
for x = 1:797
to_from(x,1)=ceil(roi_hyperconnected(x)/150);
to_from(x,2) = roi_hyperconnected(x)- (floor(roi_hyperconnected(x)/150)*150);
end
unique_hyperconnected = unique(to_from); %unique does not work in the correct way
unique_hyperconnected_pair = to_from(unique_hyperconnected,:);
save('hyperconnected_pairs_beta_mtd_cope_5.mat','to_from')


unique(to_from)
roi_anticouple = find(sig_mtd(sort_mat_id,sort_mat_id)<0);
for x = 1:226
to_from(x,1)=ceil(roi_anticouple(x)/150);
to_from(x,2) = roi_anticouple(x)- (floor(roi_anticouple(x)/150)*150);
end
unique_anticouple = unique(to_from);
unique_anticouple_pair = to_from(unique_anticouple,:); %returns the pairs of unique appearing connections
save('anticouple_pairs_beta_mtd_cope_5.mat','to_from')



%% Figure plots of raw results
% subplots of beta values for [plank and normal]
figure
subplot(1,3,1)
%plank_mtd = matify(mean(plank_mtd_all)')
imagesc(plank_mtd(sort_mat_id,sort_mat_id),clims)
hold on
xline(54,'-','LineWidth',5);
yline(54,'-','LineWidth',5);
xline(78,'--r','LineWidth',5);
yline(78,'--r','LineWidth',5);
yline(140,'-.','LineWidth',5);
xline(140, '-.','LineWidth',5);
title('plank walk')
subplot(1,3,2)
imagesc(normal_mtd(sort_mat_id,sort_mat_id),clims)
hold on
xline(54,'-','LineWidth',5);
yline(54,'-','LineWidth',5);
xline(78,'--r','LineWidth',5);
yline(78,'--r','LineWidth',5);
yline(140,'-.','LineWidth',5);
xline(140, '-.','LineWidth',5);
title('normal walk')
subplot(1,3,3)
%freeze_mtd = matify(mean(freeze_mtd_all)');
imagesc(freeze_mtd(sort_mat_id,sort_mat_id),clims)
hold on
xline(54,'-','LineWidth',5);
yline(54,'-','LineWidth',5);
xline(78,'--r','LineWidth',5);
yline(78,'--r','LineWidth',5);
yline(140,'-.','LineWidth',5);
xline(140, '-.','LineWidth',5);
title('freeze')
%subplot(1,4,4)
figure
subplot(1,2,1)
imagesc((plank_mtd(sort_mat_id,sort_mat_id) - normal_mtd(sort_mat_id,sort_mat_id)).*pval_mat_bin(sort_mat_id,sort_mat_id),clims)
hold on
xline(54,'-','LineWidth',5);
yline(54,'-','LineWidth',5);
xline(78,'--r','LineWidth',5);
yline(78,'--r','LineWidth',5);
yline(140,'-.','LineWidth',5);
xline(140, '-.','LineWidth',5);
title('plank - walk')
subplot(1,2,2)
imagesc((plank_mtd(sort_mat_id,sort_mat_id) - freeze_mtd(sort_mat_id,sort_mat_id)).*pval_mat_plankvsfreeze_bin(sort_mat_id,sort_mat_id),clims)
hold on
xline(54,'-','LineWidth',5);
yline(54,'-','LineWidth',5);
xline(78,'--r','LineWidth',5);
yline(78,'--r','LineWidth',5);
yline(140,'-.','LineWidth',5);
xline(140, '-.','LineWidth',5);
title('plank - freeze')




figure
imagesc(sig_cope_beta_pc(:,sort_mat_id))




figure
subplot(1,3,1)
plank_mtd_beta = matify(mean(plank_mtd_all)');
imagesc(plank_mtd_beta(sort_mat_id,sort_mat_id))
hold on
xline(54,'-','LineWidth',3);
yline(54,'-','LineWidth',3);
xline(78,'--','LineWidth',3);
yline(78,'--','LineWidth',3);
yline(140,'-.','LineWidth',3);
xline(140, '-.','LineWidth',3);
title('plank walk')
colorbar
subplot(1,3,2)
normal_mtd_beta = matify(mean(normal_mtd_all)');
imagesc(normal_mtd_beta(sort_mat_id,sort_mat_id))
hold on
xline(54,'-','LineWidth',3);
yline(54,'-','LineWidth',3);
xline(78,'--','LineWidth',3);
yline(78,'--','LineWidth',3);
yline(140,'-.','LineWidth',3);
xline(140, '-.','LineWidth',3);
title('normal walk')
colorbar
subplot(1,3,3)
freeze_mtd_beta = matify(mean(freeze_mtd_all)');
imagesc(freeze_mtd_beta(sort_mat_id,sort_mat_id))
hold on
xline(54,'-','LineWidth',3);
yline(54,'-','LineWidth',3);
xline(78,'--','LineWidth',3);
yline(78,'--','LineWidth',3);
yline(140,'-.','LineWidth',3);
xline(140, '-.','LineWidth',3);
title('freeze')
colorbar




figure
imagesc(plank_bold(:,sort_mat_id) - normal_bold(:,sort_mat_id).*pval_bold_bin(:,sort_mat_id))
hold on
xline(54,'-','LineWidth',3);
yline(54,'-','LineWidth',3);
xline(78,'--','LineWidth',3);
yline(78,'--','LineWidth',3);
xline(140, '-.','LineWidth',3);
colorbar
title('bold beta plank vs normal')
unique_bold_sig = plank_bold(:,sort_mat_id) - normal_bold(:,sort_mat_id).*pval_bold_bin(:,sort_mat_id);



[h_bold,p_bold] = ttest(cope_bold_plankvsnormal_all,zeros(60,150),'dim',1);
sig_bold = mean(cope_bold_plankvsnormal_all).*(h_bold);

figure
imagesc(sig_bold(:,sort_mat_id))
hold on
xline(54,'-','LineWidth',3);
yline(54,'-','LineWidth',3);
xline(78,'--','LineWidth',3);
yline(78,'--','LineWidth',3);
xline(140, '-.','LineWidth',3);
yline(140, '-.','LineWidth',3);
colorbar
title('significant(t-test) mean beta bold plank vs normal')


%figure for significant mean mtd beta values, significance determined by
%t-test
figure
imagesc(sig_mtd(sort_mat_id,sort_mat_id))
hold on
xline(54,'-','LineWidth',3);
yline(54,'-','LineWidth',3);
xline(78,'--','LineWidth',3);
yline(78,'--','LineWidth',3);
xline(140, '-.','LineWidth',3);
yline(140, '-.','LineWidth',3);
colorbar
title('significant(t-test) mean beta mtd plank vs normal')
%figure plots significant mean mtd beta values, with increased
%hyperconnectivity
roi_hyperconnected = find(sig_mtd(sort_mat_id,sort_mat_id)>0);
figure
imagesc(sig_mtd(sort_mat_id,sort_mat_id)>0)
hold on
xline(54,'-','LineWidth',3);
yline(54,'-','LineWidth',3);
xline(78,'--','LineWidth',3);
yline(78,'--','LineWidth',3);
xline(140, '-.','LineWidth',3);
yline(140, '-.','LineWidth',3);
colorbar
title('significant hyperconnected beta mtd values, cope plank vs normal')
%figure of anti-coupling 
roi_anticouple = find(sig_mtd(sort_mat_id,sort_mat_id)<0);
figure
imagesc(sig_mtd(sort_mat_id,sort_mat_id)<0)
hold on
xline(54,'-','LineWidth',3);
yline(54,'-','LineWidth',3);
xline(78,'--','LineWidth',3);
yline(78,'--','LineWidth',3);
xline(140, '-.','LineWidth',3);
yline(140, '-.','LineWidth',3);
colorbar
title('significant anticoupling/decoupling beta mtd values, cope plank vs normal')

%figure for mean beta values from mtd of plank and normal (difference),
%significance determined by permutation testing
%after permutation run on plank_mtd_all and normal_mtd_all
pval_mtd_bin = double(matify(pval_mtd_permute'<0.05));
plank_mtd = matify(mean(plank_mtd_all)');
normal_mtd = matify(mean(normal_mtd_all)');

figure
imagesc((plank_mtd(sort_mat_id,sort_mat_id) - normal_mtd(sort_mat_id,sort_mat_id)).*pval_mtd_bin(sort_mat_id,sort_mat_id))
hold on
xline(54,'-','LineWidth',3);
yline(54,'-','LineWidth',3);
xline(78,'--','LineWidth',3);
yline(78,'--','LineWidth',3);
xline(140, '-.','LineWidth',3);
yline(140, '-.','LineWidth',3);
colorbar
title('mean beta mtd for plank vs normal, significance permutation')



%% Entropy Analysis

for ii = 1:60
if ii < 23
abc = sprintf('%s%d%s%d%s','load(''sub-0',loader(ii,1),'_run-',loader(ii,2),'_entropy_bold_5.mat'');');
eval(abc)
entropy_bold_all(ii,:) = ent_ts';
else
abc = sprintf('%s%d%s%d%s','load(''sub-',loader(ii,1),'_run-',loader(ii,2),'_entropy_bold_5.mat'');');
eval(abc)
entropy_bold_all(ii,:) = ent_ts';
end
end

[h_entropy,p_entropy] = ttest(entropy_bold_all,zeros(60,150),'dim',1);
mean_entropy_bold = mean(entropy_bold_all).*h_entropy;

% entropy into glm

for ii = 1:60
if ii < 23
abc = sprintf('%s%d%s%d%s','load(''sub-0',loader(ii,1),'_run-',loader(ii,2),'_cope_entropy_5.mat'');');
eval(abc)
cope_entropy_all(ii,:) = sub_cope;
else
abc = sprintf('%s%d%s%d%s','load(''sub-',loader(ii,1),'_run-',loader(ii,2),'_cope_entropy_5.mat'');');
eval(abc)
cope_entropy_all(ii,:) = sub_cope;
end
end

[h_entropy,p_entropy] = ttest(cope_entropy_all,zeros(60,1),'dim',1);





%% Creating brain plots of connectivity
mat_id = [3
3
3
3
3
3
3
3
3
3
3
3
3
3
3
3
3
3
3
9
9
9
9
9
10
10
10
10
10
10
10
11
11
11
11
11
11
11
11
11
11
11
11
11
12
12
12
12
12
12
12
12
12
12
13
13
13
13
13
3
3
3
3
3
3
3
3
3
3
3
3
3
3
3
3
3
3
3
3
4
4
4
4
4
4
4
4
4
4
4
4
4
4
4
9
9
9
9
9
9
10
10
10
10
10
10
11
11
11
11
11
11
11
11
11
11
11
12
12
12
12
12
12
12
12
12
12
12
12
12
12
12
13
13
13
13
13
13
13
18
18
18
18
18
18
18
18
18
18
19];
figure
imagesc(mat_id)
[order_mat_id,sort_mat_id] = sort(mat_id);

%cluster the networks
%order_mat_id --> need to cluster motor =1, limbic = 2, cognitive = 3,
%subcort = 4
cluster_network(order_mat_id==3,1) = 1;
cluster_network(order_mat_id==4,1) = 1;
cluster_network(order_mat_id==9,1) = 2;
cluster_network(order_mat_id==10,1) = 2;
cluster_network(order_mat_id==11,1) = 3;
cluster_network(order_mat_id==12,1) = 3;
cluster_network(order_mat_id==13,1) = 3;
cluster_network(order_mat_id==18,1) = 4;
cluster_network(order_mat_id==19,1) = 4;

schaef_include = [25
26
27
28
29
30
31
32
33
34
35
36
37
38
39
40
41
42
43
109
110
111
112
113
114
115
116
117
118
119
120
121
122
123
124
125
126
127
128
129
130
131
132
133
134
135
136
137
138
139
140
141
142
143
144
145
146
147
148
224
225
226
227
228
229
230
231
232
233
234
235
236
237
238
239
240
241
242
243
244
245
246
247
248
249
250
251
252
253
254
255
256
257
258
313
314
315
316
317
318
319
320
321
322
323
324
325
326
327
328
329
330
331
332
333
334
335
336
337
338
339
340
341
342
343
344
345
346
347
348
349
350
351
352
353
354
355
356
357];
str_temp = sig_mtd(:,146)
temp = zeros(400,1);
temp(schaef_include) = str_temp;
temp(schaef_include) = str_temp(1:139);

unq_net = unique(mat_id)
for xx = 1:8
for yy = 1:8
sig_mtd_net(xx,yy) = nanmean(nanmean(sig_mtd(mat_id==(unq_net(xx)),mat_id==(unq_net(yy)))));
end
end

figure
imagesc(sig_mtd(sort_mat_id,sort_mat_id))
imagesc(order_mat_id)
sig_mtd_reordered = sig_mtd(sort_mat_id,sort_mat_id);
figure
imagesc(sig_mtd_reordered(order_mat_id==18,order_mat_id~=18))
imagesc(order_mat_id(order_mat_id~=18)')

temp = zeros(1,400);
left_limbic = sig_m:,22:25);
left_limbic_mean = mean(left_limbic')';
temp(schaef_include) = left_limbic_mean(1:139td();
temp = temp';
surf_schaef(temp,'left_limbic_average');

%creates a schaef file for regions of interest
temp_2(schaef_include) = mat_id(1:139);
regions_schaef(schaef_include) = mat_id(1:139);
regions_schaef(:,358:400) = 0;
regions_schaef = regions_schaef';

%% Compare Anxiety with Mean Cope (plank vs normal) MTD values
load('final_hvalue_sig_mtd_cope_subjects_5.mat')
sum(h_mtd)

ans =

   512


load('positive_mtd_copebeta_pairs_sorted.mat')
positive_mtd = to_from;
load('sig_mtd_beta_cope_plankvsnormal_5.mat')
load('negative_mtd_copebeta_pairs_sorted.mat')
negative_mtd = to_from;
load('cope_mtd_5_plankvsnormal_all.mat')
load('Mean_Cope_MTD_5_persubject.mat')
PAS = [];
% locations of significant correlations will be replaced with 1
[rho,pval]= corr(Mean_Cope_MTD,PAS(:,1)); %run correlation
a = find(pval<0.05); %find significant pvalues
h_val = [1,1:11175];
h_val = zeros(1,11175);
h_val(1,a) = ones;

%sig correlation values
sig_rho = rho.*(h_val');

a = find(sig_rho > 0); %positive correlations
 sig_rho_pos = zeros(11175,1);
sig_rho_pos(a,1) = sig_rho(a,:); %replaces 0 with positive correlation

b = find(sig_rho <0); %negative correlations
sig_rho_neg = zeros(11175,1);
sig_rho_neg(b,1) = sig_rho(b,:);

sig_rho_pos_flat = matify(sig_rho_pos);%flatten into 150 x 150 matrix
sig_rho_neg_flat = matify(sig_rho_neg);

figure
subplot(1,3,1)
imagesc(sig_rho_pos_flat(sort_mat_id,sort_mat_id))
hold on
xline(54,'-','LineWidth',3);
yline(54,'-','LineWidth',3);
xline(78,'--','LineWidth',3);
yline(78,'--','LineWidth',3);
xline(140, '-.','LineWidth',3);
yline(140, '-.','LineWidth',3);
colorbar
title('sig positive corr with PAS total')
subplot(1,3,2)
imagesc(sig_rho_neg_flat(sort_mat_id,sort_mat_id))
hold on
xline(54,'-','LineWidth',3);
yline(54,'-','LineWidth',3);
xline(78,'--','LineWidth',3);
yline(78,'--','LineWidth',3);
xline(140, '-.','LineWidth',3);
yline(140, '-.','LineWidth',3);
colorbar
title('sig negative corr with PAS total')
subplot(1,3,3)
imagesc(sig_mtd(sort_mat_id,sort_mat_id))
hold on
xline(54,'-','LineWidth',3);
yline(54,'-','LineWidth',3);
xline(78,'--','LineWidth',3);
yline(78,'--','LineWidth',3);
xline(140, '-.','LineWidth',3);
yline(140, '-.','LineWidth',3);
colorbar
title('sig mean mtd cope')

%find the pairs of roi that are positive correlated with anxiety
rho_roi_pairs_pos = find(sig_rho_pos_flat(sort_mat_id,sort_mat_id)>0);
for x = 1:240
to_from(x,1)=ceil(rho_roi_pairs_pos(x)/150);
to_from(x,2) = rho_roi_pairs_pos(x)- (floor(rho_roi_pairs_pos(x)/150)*150);
end
pos_rho_pairs = to_from;

rho_roi_pairs_neg = find(sig_rho_neg_flat(sort_mat_id,sort_mat_id)<0);
for x = 1:248
to_from(x,1)=ceil(rho_roi_pairs_neg(x)/150);
to_from(x,2) = rho_roi_pairs_neg(x)- (floor(rho_roi_pairs_neg(x)/150)*150);
end
neg_rho_pairs = to_from;

%need to compare the positive correlation with significant positive MTD
%value, note all these have been sorted in MAC networks
sig_mean_cope_mtd_plankvsnormal_all = mean(cope_mtd_plankvsnormal_all);
sig_mean_cope_mtd_plankvsnormal_all =  sig_mean_cope_mtd_plankvsnormal_all.*h_mtd;
a = find(sig_mean_cope_mtd_plankvsnormal_all >0);
sig_mean_cope_mtd_pos = zeros(1,11175);
sig_mean_cope_mtd_pos(1,a) = sig_mean_cope_mtd_plankvsnormal_all(1,a); %with cope mtd values
sig_mean_cope_mtd_pos(1,a) = ones;
%constraint to positive significant mtd cope values
sig_rho_pos_constraint = sig_rho_pos.*(sig_mean_cope_mtd_pos');
sig_rho_pos_constraint_flat = matify(sig_rho_pos_constraint);
rho_roi_pairs_pos_constraint = find(sig_rho_pos_constraint_flat(sort_mat_id,sort_mat_id)>0);
for x = 1:16
to_from(x,1)=ceil(rho_roi_pairs_pos_constraint(x)/150);
to_from(x,2) = rho_roi_pairs_pos_constraint(x)- (floor(rho_roi_pairs_pos_constraint(x)/150)*150);
end
rho_roi_constraint_pos = to_from;

%need to compare the negative correlation with significant negative MTD
%value
a = find(sig_mean_cope_mtd_plankvsnormal_all <0);
sig_mean_cope_mtd_neg = zeros(1,11175);
sig_mean_cope_mtd_neg(1,a) = sig_mean_cope_mtd_plankvsnormal_all(1,a); %with cope mtd values
sig_mean_cope_mtd_neg(1,a) = ones;

sig_rho_neg_constraint = sig_rho_neg.*(sig_mean_cope_mtd_neg');
sig_rho_neg_constraint_flat = matify(sig_rho_neg_constraint);
rho_roi_pairs_neg_constraint = find(sig_rho_neg_constraint_flat(sort_mat_id,sort_mat_id)<0);
for x = 1:2
to_from(x,1)=ceil(rho_roi_pairs_neg_constraint(x)/150);
to_from(x,2) = rho_roi_pairs_neg_constraint(x)- (floor(rho_roi_pairs_neg_constraint(x)/150)*150);
end
rho_roi_constraint_neg = to_from;