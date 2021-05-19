function aaaa = MRItime_process(filename, subnum, run)
%% load in processed 4D nifti file
%filename = full length name of fMRI file in nii.gz form
%subnum = subject number
%run = run, session
%change into location of denoised fmri data
home_path = sprintf('%s',subnum);

cd(home_path)
% data is fMRI scan that has been denoised using https://github.com/arielletambini/denoiser
data = MRIread(filename); %use freesurfer MRIread function
% parc is parcellation file containing binarized region of interest
parc = MRIread('/Volumes/PRJ-shiver/Shiver_Fest/Natasha/Study-2_Dynamic_Connectivity_Task-FOG/Parcellations/parc_1_schaef+roi.nii.gz'); 
nNodes = 410; % cortex = 1-400; parc1 = parc_1_schaef+roi.nii.gz

%% extract the time series

sub_ts = zeros(nNodes,data.nframes); %pre-defining a matrix that is growing with each loop makes matlab WAY faster...

for t = 1:data.nframes %for each time point
    
    temp = data.vol(:,:,:,t); %create a temporary file with the "t"-th slice of the data
    temp = temp(:); %flatten it out
    
    for j = 1:nNodes %for each region in turn
        temp1 = find(parc.vol==j); %find the points in the volume where the ROI lives
        temp2 = temp(temp1); %extract the same points from the flattened temp file
        sub_ts(j,t) = sum(temp2)/size(temp,1); %calculate the average value across the ROI points
    end
    
    sprintf('%f%s',t/data.nframes*100,' % completed') %ticker to keep you up to date on progress
    
end

save(['/Volumes/PRJ-shiver/Shiver_Fest/Natasha/Study-2_Dynamic_Connectivity_Task-FOG/TS_extract/' subnum '_' run '-ts_parc1'],'sub_ts') %main output that needs to be saved

%% calculate correlation matrix
ts_corr = corr(sub_ts'); %calculate the Pearson's correlation between regions
save(['/Volumes/PRJ-shiver/Shiver_Fest/Natasha/Study-2_Dynamic_Connectivity_Task-FOG/TS_extract/' subnum '_' run '-ts_corr_parc1'],'ts_corr') %main output that needs to be saved

%% you can break up the correlation matrix into interpretable bits if you'd like

% cort_id is the network ID for each cortical region (i.e., region 1-400)
% 1-2 = visual, 3-4 = motor, 5-6 = dorsal attention, 7-8 = ventral
% attention, 9-10 = limbic, 11-13 = control, 14-16 = default & 17 =
% temporo-parietal
cort_id = [1;1;1;1;1;1;1;1;1;1;1;1;2;2;2;2;2;2;2;2;2;2;2;2;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;5;5;5;5;5;5;5;5;5;5;5;5;5;6;6;6;6;6;6;6;6;6;6;6;6;6;7;7;7;7;7;7;7;7;7;7;7;7;7;7;7;8;8;8;8;8;8;8;8;9;9;9;9;9;10;10;10;10;10;10;10;11;11;11;11;11;11;11;11;11;11;11;11;11;12;12;12;12;12;12;12;12;12;12;13;13;13;13;13;14;14;14;14;14;14;14;14;14;14;14;14;14;14;14;14;14;14;15;15;15;15;15;15;15;15;15;15;15;15;15;15;15;15;15;15;15;15;15;16;16;16;16;16;16;16;17;17;17;17;17;17;1;1;1;1;1;1;1;1;1;1;1;1;2;2;2;2;2;2;2;2;2;2;2;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;3;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;5;5;5;5;5;5;5;5;5;5;5;5;5;5;6;6;6;6;6;6;6;6;6;6;6;6;7;7;7;7;7;7;7;7;7;7;7;7;7;7;7;7;7;7;7;8;8;8;8;8;8;8;8;8;9;9;9;9;9;9;10;10;10;10;10;10;11;11;11;11;11;11;11;11;11;11;11;12;12;12;12;12;12;12;12;12;12;12;12;12;12;12;13;13;13;13;13;13;13;14;14;14;14;14;14;14;14;14;14;14;14;14;14;14;14;15;15;15;15;15;15;15;15;15;15;15;16;16;16;16;16;16;17;17;17;17;17;17;17;17;17;17];
nNet = max(cort_id);

% you could take the mean of each like this:
ts_corr_net = zeros(nNet);

for x = 1:nNet
    for y = 1:nNet
        ts_corr_net(x,y) = mean(mean(ts_corr(cort_id==x,cort_id==y)));
    end
    save(['/Volumes/PRJ-shiver/Shiver_Fest/Natasha/Study-2_Dynamic_Connectivity_Task-FOG/TS_extract/' subnum '_' run '-ts_corr_net_parc1'],'ts_corr_net')
end

cd .. %takes back to original directory

%% if you can get the data from each subject processed through this pipeline



