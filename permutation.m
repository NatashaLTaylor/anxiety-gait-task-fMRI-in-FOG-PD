function [pval]= permutation(data1, iter)
%data1 in this case will be the mean of mtd/bold COPE for each trial
%iter number of interations (standard 5000)
% 
%vector with original values 
%MTD_COPE or BOLD_COPE with same number of rows 

%make a group ID, with the first N values = 1 (i.e., the original values) and the second N values = 2 (i.e., the zeros)


size1 = size(data1);
data2 = zeros(size1); 
stdv_data1 = std(data1);
% for a = 1:150
% data2(:,a) = normrnd(0,[stdv_data1(1,a)],60,1);
% end


data_permute = vertcat(data1,data2);
grp_size1 = size(data1,1);
grp_id(1:(grp_size1),1) =1; %first 1:61 is original values
grp_id(grp_size1 + 1:grp_size1*2,1)=2; %62:122 is zero values


for x = 1:iter
    rand_vec = rand(grp_size1*2,1);
    [~,sort_rand] = sort(rand_vec);
    grp_rand = grp_id(sort_rand);
    new_diff(x,:) = mean(data_permute(grp_rand==1,:)) - mean(data_permute(grp_rand==2,:));
    %calculate difference between 1s and 0s for each edge for 5000
    
end
 
orig_diff = mean(data_permute(grp_id==1,:)) - mean(data_permute(grp_id==2,:));
%significance of it happening
for xx = 1:150
    delta_above(xx) = sum(orig_diff(xx)>new_diff(:,xx))/iter;
    delta_below(xx) = sum(orig_diff(xx)<new_diff(:,xx))/iter;
end

pval = min(vertcat(delta_above,delta_below));
end



