function [pval]= permutation2(data1,data2, iter)
%data1 in this case will be the data for one group
%data2 will be data for another group
%iter number of interations (standard 5000)
% 
%vector with original values  

%make a group ID, with the first N values = 1 (i.e., the original values) and the second N values = 2 (i.e., the zeros)


size1 = size(data1);
data_permute = vertcat(data1,data2);
grp_size1 = size(data1,1);
grp_size2 = size(data2,1);
grp_id1(1:(grp_size1),1) =1; %first number of subject in data1 is original values
grp_id2(1:(grp_size2),1)=2;
grp_id = vertcat(grp_id1,grp_id2);
grp_size = size(grp_id,1);


for x = 1:iter
    rand_vec = rand(grp_size,1);
    [~,sort_rand] = sort(rand_vec);
    grp_rand = grp_id(sort_rand);
    new_diff(x,:) = mean(data_permute(grp_rand==1,:)) - mean(data_permute(grp_rand==2,:));
    %calculate difference between 1s and 0s for each edge for 5000
    sprintf('%d',x)
end
 
orig_diff = mean(data_permute(grp_id==1,:)) - mean(data_permute(grp_id==2,:));
%significance of it happening
a= size(data1,2);
for xx = 1:a
    delta_above(xx) = sum(orig_diff(xx)>new_diff(:,xx))/iter;
    delta_below(xx) = sum(orig_diff(xx)<new_diff(:,xx))/iter;
    sprintf('%d',xx)
end

pval = min(vertcat(delta_above,delta_below));
end



