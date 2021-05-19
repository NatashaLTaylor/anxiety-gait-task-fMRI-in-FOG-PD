function [ci,q,part,z,hc] = integration(data,gamma)
%INTEGRATION       Creates time-resolved network topological measures
%
%  [ci,q,p,z,hc] = integration(data,gammambeta);
%
%  This code takes a time-resolved connectivity matrix and 
%  estimates community structure, modularity and the cartographic profile
%  for each region within a region x region x time connectivity matrix.
%  See https://arxiv.org/abs/1511.02976 for more details.
%
%  Requirements: https://sites.google.com/site/bctnet/
%
%  Input:      data     time-series organized in 'nodes x nodes x time' matrix
%              gamma    tuning parameter for louvain algorithm (low = large modules & high = small modules)
%              beta     similarity measure used to determine clustering assignment (hungarian algorithm) - requires munkres.m & apcluster.m
%
%  Output:     ci       time-resolved community assignment
%              q        time-resolved modularity
%              p        time-resolved participation coefficient
%              z        time-resolved module-degree z-score
%              hc       cartographic profile
%              f        flexibility


    %define variables
    
    [nodes,~,time] = size(data);

    ci = zeros(nodes,time); q = zeros(time,1); part = zeros(nodes,time); z = zeros(nodes,time);

    for tt = 1:time
        [ci(:,tt),q(tt,1)] = community_louvain(data(:,:,tt),gamma,1:1:nodes,'negative_asym');
        part(:,tt) = participation_coef_sign(data(:,:,tt),ci(:,tt));
        z(:,tt) = module_degree_zscore(data(:,:,tt),ci(:,tt));
    end

    %cartographic profile

    xbins = [0:0.01:1]; ybins = [5:-.1:-5];
    hc = zeros(size(xbins,2),size(ybins,2),time);
    xNumBins = numel(xbins); yNumBins = numel(ybins);

    for tt = 1:time
        Xi = round(interp1(xbins,1:xNumBins,part(:,tt),'linear','extrap'));
        Yi = round(interp1(ybins,1:yNumBins,z(:,tt),'linear','extrap'));
        Xi = max(min(Xi,xNumBins),1);
        Yi = max(min(Yi,yNumBins),1);
        hc(:,:,tt) = accumarray([Yi(:) Xi(:)], 1, [yNumBins xNumBins]);
    end

end


