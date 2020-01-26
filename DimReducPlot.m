function [labels ,centroids]=DimReducPlot(snippets,initLabels,option)

if size(snippets,1)<size(snippets,2)
    snippets=snippets';
end

if strcmp(option{1},'tsne')
    comps = tsne(single(snippets));
elseif strcmp(option{1},'pca')
    [~,comps] = pca(single(snippets));
elseif strcmp(option{1},'ftsne')
    numDims=2;
    pcaDims=50;
    perplexity=50;
    theta=0.5;
    alg='svd';
    comps = fast_tsne(snippets, numDims, pcaDims, perplexity, theta, alg);
end

 
% dist = pdist2(comps,comps);

% [L, center_idxs] = cluster_dp(dist); %get from https://github.com/shihai1991/kernel-density-peaks

% labels = cluster_dp(comps, 3, 0, 0);

[labels ,centroids]= kmeans(comps,3,'Distance','sqeuclidean','Replicates',5);

if numel(option)>1 && strcmp(option{2},'plot')
    figure; hold on
%     gscatter(comps(:,1), comps(:,2), initLabels ,[parula ones(size(parula,1),1)*0.5],'.',8,'doleg','off');
    cmap=lines; cmap=[0,0,0;cmap];
    gscatter(comps(:,1), comps(:,2), labels ,cmap,'.',7);
end

%
% clusterCat = gmdistribution(PrComps(1,:),...
%     [abs(min(PrComps(:,1)))+abs(max(PrComps(:,1))),...
%     -0; -0, abs(min(PrComps(:,2)))+abs(max(PrComps(:,2)))]);
% mahalDis = mahal(clusterCat,PrComps(2:end,:));
% % scatter(PrComps(:,1),PrComps(:,2),50,mahalDis(:,1),'.')
% % hb = colorbar;
% % ylabel(hb,'Mahalanobis Distance to Component 1')
%
% similIdx=find(mahalDis<max(mahalDis)/10);
% % figure
% % scatter(PrComps(similIdx,1), PrComps(similIdx,2), 'k.');



