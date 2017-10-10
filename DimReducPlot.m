function labels=DimReducPlot(snippets,labels,option)

%load('040v_1002_2Hz20ms_30mW_22Ch_nopp_Ch12.mat')



if strcmp(option,'tsne')
    comps = tsne(single(snippets));
elseif strcmp(option,'pca')
    [~,comps] = pca(single(snippets));
elseif strcmp(option,'ftsne')
    numDims=2;
    pcaDims=50;
    perplexity=50;
    theta=0.5;
    alg='svd';
    comps = fast_tsne(snippets, numDims, pcaDims, perplexity, theta, alg);
end

cmap=lines; cmap=[0,0,0;cmap];
figure; scatter(comps(:,2), comps(:,2), labels ,cmap,'.',7);


dist = pdist2(comp,comp);
[L, center_idxs] = cluster_dp(dist); %get from https://github.com/shihai1991/kernel-density-peaks


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



% options.method = 'gaussian';
% options.percent = 2;
% options.sigma = 20;
% options.trim_halo = false;
% options.n_clusters = Inf;
% options.show_plot = false;
% options.M = 10;


