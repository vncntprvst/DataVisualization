function plot_rasters_psth(rasterData,alignSpecs,figureStyle,legends)

    %# method 1
    figure('Position',[1050 120 750 790]);
    subplot(4,1,1:3)
    colormap(figureStyle{1});
    imagesc(rasterData); %     imagesc(zscore(spikeRasters,[])); % spikeArray
    % imagesc(MeanChan);
    ylabel(legends{2},'FontWeight','bold','FontSize',12);
    % draw alignment bar
    currylim=get(gca,'YLim');
%     currxlim=get(gca,'XLim');
    set(gca,'XTick',figureStyle{2},'TickDir','out');
    set(gca,'XTickLabel',figureStyle{3});
    %opto stim patch
    patch([repmat(alignSpecs{1},1,2) repmat(alignSpecs{1}+alignSpecs{3},1,2)], ...
        [[0 currylim(2)] fliplr([0 currylim(2)])], ...
        [0 0 0 0],alignSpecs{4},'EdgeColor','none','FaceAlpha',0.3);
    box off; axis tight
    title(legends{3});
    %     hcb = colorbar('eastoutside');
    %     hcb.Label.String = legends{4};
    
    %#method 2
    %     [indy, indx] = ind2sub(size(spikeRasters),find(spikeRasters)); %find row and column coordinates of spikes
    %     figure;
    %     plot([indx';indx'],[indy';indy'+1],'color','k','LineStyle','-'); % plot rasters
    
    %% plot psth
    subplot(4,1,4)
    barPlot=bar(mean(rasterData));
%     barPlot.FaceColor=[0.1 0.4 0.8];
    barPlot.EdgeColor='none';
    barPlot.BarWidth=1;
    set(gca,'XTick',figureStyle{2},'TickDir','out');
    set(gca,'XTickLabel',figureStyle{3});
    box off; axis tight
%     set(gca,'ylim',[0 0.5]);
    xlabel(legends{1},'FontWeight','bold','FontSize',12);
end


