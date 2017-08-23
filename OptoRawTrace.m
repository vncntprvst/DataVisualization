function OptoRawTrace(recTrace,spikeTimes,msConv,TTLtimes,axisHandle)


for cellNum=1:length(spikeTimes)
    if ~exist('axisHandle','var') | isempty(axisHandle)
        figure('Position',[1092 149 708 761]); hold on
    end
    colormap(parula); cmap=colormap;
    %plot raw trace
    plot(recTrace.data,'color','k');
    %plot spike id labels
    %     spkLabelYLoc=ones(1,size(spikeTimes{cellNum},2))*(min(get(gca,'ylim'))/4*3);
    %     plot(single(spikeTimes{cellNum})-(recTrace.location-recTrace.excerptSize),...
    %         spkLabelYLoc,'Color',[cmap(cellNum,:),0.4],...
    %         'linestyle','none','Marker','^');
    
    for TTLNum=1:length(TTLtimes)
        patch([TTLtimes(TTLNum), TTLtimes(TTLNum),...
            TTLtimes(TTLNum)+2*msConv, TTLtimes(TTLNum)+2*msConv], ...
            [get(gca,'ylim') fliplr(get(gca,'ylim'))], ...
            [0 0 0 0],[0.3 0.75 0.93],'EdgeColor','none','FaceAlpha',0.5);
    end
    if length(TTLtimes)==1
        set(gca,'xlim',[TTLtimes-50*msConv TTLtimes+50*msConv])
        set(gca,'xtick',TTLtimes-50*msConv:10*msConv:TTLtimes+50*msConv,'xticklabel',-50:10:50);
        set(gca,'ytick',[],'yticklabel',[],'TickDir','out');
        
    else
        %     set(gca,'xtick',recTrace.xTicks,'xticklabel',recTrace.xTicklabels,...
        %         'ytick',[],'yticklabel',[],'TickDir','out');
        set(gca,'xtick',linspace(0,2*recTrace.excerptSize,50),'xticklabel',linspace(0,2*recTrace.excerptSize,50)/double(msConv),...
            'ytick',[],'yticklabel',[],'TickDir','out');
        
        
    end
    box off;
    xlabel('Time (ms)');
    set(gca,'Color','white','FontSize',12,'FontName','calibri');
end