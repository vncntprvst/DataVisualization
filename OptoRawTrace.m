function OptoRawTrace(recTrace,spikeTimes,msConv,TTLtimes,pulseDur,option,axisHandle)

if ~exist('axisHandle','var') | isempty(axisHandle)
    figure('Position',[1092 149 708 761]); hold on
end
cmap=parula;
%plot raw trace
plot(recTrace.data,'color','k');

if ~isempty(TTLtimes)
    for TTLNum=1:length(TTLtimes)
        patch([TTLtimes(TTLNum), TTLtimes(TTLNum),...
            TTLtimes(TTLNum)+pulseDur*msConv, TTLtimes(TTLNum)+pulseDur*msConv], ...
            [get(gca,'ylim') fliplr(get(gca,'ylim'))], ...
            [0 0 0 0],[0.3 0.75 0.93],'EdgeColor','none','FaceAlpha',0.7);
    end
    if strcmp(option,'center') %length(TTLtimes)==1
        set(gca,'xlim',[TTLtimes(end)-50*msConv TTLtimes(end)+200*msConv])
        set(gca,'xtick',TTLtimes(end)-50*msConv:10*msConv:TTLtimes(end)+200*msConv,...
            'xticklabel',-50:10:200);
        set(gca,'ytick',[],'yticklabel',[],'TickDir','out');
    else
        %     set(gca,'xtick',recTrace.xTicks,'xticklabel',recTrace.xTicklabels,...
        %         'ytick',[],'yticklabel',[],'TickDir','out');
        set(gca,'xtick',linspace(0,2*recTrace.excerptSize,9),...
            'xticklabel',round(linspace(0,2*recTrace.excerptSize,9)/double(msConv)),...
            'ytick',[],'yticklabel',[],'TickDir','out');
    end
end
box off;
xlabel('Time (ms)');
set(gca,'Color','white','FontSize',18,'FontName','Helvetica');

for cellNum=1:length(spikeTimes)
    if ~isempty(spikeTimes{cellNum}) & ~isnan(spikeTimes{cellNum})
        %plot spike id labels
        spkLabelYLoc=ones(1,size(spikeTimes{cellNum},2))*(min(get(gca,'ylim'))/4*3);
        plot(double(spikeTimes{cellNum})-double(recTrace.location-recTrace.excerptSize),...
            spkLabelYLoc,'Color',[cmap(cellNum,:),0.4],...
            'linestyle','none','Marker','^');
        
        % plot unit markers
        %                         plot(spkTimes{unitP}-(handles.rawDataInfo.excerptLocation-handles.rawDataInfo.excerptSize),...
        %                         rasterHeight,'Color','k',...
        %                         'linestyle','none','Marker','*');
        
        
    end
end