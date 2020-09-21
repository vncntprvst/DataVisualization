function OptoJitter(spikeData,TTLtimes,keepCell,duration,axisHandle)
%% plots spike latency standard deviation
% Input: alignedRasters - cell array with stim aligned rasters
%         alignTime - time used to align responses
%         pulseDur - pulse duration 

for cellNum=1:length(keepCell)
    %find(mean(meanChan)==max(mean(meanChan)));
    if ~exist('axisHandle','var') || isempty(axisHandle)
        figure('Position',[1092 149 708 761]); hold on
    end
    
    times=spikeData.times(spikeData.unitID==keepCell(cellNum),:);
    
    %get wich spike times occur during TTL
    spikeLatency=nan(numel(TTLtimes),1);
    for TTLNum=1:length(TTLtimes)
        try
            spikeLatency(TTLNum)=times(find(times>=TTLtimes(TTLNum),1))-TTLtimes(TTLNum);
        catch
            continue
        end
    end

    hold on
    errorbar(1,mean(spikeLatency),std(spikeLatency),'color','b')
    plot(1,mean(spikeLatency),'marker','o','MarkerFaceColor','b',...
        'MarkerEdgeColor','none','MarkerSize',8);%/sqrt(numel(spikeLatency))*1.96
    errorbar(2,mean(diff(times)),std(double(diff(times))),'color','k')
    plot(2,mean(diff(times)),'marker','o','MarkerFaceColor','k',...
        'MarkerEdgeColor','none','MarkerSize',8);%/sqrt(numel(spikeLatency))*1.96
    set(gca,'xlim',[0 3],'ylim',[0 max(get(gca,'ylim'))],...
        'xTick',[1,2],'xTickLabel',{'latency','mean ISI'},...
        'FontSize',10,'FontName','Helvetica','TickDir','out');
    box off;
    xlabel({'mean response latency vs';' mean ISI (+/- SD)'})
    ylabel('Time (s)');
    set(gca,'Color','white','FontSize',12,'FontName','Helvetica');

end
