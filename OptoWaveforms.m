function onSpikes=OptoWaveforms(spikeData,TTLtimes,keepCell,duration,axisHandle)

for cellNum=1:length(keepCell)
    %find(mean(meanChan)==max(mean(meanChan)));
    if ~exist('axisHandle','var') || isempty(axisHandle)
        figure('Position',[1092 149 708 761]); hold on
    end
    colormap(parula); %cmap=colormap;
    %% Spike waveform
    waveforms=double(spikeData.waveforms(spikeData.unitID==keepCell(cellNum),:));
    times=single(spikeData.times(spikeData.unitID==keepCell(cellNum),:));
    
    %get wich spike times occur during TTL
    pulseIdx=false(size(times,1),size(TTLtimes,1));
%     figure; hold on; plot(times,ones(length(times),1),'*');...
%     plot(TTLtimes,ones(length(TTLtimes),1),'d')
    for TTLNum=1:length(TTLtimes)
        pulseIdx(:,TTLNum)=times>=TTLtimes(TTLNum) & times<=TTLtimes(TTLNum)+max([duration 0.010]); %keep min window of 10ms 
    end
    onSpikes=any(pulseIdx,2);
    
    %off-pulse waveforms
    offpulseSDFploth=plot(mean(waveforms(~onSpikes &...
        times>=TTLtimes(1) & times<=TTLtimes(end),:)),'linewidth',2,'color','k'); %cmap(cellNum,:)
    
    wfSEM=std(waveforms(~onSpikes,:))/ sqrt(size(waveforms(~onSpikes,:),2)); %standard error of the mean
    wfSEM = wfSEM * 1.96; % 95% of the data will fall within 1.96 standard deviations of a normal distribution
    patch([1:length(wfSEM),fliplr(1:length(wfSEM))],...
        [mean(waveforms(~onSpikes,:))-wfSEM,fliplr(mean(waveforms(~onSpikes,:))+wfSEM)],...
        'k','EdgeColor','none','FaceAlpha',0.2); %cmap(cellNum,:)
    
    % on-pulse waveforms
    onpulseSDFploth=plot(mean(waveforms(onSpikes,:)),'linewidth',2,'color',[0.3 0.75 0.93]);
    
    wfSEM=std(waveforms(onSpikes,:))/ sqrt(size(waveforms(onSpikes,:),2)); %standard error of the mean
    wfSEM = wfSEM * 1.96; % 95% of the data will fall within 1.96 standard deviations of a normal distribution
    patch([1:length(wfSEM),fliplr(1:length(wfSEM))],...
        [mean(waveforms(onSpikes,:))-wfSEM,fliplr(mean(waveforms(onSpikes,:))+wfSEM)],...
        [0.3 0.75 0.93],'EdgeColor','none','FaceAlpha',0.5);
    
    set(gca,'xTick',[0,15,30,45])
    figXlabels=get(gca,'xTickLabel');
    set(gca,'xTickLabel',cellfun(@(xlabl) num2str(str2double(xlabl)/30),...
        figXlabels,'UniformOutput',false),'FontSize',10,'FontName','Helvetica','TickDir','out');
    axis('tight');box off;
    xlabel('Time (ms)')
    ylabel('Voltage (\muV)');
    set(gca,'Color','white','FontSize',10,'FontName','Calibri');
%     legend([onpulseSDFploth,offpulseSDFploth],{'Pulse-evoked spikes','Spontaneous spikes'},'FontSize',8,'location','southeast');
%     legend('boxoff')
    
end