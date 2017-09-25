function OptoWaveforms(spikeData,TTLtimes,keepCell,axisHandle)

for cellNum=1:length(keepCell)
    %find(mean(meanChan)==max(mean(meanChan)));
    if ~exist('axisHandle','var') || isempty(axisHandle)
        figure('Position',[1092 149 708 761]); hold on
    end
    colormap(parula); %cmap=colormap;
    %% Spike waveform
    waveForms=double(spikeData.waveForms(:,spikeData.unitsIdx==keepCell(cellNum))');
    spikeTimes=spikeData.spikeTimes(spikeData.unitsIdx==keepCell(cellNum),:)/(spikeData.samplingRate/1000);
    
    %get wich spike time occur during TTL
    pulseIdx=false(size(spikeTimes,1),size(TTLtimes,1));
    %     figure; hold on; plot(spikeTimes,'*'); plot(TTLtimes,'d')
    for TTLNum=1:size(TTLtimes,1)
        pulseIdx(:,TTLNum)=spikeTimes>TTLtimes(TTLNum) & spikeTimes<TTLtimes(TTLNum)+5;
    end
    onSpikes=logical(sum(pulseIdx,2));
    
    %off-pulse waveforms
    offpulseSDFploth=plot(mean(waveForms(~onSpikes,:)),'linewidth',2,'color','k'); %cmap(cellNum,:)
    
    wfSEM=std(waveForms(~onSpikes,:))/ sqrt(size(waveForms(~onSpikes,:),2)); %standard error of the mean
    wfSEM = wfSEM * 1.96; % 95% of the data will fall within 1.96 standard deviations of a normal distribution
    patch([1:length(wfSEM),fliplr(1:length(wfSEM))],...
        [mean(waveForms(~onSpikes,:))-wfSEM,fliplr(mean(waveForms(~onSpikes,:))+wfSEM)],...
        'k','EdgeColor','none','FaceAlpha',0.2); %cmap(cellNum,:)
    
    % on-pulse waveforms
    onpulseSDFploth=plot(mean(waveForms(onSpikes,:)),'linewidth',2,'color',[0.3 0.75 0.93]);
    
    wfSEM=std(waveForms(onSpikes,:))/ sqrt(size(waveForms(onSpikes,:),2)); %standard error of the mean
    wfSEM = wfSEM * 1.96; % 95% of the data will fall within 1.96 standard deviations of a normal distribution
    patch([1:length(wfSEM),fliplr(1:length(wfSEM))],...
        [mean(waveForms(onSpikes,:))-wfSEM,fliplr(mean(waveForms(onSpikes,:))+wfSEM)],...
        [0.3 0.75 0.93],'EdgeColor','none','FaceAlpha',0.5);
    
    
    set(gca,'xTick',[0,15,30,45])
    figXlabels=get(gca,'xTickLabel');
    set(gca,'xTickLabel',cellfun(@(xlabl) num2str(str2double(xlabl)/30),...
        figXlabels,'UniformOutput',false),'FontSize',10,'FontName','calibri','TickDir','out');
    axis('tight');box off;
    xlabel('Time (ms)')
    ylabel('Voltage (\muV)');
    set(gca,'Color','white','FontSize',18,'FontName','calibri');
    legend([onpulseSDFploth,offpulseSDFploth],{'Pulse-evoked spikes','Spontaneous spikes'},'FontSize',12,'location','southeast');
    legend('boxoff')
    
end