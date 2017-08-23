function OptoISI(spikeData,TTLtimes,keepCell,axisHandle)

for cellNum=1:length(keepCell)
    %find(mean(meanChan)==max(mean(meanChan)));
    if ~exist('axisHandle','var') || isempty(axisHandle)
        figure('Position',[1092 149 708 761]); hold on
    end
    colormap(parula); %cmap=colormap;
    %% Spike times
    spikeTimes=spikeData.spikeTimes(spikeData.unitsIdx==keepCell(cellNum),:)/(spikeData.samplingRate/1000);
    
    %get wich spike time occur during TTL
    pulseIdx=false(size(spikeTimes,1),size(TTLtimes,1));
    %     figure; hold on; plot(spikeTimes,'*'); plot(TTLtimes,'d')
    for TTLNum=1:size(TTLtimes,1)
        pulseIdx(:,TTLNum)=spikeTimes>TTLtimes(TTLNum) & spikeTimes<TTLtimes(TTLNum)+5;
    end
    onSpikes=logical(sum(pulseIdx,2));
    
    unitST_onPulse=spikeTimes(onSpikes);
    unitST_offPulse=spikeTimes(~onSpikes);
    % compute interspike interval
    ISI_onPulse=diff(unitST_onPulse);
    ISI_offPulse=diff(unitST_offPulse);
    
    ISI_offPulsehist=histogram(double(ISI_offPulse),0:10:max(ISI_offPulse)+1);  %,'Normalization','probability'
    ISI_offPulsehist.FaceColor = 'w';
    ISI_offPulsehist.EdgeColor = 'k';
    ISI_onPulsehist=histogram(double(ISI_onPulse),0:10:max(ISI_onPulse)+1);  %,'Normalization','probability'
    ISI_onPulsehist.FaceColor = [0.3 0.75 0.93];
    ISI_onPulsehist.EdgeColor = 'k';
    xlabel('Inter-spike Interval distribution (ms)')
    axis('tight');box off;
    legend([ISI_onPulsehist,ISI_offPulsehist],{'Pulse On','Pulse Off'},'location','northeast')
    set(gca,'xlim',[0 1050],'Color','white','FontSize',10,'FontName','calibri','TickDir','out');
    hold off
%     title(['Neuron ' num2str(keepCell(cellNum)) ' ISI'])  
    set(gca,'Color','white','FontSize',12,'FontName','calibri');
    
end


