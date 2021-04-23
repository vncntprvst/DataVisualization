function OptoISI(spikeData,TTLtimes,keepCell,pulseDur,axisHandle)

for cellNum=1:length(keepCell)
    %find(mean(meanChan)==max(mean(meanChan)));
    if ~exist('axisHandle','var') || isempty(axisHandle)
        figure('Position',[1092 149 708 761]); hold on
    end
    colormap(parula); %cmap=colormap;
    %% Spike times
    spikeTimes=spikeData.times(spikeData.unitID==keepCell(cellNum),:);
    
    %get wich spike time occur during TTL
    pulseIdx=false(size(spikeTimes,1),size(TTLtimes,1));
    %     figure; hold on; plot(spikeTimes,'*'); plot(TTLtimes,'d')
    for TTLNum=1:size(TTLtimes,1)
        pulseIdx(:,TTLNum)=spikeTimes>TTLtimes(TTLNum) & spikeTimes<TTLtimes(TTLNum)+pulseDur;
    end
    onSpikes=any(pulseIdx,2);
    
    unitST_onPulse=spikeTimes(onSpikes);
    unitST_offPulse=spikeTimes(~onSpikes);
    % compute interspike interval
    ISI_onPulse=diff(unitST_onPulse)*1000;
    ISI_offPulse=diff(unitST_offPulse)*1000;
    
    
    ISI_offPulsehist=histogram(double(ISI_offPulse),logspace(0, 4, 50),...
        'Normalization','probability','LineWidth',0.5);
    ISI_offPulsehist.FaceColor = [1.0000    0.6784    0.0980]; %[1.0000    0.8941    0.7686];
    ISI_offPulsehist.EdgeColor = 'k';
    ISI_onPulsehist=histogram(double(ISI_onPulse),logspace(0, 4, 50),...
        'Normalization','probability','LineWidth',0.5);
    ISI_onPulsehist.FaceColor = [0.3 0.75 0.93];
    ISI_onPulsehist.EdgeColor = 'k';
    
        
    xlabel('ISI distribution (ms)')
    ylabel('Spike probability')
    axis('tight'); box off; grid('on'); set(gca,'xscale','log','GridAlpha',0.25,'MinorGridAlpha',1);
%     legend([ISI_onPulsehist,ISI_offPulsehist],{'Pulse On','Pulse Off'},'location','northeast')
%     legend('boxoff')
    set(gca,'xlim',[0 10^3],'Color','white','FontSize',10,'FontName','calibri','TickDir','out');
    set(gca,'XTick',[1,10,100,1000],'XTickLabel',[1,10,100,1000]);
    hold off
    %     title(['Neuron ' num2str(keepCell(cellNum)) ' ISI'])
    
end


