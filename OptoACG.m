function OptoACG(spikeData,TTLtimes,keepCell,pulseDur,axisHandle)

for cellNum=1:length(keepCell)
    %find(mean(meanChan)==max(mean(meanChan)));
    if ~exist('axisHandle','var') || isempty(axisHandle)
        figure('Position',[1092 149 708 761]); hold on
    end
    colormap(parula); cmap=colormap;
    %% Spike times
    spikeTimes=spikeData.times(spikeData.unitID==keepCell(cellNum),:);
    
    %get wich spike time occur during TTL
    pulseIdx=false(size(spikeTimes,1),size(TTLtimes,1));
    %     figure; hold on; plot(spikeTimes,'*'); plot(TTLtimes,'d')
    for TTLNum=1:size(TTLtimes,1)
        pulseIdx(:,TTLNum)=spikeTimes>TTLtimes(TTLNum) & spikeTimes<TTLtimes(TTLNum)+pulseDur;
    end
    onSpikes=logical(sum(pulseIdx,2));
    
    
    binSize=1/2;
    unitST_onPulse=int32(spikeTimes(onSpikes)*1000/binSize);
    unitST_offPulse=int32(spikeTimes(~onSpikes)*1000/binSize);
    % bin spikes
    if ~isempty(unitST_onPulse)
        spikeTimeIdx_onPulse=zeros(1,unitST_onPulse(end));
        spikeTimeIdx_onPulse(unitST_onPulse)=1;
        numBin=ceil(length(spikeTimeIdx_onPulse)/binSize);
        binUnits_onPulse = histcounts(double(unitST_onPulse), linspace(0,length(spikeTimeIdx_onPulse),numBin));
        binUnits_onPulse(binUnits_onPulse>1)=1; %no more than 1 spike per ms
        % compute autocorrelogram
        [ACG_onPulse,lags_onPulse]=xcorr(double(binUnits_onPulse),200,'unbiased'); %'coeff'
        ACG_onPulse(lags_onPulse==0)=0;
    end
    spikeTimeIdx_offPulse=zeros(1,unitST_offPulse(end));
    spikeTimeIdx_offPulse(unitST_offPulse)=1;
    numBin=ceil(length(spikeTimeIdx_offPulse)/binSize);
    binUnits_offPulse = histcounts(double(unitST_offPulse), linspace(0,length(spikeTimeIdx_offPulse),numBin));
    binUnits_offPulse(binUnits_offPulse>1)=1; %no more than 1 spike per ms
    % compute autocorrelogram
    [ACG_offPulse,lags_offPulse]=xcorr(double(binUnits_offPulse),200,'unbiased'); %'coeff'
    ACG_offPulse(lags_offPulse==0)=0;
    
    ACGh_offPulse=bar(lags_offPulse,ACG_offPulse,'BarWidth', 1.6);
    ACGh_offPulse.FaceColor = [1.0000    0.6784    0.0980]; %0.8941    0.7686
    ACGh_offPulse.EdgeColor = 'none';
    if ~isempty(unitST_onPulse)
        ACGh_onPulse=bar(lags_onPulse,ACG_onPulse,'BarWidth', 1.6);
        ACGh_onPulse.FaceColor = [0.3 0.75 0.93];
        ACGh_onPulse.EdgeColor = 'none';
    end
    
    axis('tight');box off; grid('on');
    xlabel('ACG (ms)');
    %     legend([ACGh_onPulse,ACGh_offPulse],{'Pulse On','Pulse Off'},'location','northeast')
    %     legend('boxoff')
    set(gca,'xlim',[-50 50]/binSize,... %'ylim',[0 max([max(get(gca,'ylim')) 10^1])]
        'xtick',-50/binSize:20:50/binSize,'xticklabel',-50:20*binSize:50);
    set(gca,'Color','white','FontSize',10,'FontName','calibri','TickDir','out');
    hold off
end