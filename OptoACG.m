function OptoACG(spikeData,TTLtimes,keepCell,axisHandle)

for cellNum=1:length(keepCell)
    %find(mean(meanChan)==max(mean(meanChan)));
    if ~exist('axisHandle','var') || isempty(axisHandle)
        figure('Position',[1092 149 708 761]); hold on
    end
    colormap(parula); cmap=colormap;
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
    % bin spikes
    spikeTimeIdx_onPulse=zeros(1,unitST_onPulse(end));
    spikeTimeIdx_onPulse(unitST_onPulse)=1;
    binSize=1;
    numBin=ceil(length(spikeTimeIdx_onPulse)/binSize);
    binUnits_onPulse = histcounts(double(unitST_onPulse), linspace(0,length(spikeTimeIdx_onPulse),numBin));
    binUnits_onPulse(binUnits_onPulse>1)=1; %no more than 1 spike per ms
    
    spikeTimeIdx_offPulse=zeros(1,unitST_offPulse(end));
    spikeTimeIdx_offPulse(unitST_offPulse)=1;
    numBin=ceil(length(spikeTimeIdx_offPulse)/binSize);
    binUnits_offPulse = histcounts(double(unitST_offPulse), linspace(0,length(spikeTimeIdx_offPulse),numBin));
    binUnits_offPulse(binUnits_offPulse>1)=1; %no more than 1 spike per ms
    
    % compute autocorrelogram
    [ACG_onPulse,lags_onPulse]=xcorr(double(binUnits_onPulse),'coeff'); %'coeff'
    ACG_onPulse(lags_onPulse==0)=0;
    [ACG_offPulse,lags_offPulse]=xcorr(double(binUnits_offPulse),'coeff'); %'coeff'
    ACG_offPulse(lags_offPulse==0)=0;
    
    ACGh_onPulse=bar(lags_onPulse,ACG_onPulse);
    ACGh_offPulse=bar(lags_offPulse,ACG_offPulse);
    ACGh_onPulse.FaceColor = [0.3 0.75 0.93];
    ACGh_onPulse.EdgeColor = 'none';
    ACGh_offPulse.FaceColor = cmap(cellNum,:);
    ACGh_offPulse.EdgeColor = 'none';
    axis('tight');box off;
    xlabel('Autocorrelogram (5 ms bins)');
    legend([ACGh_onPulse,ACGh_offPulse],{'Pulse On','Pulse Off'},'location','northeast')
    set(gca,'xlim',[-550 550],'Color','white','FontSize',10,'FontName','calibri','TickDir','out');
    set(gca,'Color','white','FontSize',12,'FontName','calibri');
    
end


