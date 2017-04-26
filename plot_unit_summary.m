%% load data file
channelNum=22;
[spikeData,fileInfo]=GetSpikeData(channelNum); %
Trials = LoadTTL([fileInfo.Filename fileInfo.FileExt]);

%% get spike times and convert to binary array
for clusNum=1:size(spikeData,2)
    %% convert to 1 millisecond bins and plot excerpt
    binSize=1;
    numBin=ceil((max(spikeData(clusNum).spikeTimes)+1)/...
        (Trials.sampleRate/1000)/binSize);
    
    [spikeCount,spikeTime]=histcounts(double(spikeData(clusNum).spikeTimes)/...
        double(Trials.sampleRate/1000), numBin);
    
    %     figure; bar(spikeTime(1:6000),spikeCount(1:6000),'hist')
    
    %% spike density function
    spikeArray = zeros(1,ceil(max(spikeTime))+1);
    spikeArray(ceil(spikeTime(1:end-1)))=spikeCount;
    sigma=1;
    convSpikeTime = [zeros(1,sigma*3) fullgauss_filtconv(spikeArray,sigma,0)].*1000;
    %     hold on
    % plot([zeros(1,sigma*3) convspikeTime zeros(1,sigma*3)])
    %     plot( convSpikeTime(1:6000-sigma*3))
    
    %% create rasters aligned to TTL
    %define parameters
    preAlignWindow=20;
    postAlignWindow=59;
    TTLtimes=Trials.start/(Trials.sampleRate/1000);
    spikeRasters=nan(numel(Trials.start),preAlignWindow+postAlignWindow+1);
    for trialNum=1:numel(Trials.start)
        try
            spikeRasters(trialNum,:)=spikeArray(...
                TTLtimes(trialNum)-preAlignWindow:...
                TTLtimes(trialNum)+postAlignWindow);
            %smoothed:
            %             spikeRasters(trialNum,:)=convSpikeTime(...
            %                 TTLtimes(trialNum)-preAlignWindow:...
            %                 TTLtimes(trialNum)+postAlignWindow);
        catch
            continue
        end
    end
    spikeRasters=spikeRasters(~isnan(sum(spikeRasters,2)),:);
    
    %% plot rasters aligned to TTL
    legends{1}= 'Time (ms)'; %xlabel
    legends{2}='Pulse #'; %ylabel
    legends{3}={['Channel ' num2str(channelNum) ...
        ', Neuron # ' num2str(clusNum)],...
        'Neural response to laser stimulation'}; % Title
    alignSpecs={5,[0.3 0.75 0.93]};
    figureStyle={'copper',...
        [1 preAlignWindow preAlignWindow+10 preAlignWindow+postAlignWindow+1],...
        [-preAlignWindow 0 10 postAlignWindow+1],...
        };
    SpikeRasters(spikeRasters,preAlignWindow,postAlignWindow,legends);
end


