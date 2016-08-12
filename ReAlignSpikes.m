function Spikes=ReAlignSpikes(CurrentWaveforms,RawTraces,SpikeTimes,UnitIDs)
for ElNum=1:size(RawTraces,1)
    %         if find(size(UnitIDs)==max(size(UnitIDs)))==1
    %             UnitIDs=UnitIDs';
    %         end
    %flip waveforms
    flipWF=0;
    if find(size(CurrentWaveforms)==max(size(CurrentWaveforms)))==2
        flipWF=1;
        CurrentWaveforms=CurrentWaveforms';
    end
    uniqIDs=unique(UnitIDs);
    if isempty(CurrentWaveforms)
    try
        % Spikes.Offline_Threshold.Units{ElNum,1}=zeros(1,numel(find(Spikes.Offline_Threshold.data{Channels(ElNum)})));
        % Spikes.Offline_Threshold.SpikeTimes{ElNum,1}=find(Spikes.Offline_Threshold.data{Channels(ElNum)});
        Spikes.Waveforms{ElNum,1}=ExtractChunks(RawTraces(ElNum,:),...
            SpikeTimes,40,'tshifted'); %'tzero' 'tmiddle' 'tshifted'
        % 0.25 bit per uV, so divide by 4 - adjust according to
        % recording system
        Spikes.Waveforms{ElNum,1}=Spikes.Waveforms{ElNum,1}./4;
        
    catch
    end
    else
        Spikes.Waveforms{ElNum,1}=CurrentWaveforms;
    end
    try
%     figure; hold on
    wfPeak=nan(length(uniqIDs),1);
    for unit=1:length(uniqIDs)
        meanWF=mean(Spikes.Waveforms{ElNum,1}(UnitIDs==uniqIDs(unit),:));
%         plot(meanWF);
        thWF=meanWF<-std(meanWF);
        minimas = find(diff(thWF) >= 1,1)+1;
        deriveWF=diff(meanWF);      %  plot(deriveWF)
        wfPeak(unit)=minimas-find(deriveWF(minimas:-1:1)>0,1)+2;
    end
    for unit=1:length(uniqIDs)
        if wfPeak(unit)>min(wfPeak)
            SpikeTimes(UnitIDs==uniqIDs(unit))=SpikeTimes(UnitIDs==uniqIDs(unit))+(wfPeak(unit)-min(wfPeak));
        end
    end
        % Spikes.Offline_Threshold.Units{ElNum,1}=zeros(1,numel(find(Spikes.Offline_Threshold.data{Channels(ElNum)})));
        % Spikes.Offline_Threshold.SpikeTimes{ElNum,1}=find(Spikes.Offline_Threshold.data{Channels(ElNum)});
        Spikes.Waveforms{ElNum,1}=ExtractChunks(RawTraces(ElNum,:),...
            SpikeTimes,40,'tshifted'); %'tzero' 'tmiddle' 'tshifted'
        % 0.25 bit per uV, so divide by 4 - adjust according to
        % recording system
        Spikes.Waveforms{ElNum,1}=Spikes.Waveforms{ElNum,1}./4;
        if flipWF==1
            Spikes.Waveforms{ElNum,1}={Spikes.Waveforms{ElNum,1}'};
        end
        Spikes.SpikeTimes={SpikeTimes};
    catch
    end

%     figure; hold on
%     for unit=1:length(unique(UnitIDs))
%         plot(mean(Spikes.Waveforms{ElNum,1}(UnitIDs==uniqIDs(unit),:)))
%     end
    
end


