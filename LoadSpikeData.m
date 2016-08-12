function Spikes=LoadSpikeData(fName,Channels,samplingRate,rawData)
if logical(regexp(fName,'Ch\d+.')) %from Spike2
    load(fName)
    Spikes.Units{Channels,1}=nw_401.codes(:,1);
    Spikes.SpikeTimes{Channels,1}=uint32(nw_401.times*samplingRate);
    Spikes.Waveforms{Channels,1}=nw_401.values;
    Spikes.samplingRate(Channels,1)=samplingRate;
elseif strfind(fName,'.hdf5') %Intan / Open-Ephys
    fName=regexp(fName,'\w+(?=\.\w+\.)','match','once');
    for chNum=1:Channels
        try
            Spikes.Units{chNum,1}=h5read([fName '.clusters.hdf5'],['/clusters_' num2str(chNum-1)]);
            Spikes.SpikeTimes{chNum,1}=h5read([fName '.clusters.hdf5'],['/times_' num2str(chNum-1)]);
            Spikes.Waveforms{chNum,1}=ExtractChunks(rawData(chNum,:),...
                Spikes.SpikeTimes{chNum,1},40,'tshifted'); %'tzero' 'tmiddle' 'tshifted'
            % 0.25 bit per uV, so divide by 4
            Spikes.Waveforms{chNum,1}=Spikes.Waveforms{chNum,1}./4;
            Spikes.samplingRate(chNum,1)=samplingRate;
        catch
        end
    end
elseif strfind(fName,'.mat') %Matlab export - all units unsorted by default
    load(fName);
    for chNum=1:numel(Channels)
        try
            Spikes.Offline_Threshold.Units{chNum,1}=zeros(1,numel(find(Spikes.Offline_Threshold.data{Channels(chNum)})));
            Spikes.Offline_Threshold.SpikeTimes{chNum,1}=find(Spikes.Offline_Threshold.data{Channels(chNum)});
            Spikes.Offline_Threshold.Waveforms{chNum,1}=ExtractChunks(rawData(chNum,:),...
                Spikes.Offline_Threshold.SpikeTimes{chNum,1},40,'tshifted'); %'tzero' 'tmiddle' 'tshifted'
            % 0.25 bit per uV, so divide by 4 - adjust according to
            % recording system
            Spikes.Offline_Threshold.Waveforms{chNum,1}=Spikes.Waveforms{chNum,1}./4;
        catch
        end
    end
end
