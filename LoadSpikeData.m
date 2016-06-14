function Spikes=LoadSpikeData(fName,Channels,samplingRate)
if logical(regexp(fName,'Ch\d+.'))
    load(fName)
    Spikes.Offline_Sorting.Units{Channels,1}=nw_401.codes(:,1);
    Spikes.Offline_Sorting.SpikeTimes{Channels,1}=uint32(nw_401.times*samplingRate);
    Spikes.Offline_Sorting.Waveforms{Channels,1}=nw_401.values;
    Spikes.Offline_Sorting.samplingRate(Channels,1)=samplingRate;
elseif strfind(fName,'.hdf5')
    fName=regexp(fName,'\w+(?=\.\w+\.)','match','once');
    for chNum=1:Channels
        try
            Spikes.Offline_SpkSort.Units{chNum,1}=h5read([fName '.clusters.hdf5'],['/clusters_' num2str(chNum-1)]);
            Spikes.Offline_SpkSort.SpikeTimes{3,1}=h5read([fName '.clusters.hdf5'],['/times_' num2str(chNum-1)]);
%             Spikes.Offline_SpkSort.Waveforms{3,1}=h5read([fName '.clusters.hdf5'],['/data_' num2str(chNum-1)]);
            %     Spikes.Offline_SpkSort.Waveforms=h5read([fName '.templates.hdf5'],'/temp_data');
            %     Spikes.Offline_SpkSort.templates{10,1}.spiketimes=h5read([fName '.result.hdf5'],'/spiketimes/temp_10');
            %     Spikes.Offline_SpkSort.templates{10,1}.amplitudes=h5read([fName '.result.hdf5'],'/amplitudes/temp_10');
        catch
        end
    end
end
