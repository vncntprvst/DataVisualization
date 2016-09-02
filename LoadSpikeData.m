function Spikes=LoadSpikeData(fName,rawData,electrodes,samplingRate,bitResolution)
if logical(regexp(fName,'Ch\d+.')) %from Spike2
    load(fName)
    Spikes.Units{electrodes,1}=nw_401.codes(:,1);
    Spikes.SpikeTimes{electrodes,1}=uint32(nw_401.times*samplingRate);
    Spikes.Waveforms{electrodes,1}=nw_401.values;
    Spikes.samplingRate(electrodes,1)=samplingRate;
elseif strfind(fName,'.hdf5') % Spyking Circus
    fName=regexp(fName,'\w+(?=\.\w+\.)','match','once');
    try
        templateToEl=h5read([fName '.templates.hdf5'],'/electrodes');
    catch
        templateToEl=h5read([fName '.clusters.hdf5'],'/electrodes');
    end
    for elNum=1:electrodes
        try
            %Clusters data (including non-clustered spikes)
%             Spikes.Units{elNum,1}=h5read([fName '.clusters.hdf5'],['/clusters_' num2str(elNum-1)]);
%             Spikes.SpikeTimes{elNum,1}=h5read([fName '.clusters.hdf5'],['/times_' num2str(elNum-1)]);
            
            %Results, after fitting templates
            thisElTemplates=find(templateToEl==elNum-1)-1;
            [spktimes,units]=deal(cell(size(thisElTemplates,1),1));
            for templt=1:size(thisElTemplates,1)
                 spktimes{templt}=h5read([fName '.result.hdf5'],['/spiketimes/temp_' num2str(thisElTemplates(templt))]);
                 units{templt}=ones(size(spktimes{templt},1),1)*templt;
            end
            % concatenate values
            Spikes.Units{elNum,1}=vertcat(units{:});
            Spikes.SpikeTimes{elNum,1}=vertcat(spktimes{:});
            % sort times, and adjust unit orders
            [Spikes.SpikeTimes{elNum,1},timeIdx]=sort(Spikes.SpikeTimes{elNum,1});
            Spikes.Units{elNum,1}=Spikes.Units{elNum,1}(timeIdx);
            % extract spike waveforms foo=Spikes.Waveforms{elNum,1}(Spikes.Units{elNum,1}==7,:);
            if isa(rawData,'memmapfile') % reading electrode data from .dat file
                Spikes.Waveforms{elNum,1}=ExtractChunks(rawData.Data(elNum:electrodes:max(size(rawData.Data))),...
                    Spikes.SpikeTimes{elNum,1},50,'tshifted'); %'tzero' 'tmiddle' 'tshifted'
            else
                Spikes.Waveforms{elNum,1}=ExtractChunks(rawData(elNum,:),...
                    Spikes.SpikeTimes{elNum,1},50,'tshifted'); %'tzero' 'tmiddle' 'tshifted'
            end
            % scale to resolution
            Spikes.Waveforms{elNum,1}=Spikes.Waveforms{elNum,1}.*bitResolution;
            Spikes.samplingRate(elNum,1)=samplingRate;
        catch
        end
    end
elseif strfind(fName,'.mat') %Matlab export - all units unsorted by default
    load(fName);
    for elNum=1:numel(electrodes)
        try
            Spikes.Offline_Threshold.Units{elNum,1}=zeros(1,numel(find(Spikes.Offline_Threshold.data{electrodes(elNum)})));
            Spikes.Offline_Threshold.SpikeTimes{elNum,1}=find(Spikes.Offline_Threshold.data{electrodes(elNum)});
            Spikes.Offline_Threshold.Waveforms{elNum,1}=ExtractChunks(rawData(elNum,:),...
                Spikes.Offline_Threshold.SpikeTimes{elNum,1},40,'tshifted'); %'tzero' 'tmiddle' 'tshifted'
            % 0.25 bit per uV, so divide by 4 - adjust according to
            % recording system
            Spikes.Offline_Threshold.Waveforms{elNum,1}=Spikes.Waveforms{elNum,1}./4;
        catch
        end
    end
end
