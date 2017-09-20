function Spikes=LoadSpikeData(fName,rawData,electrodes,samplingRate,bitResolution)
if logical(regexp(fName,'Ch\d+.')) %from Spike2
    load(fName)
    Spikes.Units{electrodes,1}=nw_401.codes(:,1);
    Spikes.SpikeTimes{electrodes,1}=uint32(nw_401.times*samplingRate);
    Spikes.Waveforms{electrodes,1}=nw_401.values;
    Spikes.samplingRate(electrodes,1)=samplingRate;
elseif strfind(fName,'.hdf5') % Spyking Circus
    fName=regexp(fName,'\w+(?=\.\w+\.)','match','once');
    templateToEl=h5read([fName '.clusters.hdf5'],'/electrodes'); % this are the *prefered* electrodes for all K templates
    for elNum=1:electrodes
        try
            %Clusters data (including non-clustered spikes)
            %             Spikes.Units{elNum,1}=h5read([fName '.clusters.hdf5'],['/clusters_' num2str(elNum-1)]);
            %             Spikes.SpikeTimes{elNum,1}=h5read([fName '.clusters.hdf5'],['/times_' num2str(elNum-1)]);
            
            %Results, after fitting templates
            thisElTemplates=find(templateToEl==elNum-1)-1;
            [spktimes,units]=deal(cell(size(thisElTemplates,1)+1,1));
            for templt=1:size(thisElTemplates,1)
                spktimes{templt}=h5read([fName '.result.hdf5'],['/spiketimes/temp_' num2str(thisElTemplates(templt))]);
                units{templt}=ones(size(spktimes{templt},1),1)*templt;               
            end
            % collect non-fitted ("garbage") spikes, with unit ID 0
            try
            spktimes{templt+1}=h5read([fName '.result.hdf5'],['/gspikes/elec_' num2str(elNum-1)]);
            units{templt+1}=zeros(size(spktimes{templt+1},1),1);
            catch 
                % no "garbage" spikes
            end
            % concatenate values
            Spikes.Units{elNum,1}=vertcat(units{:});
            Spikes.SpikeTimes{elNum,1}=vertcat(spktimes{:});
            % sort times, and adjust unit orders
            [Spikes.SpikeTimes{elNum,1},timeIdx]=sort(Spikes.SpikeTimes{elNum,1});
            Spikes.Units{elNum,1}=Spikes.Units{elNum,1}(timeIdx);
            % extract spike waveforms 
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
            
            % plots
%             foo=rawData.Data(elNum:electrodes:max(size(rawData.Data)));
%             figure; hold on
%             plot(foo(round(size(foo,1)/2)-samplingRate:round(size(foo,1)/2)+samplingRate));
%             axis('tight');box off;
%             text(100,100,num2str(round(size(foo,1)/2)))
%             text(100,50,'PrV 77 115 El 11');
%             allunits= Spikes.Units{elNum,1};
%             allspktimes=Spikes.SpikeTimes{elNum,1};
%             spkTimes=allspktimes(allspktimes>=round(size(foo,1)/2)-samplingRate &...
%                 allspktimes<round(size(foo,1)/2)+samplingRate & allunits==1);
%             rasterHeight=ones(1,size(spkTimes,2))*(min(get(gca,'ylim'))/4*3);
%             plot(spkTimes-(round(size(foo,1)/2)-samplingRate),...
%                 rasterHeight,'Color','r',...
%                 'linestyle','none','Marker','^');
            
        catch
        end
    end
elseif strfind(fName,'.mat')
    load(fName);
    if strfind(fName,'rez') || strfind(fName,'_KS') %Kilosort
        spikeTimes = uint64(rez.st3(:,1));
        spikeTemplates = uint32(rez.st3(:,2));
        templates=abs(rez.Wraw);
        templateToEl=zeros(max(unique(spikeTemplates)),1);
        for templNum=1:max(unique(spikeTemplates))
            thatTemplate=squeeze(templates(:,:,templNum));
            [elecRow,~] = ind2sub(size(thatTemplate),find(thatTemplate==max(max(thatTemplate))));
            if size(elecRow,1)>1
                if length(unique(elecRow))>1 %weird
                    %                     then look for next biggest value?
                    return
                else
                    elecRow=unique(elecRow);
                end
            end
            templateToEl(templNum)=elecRow;
        end
        for elNum=1:electrodes
            try
                %Results, after fitting templates
                thisElTemplates=find(templateToEl==elNum);
                units=false(size(spikeTemplates,1),1);
                for templt=1:size(thisElTemplates,1)
                    units=units | spikeTemplates==thisElTemplates(templt);
                end
                Spikes.Units{elNum,1}=spikeTemplates(units);
                Spikes.SpikeTimes{elNum,1}=spikeTimes(units);
                % extract spike waveforms  rawData = memmapfile('example.dat','Format','int16');
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
        
    else
        %Matlab export - all units unsorted by default
        
        
        for elNum=1:numel(electrodes)
            try
                Spikes.Units{elNum,1}=zeros(1,numel(find(Spikes.data{electrodes(elNum)})));
                Spikes.SpikeTimes{elNum,1}=find(Spikes.data{electrodes(elNum)});
                Spikes.Waveforms{elNum,1}=ExtractChunks(rawData(elNum,:),...
                    Spikes.SpikeTimes{elNum,1},40,'tshifted'); %'tzero' 'tmiddle' 'tshifted'
                % 0.25 bit per uV, so divide by 4 - adjust according to
                % recording system
                Spikes.Waveforms{elNum,1}=Spikes.Waveforms{elNum,1}./4;
            catch
            end
        end
    end
end
