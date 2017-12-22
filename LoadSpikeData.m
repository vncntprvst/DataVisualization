function Spikes=LoadSpikeData(fName,rawData,electrodes,samplingRate,bitResolution)
if logical(regexp(fName,'Ch\d+.')) %from Spike2
    load(fName)
    Spikes.Units{electrodes,1}=nw_401.codes(:,1);
    Spikes.SpikeTimes{electrodes,1}=uint32(nw_401.times*samplingRate);
    Spikes.Waveforms{electrodes,1}=nw_401.values;
    Spikes.samplingRate(electrodes,1)=samplingRate;
elseif contains(fName,'.hdf5') % Spyking Circus
    fName=regexp(fName,'\S+?(?=\.\w+\.\w+$)','match','once');
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
elseif contains(fName,'rez.mat') || contains(fName,'_KS') %Kilosort
    load(fName);
    
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
    
elseif contains(fName,'.csv') || contains(fName,'_jrc.mat') %from JRClust
    
    %% locate the _jrc file
    dirListing=dir;
    S0struct=dirListing(~cellfun('isempty',cellfun(@(x) strfind(x,'_jrc.mat'),...
        {dirListing.name},'UniformOutput',false))).name;
    
    % dimm_spk Dimensions for spike waveforms (stored in_spkwav.bin file)
    % viTime_spk Spike timing in ADC sample unit
    % cviSpk_site Cell of spike index (for _spk prefix) per site
    % miClu_log
    % P Parameter struct used for automated clustering
    % S_clu Cluster-specific information
    load(S0struct, 'dimm_spk','viTime_spk','cviSpk_site','miClu_log','P','S_clu')
    
    %% import info from cvs file export
%     clusterInfo = ImportJRClusSortInfo(fName);
    
    %% if we want to attribute each cluster to a specific electrode:
    %     allClusters=unique(clusterInfo.clusterNum);
    %     for clusNum=1:length(allClusters)
    %         bestSite=mode(clusterInfo.bestSite(clusterInfo.clusterNum==allClusters(clusNum)));
    %         clusterInfo.bestSite(clusterInfo.clusterNum==allClusters(clusNum))=bestSite;
    %     end
    
    %     Spikes.Units=clusterInfo.clusterNum;
    %     Spikes.SpikeTimes=clusterInfo.bestSite;
    
    
    %% get filtered waveforms
    vcFile=dirListing(~cellfun('isempty',cellfun(@(x) strfind(x,'_spkwav'),...
        {dirListing.name},'UniformOutput',false))).name;
    vcDataType = 'int16';
    fid=fopen(vcFile, 'r');
    % mnWav = fread_workingresize(fid, dimm, vcDataType);
    mnWav = fread(fid, prod(dimm_spk), ['*', vcDataType]);
    if numel(mnWav) == prod(dimm_spk)
        mnWav = reshape(mnWav, dimm_spk);
    else
        dimm2 = floor(numel(mnWav) / dimm_spk(1));
        if dimm2 >= 1
            mnWav = reshape(mnWav, dimm_spk(1), dimm2);
        else
            mnWav = [];
        end
    end
    if ~isempty(vcFile), fclose(fid); end
    %% degenerate. keeping largest waveforms
    %     keepSite=squeeze(prod(abs(mnWav)));[keepSite,~]=find(keepSite==max(keepSite));
    %     waveForms=nan(size(mnWav,1),size(mnWav,3));
    %     for spktTimeIdx=1:size(mnWav,3)
    %         waveForms(:,spktTimeIdx)=squeeze(mnWav(:,keepSite(spktTimeIdx),spktTimeIdx));
    %     end
    
    for elNum=1:electrodes
        try
            units=cviSpk_site{elNum}; % if data from csv file:  clusterInfo.bestSite==elNum;
            units=units(miClu_log(units,1)>=0);
            Spikes.Units{elNum,1}=miClu_log(units,1); %         clusterInfo.clusterNum(units);
            Spikes.SpikeTimes{elNum,1}=viTime_spk(units) ; %    clusterInfo.timeStamps(units)*samplingRate;
            Spikes.Waveforms{elNum,1}=squeeze(mnWav(:,1,units));
            
            %% proof that the first trace in mnWav's 2nd dimension is always from the center site:
            %             miSites_clu = P.miSites(:, S_clu.viSite_clu); % which sites correspond to mnWav's 2nd dimension
            %             rndTimeStamp=922;
            %             figure; hold on;
            %             for wfNum=1:9
            %                 plot(mnWav(:,wfNum,rndTimeStamp));
            %             end
            %             plot(mnWav(:,miSites_clu(:,miClu_log(rndTimeStamp,1))==S_clu.viSite_clu(miClu_log(rndTimeStamp,1)),rndTimeStamp),'ko')
            
            %% some more exploration
            %             mode(clusterInfo.clusterNum(units))
            %             foo=mnWav(:,:,units);
            %             figure; plot(mean(squeeze(foo(:,1,:)),2))
            %
            %             foo=mnWav(:,:,clusterInfo.clusterNum==1);
            %             subsampleIdx=round(linspace(1,24000,20));
            %             figure; hold on;
            %             for timestamps=1:20
            %                 plot(foo(:,1,subsampleIdx(timestamps)));
            %             end
            %             plot(mean(squeeze(mnWav(:,1,:)),2),'k','linewidth',1.5);
            %
            %             figure; hold on;
            %             for avwf=1:9
            %                 plot(squeeze(mnWav(:,avwf,2)));
            %             end
            %             plot(squeeze(mnWav(:,1,2)),'ko');
            %
            %             faa=Spikes.Waveforms{elNum,1};
            %             figure; hold on;
            %             for timestamps=1:20
            %                 plot(faa(timestamps,:)');
            %             end
            %             plot(mean(squeeze(mnWav(:,1,:)),2),'k','linewidth',1.5);
            
            %% alternative spike extraction
            % extract spike waveforms  rawData = memmapfile('example.dat','Format','int16');
            %             if isa(rawData,'memmapfile') % reading electrode data from .dat file
            %                 Spikes.Waveforms{elNum,1}=ExtractChunks(rawData.Data(elNum:electrodes:max(size(rawData.Data))),...
            %                     Spikes.SpikeTimes{elNum,1},50,'tshifted'); %'tzero' 'tmiddle' 'tshifted'
            %             else
            %                 Spikes.Waveforms{elNum,1}=ExtractChunks(rawData(elNum,:),...
            %                     Spikes.SpikeTimes{elNum,1},50,'tshifted'); %'tzero' 'tmiddle' 'tshifted'
            %             end
            
            %% scale to resolution
            Spikes.Waveforms{elNum,1}=Spikes.Waveforms{elNum,1}.*bitResolution;
            Spikes.samplingRate(elNum,1)=samplingRate;
        catch
            [Spikes.Units{elNum,1},Spikes.SpikeTimes{elNum,1}]=deal([]);
        end
    end
elseif contains(fName,'.mat') % just Matlab processing
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
