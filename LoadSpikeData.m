function spikes=LoadSpikeData(argin_fName,traces) %electrodes,samplingRate,bitResolution
% fName='vIRt22_2018_10_16_20_36_04_5600_50ms1Hz10mW_1_1_export.result.hdf5';
%% from Spike2
if logical(regexp(argin_fName,'Ch\d+.'))
    load(argin_fName)
    spikes.Units{electrodes,1}=nw_401.codes(:,1);
    spikes.SpikeTimes{electrodes,1}=uint32(nw_401.times*samplingRate);
    spikes.Waveforms{electrodes,1}=nw_401.values;
    spikes.samplingRate(electrodes,1)=samplingRate;
    %% Spyking Circus
elseif contains(argin_fName,'.hdf5')
    fName=regexp(argin_fName,'\S+?(?=\.\w+\.\w+$)','match','once');
    postFix='';
    if isempty(fName)
            fName=regexp(argin_fName,'\S+?(?=\.\w+\-\w+\.\w+$)','match','once'); %in case loading merged files
            if ~isempty(fName)
                postFix='-merged';
            end
    end
    % find templates and preferred electrodes
    templateToEl=h5read([fName '.clusters' postFix '.hdf5'],'/electrodes'); % this are the *preferred* electrodes for all K templates
    numTemplates=length(templateToEl); % template has equivalent meaning to cluster
    
    % get spike times, amplitudes
    resultFile = [fName '.result' postFix '.hdf5'];
    for templateNum=1:numTemplates
        spikeTimes{templateNum,1}=double(h5read(resultFile, ['/spiketimes/temp_' num2str(templateNum-1)]));
        spikeAmplitudes{templateNum,1}=double(h5read(resultFile, ['/amplitudes/temp_' ...
            num2str(templateNum-1)])); %
        spikeAmplitudes{templateNum,1}=spikeAmplitudes{templateNum,1}(1,:)';
        templatePrefElectrode{templateNum,1}=ones(size(spikeTimes{templateNum,1},1),1)*double(templateToEl(templateNum));
        unitID{templateNum,1}=ones(size(spikeTimes{templateNum,1},1),1)*templateNum;
    end
    
    % collect non-fitted ("garbage") spikes, with unit ID 0. Those are listed by electrode
    [spikeTimes{templateNum+1},templatePrefElectrode{templateNum+1}]=deal([]);
    for electrodeNum=unique(templateToEl)'
        try
            gbSpikeTimes=h5read([fName '.result' postFix '.hdf5'],['/gspikes/elec_' num2str(electrodeNum)]);
            spikeTimes{templateNum+1}=[spikeTimes{templateNum+1};gbSpikeTimes];
            templatePrefElectrode{templateNum+1}=[templatePrefElectrode{templateNum+1};...
                ones(size(gbSpikeTimes,1),1)*double(electrodeNum)];
        catch
            % no "garbage" spikes
        end
    end
    unitID{templateNum+1}=zeros(size(spikeTimes{templateNum+1},1),1);
%     numTemplates=size(spikeTimes,1);
    % concatenate values
    spikes.unitID=uint32(vertcat(unitID{:}));
    spikes.times=vertcat(spikeTimes{:});
    spikes.amplitude=[vertcat(spikeAmplitudes{:});zeros(size(spikeTimes{end},1),1)];
    spikes.preferredElectrode=uint32(vertcat(templatePrefElectrode{:}));
    % sort times, and adjust unit orders
    [spikes.times,timeIdx]=sort(spikes.times);
    spikes.unitID=spikes.unitID(timeIdx);
    spikes.amplitude=spikes.amplitude(timeIdx);
    spikes.preferredElectrode=spikes.preferredElectrode(timeIdx);
    
    % extract spike waveforms by electrode
    %     traces=load(['../' fName '.mat']);
%     traces = memmapfile(['../' fName '.dat'],'Format','int16');
    % gert number of electrodes
    clustersData=h5info([fName '.clusters' postFix '.hdf5']);
    clustersDatasetsNames={clustersData.Datasets.Name};
    electrodesId=clustersDatasetsNames(cellfun(@(x) contains(x,'data'),...
        clustersDatasetsNames));
    electrodesId=cellfun(@(x) str2double(regexp(x,'(?<=data_)\w+','match','once')),...
        electrodesId);
%         unique(spikes.preferredElectrode(spikes.unitID==templateNum))
    
    if exist('traces','var')
        spikes.waveforms=NaN(size(spikes.times,1),50);
        for electrodeNum=electrodesId
%             =templateToEl(templateNum)+1;
            if isa(traces,'memmapfile') % reading electrode data from .dat file
                spikes.waveforms(spikes.preferredElectrode==electrodeNum,:)=...
                ExtractChunks(traces.Data(electrodeNum+1:numel(electrodesId):max(size(traces.Data))),...
                    spikes.times(spikes.preferredElectrode==electrodeNum),50,'tshifted'); %'tzero' 'tmiddle' 'tshifted'   
            else              
                spikes.waveforms(spikes.preferredElectrode==electrodeNum,:)=...
                    ExtractChunks(traces(electrodeNum+1,:),...
                    spikes.times(spikes.preferredElectrode==electrodeNum),50,'tshifted'); %'tzero' 'tmiddle' 'tshifted'
            end
            % scale to resolution
%             spikes.waveforms{elNum,1}=spikes.Waveforms{elNum,1}.*bitResolution;
        end
    else
        spikes.waveforms=[];
    end
    
    %     spikes.samplingRate=samplingRate;
    
    % plots
    %             foo=traces.Data(elNum:electrodes:max(size(traces.Data)));
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
    
    % Compute ISI
    % isis = diff(spikeTimes{templateNum,1}); hold on
%     isis = double(diff(spikes.spikeTimes(spikes.unitID==templateNum)));
%     hist(isis)
    
    % Display the amplitude
%     figure
    % plot(spikeTimes{templateNum,1}, spikeAmplitudes{templateNum,1}, '.')
%     plot(spikes.spikeTimes(spikes.unitID==templateNum),spikes.amplitude(spikes.unitID==templateNum), '.')
elseif contains(argin_fName,'rez.mat') || contains(argin_fName,'_KS') %Kilosort
        load(argin_fName);
    %        
        spikes.times=uint64(rez.st3(:,1));
        spikes.unitID=uint32(rez.st3(:,2));
        spikes.waveforms=[];
        spikes.samplingRate=30000; 
        unitIDs=unique(spikes.unitID);
        templates=abs(rez.Wraw);
        templateToEl=zeros(numel(unitIDs),1);
        for templNum=1:numel(unitIDs)
            thatTemplate=squeeze(templates(:,:,unitIDs(templNum)));
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
        spikes.preferredElectrode=nan(numel(spikes.times),1);
        for unitNum=1:numel(unitIDs)
            spikes.preferredElectrode(unitIDs(unitNum)==spikes.unitID)=templateToEl(unitNum)
        end      
        
        [spikes.times,timeIdx]=sort(spikes.times);
        spikes.unitID=spikes.unitID(timeIdx);
        spikes.preferredElectrode=spikes.preferredElectrode(timeIdx);
%         spikes.waveforms=spikes.waveforms(timeIdx,:);       

%         for elNum=1:electrodes
%             try
%                 %Results, after fitting templates
%                 thisElTemplates=find(templateToEl==elNum);
%                 units=false(size(spikeTemplates,1),1);
%                 for templt=1:size(thisElTemplates,1)
%                     units=units | spikeTemplates==thisElTemplates(templt);
%                 end
%                 Spikes.Units{elNum,1}=spikeTemplates(units);
%                 Spikes.SpikeTimes{elNum,1}=spikeTimes(units);
%                 % extract spike waveforms  traces = memmapfile('example.dat','Format','int16');
%                 if isa(traces,'memmapfile') % reading electrode data from .dat file
%                     Spikes.Waveforms{elNum,1}=ExtractChunks(traces.Data(elNum:electrodes:max(size(traces.Data))),...
%                         Spikes.SpikeTimes{elNum,1},50,'tshifted'); %'tzero' 'tmiddle' 'tshifted'
%                 else
%                     Spikes.Waveforms{elNum,1}=ExtractChunks(traces(elNum,:),...
%                         Spikes.SpikeTimes{elNum,1},50,'tshifted'); %'tzero' 'tmiddle' 'tshifted'
%                 end
%                 % scale to resolution
%                 Spikes.Waveforms{elNum,1}=Spikes.Waveforms{elNum,1}.*bitResolution;
%                 Spikes.samplingRate(elNum,1)=samplingRate;
%             catch
%             end
%         end
    
elseif contains(argin_fName,'.csv') || ... 
        contains(argin_fName,'_jrc') || ... 
        contains(argin_fName,'_res') %from JRClust%from JRClust
    %
    %     %% locate the _jrc file
        
%         S0struct=dirListing(~cellfun('isempty',cellfun(@(x) strfind(x,'_jrc.mat'),...
%             {dirListing.name},'UniformOutput',false))).name;
    
    %     % dimm_spk Dimensions for spike waveforms (stored in_spkwav.bin file)
    %     % viTime_spk Spike timing in ADC sample unit
    %     % cviSpk_site Cell of spike index (for _spk prefix) per site
    %     % miClu_log
    %     % P Parameter struct used for automated clustering
    %     % S_clu Cluster-specific information
  
            % from KiloSort spikes.times=readNPY('spike_times.npy');
      
            try % JRC v3 and v4:


                load(argin_fName,'spikeTimes','spikeSites','spikeClusters','filtShape')
%                 evtWindow = [-0.25, 0.75]; %evtWindowRaw = [-0.5, 1.5]; nSiteDir = 4;
%                 waveformsFid=fopen('vIRt32_2019_04_24_16_48_53_5185_1_1_export_filt.jrc');
%                 waveforms=fread(waveformsFid,...
%                     [sum(abs(evtWindow))*30,nSiteDir,size(spikeClusters,1)],'int16');
%                 fclose(waveformsFid);
%                 figure; plot(mean(waveforms(1:4:120,spikeClusters==4)'))
                spikes.unitID=spikeClusters;
                spikes.times=spikeTimes;
                spikes.preferredElectrode=spikeSites; %Site with the peak spike amplitude %cviSpk_site Cell of the spike indices per site
                spikes.templatesIdx=[];
                spikes.templates=[];
                spikes.waveforms=[];                        
                spikes.bitResolution=[];
                spikes.samplingRate=[];
            catch
                try
                    % v2 updated structure:
                    load(argin_fName,'miClu_log','P','S_clu','dimm_spk',...
                        'viSite_spk','viTime_spk');%'cviSpk_site'
                    
                    spikes.unitID=S_clu.viClu;
                    spikes.times=viTime_spk;
                    spikes.preferredElectrode=viSite_spk; %Site with the peak spike amplitude %cviSpk_site Cell of the spike indices per site
                    spikes.templatesIdx=S_clu.viSite_clu;
                    spikes.templates=S_clu.tmrWav_spk_clu;
                    spikes.waveforms=[];
                    spikes.bitResolution=P.uV_per_bit;
                    spikes.samplingRate=P.sRateHz;
                    
                catch
                    % old structure
                    load(argin_fName,'S_clu','spikeTimes','spikeSites','P');
                    
                    spikes.unitID=S_clu.spikeClusters;
                    spikes.times=spikeTimes;
                    spikes.preferredElectrode=spikeSites;
%                     spikes.templatesIdx=S_clu.clusterTemplates;
%                     spikes.templates=S_clu.spikeTemplates;
                    spikes.waveforms=S_clu.tmrWav_spk_clu; %mean waveform
                    spikes.bitResolution=P.uV_per_bit;
                    spikes.samplingRate=P.sampleRateHz;
                    
                    % get filtered waveforms
                    dirListing=dir;
                    spikeWaveFormsFile=cellfun(@(x) strfind(x,'_spkwav'),...
                        {dirListing.name},'UniformOutput',false);
                    if ~isempty(vertcat(spikeWaveFormsFile{:}))
                        vcFile=dirListing(~cellfun('isempty',cellfun(@(x) strfind(x,'_spkwav'),...
                            {dirListing.name},'UniformOutput',false))).name;
                        vcDataType = 'int16';
                        fid=fopen(vcFile, 'r');
                        % mnWav = fread_workingresize(fid, dimm, vcDataType);
                        mnWav = fread(fid, prod(dimm_spk), ['*', vcDataType]); %(nSamples_spk x nSites_spk x nSpikes: int16)
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
                        spikes.waveforms=mnWav;
                        spikes.waveforms=permute(spikes.waveforms,[3 1 2]);
                        spikes.waveforms=squeeze(spikes.waveforms(:,:,1)); %keep best waveform only
                        if ~isempty(vcFile), fclose(fid); end
                    end
                end
            end
            %extract spike waveform
%             filtWFfile=[regexp(argin_fName,'\w+(?=_res)','match','once') '_filt.jrc'];
%             if exist('filtShape','var') & exist(fullfile(cd,filtWFfile),'file')
%                 fid = fopen(filtWFfile, 'r');
%                 spikes.waveforms= reshape(fread(fid, inf, '*int16'), filtShape);
%                 fclose(fid);
%                 spikes.waveforms = permute(spikes.waveforms,[3 1 2]);
% %                 unitIDs=unique(spikes.unitID);
% %                 for unitNum=1:numel(unitIDs)
% %                     numel(unique(spikeSites(spikes.unitID==unitIDs(unitNum))))
% %                 end
%             else
            if (isempty(spikes.waveforms) || size(spikes.waveforms,1) <  size(spikes.unitID,1))...
                    && exist('traces','var')
                spikes.waveforms=NaN(size(spikes.times,1),50);
                electrodesId=unique(spikes.preferredElectrode);
                for electrodeNum=1:numel(electrodesId)
                    if isa(traces,'memmapfile') % reading electrode data from .dat file
                        spikes.waveforms(spikes.preferredElectrode==electrodeNum,:)=...
                            ExtractChunks(traces.Data(electrodeNum:numel(electrodesId):max(size(traces.Data))),...
                            spikes.times(spikes.preferredElectrode==electrodeNum),50,'tshifted'); %'tzero' 'tmiddle' 'tshifted'
                    else
                        spikes.waveforms(spikes.preferredElectrode==electrodeNum,:)=...
                            ExtractChunks(traces(electrodeNum,:),...
                            spikes.times(spikes.preferredElectrode==electrodeNum),50,'tshifted'); %'tzero' 'tmiddle' 'tshifted'
                    end
                    % scale to resolution
                    %             spikes.waveforms{elNum,1}=spikes.Waveforms{elNum,1}.*bitResolution;
                end
            end
% figure; hold on
% plot(mean(spikes.waveforms(spikeClusters==4,:))/bitResolution); % bitResolution=0.25;
% plot(mean(waveforms(spikeClusters==4,:,1)))
% refCh=mode(spikes.preferredElectrode(spikeClusters==4));
% spikeTimes=spikes.times(spikeClusters==4);
% spikeSites=spikes.preferredElectrode(spikeClusters==4);
% unitWF=waveforms(spikeClusters==4,:,:);
% unitWF_t=spikes.waveforms(spikeClusters==4,:);
% figure; hold on;
% for chNum=1:16
%     plot(traces(chNum,1:6000)+(chNum-1)*max(max(traces(:,1:6000)))*2,'k')
% end
% for spikeNum=1:15
%     for spkchNum=1:9
%         plot(spikeTimes(spikeNum)-10:spikeTimes(spikeNum)+21,unitWF(spikeNum,:,spkchNum)+int16((refCh-1)*max(max(traces(:,1:6000)))*2))
%     end
% %     plot(spikeTimes(spikeNum)-10:spikeTimes(spikeNum)+21,unitWF_t(spikeNum,:)+double((refCh-1)*max(max(traces(:,1:6000)))*2),'b')
% end
 %
    %     %% import info from cvs file export
    %     %     clusterInfo = ImportJRClusSortInfo(fName);
    %
    %     %% if we want to attribute each cluster to a specific electrode:
    %     %     allClusters=unique(clusterInfo.clusterNum);
    %     %     for clusNum=1:length(allClusters)
    %     %         bestSite=mode(clusterInfo.bestSite(clusterInfo.clusterNum==allClusters(clusNum)));
    %     %         clusterInfo.bestSite(clusterInfo.clusterNum==allClusters(clusNum))=bestSite;
    %     %     end
    %
    %     %     Spikes.Units=clusterInfo.clusterNum;
    %     %     Spikes.SpikeTimes=clusterInfo.bestSite;
    %
    %
    %     
    %     %% degenerate. keeping largest waveforms
    %     %     keepSite=squeeze(prod(abs(mnWav)));[keepSite,~]=find(keepSite==max(keepSite));
    %     %     waveForms=nan(size(mnWav,1),size(mnWav,3));
    %     %     for spktTimeIdx=1:size(mnWav,3)
    %     %         waveForms(:,spktTimeIdx)=squeeze(mnWav(:,keepSite(spktTimeIdx),spktTimeIdx));
    %     %     end
    %
    %     for elNum=1:electrodes
    %         try
    %             units=cviSpk_site{elNum}; % if data from csv file:  clusterInfo.bestSite==elNum;
    %             units=units(miClu_log(units,1)>=0);
    %             Spikes.Units{elNum,1}=miClu_log(units,1); %         clusterInfo.clusterNum(units);
    %             Spikes.SpikeTimes{elNum,1}=viTime_spk(units) ; %    clusterInfo.timeStamps(units)*samplingRate;
    %             Spikes.Waveforms{elNum,1}=squeeze(mnWav(:,1,units));
    %
    %             %% proof that the first trace in mnWav's 2nd dimension is always from the center site:
    %             %             miSites_clu = P.miSites(:, S_clu.viSite_clu); % which sites correspond to mnWav's 2nd dimension
    %             %             rndTimeStamp=922;
    %             %             figure; hold on;
    %             %             for wfNum=1:9
    %             %                 plot(mnWav(:,wfNum,rndTimeStamp));
    %             %             end
    %             %             plot(mnWav(:,miSites_clu(:,miClu_log(rndTimeStamp,1))==S_clu.viSite_clu(miClu_log(rndTimeStamp,1)),rndTimeStamp),'ko')
    %
    %             %% some more exploration
    %             %             mode(clusterInfo.clusterNum(units))
    %             %             foo=mnWav(:,:,units);
    %             %             figure; plot(mean(squeeze(foo(:,1,:)),2))
    %             %
    %             %             foo=mnWav(:,:,clusterInfo.clusterNum==1);
    %             %             subsampleIdx=round(linspace(1,24000,20));
    %             %             figure; hold on;
    %             %             for timestamps=1:20
    %             %                 plot(foo(:,1,subsampleIdx(timestamps)));
    %             %             end
    %             %             plot(mean(squeeze(mnWav(:,1,:)),2),'k','linewidth',1.5);
    %             %
    %             %             figure; hold on;
    %             %             for avwf=1:9
    %             %                 plot(squeeze(mnWav(:,avwf,2)));
    %             %             end
    %             %             plot(squeeze(mnWav(:,1,2)),'ko');
    %             %
    %             %             faa=Spikes.Waveforms{elNum,1};
    %             %             figure; hold on;
    %             %             for timestamps=1:20
    %             %                 plot(faa(timestamps,:)');
    %             %             end
    %             %             plot(mean(squeeze(mnWav(:,1,:)),2),'k','linewidth',1.5);
    %
    %             %% alternative spike extraction
    %             % extract spike waveforms  traces = memmapfile('example.dat','Format','int16');
    %             %             if isa(traces,'memmapfile') % reading electrode data from .dat file
    %             %                 Spikes.Waveforms{elNum,1}=ExtractChunks(traces.Data(elNum:electrodes:max(size(traces.Data))),...
    %             %                     Spikes.SpikeTimes{elNum,1},50,'tshifted'); %'tzero' 'tmiddle' 'tshifted'
    %             %             else
    %             %                 Spikes.Waveforms{elNum,1}=ExtractChunks(traces(elNum,:),...
    %             %                     Spikes.SpikeTimes{elNum,1},50,'tshifted'); %'tzero' 'tmiddle' 'tshifted'
    %             %             end
    %
    %             %% scale to resolution
    %             Spikes.Waveforms{elNum,1}=Spikes.Waveforms{elNum,1}.*bitResolution;
    %             Spikes.samplingRate(elNum,1)=samplingRate;
    %         catch
    %             [Spikes.Units{elNum,1},Spikes.SpikeTimes{elNum,1}]=deal([]);
    %         end
    %     end
elseif contains(argin_fName,'.mat') % Matlab processing / export
    spikes=load(argin_fName);
    if isfield(spikes,'metadata') %was exported from OE npy files
        spikes.times=spikes.spikeTimes;
        spikes.waveforms=spikes.waveForms;
        spikes.samplingRate=30000;
        spikes.unitID=spikes.clusters;
        spikes.preferredElectrode=spikes.electrodes;
        spikes = rmfield(spikes,{'spikeTimes','waveForms','clusters','electrodes','clusters','metadata'});
    else
        numUnits=numel(spikes.Offline_Sorting.Units);
        spikes.unitID=vertcat(spikes.Offline_Sorting.Units{:});
        unitIds=unique(spikes.unitID);
        spikes.preferredElectrode=...
            cellfun(@(x,y) ones(numel(x),1)*y, spikes.Offline_Sorting.Units,...
            mat2cell([1:numUnits]',ones(numUnits,1)),'UniformOutput',false);
        spikes.preferredElectrode=vertcat(spikes.preferredElectrode{:});
        for unitNUm=1:numel(unitIds)
            unitIdx=spikes.unitID==unitIds(unitNUm);
            spikes.preferredElectrode(unitIdx)=mode(spikes.preferredElectrode(unitIdx));
        end
        spikes.times=vertcat(spikes.Offline_Sorting.SpikeTimes{:});
        spikes.waveforms=vertcat(spikes.Offline_Sorting.Waveforms{:});
        spikes.samplingRate=spikes.Offline_Sorting.samplingRate;
        [spikes.times,timeIdx]=sort(spikes.times);
        spikes.unitID=spikes.unitID(timeIdx);
        spikes.waveforms=spikes.waveforms(timeIdx,:);
        spikes.preferredElectrode=spikes.preferredElectrode(timeIdx,:);
    end
    %     %Matlab export - all units unsorted by default
    %     for elNum=1:numel(electrodes)
    %         try
    %             Spikes.Units{elNum,1}=zeros(1,numel(find(Spikes.data{electrodes(elNum)})));
    %             Spikes.SpikeTimes{elNum,1}=find(Spikes.data{electrodes(elNum)});
    %             Spikes.Waveforms{elNum,1}=ExtractChunks(traces(elNum,:),...
    %                 Spikes.SpikeTimes{elNum,1},40,'tshifted'); %'tzero' 'tmiddle' 'tshifted'
    %             % 0.25 bit per uV, so divide by 4 - adjust according to
    %             % recording system
    %             Spikes.Waveforms{elNum,1}=Spikes.Waveforms{elNum,1}./4;
    %         catch
    %         end
    %     end
    %
end