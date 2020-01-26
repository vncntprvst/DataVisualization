function PhotoTagPlots(ephysData,pulses)

% %% From SpikeVisualizationGUI export
% fileName= 'vIRt41_0919_5400_FNOptoStim_10Hz_10ms_10mW_nopp_16Ch'
% channelNum=5;
%
% spikeData=load([fileName '_Ch' num2str(channelNum) '.mat'],'waveForms','spikeTimes','unitsIdx','samplingRate','selectedUnits');
% load([fileName '_Ch' num2str(channelNum) '.mat'],'TTLs');
% load([fileName '_Ch' num2str(channelNum) '.mat'], 'traceExcerpt');
% traceData=load([fileName '_Ch' num2str(channelNum) '.mat'], 'allTraces','traceInfo');
% try % if recording info file is there
%     load([fileName '_info.mat'])
% catch
%     %?
% end
%
% % read TTL dat file
% %     cd(sessionDir);
% %     TTLFileName=[regexp(recName,'\S+?(?=_export)','match','once') '_TTLs.dat'];
% %     fid = fopen(TTLFileName, 'r');
% %     TTLs = fread(fid,[2,Inf],'int32');
% %     fclose(fid);
%
% % TTL times should be sync'ed to recoding start already ...
% TTLs.start(:,1)=TTLs.start(:,1)-double(rec_info.recordingStartTime);
% TTLs.start(:,2)=TTLs.start(:,1)/double(TTLs.samplingRate{1}/1000);
% TTLs.end(:,1)=TTLs.end(:,1)-double(rec_info.recordingStartTime);
% TTLs.end(:,2)=TTLs.end(:,1)/double(TTLs.samplingRate{1}/1000);
%
% %% From JRClust csv export
% dirListing=dir; dirName=cd;
% infoFileName=dirListing(~cellfun('isempty',cellfun(@(x) strfind(x,'_info.'),...
%     {dirListing.name},'UniformOutput',false))).name;
% recInfo=load(fullfile(dirName,infoFileName),'rec_info');recInfo=recInfo.rec_info;
%
% JRclustData=load([fileName '.csv']); %_JR
% spikeData.selectedUnits=[2]; % 3 4];
% unitIdx=JRclustData(:,2)==2; % | JRclustData(:,2)==3 | JRclustData(:,2)==4;
% spikeData.unitsIdx=int8(JRclustData(unitIdx,2));
% spikeData.spikeTimes=uint32(JRclustData(unitIdx,1)*recInfo.samplingRate);
%
% traces = memmapfile(fullfile(dirName,[fileName '.dat']),'Format','int16');
% waveForms=cell(recInfo.numRecChan,1);
% for chNum=1:recInfo.numRecChan
%     waveForms{chNum,1}=ExtractChunks(traces.Data(chNum:recInfo.numRecChan:max(size(traces.Data))),...
%         JRclustData(JRclustData(JRclustData(:,3)==chNum,2)==spikeData.selectedUnits,1)*recInfo.samplingRate,...
%         50,'tshifted'); %'tzero' 'tmiddle' 'tshifted'
% end
%
% spikeData.waveForms=int16(cell2mat(waveForms(~cellfun('isempty',waveForms)))');
% % %
% % % spikes=LoadSpikeData([fileName '_JR.csv'],[],recInfo.numRecChan,recInfo.samplingRate,recInfo.bitResolution);
% % %
% % % spikeData.selectedUnits=2;
% % % units=cellfun(@(chUnits) chUnits(chUnits==spikeData.selectedUnits),spikes.Units,'UniformOutput',false);
% % % spikeTimes=cellfun(@(chUnits,chSpkTimes) chSpkTimes(chUnits==spikeData.selectedUnits),...
% % %         spikes.Units,spikes.SpikeTimes,'UniformOutput',false);
% % % waveForms=cellfun(@(chUnits,chWf) chWf(:,chUnits==spikeData.selectedUnits),spikes.Units,...
% % %         spikes.Waveforms,'UniformOutput',false);
% % % spikeTimes=cell2mat(spikeTimes(~cellfun('isempty',units)));
% % % waveForms=cell2mat(waveForms(~cellfun('isempty',units))');
% % % units=cell2mat(units(~cellfun('isempty',units)));
%
% %% Just use LoadSpikeData
% traces = memmapfile(fullfile(cd,[fileName '_traces.bin']),'Format','int16');
% spikeData=LoadSpikeData([fileName '_jrc.mat' ],traces); %electrodes,samplingRate,bitResolution;
% % if field name change needed:
% spikeData.unitsIdx=spikeData.unitID;
% spikeData.spikeTimes=spikeData.times;
% spikeData.waveForms=spikeData.waveforms';
% rmfield(spikeData,{'unitID','times','waveforms'});

% % if TTLs arent exported, do it from DataExportGUI, separate channel
% load('vIRt41_0919_5400_FNOptoStim_10Hz_10ms_10mW_nopp_trials.mat', 'trials')
% TTLs=trials;

%% variables
fileName=ephysData.recInfo.sessionName; %'vIRt44_1210_5450';
TTLs.start=pulses.TTLTimes(1,:); TTLs.end=pulses.TTLTimes(2,:);
pulseDur=min(mode(TTLs.end-TTLs.start));
IPI=mode(diff(TTLs.start));
delay=5;
preAlignWindow=50;
postAlignWindow=200;
SRR=ephysData.recInfo.SRratio;
traceExcerpt.excerptSize=1000*SRR;

% spikeData.selectedUnits=[7,18,24]-1;
if islogical(ephysData.selectedUnits) %logical array
    ephysData.selectedUnits=find(ephysData.selectedUnits);
end

%% compute rasters
spikeRasters=EphysFun.MakeRasters(ephysData.spikes.times,ephysData.spikes.unitID,...
    ephysData.spikes.samplingRate,int32(size(ephysData.traces,2)/ephysData.spikes.samplingRate*1000));
spikeRasters=spikeRasters(ephysData.selectedUnits,:);
alignedRasters=EphysFun.AlignRasters(spikeRasters,TTLs.start,preAlignWindow,postAlignWindow);

%% compute spike density functions
% spikeRate=EphysFun.MakeSDF(spikeRasters);

%% Figures
% some issue with ttl times from npy -> see CH29 from 'vIRt22_2018-10-16_18-43-54_5100_50ms1Hz5mW_nopp' KS

for cellNum=1:size(ephysData.selectedUnits,1)
    % keep one cell
    % cellNum=2;
    
    figure('Position',[296 149 1504 761],'name',...
        [fileName ' Unit' num2str(ephysData.selectedUnits(cellNum))] ); %Ch' num2str(spikeData.selectedUnits(cellNum))
    
    %% waveforms
    subplot(3,3,[1,4]); hold on
    OptoWaveforms(ephysData.spikes,TTLs.start,ephysData.selectedUnits(cellNum),delay,gca)
    
    %% rasters
    subplot(3,3,[2,5]);
    if ~iscell(alignedRasters); alignedRasters={alignedRasters}; end
    OptoRasters(alignedRasters(cellNum),preAlignWindow,pulseDur,gca);
    % title(['Channel ' num2str(channelNum) ', Neuron ' num2str(spikeData.selectedUnits(cellNum))],'FontName','Cambria');
    
    %% SDF
    subplot(3,3,[3,6])
    OptoSDF(alignedRasters(cellNum),preAlignWindow,pulseDur,IPI,gca)
    
    % %% ISI
    % subplot(3,3,4); hold on
    % OptoISI(spikeData,TTLtimes,spikeData.selectedUnits(cellNum),gca)
    %
    % %% ACG
    % subplot(3,3,7); hold on
    % OptoACG(spikeData,TTLtimes,spikeData.selectedUnits(cellNum),gca)
    
    %% raw trace
    subplot(3,3,7:9); hold on
    if ~isfield(ephysData.recInfo,'SRratio')
        SRR=double(ephysData.spikes.samplingRate/1000);
    end
    % excerptTTLtimes=double(TTLtimes(TTLtimes>(traceExcerpt.location-traceExcerpt.excerptSize)/spikeData.recInfo.SRratio &...
    %     TTLtimes<(traceExcerpt.location+traceExcerpt.excerptSize)/spikeData.recInfo.SRratio)-...
    %     (traceExcerpt.location-traceExcerpt.excerptSize)/spikeData.recInfo.SRratio)*spikeData.recInfo.SRratio;
    % if ~isempty(excerptTTLtimes)
    % %     excerptTTLtimes=excerptTTLtimes(end); %if wants to keep only one pulse
    % else % check further out in the trace
    traceExcerpt.location=TTLs.start(1)*SRR;
    %     mod(winIdxStart,traceData.traceInfo.numChan)
    if exist('traceData','var') && isa(traceData,'memmapfile')
        winIdxStart=(traceExcerpt.location-traceExcerpt.excerptSize)*traceData.traceInfo.numChan+1;
        winSize=2; %default 1 pulse
        winIdxEnd=winIdxStart+(winSize*2*traceExcerpt.excerptSize*traceData.traceInfo.numChan);
    else
        winIdxStart=(traceExcerpt.location-traceExcerpt.excerptSize); %*traceData.traceInfo.numChan+1;
        winIdxEnd=traceExcerpt.location+traceExcerpt.excerptSize;
    end
    excerptWindow=int32(winIdxStart:winIdxEnd-1)-SRR;
    %     size(excerptWindow,2)>(2*traceExcerpt.excerptSize*traceData.traceInfo.numChan)
    if exist('traceData','var') && isa(traceData,'memmapfile')
        traceExcerpt.data=traceData.allTraces.Data(excerptWindow);
        traceExcerpt.data=reshape(traceExcerpt.data,[traceData.traceInfo.numChan traceExcerpt.excerptSize*2*winSize]);
        preprocOption={'CAR','all'};
        traceExcerpt.data=PreProcData(traceExcerpt.data,30000,preprocOption);
        traceExcerpt.data=traceExcerpt.data(channelNum,:);%     figure; plot(dataExcerpt(11,:))
    else
        prefElec=double(ephysData.spikes.preferredElectrode(ismember(...
            ephysData.spikes.unitID,ephysData.selectedUnits(cellNum))));
        keepTrace=mode(prefElec);  
        %Sometimes not the best trace. Find a way plot most relevant trace
%         [traceFreq,uniqueTraces]=hist(prefElec,unique(prefElec));
%         keepTrace=uniqueTraces(4);
        traceExcerpt.data=ephysData.traces(keepTrace,excerptWindow); 
%         figure; plot(traceExcerpt.data)
%         figure; plot(ephysData.traces(keepTrace,:))
    end
    
    excerptTTLtimes=double(TTLs.start(TTLs.start>(traceExcerpt.location-...
        traceExcerpt.excerptSize)/SRR &...
        TTLs.start<(traceExcerpt.location+traceExcerpt.excerptSize)/SRR)-...
        (traceExcerpt.location-traceExcerpt.excerptSize)/...
        SRR)*SRR;
    
    excerptSpikeTimes={NaN};
%     figure;
    OptoRawTrace(traceExcerpt,excerptSpikeTimes,...
        SRR,excerptTTLtimes,pulseDur,'',gca)
end
end

