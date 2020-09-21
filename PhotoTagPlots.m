function PhotoTagPlots(ephysData,pulses)

%% variables
fileName=ephysData.recInfo.sessionName; %'vIRt44_1210_5450';
TTLs.start=pulses.TTLTimes; %(1,:); TTLs.end=pulses.TTLTimes(2,:);
pulseDur=pulses.duration; %  min(mode(TTLs.end-TTLs.start));
IPI=mode(diff(TTLs.start));
delay=0.005;
preAlignWindow=0.050;
postAlignWindow=0.20;
SRR=ephysData.recInfo.SRratio;
traceExcerpt.excerptSize=SRR;

% spikeData.selectedUnits=[7,18,24]-1;
if islogical(ephysData.selectedUnits) %logical array
    ephysData.selectedUnits=find(ephysData.selectedUnits);
end

%% compute rasters
spikeRasters=EphysFun.MakeRasters(ephysData.spikes.times,ephysData.spikes.unitID,...
    1,int32(size(ephysData.traces,2)/ephysData.spikes.samplingRate*1000)); %ephysData.spikes.samplingRate
spikeRasters=spikeRasters(ephysData.selectedUnits,:);
alignedRasters=EphysFun.AlignRasters(spikeRasters,TTLs.start,preAlignWindow,postAlignWindow,1000);

%% compute spike density functions
% spikeRate=EphysFun.MakeSDF(spikeRasters);

%% Figures
% some issue with ttl times from npy -> see CH29 from 'vIRt22_2018-10-16_18-43-54_5100_50ms1Hz5mW_nopp' KS

%if need to load ephys data:
spikeSortingDir=[ephysData.recInfo.dirName filesep 'SpikeSorting' filesep ephysData.recInfo.sessionName];
LoadSpikeData(fullfile(spikeSortingDir, [ephysData.recInfo.sessionName '_export_res.mat'])) ;

for cellNum=1:size(ephysData.selectedUnits,1)
    % keep one cell
    % cellNum=2;
    
    figure('Position',[639 154 923 822],'name',...
        [fileName ' Unit' num2str(ephysData.selectedUnits(cellNum))] ); %Ch' num2str(spikeData.selectedUnits(cellNum))
        
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
        %Sometimes not the best trace. Find a way plot most relevant trace
        prefElec=double(ephysData.spikes.preferredElectrode(ismember(...
                ephysData.spikes.unitID,ephysData.selectedUnits(cellNum))));
%         try
%             [traceFreq,uniqueTraces]=hist(prefElec,unique(prefElec));
%             keepTrace=uniqueTraces(end);
%             traceExcerpt.data=ephysData.traces(keepTrace,excerptWindow);
%         catch
            keepTrace=mode(prefElec);
            traceExcerpt.data=ephysData.traces(keepTrace,excerptWindow);
%         end
%         figure; plot(traceExcerpt.data)
%         figure; plot(ephysData.traces(keepTrace,:))
    end
    
    excerptTTLtimes=double(TTLs.start(TTLs.start>(traceExcerpt.location-...
        traceExcerpt.excerptSize)/SRR &...
        TTLs.start<(traceExcerpt.location+traceExcerpt.excerptSize)/SRR)-...
        (traceExcerpt.location-traceExcerpt.excerptSize)/...
        SRR)*SRR;
    
    excerptSpikeTimes={NaN};
%     figure; plot(traceExcerpt.data)
    OptoRawTrace(traceExcerpt,excerptSpikeTimes,...
        SRR,excerptTTLtimes,pulseDur,'',gca)
    
    %% waveforms
    
    spikesTimes=ephysData.spikes.times(ephysData.spikes.unitID==ephysData.selectedUnits(cellNum));
    waveForms=NaN(size(spikesTimes,1),50);
%         electrodesId=unique(spikes.preferredElectrode);
    waveForms=ExtractChunks(ephysData.traces(keepTrace,:),... %foo = PreProcData(foo,30000,{'bandpass',[300 3000]});
            spikesTimes*ephysData.recInfo.samplingRate,50,'tshifted'); %'tzero' 'tmiddle' 'tshifted'
    % scale to resolution
    waveForms=waveForms.*ephysData.recInfo.bitResolution;
    ephysData.spikes.waveforms(ephysData.spikes.unitID==ephysData.selectedUnits(cellNum),:)=waveForms;
    subplot(3,3,[1,4]); hold on
    OptoWaveforms(ephysData.spikes,TTLs.start,ephysData.selectedUnits(cellNum),delay,gca)
    
    %% rasters
    subplot(3,3,[2]);
    if ~iscell(alignedRasters); alignedRasters={alignedRasters}; end
    OptoRasters(alignedRasters(cellNum),preAlignWindow*1000,pulseDur,IPI,gca);
    % title(['Channel ' num2str(channelNum) ', Neuron ' num2str(spikeData.selectedUnits(cellNum))],'FontName','Cambria');
    
    %% Jitter
    subplot(3,3,[5]);
    OptoJitter(ephysData.spikes,TTLs.start,ephysData.selectedUnits(cellNum),delay,gca)
    
    %% SDF
    subplot(3,3,[3,6])
    OptoSDF(alignedRasters(cellNum),preAlignWindow*1000,pulseDur*1000,IPI*1000,gca)
    
    % %% ISI
    % subplot(3,3,4); hold on
    % OptoISI(spikeData,TTLtimes,spikeData.selectedUnits(cellNum),gca)
    %
    % %% ACG
    % subplot(3,3,7); hold on
    % OptoACG(spikeData,TTLtimes,spikeData.selectedUnits(cellNum),gca)
    

end
end

