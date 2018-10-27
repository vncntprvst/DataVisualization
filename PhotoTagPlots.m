fileName='vIRt22_2018-10-16_19-11-34_5100_50ms1Hz10mW_nopp' %_Ch29
% 'vIRt22_2018-10-16_18-43-54_5100_50ms1Hz5mW_nopp' %_Ch29
%'SpVi16_0403_WR_4850_LS1Hz2ms100mW_nopp' % _Ch7'% 'SpVi12_1107_WR_Texture_LS500mH_24Ch_nopp' %_Ch6
channelNum=29; %11 %24;
%% 'SpVi12_133_2Hz2ms_7mW_nopp'
%'SpVi12_1206_WR_LS_500mHz_2ms_2_nopp'
% 'SpVi12_1107_WR_Texture_LS500mH_24Ch_nopp'
%'SpVi12_1107_WR_Texture_LS500mH_24Ch_nopp'
%'SpVi12_133_2Hz2ms_7mW_nopp'
%'SpVi12_1103_WR_LS_500mHz2ms70m_nopp'
%'SpVi12_1107_WR_MS_LS500mHz2ms6_nopp';
%'SpVi12_1107_WR_MS_LS500mHz2ms6_24Ch_nopp'; 
%'SpVi12_1024_KX_MLStim_26Ch_nopp'; %'039v_0925_2Hz20ms_20mW_28Ch_nopp'; 
% '039v_0927_2Hz20ms_20mW_28Ch_nopp'; % 'SpVi12_133_2Hz2ms_10mW_nopp';
% SpVi12_133_2Hz2ms_10mW_nopp_Ch %SpVi12_133_2Hz2ms_10mW_nopp_Ch
%% From SpikeVisualizationGUI export
spikeData=load([fileName '_Ch' num2str(channelNum) '.mat'],'waveForms','spikeTimes','unitsIdx','samplingRate','selectedUnits');
load([fileName '_Ch' num2str(channelNum) '.mat'],'TTLs');
load([fileName '_Ch' num2str(channelNum) '.mat'], 'traceExcerpt');
traceData=load([fileName '_Ch' num2str(channelNum) '.mat'], 'allTraces','traceInfo');

% read TTL dat file 
%     cd(sessionDir);
%     TTLFileName=[regexp(recName,'\S+?(?=_export)','match','once') '_TTLs.dat'];
%     fid = fopen(TTLFileName, 'r');
%     TTLs = fread(fid,[2,Inf],'int32');
%     fclose(fid);

% TTL times should be sync'ed to recoding start already ...
TTLs.start(:,1)=TTLs.start(:,1)-double(rec_info.recordingStartTime);
TTLs.start(:,2)=TTLs.start(:,1)/double(TTLs.samplingRate{1}/1000); 
TTLs.end(:,1)=TTLs.end(:,1)-double(rec_info.recordingStartTime);
TTLs.end(:,2)=TTLs.end(:,1)/double(TTLs.samplingRate{1}/1000); 

%% From JRClust csv export
% dirListing=dir; dirName=cd;
% infoFileName=dirListing(~cellfun('isempty',cellfun(@(x) strfind(x,'_info.'),...
%     {dirListing.name},'UniformOutput',false))).name;
% recInfo=load(fullfile(dirName,infoFileName),'rec_info');recInfo=recInfo.rec_info;
% 
% JRclustData=load([fileName '_JR.csv']);
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
% %
% % spikes=LoadSpikeData([fileName '_JR.csv'],[],recInfo.numRecChan,recInfo.samplingRate,recInfo.bitResolution);
% % 
% % spikeData.selectedUnits=2;
% % units=cellfun(@(chUnits) chUnits(chUnits==spikeData.selectedUnits),spikes.Units,'UniformOutput',false);
% % spikeTimes=cellfun(@(chUnits,chSpkTimes) chSpkTimes(chUnits==spikeData.selectedUnits),...
% %         spikes.Units,spikes.SpikeTimes,'UniformOutput',false);
% % waveForms=cellfun(@(chUnits,chWf) chWf(:,chUnits==spikeData.selectedUnits),spikes.Units,...
% %         spikes.Waveforms,'UniformOutput',false);
% % spikeTimes=cell2mat(spikeTimes(~cellfun('isempty',units)));
% % waveForms=cell2mat(waveForms(~cellfun('isempty',units))');
% % units=cell2mat(units(~cellfun('isempty',units)));

%% get spike times and convert to binary array
for clusNum=1:size(spikeData.selectedUnits,1)
    %% convert to 1 millisecond bins and plot excerpt
    binSize=1;
    numBin=ceil((max(spikeData.spikeTimes(spikeData.unitsIdx==spikeData.selectedUnits(clusNum)))+1)/...
        (spikeData.samplingRate/1000)/binSize);
    
    [spikeCount,spikeTime]=histcounts(double(spikeData.spikeTimes(...
        spikeData.unitsIdx==spikeData.selectedUnits(clusNum)))/...
        double(spikeData.samplingRate/1000), numBin);
    
    %     foo=spikeData.spikeTimes(spikeData.unitsIdx==spikeData.selectedUnits(clusNum))/30;
    %         figure; bar(spikeTime(1:6000),spikeCount(1:6000),'hist')
    
    
    %% spike density function
    spikeArray = zeros(1,ceil(max(spikeTime))+1);
    spikeArray(ceil(spikeTime(1:end-1)))=spikeCount;
%         sigma=1;
%         convSpikeTime = [zeros(1,sigma*3) fullgauss_filtconv(spikeArray,sigma,0)].*1000;
%         hold on
%     plot([zeros(1,sigma*3) convSpikeTime zeros(1,sigma*3)])
%         plot( convSpikeTime(1:6000-sigma*3))
    
    %% create rasters aligned to TTL
    %define parameters
    preAlignWindow=100;
    postAlignWindow=359;
    TTLtimes=uint32(TTLs.start(:,2)); %deal(TTLs.start(:,1)/double(TTLs.samplingRate{1}/1000)); % uint32(TTLs.start(:,2));
    raster=nan(numel(TTLs.start(:,1)),preAlignWindow+postAlignWindow+1);
    for trialNum=1:numel(TTLs.start(:,1))
        try
            raster(trialNum,:)=spikeArray(...
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
    spikeRasters{clusNum}=raster(~isnan(sum(raster,2)),:);
    %     figure; imagesc(raster)
end

pulseDur=min(mode(TTLs.end-TTLs.start));
IPI=mode(diff(TTLs.start(:,2)))+pulseDur;

%% Figures
% some issue with ttl times from npy -> see CH29 from 'vIRt22_2018-10-16_18-43-54_5100_50ms1Hz5mW_nopp' KS 
for cellNum=1:size(spikeData.selectedUnits,1)
% keep one cell 
% cellNum=2;

figure('Position',[296 149 1504 761],'name',...
    [fileName ' Ch' num2str(channelNum) ' Unit' num2str(spikeData.selectedUnits(cellNum))] );

% waveforms
subplot(3,3,[1,4]); hold on
delay=5;
OptoWaveforms(spikeData,TTLtimes,spikeData.selectedUnits(cellNum),delay,gca)

% rasters
subplot(3,3,[2,5]);
% delay=5;
OptoRasters(spikeRasters(cellNum),preAlignWindow,pulseDur,IPI,gca);
% title(['Channel ' num2str(channelNum) ', Neuron ' num2str(spikeData.selectedUnits(cellNum))],'FontName','Cambria');

% SDF
subplot(3,3,[3,6])
OptoSDF(spikeRasters(cellNum),preAlignWindow,pulseDur,IPI,gca)

% % ISI
% subplot(3,3,4); hold on
% OptoISI(spikeData,TTLtimes,spikeData.selectedUnits(cellNum),gca)
% 
% %ACG
% subplot(3,3,7); hold on
% OptoACG(spikeData,TTLtimes,spikeData.selectedUnits(cellNum),gca)

% raw trace
subplot(3,3,7:9); hold on

msConv=double(spikeData.samplingRate/1000);
excerptTTLtimes=double(TTLtimes(TTLtimes>(traceExcerpt.location-traceExcerpt.excerptSize)/msConv &...
    TTLtimes<(traceExcerpt.location+traceExcerpt.excerptSize)/msConv)-...
    (traceExcerpt.location-traceExcerpt.excerptSize)/msConv)*msConv;
if ~isempty(excerptTTLtimes)
%     excerptTTLtimes=excerptTTLtimes(end); %if wants to keep only one pulse
else % check further out in the trace
    traceExcerpt.location=TTLtimes(1)*msConv;
    winIdxStart=(traceExcerpt.location-traceExcerpt.excerptSize)*traceData.traceInfo.numChan+1;
%     mod(winIdxStart,traceData.traceInfo.numChan)
    winSize=2; %default 1 pulse 
    winIdxEnd=winIdxStart+(winSize*2*traceExcerpt.excerptSize*traceData.traceInfo.numChan);
    excerptWindow=winIdxStart:winIdxEnd-1;
%     size(excerptWindow,2)>(2*traceExcerpt.excerptSize*traceData.traceInfo.numChan) 
    
    traceExcerpt.data=traceData.allTraces.Data(excerptWindow);
    traceExcerpt.data=reshape(traceExcerpt.data,[traceData.traceInfo.numChan traceExcerpt.excerptSize*2*winSize]);
    
    preprocOption={'CAR','all'};
    traceExcerpt.data=PreProcData(traceExcerpt.data,30000,preprocOption);
%     figure; plot(dataExcerpt(11,:))
    traceExcerpt.data=traceExcerpt.data(channelNum,:);
    excerptTTLtimes=double(TTLtimes(TTLtimes>(traceExcerpt.location-traceExcerpt.excerptSize)/msConv &...
    TTLtimes<(traceExcerpt.location+traceExcerpt.excerptSize*winSize)/msConv)-...
    (traceExcerpt.location-traceExcerpt.excerptSize)/msConv)*msConv;
    traceExcerpt.spkTimes{cellNum}=NaN;
end %plot anyway
OptoRawTrace(traceExcerpt,traceExcerpt.spkTimes(cellNum),msConv,excerptTTLtimes,pulseDur,'center',gca)
end

