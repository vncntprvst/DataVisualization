fileName= 'vIRt20_0704_5119+_CutWhiskers_LS50mW10Hz5ms_nopp_Ch51'
%'vIRt16_0403_WR_4975_pole_nopp_SelectedData_Multi'; %
% 'SpVi16_0403_WR_5225_pole_nopp_SelectedData_Multi'; 
%'SpVi16_0403_WR_5225_MS_nopp_SelectedData_Multi'; 
%'SpVi16_0403_WR_4850_LS1Hz2ms100mW_nopp_SelectedData_Multi';
%'SpVi16_0403_WR_4850_LS1Hz2ms100mW_nopp_SelectedData_Ch14U33';
% fileName='SpVi16_0403_WR_4850_LS1Hz2ms100mW_nopp_SelectedData_Ch3';
recName='vIRt20_0704_5119_CutWhiskers_LS'; %'vIRt16_0403_WR_4975_pole'; %'SpVi16_0403_WR_5225_MS'; %'SpVi16_0403_WR_4850_LS1Hz2ms100mW';
%% From SpikeVisualizationGUI export
spikeData=load([fileName '.mat'],'waveForms','spikeTimes','unitsIdx','samplingRate','selectedUnits');
load([fileName '.mat'],'TTLs');
load([fileName '.mat'], 'traceExcerpt');
traceData=load([fileName '.mat'], 'allTraces','traceInfo');

% TTLs.TTLtimes=35000;
% TTLs.start=35000;
% TTLs.end=38000;

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
    
    [spikeCount,spikeTime]=histcounts(double(spikeData.spikeTimes(spikeData.unitsIdx==spikeData.selectedUnits(clusNum)))/...
        double(spikeData.samplingRate/1000), numBin);
    
    %     foo=spikeData.spikeTimes(spikeData.unitsIdx==spikeData.selectedUnits(clusNum))/30;
    %         figure; bar(spikeTime(1:6000),spikeCount(1:6000),'hist')
    
    
    %% spike density function
    spikeArray = zeros(1,ceil(max(spikeTime))+1);
    spikeArray(ceil(spikeTime(1:end-1)))=spikeCount;
    %     sigma=1;
    %     convSpikeTime = [zeros(1,sigma*3) fullgauss_filtconv(spikeArray,sigma,0)].*1000;
    %     hold on
    % plot([zeros(1,sigma*3) convSpikeTime zeros(1,sigma*3)])
    %     plot( convSpikeTime(1:6000-sigma*3))
    
    %% create rasters aligned to TTL
    %define parameters
    preAlignWindow=500;
    postAlignWindow=4000;
    TTLtimes=uint32(TTLs.TTLtimes)/(TTLs.samplingRate/1000);
    raster=nan(numel(TTLs.TTLtimes),preAlignWindow+postAlignWindow+1);
    for trialNum=1:numel(TTLs.TTLtimes)
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

pulseDur=mode(TTLs.end-TTLs.start);
IPI=mode(TTLs.interval)+pulseDur;
%% Population plots
conv_sigma=1;shiftVal=conv_sigma*3;
% allSDF=cellfun(@(x) conv_raster(x,conv_sigma,1), spikeRasters,'UniformOutput' , false);
% allSDF=vertcat(allSDF{:});
allSDF=vertcat(spikeRasters{:});
figure;
colormap(hot) %flipud(gray));
imagesc(allSDF); %
% caxis([0 200])
xlabel('Time (ms)');
ylabel('Neuron#','FontSize',12); %'FontWeight','bold'
% draw alignment bar
currylim=get(gca,'YLim');
%     currxlim=get(gca,'XLim');%midl=round(currxlim(2)/20)*10;
%     set(gca,'XTick',preAlignWindow:50:max(get(gca,'xlim')));
%     set(gca,'XTickLabel',0:50:max(get(gca,'xlim'))-preAlignWindow,'FontSize',10,'FontName','calibri','TickDir','out');
set(gca,'XLim',[0.5 preAlignWindow+660.5],'XTick',[0:10:preAlignWindow+660]-shiftVal);
set(gca,'XTickLabel',(0:10:preAlignWindow+660)-preAlignWindow,'FontSize',10,'FontName','calibri','TickDir','out');

%opto stim patch
patch([preAlignWindow-shiftVal preAlignWindow-shiftVal preAlignWindow+pulseDur preAlignWindow+pulseDur], ...
    [[0 currylim(2)] fliplr([0 currylim(2)])], ...
    [0 0 0 0],[0.3 0.75 0.93],'EdgeColor','none','FaceAlpha',0.5);
set(gca,'Color','white','FontSize',18,'FontName','Helvetica');

%% Individual Plots
for cellNum=1:size(spikeData.selectedUnits,1)
% keep one cell 
% cellNum=2;

figure('Position',[296 149 1504 761],'name',...
    [recName ' Unit' num2str(spikeData.selectedUnits(cellNum))] );

% waveforms
subplot(3,3,[1,4]); hold on
delay=5;
OptoWaveforms(spikeData,TTLtimes,spikeData.selectedUnits(cellNum),delay,gca)

% rasters
subplot(3,3,[2,5]);
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
OptoRawTrace(traceExcerpt,traceExcerpt.spkTimes(cellNum),msConv,excerptTTLtimes,'',gca)
end

