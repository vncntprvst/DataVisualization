fileName='SpVi12_1107_WR_MS_LS500mHz2ms6_nopp';
%'SpVi12_1107_WR_MS_LS500mHz2ms6_24Ch_nopp'; 
%'SpVi12_1024_KX_MLStim_26Ch_nopp'; %'039v_0925_2Hz20ms_20mW_28Ch_nopp'; 
% '039v_0927_2Hz20ms_20mW_28Ch_nopp'; % 'SpVi12_133_2Hz2ms_10mW_nopp';
channelNum=22;
% SpVi12_133_2Hz2ms_10mW_nopp_Ch %SpVi12_133_2Hz2ms_10mW_nopp_Ch
spikeData=load([fileName '_Ch' num2str(channelNum) '.mat'],'waveForms','spikeTimes','unitsIdx','samplingRate','selectedUnits');
load([fileName '_Ch' num2str(channelNum) '.mat'],'TTLs');
load([fileName '_Ch' num2str(channelNum) '.mat'], 'dataExcerpt');

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
    preAlignWindow=20;
    postAlignWindow=59;
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

%% Figures
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
OptoRasters(spikeRasters(cellNum),preAlignWindow,pulseDur,gca);
% title(['Channel ' num2str(channelNum) ', Neuron ' num2str(spikeData.selectedUnits(cellNum))],'FontName','Cambria');

% SDF
subplot(3,3,[3,6])
OptoSDF(spikeRasters(cellNum),preAlignWindow,pulseDur,gca)

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
excerptTTLtimes=double(TTLtimes(TTLtimes>(dataExcerpt.location-dataExcerpt.excerptSize)/msConv &...
    TTLtimes<(dataExcerpt.location+dataExcerpt.excerptSize)/msConv)-...
    (dataExcerpt.location-dataExcerpt.excerptSize)/msConv)*msConv;
if ~isempty(excerptTTLtimes)
%     excerptTTLtimes=excerptTTLtimes(end); %if wants to keep only one pulse
end %plot anyway
OptoRawTrace(dataExcerpt,dataExcerpt.spkTimes(cellNum),msConv,excerptTTLtimes,'',gca)
end

