% Opto rasters

%% plot raster showing all channels
figure('Position',[1050 120 750 790]);
subplot(1,1)
colormap bone;
MeanChan=cellfun(@(x) conv_raster(x),Rasters.units(:,1),'UniformOutput',false);
MeanChan=cell2mat(MeanChan);
% subplot(1,1)
imagesc(zscore(MeanChan,[])); %
% imagesc(MeanChan);
xlabel('Time (ms)');
ylabel('Channels','FontWeight','bold','FontSize',12);
% draw alignment bar
currylim=get(gca,'YLim');
currxlim=get(gca,'XLim');midl=round(currxlim/20)*10;
set(gca,'XTick',[midl-preAlignWindow/2 midl midl+postAlignWindow/2]);
set(gca,'XTickLabel',[-preAlignWindow/2 0 postAlignWindow/2]);
%opto stim patch
patch([repmat(midl,1) repmat(midl+syncTrials.end(1)-syncTrials.start(1),1)], ...
    [[0 currylim] fliplr([0 currylim])], ...
    [0 0 0 0],[0.3 0.75 0.93],'EdgeColor','none','FaceAlpha',0.1);

% patch([repmat(midl-3,1) repmat(midl+3,1)], ...
%     [[0 currylim] fliplr([0 currylim])], ...
%     [0 0 0 0],[0.8 0 0],'EdgeColor','none','FaceAlpha',0.8);
title('Neural response to 100% stimulation intensity, aligned to stimulation onset');
hcb = colorbar('southoutside');
hcb.Label.String = 'z-scored firing rate';

MeanChan=cellfun(@(x) conv_raster(x),Rasters.units(:),'UniformOutput',false);
MeanChan=cell2mat(MeanChan);
subplot(1)
imagesc(zscore(MeanChan,[]));
% imagesc(MeanChan);
xlabel('Time');
ylabel('Channel','FontWeight','bold','FontSize',12);
% draw alignment bar
currylim=get(gca,'YLim');
currxlim=get(gca,'XLim');midl=round(currxlim/2);
set(gca,'XTick',[midl-500 midl midl+500]);
set(gca,'XTickLabel',[-500 0 500]);
patch([repmat(midl-3,1) repmat(midl+3,1)], ...
    [[0 currylim] fliplr([0 currylim])], ...
    [0 0 0 0],[0.8 0 0],'EdgeColor','none','FaceAlpha',0.8);
title('Neural response to 80% stimulation intensity, aligned to stimulation onset');
hcb = colorbar('southoutside');
hcb.Label.String = 'z-scored firing rate';

%% plot sdf
BestChan=find(mean(MeanChan)==max(mean(MeanChan)));
start=1;
stop=size(Rasters.units{BestChan,1});
conv_sigma=1;
alignmtt=preAlignWindow;
xTickSteps=round(preAlignWindow/50)*10;
[sdf{1}, ~, rastsem{1}]=conv_raster(Rasters.units{BestChan,1},conv_sigma,start,stop);
[sdf{2}, ~, rastsem{2}]=conv_raster(Rasters.units{BestChan},conv_sigma,start,stop);
% figure('Position',[1469 542 417 417]);
subplot(1)
colormap default;
cmap = colormap(gcf);
hold on;

%plot sem
startAlignPloth=gca; box off; %subplot(1,1);hold on; box off;
patch([1:length(sdf{1}),fliplr(1:length(sdf{1}))],[sdf{1}-rastsem{1},fliplr(sdf{1}+rastsem{1})],...
    [0.16 0.38 0.27],'EdgeColor','none','FaceAlpha',0.2);
% endAlignPloth=subplot(1);hold on; box off;
% patch([1:length(sdf{2}),fliplr(1:length(sdf{2}))],[sdf{2}-rastsem{2},fliplr(sdf{2}+rastsem{2})],cmap(22,:),'EdgeColor','none','FaceAlpha',0.1);
%plot sdfs
FRploth=plot(startAlignPloth,sdf{1},'Color',[0.16 0.38 0.27],'LineWidth',1.8);

% set(startAlignPloth,'XTick',xTickSteps-(start+3*conv_sigma):xTickSteps:(stop-start-6*conv_sigma));
set(startAlignPloth,'XTick',xTickSteps-(start+3*conv_sigma):xTickSteps:(stop-start-6*conv_sigma));
set(startAlignPloth,'XTickLabel',-(alignmtt-xTickSteps):xTickSteps:stop-(alignmtt+xTickSteps));
axis(startAlignPloth,'tight');
set(startAlignPloth,'Color','white','TickDir','out','FontName','Cambria','FontSize',10);
hxlabel=xlabel(startAlignPloth,'Time (ms)','FontName','Cambria','FontSize',12);
hylabel=ylabel(startAlignPloth,'Firing rate (spikes/s)','FontName','Cambria','FontSize',12);

% plot(endAlignPloth,sdf{2},'Color',cmap(22,:),'LineWidth',1.8);
%
% set(endAlignPloth,'XTick',xTickSteps-(start+3*conv_sigma):xTickSteps:(stop-start-6*conv_sigma));
% set(endAlignPloth,'XTickLabel',-(alignmtt-xTickSteps):xTickSteps:stop-(alignmtt+xTickSteps));
% axis(endAlignPloth,'tight');
% set(endAlignPloth,'Color','white','TickDir','out','FontName','Cambria','FontSize',10);
% hxlabel=xlabel(endAlignPloth,'Time (ms)','FontName','Cambria','FontSize',12);
% hylabel=ylabel(endAlignPloth,'Firing rate (spikes/s)','FontName','Cambria','FontSize',12);

% draw alignment bar
currylim=get(startAlignPloth,'YLim');
axes(startAlignPloth)
% opto stim bar
OptoStimh=patch([repmat(alignmtt-(start+3*conv_sigma),1) repmat(alignmtt-(start+3*conv_sigma)+syncTrials.end(1)-syncTrials.start(1),1)], ...
    [[0 currylim] fliplr([0 currylim])], ...
    [0 0 0 0],[0.3 0.75 0.93],'EdgeColor','none','FaceAlpha',0.5);
% "regular" alignement bar
% patch([repmat((alignmtt-(start+3*conv_sigma))-2,1) repmat((alignmtt-(start+3*conv_sigma))+2,1)], ...
%     [[0 currylim] fliplr([0 currylim])], ...
%     [0 0 0 0],[1 0 0],'EdgeColor','none','FaceAlpha',0.5);
% axes(endAlignPloth)
% patch([repmat((alignmtt-(start+3*conv_sigma))-2,1) repmat((alignmtt-(start+3*conv_sigma))+2,1)], ...
%     [[0 currylim] fliplr([0 currylim])], ...
%     [0 0 0 0],[1 0 0],'EdgeColor','none','FaceAlpha',0.5);

%legend
legend([FRploth,OptoStimh],{'Average firing rate','Optical stimulation'});
legend('boxoff')
% text(xTickSteps,currylim-20,['Channel ' num2str(BestChan)],'FontName','Cambria');
title(['Channel ' num2str(BestChan)],'FontName','Cambria');


%% plot on-pulse off-pulse waveforms, ISI and ACG

pulseIdx=false(size(unitST,1),size(syncTrials.start,1));
%get wich spike time occur during TTL
for TTLNum=1:size(syncTrials.start,1)
    pulseIdx(:,TTLNum)=unitST>syncTrials.start(TTLNum) & unitST<syncTrials.end(TTLNum);
end
onSpikes=logical(sum(pulseIdx));

% plot waveforms
figure;hold on
plot(mean(spikeData.(['Clus' num2str(clusNum)]).Waveforms(~onSpikes,:)),'linewidth','color',cmap(clusNum,:))
plot(mean(spikeData.(['Clus' num2str(clusNum)]).Waveforms(onSpikes,:)),'linewidth','color',[0.3 0.75 0.93])
legend('Pulse On','Pulse Off','location','northeast')
set(gca,'xtick',linspace(0,size(spikeData.(['Clus' num2str(clusNum)]).Waveforms),5),...
    'xticklabel',round(linspace(-round(size(spikeData.(['Clus' num2str(clusNum)]).Waveforms)/2),...
    round(size(spikeData.(['Clus' num2str(clusNum)]).Waveforms)/2),5)/30),'TickDir','out');
box off; %axis('tight');
xlabel('Time (ms)')
ylabel('Voltage (mV)')
title(['Neuron ' num2str(clusNum) ' waveform on/off pulse'])
set(gca,'Color','white','FontSize',12,'FontName','calibri');

% plot ISI and ACG
figure;
unitST_onPulse=unitST(onSpikes);
unitST_offPulse=unitST(~onSpikes);
% compute interspike interval
ISI_onPulse=diff(unitST_onPulse);
ISI_offPulse=diff(unitST_offPulse);
subplot(2,1,1); hold on;
ISI_offPulsehist=histogram(double(ISI_offPulse),0:10:max(ISI_offPulse)+1);  %,'Normalization','probability'
ISI_offPulsehist.FaceColor = cmap(clusNum,:);
ISI_offPulsehist.EdgeColor = 'k';
ISI_onPulsehist=histogram(double(ISI_onPulse),0:10:max(ISI_onPulse)+1);  %,'Normalization','probability'
ISI_onPulsehist.FaceColor = [0.3 0.75 0.93];
ISI_onPulsehist.EdgeColor = 'k';
xlabel('Inter-spike Interval distribution (ms)')
axis('tight');box off;
legend('Pulse On','Pulse Off','location','northeast')
set(gca,'xlim',[0 500],'Color','white','FontSize',10,'FontName','calibri','TickDir','out');
hold off
title(['Neuron ' num2str(clusNum) ' ISI and ACG'])

% plot ACG
spikeTimeIdx_onPulse=zeros(1,unitST_onPulse(end));
spikeTimeIdx_onPulse(unitST_onPulse)=1;
binSize=1;
numBin=ceil(size(spikeTimeIdx_onPulse)/binSize);
binUnits_onPulse = histcounts(double(unitST_onPulse), linspace(0,size(spikeTimeIdx_onPulse),numBin));
binUnits_onPulse(binUnits_onPulse>1)=1; %no more than 1 spike per ms

spikeTimeIdx_offPulse=zeros(1,unitST_offPulse(end));
spikeTimeIdx_offPulse(unitST_offPulse)=1;
numBin=ceil(size(spikeTimeIdx_offPulse)/binSize);
binUnits_offPulse = histcounts(double(unitST_offPulse), linspace(0,size(spikeTimeIdx_offPulse),numBin));
binUnits_offPulse(binUnits_offPulse>1)=1; %no more than 1 spike per ms

% compute autocorrelogram
[ACG_onPulse,lags_onPulse]=xcorr(double(binUnits_onPulse)00,'coeff'); %'coeff'
ACG_onPulse(lags_onPulse==0)=0;
[ACG_offPulse,lags_offPulse]=xcorr(double(binUnits_offPulse)00,'coeff'); %'coeff'
ACG_offPulse(lags_offPulse==0)=0;

subplot(2,1); hold on
ACGh_onPulse=bar(lags_onPulse,ACG_onPulse);
ACGh_offPulse=bar(lags_offPulse,ACG_offPulse);
ACGh_onPulse.FaceColor = [0.3 0.75 0.93];
ACGh_onPulse.EdgeColor = 'none';
ACGh_offPulse.FaceColor = cmap(clusNum,:);
ACGh_offPulse.EdgeColor = 'none';
axis('tight');box off;
xlabel('Autocorrelogram (5 ms bins)');
legend('Pulse On','Pulse Off','location','northeast');
set(gca,'xlim',[-300 300],'Color','white','FontSize',10,'FontName','calibri','TickDir','out');
hold off

%% plot spikes and trials
figure; hold on
% plot(int16(Spikes.Offline_Threshold.data{3, 1})*max(rawData(3,:)),'ko')
% plot(Spikes.Offline_Sorting.data{3, 1},ones(1,size(Spikes.Offline_Sorting.data{3, 1},1))*0.5,'sr')
plot(spikeData.Clus1.SpikeTimes,...
    int16(ones(1,size(spikeData.Clus1.SpikeTimes,1)))*max(rawData(3,:)),...
    'linestyle','none','marker','o','MarkerSize',5,'MarkerEdgeColor',cmap(1,:),'MarkerFaceColor','none')
plot(spikeData.Clus2.SpikeTimes,...
    int16(ones(1,size(spikeData.Clus2.SpikeTimes,1)))*max(rawData(3,:))+50,...
    'linestyle','none','marker','o','MarkerSize',5,'MarkerEdgeColor',cmap(2,:),'MarkerFaceColor','none')
plot(rawData(3,1:129000));
plot(spikeData.Clus3.SpikeTimes,...
    int16(ones(1,size(spikeData.Clus3.SpikeTimes,1)))*max(rawData(3,:))+100,...
    'linestyle','none','marker','o','MarkerSize',5,'MarkerEdgeColor',cmap(3,:),'MarkerFaceColor','none')

yLims=get(gca,'ylim');
for TTLNum=1:size(syncTrials.start,1)
    patch([syncTrials.start(TTLNum)*30:syncTrials.end(TTLNum)*30,...
        fliplr(syncTrials.start(TTLNum)*30:syncTrials.end(TTLNum)*30)],...
        [ones(1,syncTrials.end(TTLNum)*30-syncTrials.start(TTLNum)*30+1)*yLims(1),...
        ones(1,syncTrials.end(TTLNum)*30-syncTrials.start(TTLNum)*30+1)*yLims],...
        [0.3 0.75 0.93],'EdgeColor','none','FaceAlpha',0.5);
    patch([syncTrials.start(TTLNum)*30:syncTrials.end(TTLNum)*30,...
        fliplr(syncTrials.start(TTLNum)*30:syncTrials.end(TTLNum)*30)],...
        [ones(1,syncTrials.end(TTLNum)*30-syncTrials.start(TTLNum)*30+1)*yLims(1),...
        ones(1,syncTrials.end(TTLNum)*30-syncTrials.start(TTLNum)*30+1)*(yLims(1)+100)],...
        [0 0 0],'EdgeColor','none','FaceAlpha',0.5);
end
tickNum=round((size(rawData)/rec_info.samplingRate))*100; % ~ every 10sec
set(gca,'xtick',linspace(0,size(rawData),tickNum),...
    'xticklabel',(linspace(0,size(rawData)/rec_info.samplingRate,tickNum)),'TickDir','out');
% set(gca,'ytick',linspace(0,double(max(abs(max(rawData)))*(ChN-1))-...
%         mean(mean(rawData(ChN,:))),size(rawData,1)),'yticklabel',...
%     cellfun(axis_name, num2cell(1:size(rawData,1)), 'UniformOutput', false))
axis('tight');box off;
xlabel('Time (10 ms.)')
ylabel('Voltage (uV)')
set(gca,'Color','white','FontSize',12,'FontName','calibri');

set(gca,'xlim',[127000 129000]);
