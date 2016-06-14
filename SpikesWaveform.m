%% Get spike times and waveforms 

%% Get file path
[fileName,dirName] = uigetfile({'*.mat; *.hdf5','Processed data';'*.dat','Flat data';...
    '*.*','All Files' },'Exported data','C:\Data\export');
cd(dirName);

%Select channels
% KeepChans=Spikes.Offline_Threshold.channel;
KeepChans=3;

% load file data
if strfind(fileName,'.mat')
    load(fileName);
    %load other data
    fileName=regexp(fileName,'.+(?=_\w+.\w+$)','match');
    load([fileName{:} '_raw.mat']);
    load([fileName{:} '_trials.mat']);
    load([fileName{:} '_spikes.mat']);
    load([fileName{:} '_info.mat']);
    
    % from Spike2
    Spikes.Offline_SpkSort.Units{KeepChans,1}=nw_401.codes(:,1);
    Spikes.Offline_SpkSort.SpikeTimes{3,1}=uint32(nw_401.times*rec_info.samplingRate);
    Spikes.Offline_SpkSort.Waveforms{3,1}=nw_401.values;
    
    ResortedSpikes=load('C:\Data\export\PrV75_61_BR_8Ch\PrV75_61_opt_stim_2003_spikesResorted.mat');
    
%     try
%         Behavior=processBehaviorData;
%     catch
%         Behavior=[];
%     end
elseif strfind(fileName,'.hdf5')
    fileName=regexp(fileName,'\w+(?=\.\w+\.)','match','once');
    Spikes.Offline_SpkSort.Units{3,1}=h5read([fileName '.clusters.hdf5'],'/clusters_2');
    Spikes.Offline_SpkSort.SpikeTimes{3,1}=h5read([fileName '.clusters.hdf5'],'/times_2');
    Spikes.Offline_SpkSort.Waveforms{3,1}=h5read([fileName '.clusters.hdf5'],'/data_2');
    Spikes.Offline_SpkSort.Waveforms=h5read([fileName '.templates.hdf5'],'/temp_data');
    Spikes.Offline_SpkSort.templates{10,1}.spiketimes=h5read([fileName '.result.hdf5'],'/spiketimes/temp_10');
    Spikes.Offline_SpkSort.templates{10,1}.amplitudes=h5read([fileName '.result.hdf5'],'/amplitudes/temp_10');
end

%get clusters spiketimes
for clusNum=min(unique(ResortedSpikes.Spikes.inGUI.Units{3,1})):max(unique(ResortedSpikes.Spikes.inGUI.Units{3,1}))   
%     Ch3.(['Clus' num2str(clusNum)]).SpikeTimes=Spikes.Offline_SpkSort.SpikeTimes{KeepChans,1}(Spikes.Offline_SpkSort.Units{KeepChans,1}==clusNum);
%     Ch3.(['Clus' num2str(clusNum)]).Waveforms=Spikes.Offline_SpkSort.Waveforms{KeepChans,1}(Spikes.Offline_SpkSort.Units{KeepChans,1}==clusNum,:);
    Ch3.(['Clus' num2str(clusNum+1)]).SpikeTimes=ResortedSpikes.Spikes.inGUI.SpikeTimes{KeepChans,1}(ResortedSpikes.Spikes.inGUI.Units{KeepChans,1}==clusNum);
    Ch3.(['Clus' num2str(clusNum+1)]).Waveforms=ResortedSpikes.Spikes.inGUI.Waveforms{KeepChans,1}(:,ResortedSpikes.Spikes.inGUI.Units{KeepChans,1}==clusNum)';
end
% Ch3.Clus1.Waveforms=ExtractChunks(rawData(3,:),Ch3.Clus1.SpikeTimes,82,'tzero');

%plot waveforms
figure;hold on
colormap lines; cmap=colormap;
% for clusNum=min(unique(Spikes.Offline_SpkSort.Units{3,1})):max(unique(Spikes.Offline_SpkSort.Units{3,1}))   
for clusNum=1:3
    plot(mean(Ch3.(['Clus' num2str(clusNum)]).Waveforms),'linewidth',2)
end
% legend({'0','1','2','3'})
legend('Neuron 1','Neuron 2','Neuron 3','location','southeast')
set(gca,'xtick',linspace(0,size(Ch3.Clus1.Waveforms,2),5),...
    'xticklabel',round(linspace(-round(size(Ch3.Clus1.Waveforms,2)/2),round(size(Ch3.Clus1.Waveforms,2)/2),5)/30,2),'TickDir','out');
axis('tight');box off;
xlabel('Time (ms)')
ylabel('Voltage (uV) x10 amp')
set(gca,'Color','white','FontSize',12,'FontName','calibri');

% by subplots
figure;
for clusNum=min(unique(Spikes.Offline_SpkSort.Units{3,1})):max(unique(Spikes.Offline_SpkSort.Units{3,1}))   
    subplot(2,3,double(clusNum+1))
    plot(mean(Ch3.(['Clus' num2str(clusNum)]).Waveforms),'color',cmap(double(clusNum+1),:))
end

% plot(mean(Ch3.Clus1.Waveforms)+std(double(Ch3.Clus1.Waveforms)))
% plot(mean(Ch3.Clus1.Waveforms)-std(double(Ch3.Clus1.Waveforms)))

%
Trials.start=Trials.start-Trials.startClockTime;
Trials.end=Trials.end-Trials.startClockTime;

% plot spikes and trials
figure; hold on
% plot(int16(Spikes.Offline_Threshold.data{3, 1})*max(rawData(3,:)),'ko')
% plot(Spikes.Offline_SpkSort.data{3, 1},ones(1,size(Spikes.Offline_SpkSort.data{3, 1},1))*0.5,'sr')
plot(Ch3.Clus1.SpikeTimes,...
    int16(ones(1,size(Ch3.Clus1.SpikeTimes,1)))*max(rawData(3,:)),...
    'linestyle','none','marker','o','MarkerSize',5,'MarkerEdgeColor',cmap(1,:),'MarkerFaceColor','none')
plot(Ch3.Clus2.SpikeTimes,...
    int16(ones(1,size(Ch3.Clus2.SpikeTimes,1)))*max(rawData(3,:))+50,...
    'linestyle','none','marker','o','MarkerSize',5,'MarkerEdgeColor',cmap(2,:),'MarkerFaceColor','none')
plot(rawData(3,:));
plot(Ch3.Clus3.SpikeTimes,...
    int16(ones(1,size(Ch3.Clus3.SpikeTimes,1)))*max(rawData(3,:))+100,...
    'linestyle','none','marker','o','MarkerSize',5,'MarkerEdgeColor',cmap(3,:),'MarkerFaceColor','none')

yLims=get(gca,'ylim');
for TTLNum=1:size(Trials.start,1)
patch([Trials.start(TTLNum,2)*30:Trials.end(TTLNum,2)*30,...
    fliplr(Trials.start(TTLNum,2)*30:Trials.end(TTLNum,2)*30)],...
    [ones(1,Trials.end(TTLNum,2)*30-Trials.start(TTLNum,2)*30+1)*yLims(1),...
    ones(1,Trials.end(TTLNum,2)*30-Trials.start(TTLNum,2)*30+1)*yLims(2)],...
    [0.3 0.75 0.93],'EdgeColor','none','FaceAlpha',0.5);
patch([Trials.start(TTLNum,2)*30:Trials.end(TTLNum,2)*30,...
    fliplr(Trials.start(TTLNum,2)*30:Trials.end(TTLNum,2)*30)],...
    [ones(1,Trials.end(TTLNum,2)*30-Trials.start(TTLNum,2)*30+1)*yLims(1),...
    ones(1,Trials.end(TTLNum,2)*30-Trials.start(TTLNum,2)*30+1)*(yLims(1)+100)],...
    [0 0 0],'EdgeColor','none','FaceAlpha',0.5);
end
tickNum=round((size(rawData,2)/rec_info.samplingRate)/60)*6; % ~ every 10sec
set(gca,'xtick',linspace(0,size(rawData,2),tickNum),...
    'xticklabel',round(linspace(0,size(rawData,2)/rec_info.samplingRate,tickNum)),'TickDir','out');
% set(gca,'ytick',linspace(0,double(max(abs(max(rawData)))*(ChN-1))-...
%         mean(mean(rawData(ChN,:))),size(rawData,1)),'yticklabel',...
%     cellfun(axis_name, num2cell(1:size(rawData,1)), 'UniformOutput', false))
axis('tight');box off;
xlabel('Time (s.)')
ylabel('Voltage (uV)')
set(gca,'Color','white','FontSize',12,'FontName','calibri');

%Keep only verified trials
[ephysCommonTrials, behaviorCommonTrials]=match_trials(Trials,Behavior);

%% gather data from selected channels
Rasters.channels=cell(length(KeepChans),2);
Rasters.epochnames={'BeginTrial','EndTrial'};
preAlignWindow(1)=0.08; %80ms
postAlignWindow(1)=0.08; %80ms
for chan=1:length(KeepChans)
    preAlignWindow(2)=round(uint64(Spikes.Offline_Threshold.samplingRate(chan,2))/(1/preAlignWindow(1)));
    postAlignWindow(2)=round(uint64(Spikes.Offline_Threshold.samplingRate(chan,2))/(1/postAlignWindow(1)));

    %     downSamplingRatio=uint64(Spikes.Offline_Threshold.samplingRate(chan,1)/Spikes.Offline_Threshold.samplingRate(chan,2));
    [Rasters.channels{chan,1},Rasters.channels{chan,2}]=deal(zeros(size(Trials.start,1),preAlignWindow(2)+postAlignWindow(2)+1));
    Spkt=Spikes.Offline_Threshold.data{KeepChans(chan),2}(1,:);
    if Trials.end(end)>size(Spkt,2)
        continue
    end
    for trialnb=1:size(Trials.start,1)
        
        %Collect spikes from 1st epoch (begining of trial)
        RastSWin=Trials.start(trialnb,2)-preAlignWindow(2); 
        RastEWin=Trials.start(trialnb,2)+postAlignWindow(2); 
        SpikeTimes=Spkt(RastSWin:RastEWin);
        Rasters.channels{chan,1}(trialnb,:)=SpikeTimes;
        %Collect spikes from 2nd epoch (end of trial)
        RastSWin=Trials.end(trialnb,2)-preAlignWindow(2); % 1 sec before
        RastEWin=Trials.end(trialnb,2)+postAlignWindow(2); % 1/2 sec afer
        SpikeTimes=Spkt(RastSWin:RastEWin);
        Rasters.channels{chan,2}(trialnb,:)=SpikeTimes;

        %         %Collect spikes from 1st epoch (begining of trial)
        %         RastSWin=Trials.start(trialnb,2)-uint64(SamplingRate); % 1 sec before
        %         RastEWin=Trials.start(trialnb,2)+uint64(SamplingRate); % 1 sec afer
        %         SpikeTimes=round((Spkt(Spkt>RastSWin & Spkt<RastEWin)-RastSWin)/uint64(SamplingRate/1000));
        %         SpikeTimes(SpikeTimes==0)=1; %no 0 indices
        %         SpikeTimes=unique(SpikeTimes);
        %         Rasters.channels{chan,1}(trialnb,SpikeTimes)=1;
        %         %Collect spikes from 2nd epoch (end of trial)
        %         RastSWin=Trials.end(trialnb,2)-uint64(SamplingRate); % 1 sec before
        %         RastEWin=Trials.end(trialnb,2)+uint64(SamplingRate); % 1 sec afer
        %         SpikeTimes=round((Spkt(Spkt>RastSWin & Spkt<RastEWin)-RastSWin)/uint64(SamplingRate/1000));
        %         SpikeTimes(SpikeTimes==0)=1; %no 0 indices
        %         SpikeTimes=unique(SpikeTimes);
        %         Rasters.channels{chan,2}(trialnb,SpikeTimes)=1;
    end
end
%% plot raster showing all channels
figure('Position',[1050 120 750 790]);
subplot(1,2,1)
colormap bone;
MeanChan=cellfun(@(x) conv_raster(x),Rasters.channels(:,1),'UniformOutput',false);
MeanChan=cell2mat(MeanChan);
% subplot(1,2,1)
imagesc(zscore(MeanChan,[],2)); %
% imagesc(MeanChan);
xlabel('Time (ms)');
ylabel('Channels','FontWeight','bold','FontSize',12);
% draw alignment bar
currylim=get(gca,'YLim');
currxlim=get(gca,'XLim');midl=round(currxlim(2)/20)*10;
set(gca,'XTick',[midl-preAlignWindow(2)/2 midl midl+postAlignWindow(2)/2]);
set(gca,'XTickLabel',[-preAlignWindow(2)/2 0 postAlignWindow(2)/2]);
%opto stim patch
patch([repmat(midl,1,2) repmat(midl+Trials.end(1,2)-Trials.start(1,2),1,2)], ...
    [[0 currylim(2)] fliplr([0 currylim(2)])], ...
    [0 0 0 0],[0.3 0.75 0.93],'EdgeColor','none','FaceAlpha',0.1);

% patch([repmat(midl-3,1,2) repmat(midl+3,1,2)], ...
%     [[0 currylim(2)] fliplr([0 currylim(2)])], ...
%     [0 0 0 0],[0.8 0 0],'EdgeColor','none','FaceAlpha',0.8);
title('Neural response to 100% stimulation intensity, aligned to stimulation onset');
hcb = colorbar('southoutside');
hcb.Label.String = 'z-scored firing rate';

MeanChan=cellfun(@(x) conv_raster(x),Rasters.channels(:,2),'UniformOutput',false);
MeanChan=cell2mat(MeanChan);
subplot(1,2,2)
imagesc(zscore(MeanChan,[],2));
% imagesc(MeanChan);
xlabel('Time');
ylabel('Channel','FontWeight','bold','FontSize',12);
% draw alignment bar
currylim=get(gca,'YLim');
currxlim=get(gca,'XLim');midl=round(currxlim(2)/2);
set(gca,'XTick',[midl-500 midl midl+500]);
set(gca,'XTickLabel',[-500 0 500]);
patch([repmat(midl-3,1,2) repmat(midl+3,1,2)], ...
    [[0 currylim(2)] fliplr([0 currylim(2)])], ...
    [0 0 0 0],[0.8 0 0],'EdgeColor','none','FaceAlpha',0.8);
title('Neural response to 80% stimulation intensity, aligned to stimulation onset');
hcb = colorbar('southoutside');
hcb.Label.String = 'z-scored firing rate';

%% plot sdf
BestChan=find(mean(MeanChan,2)==max(mean(MeanChan,2)));
start=1;
stop=size(Rasters.channels{BestChan,1},2);
conv_sigma=1;
alignmtt=preAlignWindow(2);
xTickSteps=round(preAlignWindow(2)/50)*10;
[sdf{1}, ~, rastsem{1}]=conv_raster(Rasters.channels{BestChan,1},conv_sigma,start,stop);
[sdf{2}, ~, rastsem{2}]=conv_raster(Rasters.channels{BestChan,2},conv_sigma,start,stop);
% figure('Position',[1469 542 417 417]);
subplot(1,2,2)
colormap default;
cmap = colormap(gcf);
hold on;

%plot sem
startAlignPloth=gca; box off; %subplot(1,2,1);hold on; box off;
patch([1:length(sdf{1}),fliplr(1:length(sdf{1}))],[sdf{1}-rastsem{1},fliplr(sdf{1}+rastsem{1})],...
    [0.16 0.38 0.27],'EdgeColor','none','FaceAlpha',0.2);
% endAlignPloth=subplot(1,2,2);hold on; box off;
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
OptoStimh=patch([repmat(alignmtt-(start+3*conv_sigma),1,2) repmat(alignmtt-(start+3*conv_sigma)+Trials.end(1,2)-Trials.start(1,2),1,2)], ...
    [[0 currylim(2)] fliplr([0 currylim(2)])], ...
    [0 0 0 0],[0.3 0.75 0.93],'EdgeColor','none','FaceAlpha',0.5);
% "regular" alignement bar
% patch([repmat((alignmtt-(start+3*conv_sigma))-2,1,2) repmat((alignmtt-(start+3*conv_sigma))+2,1,2)], ...
%     [[0 currylim(2)] fliplr([0 currylim(2)])], ...
%     [0 0 0 0],[1 0 0],'EdgeColor','none','FaceAlpha',0.5);
% axes(endAlignPloth)
% patch([repmat((alignmtt-(start+3*conv_sigma))-2,1,2) repmat((alignmtt-(start+3*conv_sigma))+2,1,2)], ...
%     [[0 currylim(2)] fliplr([0 currylim(2)])], ...
%     [0 0 0 0],[1 0 0],'EdgeColor','none','FaceAlpha',0.5);

%legend
legend([FRploth,OptoStimh],{'Average firing rate','Optical stimulation'});
legend('boxoff')
% text(xTickSteps,currylim(2)-20,['Channel ' num2str(BestChan)],'FontName','Cambria');
title(['Channel ' num2str(BestChan)],'FontName','Cambria');