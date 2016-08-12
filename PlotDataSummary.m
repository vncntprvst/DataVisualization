%% Plot "raw data" and aligned data

% function declaration
axis_name= @(x) sprintf('Chan %.0f',x);

rec_info.samplingRate=30000; %in future load saved info

%% Get file path
[fileName,dirName] = uigetfile({'*.mat; *.hdf5','Processed data';'*.dat','Flat data';...
    '*.*','All Files' },'Exported data','C:\Data\export');
cd(dirName);

%% load file data
if strfind(fileName,'.mat')
    fileName=regexp(fileName,'.+(?=_\w+.\w+$)','match');
    load([fileName{:} '_raw.mat']);
    load([fileName{:} '_trials.mat']);
    %     load([fileName{:} '_spikes.mat']);
    load([fileName{:} '_spikesResorted.mat']);
else
    % add code for .dat?
end

%% plot raw data excerpt
TenSecSample=rawData(:,round(size(rawData,2)/2)-(5*rec_info.samplingRate):round(size(rawData,2)/2)+(5*rec_info.samplingRate));
figure; hold on;
for ChN=1:size(rawData,1)
    plot(double(TenSecSample(ChN,:))+double(max(abs(max(TenSecSample)))*(ChN-1))-...
        mean(mean(TenSecSample(ChN,:))));
end
set(gca,'xtick',linspace(0,rec_info.samplingRate*10,4),...
    'xticklabel',round(linspace(round((round(size(rawData,2)/2)-(5*rec_info.samplingRate))/rec_info.samplingRate),...
    round((round(size(rawData,2)/2)+(5*rec_info.samplingRate))/rec_info.samplingRate),4)),'TickDir','out');
set(gca,'ytick',linspace(0,double(max(abs(max(TenSecSample)))*(ChN-1))-...
    mean(mean(TenSecSample(ChN,:))),size(rawData,1)),'yticklabel',...
    cellfun(axis_name, num2cell(1:size(rawData,1)), 'UniformOutput', false))
%   set(gca,gca,'ylim',[-1000,10000],'xlim',[0,1800000])
axis('tight');box off;
xlabel('10 sec mid-recording')
ylabel('Raw signal')
set(gca,'Color','white','FontSize',12,'FontName','calibri');

%select channels class
prompt='Select channels to keep';
name='Channel selection';
numlines=1;
chStr= num2str(linspace(1,size(rawData,1),size(rawData,1))');
KeepChans=listdlg('Name',name,'PromptString',prompt,'ListString',chStr);

%close initial figure
close(gcf);

% KeepChans=[1,4,9,15,16];%1:5;
rawData=rawData(KeepChans,:);

%% plot raw data excerpt
TenSecSample=rawData(:,round(size(rawData,2)/2)-(5*rec_info.samplingRate):round(size(rawData,2)/2)+(5*rec_info.samplingRate));
figure; hold on;
for ChN=1:size(rawData,1)
    plot(double(TenSecSample(ChN,:))+double(max(abs(max(TenSecSample)))*(ChN-1))-...
        mean(mean(TenSecSample(ChN,:))));
end
set(gca,'xtick',linspace(0,rec_info.samplingRate*10,4),...
    'xticklabel',round(linspace(round((round(size(rawData,2)/2)-(5*rec_info.samplingRate))/rec_info.samplingRate),...
    round((round(size(rawData,2)/2)+(5*rec_info.samplingRate))/rec_info.samplingRate),4)),'TickDir','out');
set(gca,'ytick',linspace(0,double(max(abs(max(TenSecSample)))*(ChN-1))-...
    mean(mean(TenSecSample(ChN,:))),size(rawData,1)),'yticklabel',...
    cellfun(axis_name, num2cell(1:size(rawData,1)), 'UniformOutput', false))
%   set(gca,gca,'ylim',[-1000,10000],'xlim',[0,1800000])
axis('tight');box off;
xlabel('10 sec mid-recording')
ylabel('Raw signal')
set(gca,'Color','white','FontSize',12,'FontName','calibri');

%% get Spike times
for chNum=1:length(KeepChans)
    try
        units=unique(Spikes.HandSort.Units{KeepChans(chNum),1});units=units(units>0);
        for unitNum=1:length(units)
            numUnits(KeepChans(chNum),unitNum)=sum(Spikes.HandSort.Units{KeepChans(chNum),1}==units(unitNum));
            SpikeTimes{KeepChans(chNum),unitNum}=Spikes.HandSort.SpikeTimes{KeepChans(chNum),1}...
                (Spikes.HandSort.Units{KeepChans(chNum),1}==units(unitNum));
            Waveforms{KeepChans(chNum),unitNum}=Spikes.HandSort.Waveforms{KeepChans(chNum),1}...
                (:,Spikes.HandSort.Units{KeepChans(chNum),1}==units(unitNum));
            SpikeTimesArray{KeepChans(chNum),unitNum}=zeros(1,ceil(SpikeTimes{KeepChans(chNum),unitNum}(end)...
                /uint32(Spikes.Online_Sorting.samplingRate(KeepChans(chNum))/1000)));
            SpikeTimesArray{KeepChans(chNum),unitNum}(round(SpikeTimes{KeepChans(chNum),unitNum}/...
                uint32(Spikes.Online_Sorting.samplingRate(KeepChans(chNum))/1000)))=1;
        end
%         Waveforms{KeepChans(chNum)}=ExtractChunks(rawData(KeepChans(chNum),:),SpikeTimes{KeepChans(chNum)},60,'tmiddle'); %'tzero' 'tmiddle'
    end
end
numUnits=numUnits(KeepChans,:);
SpikeTimes=SpikeTimes(KeepChans,:);
SpikeTimesArray=SpikeTimesArray(KeepChans,:);
Waveforms=Waveforms(KeepChans,:);

%% adjust start time
if isfield(Trials,'start')
    Trials.start=Trials.start-Trials.startClockTime;
    Trials.end=Trials.end-Trials.startClockTime;
else
    [Trials.start,Trials.end]=deal([]);
end

%% plot raw data - all channels, full recording

% figure; hold on;
% cmap=colormap;
% for ChN=1:size(rawData,1)
%     plot(double(rawData(ChN,:))+double(max(abs(max(rawData)))*(ChN-1))-...
%         mean(mean(rawData(ChN,:))));
% end
% yLims=get(gca,'ylim');
% for TTLNum=1:size(Trials.start,1)
%     patch([Trials.start(TTLNum,2)*30:Trials.end(TTLNum,2)*30,...
%         fliplr(Trials.start(TTLNum,2)*30:Trials.end(TTLNum,2)*30)],...
%         [ones(1,Trials.end(TTLNum,2)*30-Trials.start(TTLNum,2)*30+1)*yLims(1),...
%         ones(1,Trials.end(TTLNum,2)*30-Trials.start(TTLNum,2)*30+1)*yLims(2)],...
%         [0.3 0.75 0.93],'EdgeColor','none','FaceAlpha',0.3);
%     patch([Trials.start(TTLNum,2)*30:Trials.end(TTLNum,2)*30,...
%     fliplr(Trials.start(TTLNum,2)*30:Trials.end(TTLNum,2)*30)],...
%     [ones(1,Trials.end(TTLNum,2)*30-Trials.start(TTLNum,2)*30+1)*yLims(1),...
%     ones(1,Trials.end(TTLNum,2)*30-Trials.start(TTLNum,2)*30+1)*(yLims(1)+100)],...
%     [0 0 0],'EdgeColor','none','FaceAlpha',0.5);
% end
% tickNum=round((size(rawData,2)/rec_info.samplingRate)/60)*6; % ~ every 10sec
% set(gca,'xtick',linspace(0,size(rawData,2),tickNum),...
%     'xticklabel',round(linspace(0,size(rawData,2)/rec_info.samplingRate,tickNum)),'TickDir','out');
% set(gca,'ytick',linspace(0,double(max(abs(max(rawData)))*(ChN-1))-...
%         mean(mean(rawData(ChN,:))),size(rawData,1)),'yticklabel',...
%     cellfun(axis_name, num2cell(1:size(rawData,1)), 'UniformOutput', false))
% axis('tight');box off;
% xlabel('Time (s.)')
% ylabel('Recording Channels')
% set(gca,'Color','white','FontSize',12,'FontName','calibri');
% cmap=colormap;

%% plot best channel with "trials"
bestChan=find(sum(numUnits,2)==max(sum(numUnits,2)));
figure; hold on;
cmap=colormap(lines);
plot(rawData(bestChan,:),'Color',cmap(1,:));
yLims=get(gca,'ylim');
for TTLNum=1:size(Trials.start,1)
    patch([Trials.start(TTLNum,2)*30:Trials.end(TTLNum,2)*30,...
        fliplr(Trials.start(TTLNum,2)*30:Trials.end(TTLNum,2)*30)],...
        [ones(1,Trials.end(TTLNum,2)*30-Trials.start(TTLNum,2)*30+1)*yLims(1),...
        ones(1,Trials.end(TTLNum,2)*30-Trials.start(TTLNum,2)*30+1)*yLims(2)],...
        [0.3 0.75 0.93],'EdgeColor','none','FaceAlpha',0.8);
end
set(gca,'xtick',linspace(0,size(rawData,2),tickNum),...
    'xticklabel',round(linspace(0,size(rawData,2)/rec_info.samplingRate,tickNum)),'TickDir','out');
% set(gca,'ytick',linspace(0,double(max(abs(max(rawData)))*(ChN-1))-...
%         mean(mean(rawData(ChN,:))),size(rawData,1)),'yticklabel',...
%     cellfun(axis_name, num2cell(1:size(rawData,1)), 'UniformOutput', false))
axis('tight');box off;
xlabel('Time (s.)')
ylabel('Voltage (uV)')
set(gca,'Color','white','FontSize',12,'FontName','calibri');

%% same, with parameters to zoom in
% figure; hold on;
% colormap lines; cmap=colormap;
% plot(rawData(3,:),'k','LineWidth',1.5);
% % plot(rawData(3,:),'Color',cmap(3,:),'LineWidth',2);
% % plot(rawData(3,:),'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',0.1);
% yLims=get(gca,'ylim');
% yLims(1)=max(rawData(3,:))/2;
% for TTLNum=1:size(Trials.start,1)
%     patch([Trials.start(TTLNum,2)*30:Trials.end(TTLNum,2)*30,...
%         fliplr(Trials.start(TTLNum,2)*30:Trials.end(TTLNum,2)*30)],...
%         [ones(1,Trials.end(TTLNum,2)*30-Trials.start(TTLNum,2)*30+1)*yLims(1),...
%         ones(1,Trials.end(TTLNum,2)*30-Trials.start(TTLNum,2)*30+1)*yLims(2)],...
%         [0.3 0.75 0.93],'EdgeColor','none','FaceAlpha',0.8);
% end
% tickNum=size(rawData,2)/3000; %every 10ms
% set(gca,'xtick',linspace(0,size(rawData,2),tickNum),...
%     'xticklabel',round(linspace(0,size(rawData,2)/(rec_info.samplingRate/1000),tickNum)),'TickDir','out');
% axis('tight');box off;
% xlabel('Time (ms)')
% ylabel('Voltage (uV)')
% set(gca,'Color','white','FontSize',12,'FontName','calibri');

%% Spikes plots
% gather data from selected channels
Rasters.channels=cell(length(KeepChans),2);
Rasters.epochnames={'BeginTrial','EndTrial'};
preAlignWindow(1)=0.2; %0.08 = 80ms
postAlignWindow(1)=0.5; 
for chan=1:length(KeepChans)
    preAlignWindow(2)=round(uint64(Spikes.HandSort.samplingRate(chan,2))/(1/preAlignWindow(1)));
    postAlignWindow(2)=round(uint64(Spikes.HandSort.samplingRate(chan,2))/(1/postAlignWindow(1)));

    %     downSamplingRatio=uint64(Spikes.HandSort.samplingRate(chan,1)/Spikes.HandSort.samplingRate(chan,2));
    [Rasters.channels{chan,1},Rasters.channels{chan,2}]=deal(zeros(size(Trials.start,1),preAlignWindow(2)+postAlignWindow(2)+1));
%     Spkt=Spikes.HandSort.data{KeepChans(chan),2}(1,:);
    mostUnits=numUnits(chan,:)==max(numUnits(chan,:));
    Spkt=SpikeTimesArray{chan,mostUnits};
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
    end
end
%% plot raster showing all channels
figure('Position',[1050 120 750 790]);
colormap bone;
MeanChan=cellfun(@(x) conv_raster(x),Rasters.channels(:,1),'UniformOutput',false);
MeanChan=cell2mat(MeanChan);
subplot(3,2,[1,3,5])
% imagesc(zscore(MeanChan,[],2));
imagesc(MeanChan);
xlabel('Time (ms)','FontWeight','bold','FontName','Cambria','FontSize',12);
ylabel('Channel','FontWeight','bold','FontName','Cambria','FontSize',12);
% draw alignment bar
currylim=get(gca,'YLim');
% currxlim=get(gca,'XLim');midl=round(currxlim(2)/2);
midl=preAlignWindow(2)-30;
set(gca,'XTick',[midl-100 midl midl+300]);
set(gca,'XTickLabel',[-100 0 300]);
set(gca,'YTick',1:size(MeanChan,1));
set(gca,'YTickLabel',num2str(KeepChans'));

patch([repmat(midl-3,1,2) repmat(midl+3,1,2)], ...
    [[0 currylim(2)] fliplr([0 currylim(2)])], ...
    [0 0 0 0],[0.8 0 0],'EdgeColor','none','FaceAlpha',0.8);
% title('Neural response to 80% stimulation intensity, aligned to stimulation onset');
title('Neural response, aligned to trial onset')
hcb = colorbar('northoutside');
hcb.Label.String = 'z-scored firing rate';

%% plot rasters
BestChan=find(mean(MeanChan,2)==max(mean(MeanChan,2)));
if size(BestChan,1)>1
    BestChan=BestChan(1);
end
start=1;
stop=size(Rasters.channels{BestChan,1},2);
crop_rasters = Rasters.channels{BestChan,1}(:,start:stop); % Isolate rasters of interest
[indy, indx] = ind2sub(size(crop_rasters),find(crop_rasters)); %find row and column coordinates of spikes
subplot(3,2,2)
plot([indx';indx'],[indy';indy'+1],'color',cmap(BestChan,:),...
    'LineStyle','-','LineWidth',1.8); % plot rasters
set(gca,'xlim',[1 length(start:stop)]);
axis(gca, 'off'); % axis tight sets the axis limits to the range of the data.


%% plot sdf
conv_sigma=10;
alignmtt=preAlignWindow(2);
xTickSteps=round(preAlignWindow(2)/50)*10;
[sdf{1}, ~, rastsem{1}]=conv_raster(Rasters.channels{BestChan,1},conv_sigma,start,stop);
% [sdf{2}, ~, rastsem{2}]=conv_raster(Rasters.channels{BestChan,2},conv_sigma,start,stop);
% figure('Position',[1469 542 417 417]);
subplot(3,2,[4,6])
% colormap default;
% cmap = colormap(gcf);
hold on;

%plot sem
startAlignPloth=gca; box off; %subplot(1,2,1);hold on; box off;
FRploth=plot(startAlignPloth,sdf{1},'Color',[0.16 0.38 0.27],'LineWidth',1.8);

set(startAlignPloth,'XTick',[midl-100 midl midl+300]);
set(startAlignPloth,'XTickLabel',[-100 0 300]);
axis(startAlignPloth,'tight');
set(startAlignPloth,'Color','white','TickDir','out','FontName','Cambria','FontSize',10);
hxlabel=xlabel(startAlignPloth,'Time (ms)','FontWeight','bold','FontName','Cambria','FontSize',12);
hylabel=ylabel(startAlignPloth,'Firing rate (spikes/s)','FontWeight','bold','FontName','Cambria','FontSize',12);

% draw alignment bar
currylim=get(startAlignPloth,'YLim');
axes(startAlignPloth)
patch([repmat((alignmtt-(start+3*conv_sigma))-2,1,2) repmat((alignmtt-(start+3*conv_sigma))+2,1,2)], ...
    [[0 currylim(2)] fliplr([0 currylim(2)])], ...
    [0 0 0 0],[1 0 0],'EdgeColor','none','FaceAlpha',0.5);
title(['Channel ' num2str(KeepChans(BestChan))],'FontName','Cambria');