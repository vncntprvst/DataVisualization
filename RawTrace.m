%% Display raw (filtered) traces at multiple magnifications
% Overlay condition (e.g., stimulation)

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
    load([fileName{:} '_spikes.mat']);
else
    % add code for .dat?
end

KeepChans=[1,4,9,15,16];%1:5;
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

%% plot raw data - all channels, full recording
if isfield(Trials,'start')
    Trials.start=Trials.start-Trials.startClockTime;
    Trials.end=Trials.end-Trials.startClockTime;
else
    [Trials.start,Trials.end]=deal([]);
end
    
figure; hold on;
cmap=colormap;
for ChN=1:size(rawData,1)
    plot(double(rawData(ChN,:))+double(max(abs(max(rawData)))*(ChN-1))-...
        mean(mean(rawData(ChN,:))));
end
yLims=get(gca,'ylim');
for TTLNum=1:size(Trials.start,1)
    patch([Trials.start(TTLNum,2)*30:Trials.end(TTLNum,2)*30,...
        fliplr(Trials.start(TTLNum,2)*30:Trials.end(TTLNum,2)*30)],...
        [ones(1,Trials.end(TTLNum,2)*30-Trials.start(TTLNum,2)*30+1)*yLims(1),...
        ones(1,Trials.end(TTLNum,2)*30-Trials.start(TTLNum,2)*30+1)*yLims(2)],...
        [0.3 0.75 0.93],'EdgeColor','none','FaceAlpha',0.3);
    patch([Trials.start(TTLNum,2)*30:Trials.end(TTLNum,2)*30,...
    fliplr(Trials.start(TTLNum,2)*30:Trials.end(TTLNum,2)*30)],...
    [ones(1,Trials.end(TTLNum,2)*30-Trials.start(TTLNum,2)*30+1)*yLims(1),...
    ones(1,Trials.end(TTLNum,2)*30-Trials.start(TTLNum,2)*30+1)*(yLims(1)+100)],...
    [0 0 0],'EdgeColor','none','FaceAlpha',0.5);
end
tickNum=round((size(rawData,2)/rec_info.samplingRate)/60)*6; % ~ every 10sec
set(gca,'xtick',linspace(0,size(rawData,2),tickNum),...
    'xticklabel',round(linspace(0,size(rawData,2)/rec_info.samplingRate,tickNum)),'TickDir','out');
set(gca,'ytick',linspace(0,double(max(abs(max(rawData)))*(ChN-1))-...
        mean(mean(rawData(ChN,:))),size(rawData,1)),'yticklabel',...
    cellfun(axis_name, num2cell(1:size(rawData,1)), 'UniformOutput', false))
axis('tight');box off;
xlabel('Time (s.)')
ylabel('Recording Channels')
set(gca,'Color','white','FontSize',12,'FontName','calibri');
cmap=colormap;

%% plot best channel with "trials"
figure; hold on;
cmap=colormap(lines);
plot(rawData(5,:),'Color',cmap(1,:));
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

%same, with parameters to zoom in
figure; hold on;
colormap lines; cmap=colormap;
plot(rawData(3,:),'k','LineWidth',1.5);
% plot(rawData(3,:),'Color',cmap(3,:),'LineWidth',2);
% plot(rawData(3,:),'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',0.1);
yLims=get(gca,'ylim');
yLims(1)=max(rawData(3,:))/2;
for TTLNum=1:size(Trials.start,1)
    patch([Trials.start(TTLNum,2)*30:Trials.end(TTLNum,2)*30,...
        fliplr(Trials.start(TTLNum,2)*30:Trials.end(TTLNum,2)*30)],...
        [ones(1,Trials.end(TTLNum,2)*30-Trials.start(TTLNum,2)*30+1)*yLims(1),...
        ones(1,Trials.end(TTLNum,2)*30-Trials.start(TTLNum,2)*30+1)*yLims(2)],...
        [0.3 0.75 0.93],'EdgeColor','none','FaceAlpha',0.8);
end
tickNum=size(rawData,2)/3000; %every 10ms
set(gca,'xtick',linspace(0,size(rawData,2),tickNum),...
    'xticklabel',round(linspace(0,size(rawData,2)/(rec_info.samplingRate/1000),tickNum)),'TickDir','out');
axis('tight');box off;
xlabel('Time (ms)')
ylabel('Voltage (uV)')
set(gca,'Color','white','FontSize',12,'FontName','calibri');
