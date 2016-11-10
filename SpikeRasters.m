function SpikeRasters

% KeepChans=10;
%% Get spike times and waveforms
[spikeData,syncTrials,fileInfo]=GetSpikeData;
%set sync trial to correct start time
syncTrials.start=syncTrials.start-syncTrials.startClockTime;
syncTrials.end=syncTrials.end-syncTrials.startClockTime;

%downsample to ms precision
samplingRate=uint64(fileInfo.samplingRate);
dsDivider=samplingRate/1000;
if strcmp(fileInfo.sys,'OEph')
% if Open Ephys, sync rate will be 30kHz -> downsample to 1kHz
syncTrials.start=syncTrials.start/dsDivider;
syncTrials.end=syncTrials.end/dsDivider;
end

%% get whisker time indices
[Sweep,HeadPos,behaviorEvents,SweepTrialIdx,exportVDir]=GetContactTime; % trial times from behavioral events are in ms
SweepTrialIdx=SweepTrialIdx+1;
% exportVDir='E:\Data\Video\PrV77_56_HSCamClips';
% first match sync trials with behavioral trials
[ephysCommonTrials, behaviorCommonTrials]=MatchTrials(syncTrials,behaviorEvents);
% very little exact agreement so far, lots of close ones, need to debug this 

%keep good trials
syncTrials.start=syncTrials.start(ephysCommonTrials);
syncTrials.end=syncTrials.end(ephysCommonTrials);

spikeData.SpikeTimes=spikeData.SpikeTimes/dsDivider;

if size(syncTrials.start,2)>size(syncTrials.start,1)
    syncTrials.start=syncTrials.start';
    syncTrials.end=syncTrials.end';
end

%% gather data from each neuron
% Rasters=cell(length(KeepChans));
% Rasters.epochnames={'BeginTrial','EndTrial'};
preAlignWindow=100; 
postAlignWindow=200; 
Rasters=struct('startTrial',[],'endTrial',[],'contactTimes',[]);
units=unique(spikeData.Units);
for unitNum=1:length(units)  
    spikeTimes=spikeData.SpikeTimes(spikeData.Units==units(unitNum));
% figure; hold on
% plot(syncTrials.start,ones(length(syncTrials.start),1)*0.6,'*')
% plot(spikeTimes,ones(length(spikeTimes),1)*0.5,'d')
% set(gca,'ylim',[0 1])

    if syncTrials.end(end)>spikeTimes(end)
        continue
    end
    [Rasters(unitNum).startTrial,Rasters(unitNum).endTrial]=...
        deal(zeros(size(syncTrials.start,1),preAlignWindow+postAlignWindow+1));
    Rasters(unitNum).sweepContactTimes=cell(size(syncTrials.start,1),1);
    Rasters(unitNum).panelContactTimes=cell(size(syncTrials.start,1),1);
    for trialnb=1:size(Sweep,2)
        %Collect spikes from 1st epoch (begining of trial)
        epochWin=[syncTrials.start(trialnb)-preAlignWindow,syncTrials.start(trialnb)+postAlignWindow];
        Rasters(unitNum).startTrial(trialnb,spikeTimes(spikeTimes>epochWin(1) & spikeTimes<epochWin(2))-epochWin(1))=1;
        
        % collect spikes for each sweep/panel contact, for the duration of
        % the contact
        if ismember(trialnb,SweepTrialIdx)
            contactTimes=Sweep(SweepTrialIdx==trialnb).Contact;
            %             contactTimes=sort(unique(contactTimes(~isnan(contactTimes))));
            if ~isempty(contactTimes)
                contactTimes=sort(unique(vertcat(contactTimes{:})));
                %             for sweepNum=1:size(contactTimes,1)
                %                 foo{sweepNum,1}=contactTimes(sweepNum,:);
                %                 foo{sweepNum,1}=foo{sweepNum}(~isnan(foo{sweepNum}));
                %             end
                for sweepNum=1:size(contactTimes,1)
                    try
                        swc_epochWin=epochWin+contactTimes(sweepNum);
                        Rasters(unitNum).sweepContactTimes{trialnb}=zeros(sweepNum,preAlignWindow+postAlignWindow+1);
                        Rasters(unitNum).sweepContactTimes{trialnb}(sweepNum,...
                            spikeTimes(spikeTimes>swc_epochWin(1) & spikeTimes<swc_epochWin(2))-swc_epochWin(1))=1;
                    catch
                        Rasters(unitNum).sweepContactTimes{trialnb}=[];
                    end
                end
            end
        end
        
        % collect spikes for each panel contact
        panelContacts=find(Sweep(trialnb).panelTouchIdx); %simple ROI activation idx
        for contactNum=1:size(panelContacts,1)
            try
                pnc_epochWin=epochWin+panelContacts(contactNum);
                Rasters(unitNum).panelContactTimes{trialnb}=zeros(contactNum,preAlignWindow+postAlignWindow+1);
                Rasters(unitNum).panelContactTimes{trialnb}(contactNum,spikeTimes(spikeTimes>pnc_epochWin(1) & spikeTimes<pnc_epochWin(2))-pnc_epochWin(1))=1;
            catch
                Rasters(unitNum).panelContactTimes{trialnb}=[];
            end
        end
        
        %Collect spikes from 2nd epoch (end of trial)
        epochWin=[syncTrials.end(trialnb)-preAlignWindow,syncTrials.end(trialnb)+postAlignWindow];
        Rasters(unitNum).endTrial(trialnb,spikeTimes(spikeTimes>epochWin(1) & spikeTimes<epochWin(2))-epochWin(1))=1;
    end
end

%% plot peri-start trial raster showing all units
figure('Position',[1050 120 750 790]);
colormap bone;
meanUnits=cellfun(@(x) conv_raster(x),{Rasters.startTrial},'UniformOutput',false);
meanUnits=cell2mat(meanUnits');
imagesc(zscore(meanUnits,[])); %
% imagesc(meanUnits);
xlabel('Time (ms)');
ylabel('Units','FontWeight','bold','FontSize',12);
% draw alignment bar
currylim=get(gca,'YLim');
currxlim=get(gca,'XLim');
midl=round(currxlim(2)/20)*10;
set(gca,'XTick',round([midl-preAlignWindow/2 midl midl+postAlignWindow/2]));
set(gca,'XTickLabel',round([-preAlignWindow/2 0 postAlignWindow/2]));
%alignment bar
midl=uint64(midl);
patch([repmat(midl-2,1,2) repmat(midl+2,1,2)], ...
    [[0 currylim(2)] fliplr([0 currylim(2)])], ...
    [0 0 0 0],[1 0 0],'EdgeColor','none','FaceAlpha',0.5);
title('peri-start trial responses');
hcb = colorbar('southoutside');
hcb.Label.String = 'z-scored firing rate';
%% plot peri-end trial raster showing all units
figure('Position',[1050 120 750 790]);
colormap bone;
meanUnits=cellfun(@(x) conv_raster(x),{Rasters.endTrial},'UniformOutput',false);
meanUnits=cell2mat(meanUnits');
imagesc(zscore(meanUnits,[])); %
% imagesc(meanUnits);
xlabel('Time (ms)');
ylabel('Units','FontWeight','bold','FontSize',12);
% draw alignment bar
currylim=get(gca,'YLim');
currxlim=get(gca,'XLim');
midl=round(currxlim(2)/20)*10;
set(gca,'XTick',round([midl-preAlignWindow/2 midl midl+postAlignWindow/2]));
set(gca,'XTickLabel',round([-preAlignWindow/2 0 postAlignWindow/2]));
%alignment bar
midl=uint64(midl);
patch([repmat(midl-2,1,2) repmat(midl+2,1,2)], ...
    [[0 currylim(2)] fliplr([0 currylim(2)])], ...
    [0 0 0 0],[1 0 0],'EdgeColor','none','FaceAlpha',0.5);
title('peri-end trial responses');
hcb = colorbar('southoutside');
hcb.Label.String = 'z-scored firing rate';
%% plot spikes aligned to sweep contact
cmap=colormap(jet);
for unitNum=1:length(units)
    touchRasters=Rasters(unitNum).sweepContactTimes;
    if ~isempty(touchRasters)
        trialType=behaviorEvents.trials.trialType(~cellfun('isempty',touchRasters));
        touchRasters=touchRasters(~cellfun('isempty',touchRasters));
        figure('Position',[1050 120 750 790]);
        colormap jet;
        subplot(3,1,1);
        [meanUnitTouch,~,meanUnitTouchSEM]=cellfun(@(x) conv_raster(x),touchRasters(trialType==1),'UniformOutput',false);
        meanUnitTouchSEM=meanUnitTouchSEM(cellfun(@(x) size(x,2)>1, meanUnitTouch));
        meanUnitTouch=meanUnitTouch(cellfun(@(x) size(x,2)>1, meanUnitTouch));
        meanUnitTouch=cell2mat(meanUnitTouch);
        imagesc(meanUnitTouch);
%         imagesc(zscore(meanUnitTouch,[]));
        xlabel('Time (ms)');
        ylabel('Units','FontWeight','bold','FontSize',12);
        % draw alignment bar
        currylim=get(gca,'YLim');
        currxlim=get(gca,'XLim');
        midl=preAlignWindow-(preAlignWindow+postAlignWindow-currxlim(2))/2;
        set(gca,'XTick',round([midl-preAlignWindow/2 midl midl+postAlignWindow/2]));
        set(gca,'XTickLabel',round([-preAlignWindow/2 0 postAlignWindow/2]));
        %alignment bar
        midl=uint64(midl);
        patch([repmat(midl-2,1,2) repmat(midl+2,1,2)], ...
            [[0 currylim(2)] fliplr([0 currylim(2)])], ...
            [0 0 0 0],[1 0 0],'EdgeColor','none','FaceAlpha',0.5);
        title('sweep-contact responses, texture');
%         hcb = colorbar('southoutside');
%         hcb.Label.String = 'z-scored firing rate';
        subplot(3,1,3); hold on 
        plot(nanmean(meanUnitTouch));
        %plot confidence intervals
        meanUnitTouchSEM=cell2mat(meanUnitTouchSEM);
        patch([1:size(meanUnitTouch,2),fliplr(1:size(meanUnitTouch,2))],...
        [nanmean(meanUnitTouch)-nanmean(meanUnitTouchSEM),...
        fliplr(nanmean(meanUnitTouch)+nanmean(meanUnitTouchSEM))],cmap(1,:),'EdgeColor','none','FaceAlpha',0.1);
        subplot(3,1,2);
        [meanUnitTouch,~,meanUnitTouchSEM]=cellfun(@(x) conv_raster(x),touchRasters(trialType==2),'UniformOutput',false);
        meanUnitTouchSEM=meanUnitTouchSEM(cellfun(@(x) size(x,2)>1, meanUnitTouch));
        meanUnitTouch=meanUnitTouch(cellfun(@(x) size(x,2)>1, meanUnitTouch));
        meanUnitTouch=cell2mat(meanUnitTouch);
        imagesc(meanUnitTouch);
%         imagesc(zscore(meanUnitTouch,[]));
        xlabel('Time (ms)');
        ylabel('Units','FontWeight','bold','FontSize',12);
        % draw alignment bar
        currylim=get(gca,'YLim');
        currxlim=get(gca,'XLim');
        midl=preAlignWindow-(preAlignWindow+postAlignWindow-currxlim(2))/2;
        set(gca,'XTick',round([midl-preAlignWindow/2 midl midl+postAlignWindow/2]));
        set(gca,'XTickLabel',round([-preAlignWindow/2 0 postAlignWindow/2]));
        %alignment bar
        midl=uint64(midl);
        patch([repmat(midl-2,1,2) repmat(midl+2,1,2)], ...
            [[0 currylim(2)] fliplr([0 currylim(2)])], ...
            [0 0 0 0],[1 0 0],'EdgeColor','none','FaceAlpha',0.5);
        title('sweep-contact responses, no texture');
%         hcb = colorbar('southoutside');
%         hcb.Label.String = 'z-scored firing rate';
        subplot(3,1,3);
        plot(mean(meanUnitTouch));
        %plot confidence intervals
        meanUnitTouchSEM=cell2mat(meanUnitTouchSEM);
        patch([1:size(meanUnitTouch,2),fliplr(1:size(meanUnitTouch,2))],...
        [nanmean(meanUnitTouch)-nanmean(meanUnitTouchSEM),...
        fliplr(nanmean(meanUnitTouch)+nanmean(meanUnitTouchSEM))],cmap(2,:),'EdgeColor','none','FaceAlpha',0.1);
        legend('texture','no texture');
    end
end

%% plot spikes aligned to panel contact
for unitNum=1:length(units)
    touchRasters=Rasters(unitNum).panelContactTimes;
    if ~isempty(touchRasters)
        trialType=behaviorEvents.trials.trialType(~cellfun('isempty',touchRasters));
        touchRasters=touchRasters(~cellfun('isempty',touchRasters));
        figure('Position',[1050 120 750 790]);
        colormap jet;
        subplot(3,1,1);
        [meanUnitTouch,~,meanUnitTouchSEM]=cellfun(@(x) conv_raster(x),touchRasters(trialType==1),'UniformOutput',false);
        meanUnitTouchSEM=meanUnitTouchSEM(cellfun(@(x) size(x,2)>1, meanUnitTouch));
        meanUnitTouch=meanUnitTouch(cellfun(@(x) size(x,2)>1, meanUnitTouch));
        meanUnitTouch=cell2mat(meanUnitTouch);
        imagesc(meanUnitTouch);
%         imagesc(zscore(meanUnitTouch,[]));
        xlabel('Time (ms)');
        ylabel('Units','FontWeight','bold','FontSize',12);
        % draw alignment bar
        currylim=get(gca,'YLim');
        currxlim=get(gca,'XLim');
        midl=preAlignWindow-(preAlignWindow+postAlignWindow-currxlim(2))/2;
        set(gca,'XTick',round([midl-preAlignWindow/2 midl midl+postAlignWindow/2]));
        set(gca,'XTickLabel',round([-preAlignWindow/2 0 postAlignWindow/2]));
        %alignment bar
        midl=uint64(midl);
        patch([repmat(midl-2,1,2) repmat(midl+2,1,2)], ...
            [[0 currylim(2)] fliplr([0 currylim(2)])], ...
            [0 0 0 0],[1 0 0],'EdgeColor','none','FaceAlpha',0.5);
        title('sweep-contact responses, texture');
%         hcb = colorbar('southoutside');
%         hcb.Label.String = 'z-scored firing rate';
        subplot(3,1,3); hold on 
        plot(nanmean(meanUnitTouch));
        %plot confidence intervals
        meanUnitTouchSEM=cell2mat(meanUnitTouchSEM);
        patch([1:size(meanUnitTouch,2),fliplr(1:size(meanUnitTouch,2))],...
        [nanmean(meanUnitTouch)-nanmean(meanUnitTouchSEM),...
        fliplr(nanmean(meanUnitTouch)+nanmean(meanUnitTouchSEM))],cmap(1,:),'EdgeColor','none','FaceAlpha',0.1);
        subplot(3,1,2);
        [meanUnitTouch,~,meanUnitTouchSEM]=cellfun(@(x) conv_raster(x),touchRasters(trialType==2),'UniformOutput',false);
        meanUnitTouchSEM=meanUnitTouchSEM(cellfun(@(x) size(x,2)>1, meanUnitTouch));
        meanUnitTouch=meanUnitTouch(cellfun(@(x) size(x,2)>1, meanUnitTouch));
        meanUnitTouch=cell2mat(meanUnitTouch);
        imagesc(meanUnitTouch);
%         imagesc(zscore(meanUnitTouch,[]));
        xlabel('Time (ms)');
        ylabel('Units','FontWeight','bold','FontSize',12);
        % draw alignment bar
        currylim=get(gca,'YLim');
        currxlim=get(gca,'XLim');
        midl=preAlignWindow-(preAlignWindow+postAlignWindow-currxlim(2))/2;
        set(gca,'XTick',round([midl-preAlignWindow/2 midl midl+postAlignWindow/2]));
        set(gca,'XTickLabel',round([-preAlignWindow/2 0 postAlignWindow/2]));
        %alignment bar
        midl=uint64(midl);
        patch([repmat(midl-2,1,2) repmat(midl+2,1,2)], ...
            [[0 currylim(2)] fliplr([0 currylim(2)])], ...
            [0 0 0 0],[1 0 0],'EdgeColor','none','FaceAlpha',0.5);
        title('sweep-contact responses, no texture');
%         hcb = colorbar('southoutside');
%         hcb.Label.String = 'z-scored firing rate';
        subplot(3,1,3);
        plot(mean(meanUnitTouch));
        %plot confidence intervals
        meanUnitTouchSEM=cell2mat(meanUnitTouchSEM);
        patch([1:size(meanUnitTouch,2),fliplr(1:size(meanUnitTouch,2))],...
        [nanmean(meanUnitTouch)-nanmean(meanUnitTouchSEM),...
        fliplr(nanmean(meanUnitTouch)+nanmean(meanUnitTouchSEM))],cmap(2,:),'EdgeColor','none','FaceAlpha',0.1);
        legend('texture','no texture');
    end
end
