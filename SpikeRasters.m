function SpikeRasters

% KeepChans=10;
%% Get spike times and waveforms
[spikeData,syncTrials]=GetSpikeData;

%% get whisker time indices
[Sweep,HeadPos,behaviorEvents]=GetContactTime; % trial times from behavioral events are in ms

% first match sync trials with behavioral trials
[ephysCommonTrials, behaviorCommonTrials]=MatchTrials(syncTrials,behaviorEvents);
% very little exact agreement so far, lots of close ones, need to debug this 

%keep good trials
syncTrials.start=syncTrials.start(ephysCommonTrials);
syncTrials.end=syncTrials.end(ephysCommonTrials);

%set sync trial to correct start time
syncTrials.start=syncTrials.start-syncTrials.startClockTime;
syncTrials.end=syncTrials.end-syncTrials.startClockTime;

%downsample to ms precision
samplingRate=30000;
dsDivider=samplingRate/1000;
% if Open Ephys, sync rate will be 30kHz -> downsample 
syncTrials.start=syncTrials.start/dsDivider;
syncTrials.end=syncTrials.end/dsDivider;
spikeData.SpikeTimes=spikeData.SpikeTimes/dsDivider;

if size(syncTrials.start,2)>size(syncTrials.start,1)
    syncTrials.start=syncTrials.start';
    syncTrials.end=syncTrials.end';
end

%% gather data from each neuron
% Rasters.units=cell(length(KeepChans));
Rasters.epochnames={'BeginTrial','EndTrial'};
preAlignWindow=500; 
postAlignWindow=499; 

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
    [Rasters.units(unitNum).start,Rasters.units(unitNum).end]=...
        deal(zeros(size(syncTrials.start,1),preAlignWindow+postAlignWindow+1));
    for trialnb=1:size(syncTrials.start,1)
        %Collect spikes from 1st epoch (begining of trial)
        epochWin=[syncTrials.start(trialnb)-preAlignWindow,syncTrials.start(trialnb)+postAlignWindow];
        Rasters.units(unitNum).start(trialnb,spikeTimes(spikeTimes>epochWin(1) & spikeTimes<epochWin(2))-epochWin(1))=1;
        %Collect spikes from 2nd epoch (end of trial)
        epochWin=[syncTrials.end(trialnb)-preAlignWindow,syncTrials.end(trialnb)+postAlignWindow];
        Rasters.units(unitNum).end(trialnb,spikeTimes(spikeTimes>epochWin(1) & spikeTimes<epochWin(2))-epochWin(1))=1;
    end
end

%% plot raster showing all channels
figure('Position',[1050 120 750 790]);
colormap bone;
MeanChan=cellfun(@(x) conv_raster(x),{Rasters.units.start},'UniformOutput',false);
MeanChan=cell2mat(MeanChan');
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
