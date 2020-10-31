% Make video with spike display 

% vIRt44_1210_5450
% Bout19
% x20
% Unit37
% 
clearvars

cd Z:\Vincent\Ephys\vIRt44\vIRt44_1210\Analysis\vIRt44_1210_5450
baseName = 'vIRt44_1210_5450';

dataDir=cd;
load(fullfile(dataDir,[baseName '_ephys.mat']));
load(fullfile(dataDir,[baseName '_whisker.mat']));
syncFile = fopen(fullfile(dataDir,[baseName '_vSyncTTLs.dat']), 'r');
vidTimes = fread(syncFile,'single');%'int32' % VideoFrameTimes: was fread(fid,[2,Inf],'double'); Adjust export
fclose(syncFile);

%% find whisking epochs (based on first trace, if multiple whisker tracked)
ampThd=18; %12; %18 %amplitude threshold 
freqThld=1; %frequency threshold
minBoutDur=1000; %500; % 1000 % minimum whisking bout duration: 1s
whiskingEpochs=WhiskingFun.FindWhiskingEpochs(...
   whisker.Amplitude(1,:),whisker.Frequency(1,:),...
   ampThd, freqThld, minBoutDur); 
whiskingEpochs(isnan(whiskingEpochs))=false; %just in case
whiskingEpochsList=bwconncomp(whiskingEpochs);
[~,wBoutDurSort]=sort(cellfun(@length,whiskingEpochsList.PixelIdxList),'descend');
whiskingEpochsList.PixelIdxListSorted=whiskingEpochsList.PixelIdxList(wBoutDurSort);

        
%% make video of whisking bouts
boutNum=19 ; %4 ; %19 7
cellNum=37; %54; %37 26 
traceIndex=whiskingEpochsList.PixelIdxList{boutNum}; %350000:352000;%whiskingEpochsList.PixelIdxList{2}
[boutFrames,frameIndex]=WhiskingBoutVideo(ephys.recInfo.likelyVideoFile,ephys.recInfo.dirName,...
    traceIndex,vidTimes,false);
%remove any drift by adding frames when fps > intended fps  (crude)
boutTimeLine=whisker.vidTimes(frameIndex);boutTimeLine=boutTimeLine-boutTimeLine(1);
extraFrameIdx=find(diff(round(boutTimeLine))>mode(diff(round(boutTimeLine))));  
for frameIdx=1:numel(extraFrameIdx)
boutFrames=[boutFrames(1:extraFrameIdx(frameIdx)+frameIdx-1),...
    boutFrames(extraFrameIdx(frameIdx)+frameIdx-1),...
    boutFrames(extraFrameIdx(frameIdx)+frameIdx:end)];
end
boutFrames=boutFrames(1:end-numel(extraFrameIdx));
% vidDims=size(boutFrames(1).cdata);
% figure('position',[1500 450  vidDims(2) vidDims(1)],'color','k');
% movie(boutFrames,1,500);

% with trace added
% add audiovidDims=size(boutFrames(1).cdata);
spikeTimes = movmean(ephys.rasters(cellNum,traceIndex),2);
spikeTimes = logical(spikeTimes(1:2:end));
figure('position',[500 450  vidDims(2) vidDims(1)],'color','k');
boutFrames=FrameByFrame_Overlay(boutFrames,[whisker.Phase(traceIndex(1:2:end));spikeTimes]); %whisker.Angle

%add audio for a given cell, given speed
slowFactor=20; %500/25

% Code perso
% wBoutAudio=WhiskingBoutAudio(ephys.rasters(cellNum,traceIndex(1:5000)),1000,20*slowFactor);
% FacePro
wBoutAudio = FacePro.MakeSpikeAudio(ephys.rasters(cellNum,traceIndex(1:5000)),...
    slowFactor*10, [1 2500], 1000);
% 
% samplingRatio=round(numel(traceIndex)/numel(boutFrames));
% wBoutAudio = zeros(200,numel(boutFrames));
% spikeTimes = round(find(ephys.rasters(cellNum,traceIndex))/samplingRatio);
% waveforms = gausswin(20); 
% 
% for spikeNum = 1 : numel(spikeTimes)
% %     wBoutAudio(spikeTimes(spikeNum)-4:spikeTimes(spikeNum)+size(waveforms,1)-5) = waveforms;
%     wBoutAudio(1:20,spikeTimes(spikeNum)) = waveforms;
% end
% wBoutAudio = int16(wBoutAudio / max(abs(wBoutAudio(:))) * double(intmax('int16')));

% figure; plot(wBoutAudio)
% figure; imagesc(wBoutAudio)

%write video
frameRate=500/slowFactor;
videoFWriter = vision.VideoFileWriter(fullfile(cd,...
    [ephys.recInfo.likelyVideoFile(1:end-4) '_Bout' num2str(boutNum)...
    'x' num2str(slowFactor) '_Unit' num2str(cellNum) '_FP_PhaseTuning.avi']));
videoFWriter.FrameRate =frameRate ;
videoFWriter.AudioInputPort = true;
videoFWriter.VideoCompressor = 'None (uncompressed)'; % 'MJPEG Compressor';

for frameNum = 1 : 2500 %numel(boutFrames)
%     videoFWriter(boutFrames(frameNum).cdata, wBoutAudio(:,(frameNum*2)-1));  %if double the number of frames
    videoFWriter(boutFrames(frameNum).cdata, wBoutAudio(:,frameNum)); %if matrix form factor
%     videoFWriter(boutFrames(frameNum).cdata, wBoutAudio(frameNum));
end
release(videoFWriter);




