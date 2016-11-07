% Spike-triggered image average

if ~exist('videoFile','var')
    % fileName='PrV77_52_HSCam2016-03-17T19_08_11'; 
[videoFile,vDir] = uigetfile({'*.avi','AVI files';'*.*','All Files' },...
    'Video files','E:\Data\Video');
% videoFile=[vDir videoFile];
end

videoFile = VideoReader([vDir videoFile]);
% 	numberOfFrames = videoFile.NumberOfFrames;
	vidHeight = videoFile.Height;
	vidWidth = videoFile.Width;
    
spikeTimes=unique(round(10*(sort(rand(10,1)))).*round(17*(sort(rand(10,1)))))+1;
% video = videoFile.read();
% foo=video(:,:,:,spikeTimes);

% currAxes = axes;
for frNum=1:length(spikeTimes)
    videoFile.CurrentTime = spikeTimes(frNum)/1000;
    vidFrame(:,:,frNum) = rgb2gray(readFrame(videoFile));
    figure; imshow(vidFrame(:,:,frNum))
%     imshow(grayImage, 'Parent', currAxes);
%     currAxes.Visible = 'off';
%     pause(1/videoFile.FrameRate);
end

avFrame=median(vidFrame,3);
figure; imshow(avFrame)