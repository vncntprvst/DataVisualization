% Spike-triggered image average
function spikeTrigVideo(spikeTimes,videoFile,vDir)
% random times
% spikeTimes=unique(round(10*(sort(rand(10,1)))).*round(17*(sort(rand(10,1)))))+1;

if ~exist('videoFile','var')
    % fileName='PrV77_52_HSCam2016-03-17T19_08_11'; 
[videoFile,vDir] = uigetfile({'*.avi','AVI files';'*.*','All Files' },...
    'Video files','E:\Data\Video');
% videoFile=[vDir videoFile];
end

videoObj = VideoReader([vDir videoFile]);
% 	numberOfFrames = videoFile.NumberOfFrames;
% 	vidHeight = videoObj.Height;
% 	vidWidth = videoObj.Width;
%     
% to read all video
% allVid = videoObj.read();
% foo=allVid(:,:,:,spikeTimes);

% % to read just specific frames
% for frNum=1:length(spikeTimes)
%     videoObj.CurrentTime = spikeTimes(frNum)/1000;
%     vidFrames(:,:,frNum) = rgb2gray(readFrame(videoObj));
% %     figure; imshow(vidFrame(:,:,frNum))
% %     imshow(grayImage, 'Parent', currAxes);
% %     currAxes.Visible = 'off';
% %     pause(1/videoFile.FrameRate);
% end

frNum = 1;
while hasFrame(videoObj)
    vidFrames(:,:,frNum) = rgb2gray(readFrame(videoObj));
    frNum = frNum+1;
end

trigFrames=vidFrames(:,:,spikeTimes(1:4));
bgFrame=median(vidFrames,3);
% figure; imshow(bgFrame)
% bgsub_trigFrames = imsubtract(trigFrames,bgFrame);

avFrame=median(trigFrames,3);
% avFrame=imsubtract(avFrame,bgFrame);
diffFrame = imabsdiff(avFrame,bgFrame);
edges = edge(diffFrame,'Sobel'); %or Prewitt , or Roberts
figure; imshowpair(diffFrame,edges,'montage')