datDir=['/mnt/data/Vincent/ephys/vIRt22/vIRt22_1016/'...
    'vIRt22_2018-10-16_17-00-06_3900_50ms1Hz5mW/'...
    'experiment1/recording1/continuous/Rhythm_FPGA-100.0/'];
datFile='continuous.dat';
traces = memmapfile(fullfile(datDir,datFile),'Format','int16');

%% recording info
rec_info.exportedChan=1:32; rec_info.samplingRate=30000;
fid  = fopen([fullfile(datDir,datFile(1:end-4)) '.params'],'r');
params=fread(fid,'*char')';
if contains(regexp(params,'(?<=filter_done\s+=\s)\w+(?=\s)','match','once'), 'True')
    tracesInfo.preproc=1;
else
    tracesInfo.preproc=0;
end
rec_info.threshold=regexp(params,'(?<=spike_thresh\s+=\s)\w+(?=\s)','match','once');
fclose(fid);

tracesInfo= struct('name','rawData',...
    'size',size(traces.Data),...
    'numChan',32,...
    'source','dat');
tracesInfo.excerptSize=rec_info.samplingRate/2;
tracesInfo.excerptLocation=round(max(tracesInfo.size)/...
    2/numel(rec_info.exportedChan));  %mid-recording as default
tracesInfo.preproc=0; % should not be filtered yet

%% plotting parameters
winIdxStart=((tracesInfo.excerptLocation-...
    double(tracesInfo.excerptSize))*numel(rec_info.exportedChan))+1;
if mod(winIdxStart,numel(rec_info.exportedChan))~=1 %set window index to correct point in data vector
    winIdxStart=winIdxStart-...
        mod(winIdxStart,numel(rec_info.exportedChan))-...
        numel(rec_info.exportedChan)+1;    % set index loc to first electrode
    %             (numel(rec_info.exportedChan) - channelNum);     % set index loc to selected electrode
end
winIdxEnd=winIdxStart+...
    (2*tracesInfo.excerptSize*numel(rec_info.exportedChan));
excerptWindow=winIdxStart:winIdxEnd-1;
if size(excerptWindow,2)>(2*tracesInfo.excerptSize*...
        numel(rec_info.exportedChan)) %for some reason
    excerptWindow=excerptWindow(1:end-(size(excerptWindow,2)-...
        (2*tracesInfo.excerptSize*numel(rec_info.exportedChan))));
end
% excerptWindow=1:960000;
dataExcerpt=traces.Data(excerptWindow);
dataExcerpt=reshape(dataExcerpt,[numel(rec_info.exportedChan)...
    tracesInfo.excerptSize*2]);


if tracesInfo.preproc==0 % raw data is presumed bandpassed filtered at this point
    preprocOption={'CAR','all'};
    dataExcerpt=PreProcData(dataExcerpt,rec_info.samplingRate,preprocOption);
end

%% plotting
channelNum=1;
dataExcerpt=int32(dataExcerpt(channelNum,:));

figure
plot(dataExcerpt,'k','linewidth',0.1); hold on;