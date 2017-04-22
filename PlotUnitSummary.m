%% load data file
[spikeData,fileInfo]=GetSpikeData(22); %
Trials = LoadTTL([fileInfo.Filename fileInfo.FileExt]);

%% get spike times and convert to binary array
for clusNum=1:size(spikeData,2)
    %% convert to 1 millisecond bins and plot excerpt
    binSize=1;
    numBin=ceil((max(spikeData(clusNum).spikeTimes)+1)/...
        (Trials.sampleRate/1000)/binSize);
    
    [spikeCount,spikeTime]=histcounts(double(spikeData(clusNum).spikeTimes)/...
        double(Trials.sampleRate/1000), numBin);

    figure; bar(spikeTime(1:6000),spikeCount(1:6000),'hist')

    %% spike density function
    raster = zeros(1,ceil(max(spikeTime))+1);
    raster(ceil(spikeTime(1:end-1)))=spikeCount;
    sigma=20;
    convSpikeTime = fullgauss_filtconv(raster,sigma,0).*1000;
    hold on 
    % plot([zeros(1,sigma*3) convspikeTime zeros(1,sigma*3)])
    plot([zeros(1,sigma*3) convSpikeTime(1:6000-sigma*3)])

    %% plot rasters aligned to laser pulse TTL
    OptoRasters
end
