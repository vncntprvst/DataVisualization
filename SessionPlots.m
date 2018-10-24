% plot results from recording session plots
spikeSortingFiles = cellfun(@(fileFormat) dir([cd filesep '**' filesep fileFormat]),...
    {'*.result.hdf5','*_jrc.mat','*.csv'},'UniformOutput', false);
spikeSortingFiles=vertcat(spikeSortingFiles{~cellfun('isempty',spikeSortingFiles)});
% do not include those files:
spikeSortingFiles=spikeSortingFiles(~cellfun(@(flnm) contains(flnm,{'DeepCut'}),...
    {spikeSortingFiles.name}));
sessionDir=cd;

dataFiles = cellfun(@(fileFormat) dir([cd filesep '**' filesep fileFormat]),...
    {'*.dat','*raw.kwd','*RAW*Ch*.nex','*.ns*'},'UniformOutput', false);
dataFiles=vertcat(dataFiles{~cellfun('isempty',dataFiles)});
% keep those files
TTLFiles=dataFiles(cellfun(@(flnm) contains(flnm,{'_TTLs'}),...
    {dataFiles.name}));
dataFiles=dataFiles(cellfun(@(flnm) contains(flnm,{'_export'}),...
    {dataFiles.name}));

% for recNum=1:size(spikeSortingFiles,1)
        recNum=3;
    recDir=spikeSortingFiles(recNum).folder;
    recName=spikeSortingFiles(recNum).name;
    cd(recDir)
    dataFileIdx=cellfun(@(datF) contains(datF,regexp(recName,'\S+?(?=\.\w+\.\w+$)','match','once')) ,...
        {dataFiles.name});
    dataFileName=dataFiles(dataFileIdx).name;
    dataFileDir=dataFiles(dataFileIdx).folder;
    
    traces = memmapfile(fullfile(dataFileDir,dataFileName),'Format','int16');
    spikes=LoadSpikeData(recName,traces);
    
    %% load TTLs
    cd(sessionDir);
    TTLFileName=[regexp(recName,'\S+?(?=_export)','match','once') '_TTLs.dat'];
    fid = fopen(TTLFileName, 'r');
    TTLs = fread(fid,[2,Inf],'int32');
    fclose(fid);
    
    %% add voltage scaling factor and sampling rate
    bitResolution=0.195; %for Open Ephys
    spikes.waveforms=spikes.waveforms.*bitResolution;
    samplingRate=30000;
    %
    spikes.unitID=double(spikes.unitID);
    % find most frequent units
    [unitFreq,unitIDs]=hist(spikes.unitID,unique(spikes.unitID));
    [unitFreq,freqIdx]=sort(unitFreq','descend');
    unitFreq=unitFreq./sum(unitFreq)*100; unitIDs=unitIDs(freqIdx);
    bestUnitsIdx=find(unitFreq>2);
    bestUnits=unitIDs(unitIDs(bestUnitsIdx)>=~0);
%     bestUnits=0;
    %     %% generate rasters
        preAlignWindow=500;    postAlignWindow=2000;
        spikeRasters=PopulationRaster(spikes.times,TTLs(1,:),bestUnits,...
            spikes.unitID,samplingRate,preAlignWindow,postAlignWindow);
    
        %% plots
        % phototagging
        pulseDur=mode(diff(TTLs));
        IPI=mode(diff(TTLs(1,:)))+pulseDur;
    
        figure;   SDFh=subplot(1,1,1);
        OptoSDF(spikeRasters,preAlignWindow,pulseDur,IPI,SDFh)
        
%         sdf=conv_raster(spikeRasters{4},conv_sigma,1);

        
        
    %
    %
    %     conv_sigma=1;shiftVal=conv_sigma*3;
    %     allSDF=vertcat(spikeRasters{:});
    %     figure;
    %     colormap(hot) %flipud(gray));
    %     imagesc(allSDF); %
    %     % caxis([0 200])
    %     xlabel('Time (ms)');
    %     ylabel('Neuron#','FontSize',12); %'FontWeight','bold'
    %     % draw alignment bar
    %     currylim=get(gca,'YLim');
    %     %     currxlim=get(gca,'XLim');%midl=round(currxlim(2)/20)*10;
    %     %     set(gca,'XTick',preAlignWindow:50:max(get(gca,'xlim')));
    %     %     set(gca,'XTickLabel',0:50:max(get(gca,'xlim'))-preAlignWindow,'FontSize',10,'FontName','calibri','TickDir','out');
    %     set(gca,'XLim',[0.5 preAlignWindow+660.5],'XTick',[0:10:preAlignWindow+660]-shiftVal);
    %     set(gca,'XTickLabel',(0:10:preAlignWindow+660)-preAlignWindow,'FontSize',10,'FontName','calibri','TickDir','out');
    %
    %     %opto stim patch
    %     patch([preAlignWindow-shiftVal preAlignWindow-shiftVal preAlignWindow+pulseDur preAlignWindow+pulseDur], ...
    %         [[0 currylim(2)] fliplr([0 currylim(2)])], ...
    %         [0 0 0 0],[0.3 0.75 0.93],'EdgeColor','none','FaceAlpha',0.5);
    %     set(gca,'Color','white','FontSize',18,'FontName','Helvetica');
    
    % end
    
    % spike summary
    figure('name',regexp(recName,'\S+?(?=\.\w+\.\w+$)','match','once'))
    colormapSeed=lines;
    cmap=[colormapSeed(1:7,:);(colormapSeed+flipud(colormap(copper)))/2;autumn];
    
    for unitNum=bestUnits'
        unitIdx=find(unitIDs==unitNum)-1;
        %spike times for that unit
        unitSpikeTimes=spikes.times(spikes.unitID==unitNum);
        if ~isempty(diff(unitSpikeTimes))
            % Plot ISI
            subplot(3,numel(bestUnits),unitIdx)
            % compute interspike interval
            ISI=diff(unitSpikeTimes)/(samplingRate/1000);
            ISIhist=histogram(double(ISI),logspace(0, 4, 50),'DisplayStyle','stairs','LineWidth',1.5);  %,'Normalization','probability'
            %     ISIhist.FaceColor = handles.cmap(unitID(unitID==selectedUnits),:);
            ISIhist.EdgeColor = cmap(unitIdx,:); %'k';
            xlabel('Interspike Interval (ms)')
            axis('tight');box off; grid('on'); set(gca,'xscale','log','GridAlpha',0.25,'MinorGridAlpha',1);
            set(gca,'xlim',[0 10^4],... %'XTick',linspace(0,40,5),'XTickLabel',linspace(0,40,5),...
                'TickDir','out','Color','white','FontSize',10,'FontName','Calibri');
            
            % Plot autocorrelogram
            subplot(3,numel(bestUnitsIdx)-1,unitIdx+numel(bestUnitsIdx)-1)
            % change spiketimes to ms timescale
            unitSpikeTimes=unitSpikeTimes/(samplingRate/1000);
            spikeTimeIdx=zeros(1,unitSpikeTimes(end));
            spikeTimeIdx(unitSpikeTimes)=1;
            binSize=1;
            numBin=ceil(size(spikeTimeIdx,2)/binSize);
            binUnits = histcounts(double(unitSpikeTimes), linspace(0,size(spikeTimeIdx,2),numBin));
            binUnits(binUnits>1)=1; %no more than 1 spike per ms
            % compute autocorrelogram
            [ACG,lags]=xcorr(double(binUnits),200,'unbiased'); %'coeff'
            ACG(lags==0)=0;
            ACGh=bar(lags,ACG);
            ACGh.FaceColor = cmap(unitIdx,:);
            ACGh.EdgeColor = 'none';
            % axis('tight');
            box off; grid('on'); %set(gca,'yscale','log','GridAlpha',0.25,'MinorGridAlpha',1);
            xlabel('Autocorrelogram') %(1 ms bins)
            set(gca,'xlim',[-20 20],...
                'ylim',[0 max(get(gca,'ylim'))],...
                'Color','white','FontSize',10,'FontName','Calibri','TickDir','out');
        end
        
        % Plot the mean waveform
        subplot(3,numel(bestUnitsIdx)-1,unitIdx+(numel(bestUnitsIdx)-1)*2)
        unitWF=single(spikes.waveforms(spikes.unitID==unitNum,:));
        if isempty(find(isnan(mean(unitWF)), 1))
            plot(mean(unitWF),'linewidth',2,'Color',[cmap(unitIdx,:),0.7]);
            wfSEM=std(unitWF)/ sqrt(size(unitWF,2)); %standard error of the mean
            wfSEM = wfSEM * 1.96; % 95% of the data will fall within 1.96 standard deviations of a normal distribution
            patch([1:length(wfSEM),fliplr(1:length(wfSEM))],...
                [mean(unitWF)-wfSEM,fliplr(mean(unitWF)+wfSEM)],...
                cmap(unitIdx,:),'EdgeColor','none','FaceAlpha',0.2);
        end
        set(gca,'XTick',linspace(0,size(unitWF,2),5),...
            'XTickLabel',round(linspace(-round(size(unitWF,2)/2),...
            round(size(unitWF,2)/2),5)/(double(samplingRate)/1000),2),'TickDir','out');
        axis('tight');box off;
        xlabel('Time (ms)');
        ylabel('Voltage (\muV)');
    end
    
%     %plot the amplitude
%         for unitNum=bestUnits'
%             figure
%             plot(spikes.times(spikes.unitID==unitNum-1),spikes.amplitude(spikes.unitID==unitNum-1), '.')
%         end
% end