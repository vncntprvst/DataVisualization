function varargout = SpikeVisualizationGUI(varargin)
% MATLAB code for SpikeVisualizationGUI.fig


% Last Modified by GUIDE v2.5 13-Jun-2016 20:24:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @SpikeVisualizationGUI_OpeningFcn, ...
    'gui_OutputFcn',  @SpikeVisualizationGUI_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


%% --- Executes just before SpikeVisualizationGUI is made visible.
function SpikeVisualizationGUI_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for SpikeVisualizationGUI
handles.output = hObject;

if ~isempty(varargin)
    handles=catstruct(handles,varargin{:}); % catstruct available here:
    % http://www.mathworks.com/matlabcentral/fileexchange/7842-catstruct
    
    if isfield(handles,'fname')
        set(handles.FileName,'string',handles.fname(1:end-4));
    else
        set(handles.FileName,'string','');
    end
else
    userinfo=UserDirInfo;
    %% get most recently changed data folder
    exportDir=regexprep(userinfo.directory,'\\\w+$','\\export');
    dataDirListing=dir(exportDir);
    %removing dots
    dataDirListing=dataDirListing(cellfun('isempty',cellfun(@(x) strfind(x,'.'),...
        {dataDirListing.name},'UniformOutput',false)));
    %removing other folders
    dataDirListing=dataDirListing(cellfun('isempty',cellfun(@(x)...
        regexp('list | all | unwanted | folders | here ',x),...
        {dataDirListing.name},'UniformOutput',false)));
    [~,fDateIdx]=sort([dataDirListing.datenum],'descend');
    recentDataFolder=[exportDir userinfo.slash dataDirListing(fDateIdx(1)).name userinfo.slash];

    % open user input window
    [handles.spikeFile,handles.exportDir] = uigetfile({'*.mat;*.hdf5','Export Formats';...
        '*.dat','Raw data';'*.*','All Files' },'Most recent data',recentDataFolder);
    if handles.spikeFile==0
    handles.spikeFile='';
    handles.exportDir=recentDataFolder;
    end
end
%define figure colormap
colormapSeed=colormap(lines);
handles.cmap=[colormapSeed(1:7,:);(colormapSeed+flipud(colormap(copper)))/2];

if isfield(handles,'spikeFile') && sum(strfind(handles.spikeFile,'.dat'))
    GetSortedSpikes_PB_Callback(hObject, handles);
end
handles=LoadSpikes(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SpikeVisualizationGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

%% Load data function
function handles=LoadSpikes(handles)
if isfield(handles,'subset')
    handles = rmfield(handles,'subset');
end
if isfield(handles,'spikeFile')
    handles = rmfield(handles,'spikeFile');
end
if isfield(handles,'Spikes')
    handles = rmfield(handles,'Spikes');
end
% function declaration
axis_name= @(x) sprintf('Chan %.0f',x);

if isfield(handles,'fname') && strcmp(handles.fname,'')
    set(handles.FileName,'string','')
else
    if ~isfield(handles,'spikeFile')
        cd(handles.exportDir);
        userinfo=UserDirInfo;
        exportDirListing=dir;
        handles.spikeFile={exportDirListing(~cellfun('isempty',cellfun(@(x) strfind(x,'_spikes.'),...
            {exportDirListing.name},'UniformOutput',false))).name};
        if size(handles.spikeFile,2)>1
            nameComp=cellfun(@(x) sum(ismember(x,handles.fname)) ,handles.spikeFile);
            if abs(diff(nameComp))<2 %that can be tricky for some files
                %select the most recent
                fileDates=datenum({exportDirListing(~cellfun('isempty',cellfun(@(x) strfind(x,'_spikes.'),...
                    {exportDirListing.name},'UniformOutput',false))).date});
                handles.spikeFile=handles.spikeFile{fileDates==max(fileDates)};
            else
                handles.spikeFile=handles.spikeFile{nameComp==max(nameComp)};
            end
        else
            handles.spikeFile=handles.spikeFile{:};
        end
    end
%     handles.exportDir='C:\Data\export\PrV75_61_optostim2_BR_6Ch_SyncCh_CAR';
%     cd(handles.exportDir);
%     handles.spikeFile='PrV75_61_optostim2_BR_6Ch_SyncCh_CAR_Ch3.mat';
    set(handles.FileName,'string',[handles.exportDir userinfo.slash handles.spikeFile])
    
    %% Load spike data
    if isfield(handles,'offlineSpikeSort')
        if logical(regexp(handles.offlineSpikeSort,'Ch\d+.')) % Spike2
            Spikes.Offline_Sorting=LoadSpikeData([handles.offlineSpikeSortDir...
                handles.offlineSpikeSort],handles.rec_info.numRecChan,30000);
        elseif logical(regexp(handles.offlineSpikeSort,'.hdf5'))
            Spikes.Offline_Sorting=LoadSpikeData([handles.offlineSpikeSortDir...
                handles.offlineSpikeSort],handles.rec_info.numRecChan);
        end
    else
        spikeData=load(handles.spikeFile);
        handles=catstruct(handles,spikeData);
        clear spikeData;
    end
   
    if isfield(handles.Spikes,'Online_Sorting') || isfield(handles.Spikes,'Offline_Sorting')
        if isfield(handles.Spikes,'Offline_Sorting') %preferred
            set(handles.Spikes_SortOff_RB,'value',1);
            set(handles.Spikes_SortOn_RB,'value',0);
            handles.Spikes.inGUI.Units=handles.Spikes.Offline_Sorting.Units;
            handles.Spikes.inGUI.SpikeTimes=handles.Spikes.Offline_Sorting.SpikeTimes;
            handles.Spikes.inGUI.Waveforms=handles.Spikes.Offline_Sorting.Waveforms;
            handles.Spikes.inGUI.samplingRate=handles.Spikes.Offline_Sorting.samplingRate;
        else
            set(handles.Spikes_SortOff_RB,'value',0);
            set(handles.Spikes_SortOn_RB,'value',1);
            handles.Spikes.inGUI.Units=handles.Spikes.Online_Sorting.Units;
            handles.Spikes.inGUI.SpikeTimes=handles.Spikes.Online_Sorting.SpikeTimes;
            handles.Spikes.inGUI.Waveforms=handles.Spikes.Online_Sorting.Waveforms;
            handles.Spikes.inGUI.samplingRate=handles.Spikes.Online_Sorting.samplingRate;
        end
    end
    
    % downsample spike times to 1 millisecond bins
%     spikeTimes=handles.Spikes.inGUI.SpikeTimes;
    for chNum=1:size(handles.Spikes.inGUI.samplingRate,1)
        handles.Spikes.inGUI.samplingRate(chNum,2)=1000;
        handles.Spikes.inGUI.SpikeTimes{chNum,2}=...
            handles.Spikes.inGUI.SpikeTimes{chNum,1}/...
            (handles.Spikes.inGUI.samplingRate(chNum,1)/...
            handles.Spikes.inGUI.samplingRate(chNum,2));
        
%         spikeTimeIdx=zeros(1,size(Spikes.Offline_Threshold.data{ChExN,1},2));
%         spikeTimeIdx(Spikes.Offline_Threshold.data{ChExN,1})=1;
%         
%         binSize=1;
%         numBin=ceil(max(spikeTimes{chNum})/...
%             (handles.Spikes.inGUI.samplingRate(chNum,1)/...
%             handles.Spikes.inGUI.samplingRate(chNum,2))/binSize);
%         % binspikeTime = histogram(double(spikeTimes), numBin); %plots directly histogram
%         [data,binEdges] = histcounts(double(spikeTimes{chNum}),...
%             linspace(0,max(double(spikeTimes{chNum})),numBin));
%         data(data>1)=1; %no more than 1 spike per ms
    end
    
    %% Set number of electrodes and units, 
    % if opening GUI, select electrode with most units and spikes
    if strcmp(get(handles.SelectElectrode_LB,'string'),'none')
        set(handles.SelectElectrode_LB,'string',num2str(handles.Spikes.Offline_Threshold.electrode'));
        if isfield(handles.Spikes,'Online_Sorting') || isfield(handles.Spikes,'Offline_Sorting')
            numUnits=cellfun(@(x) sum(length(x)*unique(x)), handles.Spikes.inGUI.Units);
            electrodeNum=find(numUnits==max(numUnits),1);
            set(handles.SelectElectrode_LB,'value',electrodeNum);
        else
            set(handles.SelectElectrode_LB,'value',1)
        end
    else
        electrodeNum=get(handles.SelectElectrode_LB,'value');
    end
        
    %% initialize variables
    unitsIdx=handles.Spikes.inGUI.Units{electrodeNum};
    waveForms=handles.Spikes.inGUI.Waveforms{electrodeNum};
        
    % how many units on that electrode?
    unitsID=unique(unitsIdx); %number of clustered units
    set(handles.SelectUnit_LB,'string',num2str(unitsID'));
    set(handles.SelectUnit_LB,'value',find(unitsID~=0));
    
    %% take out big ouliers
    WFmeanZ=mean(abs(zscore(single(waveForms'))),2);
    % if more than one, plot it and keep
    if sum(WFmeanZ>6)>1
        figure('name', 'Artifacts','position',[30   500   500   400]);
        plot(waveForms(:,WFmeanZ>6)','linewidth',2.5); hold on;
        plot(mean(waveForms,2),'linewidth',2.5);
        legend({'Artifact','Mean waveforms'});
        title('Potential artifacts removed, mean sigma > 6');
    else
        handles.Spikes.inGUI.Units{electrodeNum}(WFmeanZ>6)=-9;%artifacts
    end
    
    % here's how this can work:
    %     Spikes detected from threshold are the benchmark (unit code = 0).
    %     May want to extract waveforms, but most important is the time.
    %     Sorted units (whatever the source) "color" those units. One spike per ms max.
    handles=Plot_Unsorted_WF(handles);
    handles=Plot_Sorted_WF(handles);
    Plot_Mean_WF(handles);
    Plot_Raster_TW(handles);
    Plot_ISI(handles);
    Plot_ACG(handles);
    Plot_XCG(handles);
end

%% Plot ISI
function Plot_ISI(handles)

% get which unit to plot
if get(handles.ShowAllUnits_RB,'value')
    unitID=str2num(get(handles.SelectUnit_LB,'string'));
    selectedUnitsListIdx=find(unitID>0);
    selectedUnits=unitID(selectedUnitsListIdx);
    % keep the first one
else
    unitID=str2num(get(handles.SelectUnit_LB,'string'));
    selectedUnitsListIdx=get(handles.SelectUnit_LB,'value');
    selectedUnits=unitID(selectedUnitsListIdx);
end
if isempty(selectedUnits)
    cla(handles.ISI_Axes);
    return
end
electrodeNum=get(handles.SelectElectrode_LB,'value');
spikeTimes=handles.Spikes.inGUI.SpikeTimes{electrodeNum,1};
unitsIdx=handles.Spikes.inGUI.Units{electrodeNum};
samplingRate=handles.Spikes.inGUI.samplingRate(electrodeNum,1);

%keep the most numerous if more than one
if size(selectedUnits,1)>1
    keepU=1;
    for uidx=1:size(selectedUnits,1)
        if sum(unitsIdx==selectedUnits(uidx))>sum(unitsIdx==selectedUnits(keepU))
            keepU=uidx;
        end
    end
    selectedUnits=selectedUnits(keepU);    
end

unitST=spikeTimes(unitsIdx==selectedUnits);
% compute interspike interval
ISI=diff(unitST)/(samplingRate/1000);
axes(handles.ISI_Axes); hold on; 
cla(handles.ISI_Axes);
set(handles.ISI_Axes,'Visible','on');
ISIhist=histogram(double(ISI),0:5:max(ISI)+1);  %,'Normalization','probability'
ISIhist.FaceColor = handles.cmap(unitID(unitID==selectedUnits),:);
ISIhist.EdgeColor = 'k';
xlabel('Inter-spike Interval distribution (ms)')
axis('tight');box off;
set(gca,'xlim',[0 100],'Color','white','FontSize',10,'FontName','calibri');
hold off

%% Plot autocorrelogram
function Plot_ACG(handles)
% get which unit to plot
if get(handles.ShowAllUnits_RB,'value')
    unitID=str2num(get(handles.SelectUnit_LB,'string'));
    selectedUnitsListIdx=find(unitID>0);
    selectedUnits=unitID(selectedUnitsListIdx);
    % keep the first one
else
    unitID=str2num(get(handles.SelectUnit_LB,'string'));
    selectedUnitsListIdx=get(handles.SelectUnit_LB,'value');
    selectedUnits=unitID(selectedUnitsListIdx);
end
if isempty(selectedUnits)
    cla(handles.ACG_Axes);
    return
end

electrodeNum=get(handles.SelectElectrode_LB,'value');
spikeTimes=handles.Spikes.inGUI.SpikeTimes{electrodeNum,1};
unitsIdx=handles.Spikes.inGUI.Units{electrodeNum};
samplingRate=handles.Spikes.inGUI.samplingRate(electrodeNum,1);

%keep the most numerous if more than one
if size(selectedUnits,1)>1
    keepU=1;
    for uidx=1:size(selectedUnits,1)
        if sum(unitsIdx==selectedUnits(uidx))>sum(unitsIdx==selectedUnits(keepU))
            keepU=uidx;
        end
    end
    selectedUnits=selectedUnits(keepU);    
end
%get unit spike times
unitST=spikeTimes(unitsIdx==selectedUnits);
% change to ms timescale
unitST=unitST/(samplingRate/1000);
%get ISI
ISI=diff(unitST)/(samplingRate/1000);
%bin
spikeTimeIdx=zeros(1,unitST(end));
spikeTimeIdx(unitST)=1;
binSize=5;
numBin=ceil(size(spikeTimeIdx,2)/binSize);
binUnits = histcounts(double(unitST), linspace(0,size(spikeTimeIdx,2),numBin));
binUnits(binUnits>1)=1; %no more than 1 spike per ms
% compute autocorrelogram 
[ACG,lags]=xcorr(double(binUnits),200,'unbiased'); %'coeff'
ACG(lags==0)=0;
axes(handles.ACG_Axes); hold on; 
cla(handles.ACG_Axes);
set(handles.ACG_Axes,'Visible','on');
ACGh=bar(lags,ACG);
ACGh.FaceColor = handles.cmap(unitID(unitID==selectedUnits),:);
ACGh.EdgeColor = 'none';
axis('tight');box off;
xlabel('Autocorrelogram (5 ms bins)')
set(gca,'xlim',[-100 100],'Color','white','FontSize',10,'FontName','calibri','TickDir','out');
hold off

%% Plot cross-correlogram
function Plot_XCG(handles)


%% Plot Unsorted Spikes
function handles=Plot_Unsorted_WF(handles)
electrodeNum=get(handles.SelectElectrode_LB,'value');
waveForms=handles.Spikes.inGUI.Waveforms{electrodeNum};
unitsIdx=handles.Spikes.inGUI.Units{electrodeNum};
samplingRate=handles.Spikes.inGUI.samplingRate(electrodeNum,1);
%% Plot unsorted spikes
if sum(unitsIdx==0)>2000 %then only plot subset of waveforms
    subset=find(unitsIdx==0);
    handles.subset{1}=subset(1:round(sum(unitsIdx==0)/2000):end);
else
    handles.subset{1}=find(unitsIdx==0);
end
axes(handles.UnsortedUnits_Axes); hold on;
cla(handles.UnsortedUnits_Axes);
set(handles.UnsortedUnits_Axes,'Visible','on');
plot(waveForms(:,handles.subset{1}),'linewidth',1,'Color',[0 0 0 0.2]);
lineH=flipud(findobj(gca,'Type', 'line'));
% childH=flipud(get(gca,'Children'));
for lineTag=1:size(lineH,1)
    lineH(lineTag).Tag=num2str(handles.subset{1}(lineTag)); %Tag the unit ID
end
% foo=reshape([lineH.YData],size([lineH.YData],2)/size(lineH,1),size(lineH,1));
% faa=reshape([childH.YData],size([childH.YData],2)/size(childH,1),size(childH,1));
% lineH(80).Tag
% figure; hold on;
% plot(foo(:,80));
% plot(waveForms(:,handles.subset{1}(80)));
% [lineH.Tag]=deal(num2str(handles.subset{1}));
set(gca,'xtick',linspace(0,size(waveForms(:,handles.subset{1}),1),5),...
    'xticklabel',round(linspace(-round(size(waveForms(:,handles.subset{1}),1)/2),...
    round(size(waveForms(:,handles.subset{1}),1)/2),5)/(double(samplingRate)/1000),2),'TickDir','out');
% legend('Unclustered waveforms','location','southeast')
axis('tight');box off;
xlabel('Time (ms)')
ylabel('Voltage (0.25uV)')
set(gca,'Color','white','FontSize',10,'FontName','calibri');

%% Plot clusters
function handles=Plot_Sorted_WF(handles)
electrodeNum=get(handles.SelectElectrode_LB,'value');
waveForms=handles.Spikes.inGUI.Waveforms{electrodeNum};
unitsIdx=handles.Spikes.inGUI.Units{electrodeNum};
samplingRate=handles.Spikes.inGUI.samplingRate(electrodeNum,1);
% selected unit ids
axes(handles.SortedUnits_Axes); hold on;%colormap lines; cmap=colormap;
cla(handles.SortedUnits_Axes);
set(handles.SortedUnits_Axes,'Visible','on');
if get(handles.ShowAllUnits_RB,'value')
    unitID=str2num(get(handles.SelectUnit_LB,'string'));
    selectedUnitsListIdx=find(unitID>0);
    selectedUnits=unitID(selectedUnitsListIdx);
else
    unitID=str2num(get(handles.SelectUnit_LB,'string'));
    selectedUnitsListIdx=get(handles.SelectUnit_LB,'value');
    selectedUnits=unitID(selectedUnitsListIdx);
end
for unitP=1:length(selectedUnits)
    if sum(unitsIdx==selectedUnits(unitP))>2000 %then only plot subset of waveforms
        subset=find(unitsIdx==selectedUnits(unitP));
        handles.subset{selectedUnitsListIdx(unitP)}=subset(1:round(sum(unitsIdx==selectedUnits(unitP))/2000):end);
    else
        handles.subset{selectedUnitsListIdx(unitP)}=find(unitsIdx==selectedUnits(unitP));
    end
    plot(waveForms(:,handles.subset{selectedUnitsListIdx(unitP)}),...
        'linewidth',1,'Color',[handles.cmap(unitID(selectedUnitsListIdx(unitP)),:),0.4],...
        'Tag',num2str(selectedUnits(unitP)));
    lineH=flipud(findobj(gca,'Type', 'line'));
    passedTags=0;
    for lineTag=1:size(lineH,1)
        % check if it's the right cluster
        if strcmp(lineH(lineTag).Tag,num2str(selectedUnits(unitP)))
            %apply corresponding Tag
            lineH(lineTag).Tag=num2str(handles.subset{selectedUnitsListIdx(unitP)}(lineTag-passedTags)); %Tag the unit ID
        else
            passedTags=passedTags+1;
        end
    end
%     figure;hold on
%     plot(lineH(lineTag).YData);
%     plot(waveForms(:,handles.subset{selectedUnitsListIdx(unitP)}(lineTag-passedTags)))
end
set(gca,'xtick',linspace(0,size(waveForms(:,handles.subset{selectedUnitsListIdx(unitP)}),1),5),...
    'xticklabel',round(linspace(-round(size(waveForms(:,handles.subset{selectedUnitsListIdx(unitP)}),1)/2),...
    round(size(waveForms(:,handles.subset{selectedUnitsListIdx(unitP)}),1)/2),5)/(double(samplingRate)/1000),2),'TickDir','out');
% legend('Unclustered waveforms','location','southeast')
axis('tight');box off;
xlabel('Time (ms)')
ylabel('Voltage (0.25uV)')
set(gca,'Color','white','FontSize',10,'FontName','calibri');
hold off

%% Plot mean waveforms
function Plot_Mean_WF(handles)
electrodeNum=get(handles.SelectElectrode_LB,'value');
waveForms=handles.Spikes.inGUI.Waveforms{electrodeNum};
unitsIdx=handles.Spikes.inGUI.Units{electrodeNum};
samplingRate=handles.Spikes.inGUI.samplingRate(electrodeNum,1);
axes(handles.MeanSortedUnits_Axes); hold on;%colormap lines; 
cla(handles.MeanSortedUnits_Axes);
set(handles.MeanSortedUnits_Axes,'Visible','on');
if get(handles.ShowAllUnits_RB,'value')
    unitID=str2num(get(handles.SelectUnit_LB,'string'));
    selectedUnitsListIdx=find(unitID>0);
    selectedUnits=unitID(selectedUnitsListIdx);
else
    unitID=str2num(get(handles.SelectUnit_LB,'string'));
    selectedUnitsListIdx=get(handles.SelectUnit_LB,'value');
    selectedUnits=unitID(selectedUnitsListIdx);
end
for unitP=1:length(selectedUnits)
    selectWF=single(waveForms(:,unitsIdx==selectedUnits(unitP))');
    if ~isnan(mean(selectWF))
        plot(mean(selectWF),'linewidth',2,'Color',[handles.cmap(unitID(selectedUnitsListIdx(unitP)),:),0.7]);
        wfSEM=std(selectWF)/ sqrt(size(selectWF,2)); %standard error of the mean
        wfSEM = wfSEM * 1.96; % 95% of the data will fall within 1.96 standard deviations of a normal distribution
        patch([1:length(wfSEM),fliplr(1:length(wfSEM))],...
            [mean(selectWF)-wfSEM,fliplr(mean(selectWF)+wfSEM)],...
            handles.cmap(unitID(selectedUnitsListIdx(unitP)),:),'EdgeColor','none','FaceAlpha',0.2);
        %duplicate mean unit over unsorted plot
        %             plot(handles.UnsortedUnits_Axes,mean(selectWF),'linewidth',2,'Color',handles.cmap(unitP,:));
        if unitP==1
            delete(findobj(handles.UnsortedUnits_Axes,'Type', 'patch'));
        end
        patch([1:length(wfSEM),fliplr(1:length(wfSEM))],...
            [mean(selectWF)-wfSEM,fliplr(mean(selectWF)+wfSEM)],...
            handles.cmap(unitID(selectedUnitsListIdx(unitP)),:),'EdgeColor','none','FaceAlpha',0.5,'Parent', handles.UnsortedUnits_Axes);
    end
end
set(gca,'xtick',linspace(0,size(waveForms(:,handles.subset{selectedUnitsListIdx(unitP)}),1),5),...
    'xticklabel',round(linspace(-round(size(waveForms(:,handles.subset{selectedUnitsListIdx(unitP)}),1)/2),...
    round(size(waveForms(:,handles.subset{selectedUnitsListIdx(unitP)}),1)/2),5)/(double(samplingRate)/1000),2),'TickDir','out');
% legend('Unclustered waveforms','location','southeast')
axis('tight');box off;
xlabel('Time (ms)')
ylabel('Voltage (0.25uV)')
set(gca,'Color','white','FontSize',10,'FontName','calibri');
hold off

function  Plot_Raster_TW(handles)
%% plot rasters
electrodeNum=get(handles.SelectElectrode_LB,'value');
spikeTimes=handles.Spikes.inGUI.SpikeTimes{electrodeNum,2};

% plot 10 sec or 2000 waveforms max

% --- Executes on mouse press over axes background.
function UnsortedUnits_Axes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to UnsortedUnits_Axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

electrodeNum=get(handles.SelectElectrode_LB,'value');

%% initialize variables
unitsIdx=find(handles.Spikes.inGUI.Units{electrodeNum}==0);
% waveForms=handles.Spikes.inGUI.Waveforms{electrodeNum};
% spikeTimes=handles.Spikes.inGUI.SpikeTimes{electrodeNum,2};
% samplingRate=handles.Spikes.inGUI.samplingRate(electrodeNum,1);

lineH=findobj(gca,'Type', 'line');

%adjust number of waveforms to those displayed.
% or pull it from handles
visibleLines=cellfun(@(x) strcmp(x,'on'), {lineH.Visible});

waveForms=fliplr(reshape([lineH(visibleLines).YData],...
    size([lineH(visibleLines).YData],2)/size(lineH(visibleLines),...
    1),size(lineH(visibleLines),1)));
waveForms=waveForms';%one waveform per row

% foo=fliplr(reshape([lineH(~visibleLines).YData],...
%     size([lineH(~visibleLines).YData],2)/size(lineH(~visibleLines),...
%     1),size(lineH(~visibleLines),1)));
% figure;hold on;plot(foo,'k');

%This is the unsorted units plot
% make cluster classes 0 (unsorted), and -1 (hidden)
% clusterClasses=ones(size(waveForms,2),1)*-1;
% clusterClasses(handles.subset{1})=0;%mark subset visible (might be all of them)
% clusterClasses=clusterClasses(unitsIdx==0);%remove those already sorted
clusterClasses=zeros(size(waveForms,1),1);

%some other plots may have been overlayed. All unsorted lines should be black
if sum(cellfun(@(x) sum(x), {lineH.Color})~=0)
    %     %add waveform data to waveforms matrix
    %     waveForms=[reshape([lineH(cellfun(@(x) sum(x), {lineH.Color})~=0).YData],...
    %         size(waveForms,1),sum(cellfun(@(x) sum(x), {lineH.Color})~=0)),waveForms];
    %     clusterClasses=[ones(sum(cellfun(@(x) sum(x), {lineH.Color})~=0),1)*-10;...
    %         clusterClasses];
    
    % change ClusterClass to -10
end

clusterClasses=InteractiveClassification(waveForms,clusterClasses,0); % viewClasses=0
% foo=handles.Spikes.inGUI.Waveforms{electrodeNum}; foo=foo';
% figure;plot(foo(unitsIdx(logical(clusterClasses)),:)');hold on
% plot(lineH(flip(logical(clusterClasses))).YData)
handles.Spikes.inGUI.Units{electrodeNum}(unitsIdx(logical(clusterClasses)))=...
    clusterClasses(logical(clusterClasses));
unitsID=unique(handles.Spikes.inGUI.Units{electrodeNum});
set(handles.SelectUnit_LB,'String',num2str(unitsID(unitsID>=0)'))
if find(clusterClasses(logical(clusterClasses))>0,1)
    handles=Plot_Sorted_WF(handles);
    Plot_Mean_WF(handles);
    Plot_ISI(handles);
    Plot_ACG(handles);
    Plot_XCG(handles);
end
%  Update handles structure
guidata(hObject, handles);

% --- Executes on mouse press over axes background.
function SortedUnits_Axes_ButtonDownFcn(hObject, eventdata, handles)
electrodeNum=get(handles.SelectElectrode_LB,'value');

%% initialize variables
unitID=str2num(get(handles.SelectUnit_LB,'string'));
% for uIdxNum=1:length(unitsID)
%     unitsIdx{uIdxNum}=find(handles.Spikes.inGUI.Units{electrodeNum}==unitsID(uIdxNum));
% end
unitsIdx=handles.Spikes.inGUI.Units{electrodeNum};

lineH=findobj(gca,'Type', 'line');
visibleLines=cellfun(@(x) strcmp(x,'on'), {lineH.Visible});

waveForms=fliplr(reshape([lineH(visibleLines).YData],...
    size([lineH(visibleLines).YData],2)/size(lineH(visibleLines),...
    1),size(lineH(visibleLines),1)));
waveForms=waveForms';%one waveform per row

linesTags=flip(cellfun(@(x) str2double(x), {lineH(visibleLines).Tag}));
% ismember(find(unitsIdx==2),linesClasses)
clusterClasses=unitsIdx(linesTags)';
% figure;hold on
% plot(waveForms(530:537,:)','r')
% plot(waveForms(1:8,:)','b')
if get(handles.ShowAllUnits_RB,'value')
    selectedUnitsListIdx=find(unitID>0);
    viewClasses=unitID(selectedUnitsListIdx);
else
    selectedUnitsListIdx=get(handles.SelectUnit_LB,'value');
    viewClasses=unitID(selectedUnitsListIdx);
end

[clusterClasses,lineSelecIdx]=InteractiveClassification(waveForms,clusterClasses,viewClasses); % viewClasses=0
% foo=handles.Spikes.inGUI.Waveforms{electrodeNum}; foo=foo';
% figure;plot(foo(unitsIdx(logical(clusterClasses)),:)');hold on
% plot(lineH(flip(logical(clusterClasses))).YData)
% 
% waveForms=handles.Spikes.inGUI.Waveforms{electrodeNum};
% figure; plot(waveForms(:,linesTags(lineSelecIdx)))
if lineSelecIdx==0
    return
end
handles.Spikes.inGUI.Units{electrodeNum}(linesTags(lineSelecIdx))=...
    clusterClasses(lineSelecIdx);
unitsID=unique(handles.Spikes.inGUI.Units{electrodeNum});
set(handles.SelectUnit_LB,'String',num2str(unitsID(unitsID>=0)'))
if find(clusterClasses(logical(clusterClasses))>0,1)
    handles=Plot_Sorted_WF(handles);
end
if find(clusterClasses(logical(clusterClasses))==0,1)
    handles=Plot_Unsorted_WF(handles);
end
Plot_Mean_WF(handles);
Plot_ISI(handles);
Plot_ACG(handles);
Plot_XCG(handles);
%  Update handles structure
guidata(hObject, handles);

%% --- Executes on selection change in SelectElectrode_LB.
function SelectElectrode_LB_Callback(hObject, eventdata, handles)
if strcmp(get(gcf,'SelectionType'),'normal')
handles=LoadSpikes(handles);
% Update handles structure
guidata(hObject, handles);
end

%% Get chanel and unit selection
%     channelMenu=get(handles.SelectElectrode_LB,'string');
%     channelSelected=get(handles.SelectElectrode_LB,'value');
%     channelSelected=channelMenu(channelSelected);
%
%     unitMenu=get(handles.SelectUnit_LB,'string');
%     unitsSelected=get(handles.SelectUnit_LB,'value');
%     unitsSelected=unitMenu(unitsSelected);


%% --- Executes on selection change in SelectUnit_LB.
function SelectUnit_LB_Callback(hObject, eventdata, handles)
set(handles.ShowAllUnits_RB,'value',0)
if strcmp(get(gcf,'SelectionType'),'normal')
    %     = cellstr(get(hObject,'String'))
    %        contents{get(hObject,'Value')}
    handles=Plot_Sorted_WF(handles);
    Plot_Mean_WF(handles);
    Plot_ISI(handles);
    Plot_ACG(handles);
    Plot_XCG(handles);
    % elseif strcmp(get(gcf,'SelectionType'),'open') % double click
    % else
    %  Update handles structure
    guidata(hObject, handles);
end

%% --- Executes on button press in ShowAllUnits_RB.
function ShowAllUnits_RB_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    handles=Plot_Sorted_WF(handles);
    Plot_Mean_WF(handles);
    Plot_ISI(handles);
    Plot_ACG(handles);
    Plot_XCG(handles);
    guidata(hObject, handles);
end

%% --- Executes on mouse press over axes background.
function MeanSortedUnits_Axes_ButtonDownFcn(hObject, eventdata, handles)

% InteractiveClassification; % viewClasses=0

% unitsIdx
%  Update handles structure
guidata(hObject, handles);

%% --- Executes on button press in PreviewTh_PB.
function PreviewTh_PB_Callback(hObject, eventdata, handles)
noppDatFile=[cell2mat(regexp(handles.spikeFile,'.+(?=\_\w+\_\w+\.)','match')) '_nopp'];
handles.PreviewTh_PB.ForegroundColor=[0.2 0.5 0.7];
[status,cmdout]=RunSpykingCircus(cd,noppDatFile,'previewspkc');
handles.PreviewTh_PB.ForegroundColor=[0.3490 0.2000 0.3294];

%% --- Executes on button press in GetSortedSpikes_PB.
function GetSortedSpikes_PB_Callback(hObject, eventdata, handles)
noppDatFile=[cell2mat(regexp(handles.spikeFile,'.+(?=\_\w+\_\w+\.)','match')) '_nopp'];
[status,cmdout]=RunSpykingCircus(cd,noppDatFile,'runspkc,exportspikes');
if handles.GetSortedSpikes_PB.ForegroundColor(1)==0.3490
    handles.GetSortedSpikes_PB.ForegroundColor=[0.2 0.5 0.7];
else
    handles.GetSortedSpikes_PB.ForegroundColor=[0.3490 0.2000 0.3294];
end

%% --- Executes on button press in RefineSort_PB.
function RefineSort_PB_Callback(hObject, eventdata, handles)
noppDatFile=[cell2mat(regexp(handles.spikeFile,'.+(?=\_\w+\_\w+\.)','match')) '_nopp'];
handles.RefineSort_PB.ForegroundColor=[0.2 0.5 0.7];
[status,cmdout]=RunSpykingCircus(cd,noppDatFile,'startGUI');
handles.RefineSort_PB.ForegroundColor=[0.3490 0.2000 0.3294];

%% --- Outputs from this function are returned to the command line.
function varargout = SpikeVisualizationGUI_OutputFcn(hObject, eventdata, handles)

% Get default command line output from handles structure
varargout{1} = handles.output;


%% --- Executes on button press in Spikes_SortOff_RB.
function Spikes_SortOff_RB_Callback(hObject, eventdata, handles)

[handles.offlineSpikeSort,handles.offlineSpikeSortDir] = uigetfile({'*.mat;*.hdf5','All Data Formats';...
    '*.*','All Files' },'Export folder',handles.exportDir);
handles=LoadSpikes(handles);
% Update handles structure
guidata(hObject, handles);

%% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)

function edit1_Callback(hObject, eventdata, handles)

%% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)


%% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)

%% --- Executes during object creation, after setting all properties.
function SelectElectrode_LB_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)

%% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

%% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)

%% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)

%% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)

%% --- Executes on button press in Spikes_Th_RB.
function Spikes_Th_RB_Callback(hObject, eventdata, handles)

%% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)

%% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit2_Callback(hObject, eventdata, handles)


%% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% --- Executes on button press in Spikes_SortOn_RB.
function Spikes_SortOn_RB_Callback(hObject, eventdata, handles)

%% --- Executes during object creation, after setting all properties.
function SelectUnit_LB_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% --- Executes on slider movement.
function TW_slider_Callback(hObject, eventdata, handles)

%% --- Executes during object creation, after setting all properties.
function TW_slider_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

%% --- Executes on button press in TWplus_PB.
function TWplus_PB_Callback(hObject, eventdata, handles)

%% --- Executes on button press in TWminus_PB.
function TWminus_PB_Callback(hObject, eventdata, handles)

%% --- Executes on button press in TWall_PB.
function TWall_PB_Callback(hObject, eventdata, handles)

%% --- Executes on button press in LoadFile_PB.
function LoadFile_PB_Callback(hObject, eventdata, handles)

%% --- Executes on button press in Reload_PB.
function Reload_PB_Callback(hObject, eventdata, handles)

%% --- Executes on button press in Save_PB.
function Save_PB_Callback(hObject, eventdata, handles)
    userinfo=UserDirInfo;
    save([handles.exportDir userinfo.slash cell2mat(regexp(handles.fname,'.+(?=\.)','match'))...
        '_spikesResorted'],'-struct','handles','Spikes','-v7.3');
