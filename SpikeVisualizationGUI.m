function varargout = SpikeVisualizationGUI(varargin)
% MATLAB code for SpikeVisualizationGUI.fig
% Last Modified by GUIDE v2.5 17-Nov-2017 16:03:13
% version 0.7 (nov 2017), tested in R2017b
% version 0.5 (sept 2016), tested in R2014b
% Vincent Prevosto
% email: vp35 at duke.edu

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
function SpikeVisualizationGUI_OpeningFcn(hObject, ~, handles, varargin)

% Choose default command line output for SpikeVisualizationGUI
handles.output = hObject;

if ~isempty(varargin)
    handles=CatStruct(handles,varargin{:}); % CatStruct available here:
    % http://www.mathworks.com/matlabcentral/fileexchange/7842-CatStruct
    
    if isfield(handles,'fname')
        set(handles.TB_FileName,'string',handles.fname(1:end-4));
    else
        set(handles.TB_FileName,'string','');
    end
    if isempty(handles.spikeFile)
        % check if spike sorting results are present in export folder
        cd(handles.exportDir);
        listExportDir=dir;
        listExportDir=listExportDir(cellfun('isempty',cellfun(@(x) strfind('.',x(end)),{listExportDir.name},'UniformOutput',false)));
        spikeFileIdx=cellfun(@(x) contains(x,'_spikes.mat'),{listExportDir.name});
        if logical(sum(spikeFileIdx))
            handles.spikeFile=listExportDir(spikeFileIdx).name;
        else %maybe spike sorting results are sub folder
            spikeSortingFolderIdx=[listExportDir.isdir];
            spikeSortingFolder=[handles.exportDir filesep listExportDir(spikeSortingFolderIdx).name];
            listSpikeSortingFolder=dir(spikeSortingFolder);
            spikeFileIdx=cellfun(@(x) contains(x,'result'),{listSpikeSortingFolder.name});
            if logical(sum(spikeFileIdx))
                handles.exportDir=spikeSortingFolder;
                handles.spikeFile=listSpikeSortingFolder(spikeFileIdx).name;
            else
                % ask user to select file to load
                [handles.spikeFile,handles.exportDir] = uigetfile({'*.mat;*.hdf5;*.csv','Export Formats';...
                    '*.dat','Raw data';'*.*','All Files' },'Select data to load',cd);
                if handles.spikeFile==0
                    handles=rmfield(handles,'spikeFile');
                    handles.exportDir=cd;
                end
            end
        end
    end
else
    try
        handles.userinfo=UserDirInfo;
    catch
        handles.userinfo=[];
        handles.userinfo.user=getenv('username');
    end
    %% get most recently changed data folder (looks for a folder named "export")
    if isfield(handles.userinfo,'directory')
        exportDir=regexprep(handles.userinfo.directory,'\\\w+$','\\export');
    else
        [exportDir,handles.userinfo.directory]=deal(cd);
    end
    dataDirListing=dir(exportDir);
    if ~isempty(dataDirListing)
        %removing dots
        dataDirListing=dataDirListing(cellfun('isempty',cellfun(@(x) strfind(x,'.'),...
            {dataDirListing.name},'UniformOutput',false)));
        %removing other folders
        dataDirListing=dataDirListing(cellfun('isempty',cellfun(@(x)...
            regexp('Behav | Video | Impedance',x),... % list | all | unwanted | folders | here
            {dataDirListing.name},'UniformOutput',false)));
        [~,fDateIdx]=sort([dataDirListing.datenum],'descend');
        recentDataFolder=[exportDir filesep dataDirListing(fDateIdx(1)).name filesep];
    else
        recentDataFolder=cd;
    end
    
    % ask user to select file to load
    [handles.spikeFile,handles.exportDir] = uigetfile({'*.mat;*.hdf5;*.csv','Export Formats';...
        '*.dat','Raw data';'*.*','All Files' },'Select data to load',recentDataFolder);
    if handles.spikeFile==0
        handles.spikeFile='';
        handles.exportDir=recentDataFolder;
    end
end
%define figure colormap
colormapSeed=lines;
handles.cmap=[colormapSeed(1:7,:);(colormapSeed+flipud(colormap(copper)))/2;autumn];

if isfield(handles,'spikeFile') && ~isempty(handles.spikeFile) && sum(strfind(handles.spikeFile,'.dat'))
    GetSortedSpikes_PB_Callback(hObject, handles);
end
handles.fileLoaded=0;
%create classification table
handles.classification = table([],[],[],categorical(),{},'VariableNames',{'SortID','Channel','UnitNumber','Classification','Comment'});

if isfield(handles,'spikeFile')
    handles=LoadSpikes(handles);
    handles=LoadRawData(handles);
end
% Update handles structure
guidata(hObject, handles);

%% --- Executes (conditional) at startup and on button press in GetSortedSpikes_PB.
function GetSortedSpikes_PB_Callback(hObject, ~, handles)
if ~isfield(handles,'datFile')
    handles.datFile=[cell2mat(regexp(handles.spikeFile,'.+(?=\_\w+\_\w+\.)','match')) '_nopp'];
end
[status,cmdout]=RunSpykingCircus(cd,handles.datFile,{'runspkc,exportspikes'});
if handles.GetSortedSpikes_PB.ForegroundColor(1)==0.3490
    handles.GetSortedSpikes_PB.ForegroundColor=[0.2 0.5 0.7];
else
    handles.GetSortedSpikes_PB.ForegroundColor=[0.3490 0.2000 0.3294];
end

%% Load data function
function handles=LoadSpikes(handles)
if isfield(handles,'subset')
    handles = rmfield(handles,'subset');
end
% if isfield(handles,'spikeFile')
%     handles = rmfield(handles,'spikeFile');
% end
if isfield(handles,'offlineSort_SpikeFile')
    handles = rmfield(handles,'offlineSort_SpikeFile');
end
% function declaration
axis_name= @(x) sprintf('Chan %.0f',x);

if isfield(handles,'fname') && strcmp(handles.fname,'')
    set(handles.TB_FileName,'string','')
else
    if ~isfield(handles,'fname') && isfield(handles,'rec_info')
        try
            handles.fname=handles.rec_info.exportname;
        catch
            handles.fname=handles.datFile;
        end
    end
    if handles.fileLoaded==0 & strfind(handles.spikeFile,'Resorted.mat')
        cd(handles.exportDir);
        spikeData=load(handles.spikeFile);
        handles = rmfield(handles,{'spikeFile','exportDir'});
        if isfield(handles,'rec_info')
            handles = rmfield(handles,'rec_info');
        end
        if isfield(handles,'Spikes')
            handles = rmfield(handles,'Spikes');
        end
        if isfield(spikeData,'classification')
            handles = rmfield(handles,'classification');
        end
        handles=CatStruct(handles,spikeData);
        clear spikeData;
    elseif handles.fileLoaded==0
        %% load spike data from matlab export, if exists
        cd(handles.exportDir);
        exportDirListing=dir;
        if isempty({exportDirListing(~cellfun('isempty',cellfun(@(x) strfind(x,'_spikes.'),...
                {exportDirListing.name},'UniformOutput',false))).name})
            handles.offlineSort_SpikeDir=handles.exportDir;
            handles.offlineSort_SpikeFile=handles.spikeFile;
            %go one folder up
            exportDirUp=handles.exportDir(1:regexp(handles.exportDir,'\\\w+\\$'));
            exportDirListing=dir(exportDirUp);
            handles.spikeFile={exportDirListing(~cellfun('isempty',cellfun(@(x) strfind(x,'_spikes.'),...
                {exportDirListing.name},'UniformOutput',false))).name};
            if ~isempty(handles.spikeFile)
                handles.exportDir=exportDirUp;
                cd(handles.exportDir);
            end
        end
        %         if ~isfield(handles,'fname') && isfield(handles,'offlineSort_SpikeFile')
        %             handles.fname=regexp(handles.offlineSort_SpikeFile,'.+?(?=\.)','match');
        %             handles.fname=handles.fname{1};
        %         end
        if iscell(handles.spikeFile) && size(handles.spikeFile,2)>1
            % ask which export file is the good one
            [whichOne,status] = listdlg('PromptString','Select correct spike file:',...
                'SelectionMode','single',...
                'ListString',handles.spikeFile);
            if status==1
                handles.spikeFile=handles.spikeFile{whichOne};
            else
                handles.spikeFile=[];
            end
            
            %             nameComp=cellfun(@(name) sum(ismember(handles.fname(1:end-4),...
            %                 name(1:size(handles.fname(1:end-4),2)))) ,handles.spikeFile);
            %             if abs(diff(nameComp))<2 %that can be tricky for some files
            %                 %select the most recent
            %                 fileDates=datenum({exportDirListing(~cellfun('isempty',cellfun(@(x) strfind(x,'_spikes.'),...
            %                     {exportDirListing.name},'UniformOutput',false))).date});
            %                 handles.spikeFile=handles.spikeFile{fileDates==max(fileDates)};
            %             else
            %                 handles.spikeFile=handles.spikeFile{nameComp==max(nameComp)};
            %             end
        else
            if iscell(handles.spikeFile)
                try
                    handles.spikeFile=handles.spikeFile{:};
                catch
                    handles.spikeFile=[];
                end
            end
        end
        %     end
        %     handles.exportDir='C:\Data\export\PrV75_61_optostim2_BR_6Ch_SyncCh_CAR';
        %     cd(handles.exportDir);
        %     handles.spikeFile='PrV75_61_optostim2_BR_6Ch_SyncCh_CAR_Ch3.mat';
        set(handles.TB_FileName,'string',[handles.exportDir handles.spikeFile])
        
        %% Load spike data
        if ~isempty(handles.spikeFile) && ~strcmp(get(handles.Spikes_PrevSorted_Menu,'Checked'),'on')
            spikeData=load(handles.spikeFile);
            if isfield(handles,'Spikes')
                % This will delete all Spikes data including any changes to HandSort
                % Might want to change that in future version
                handles = rmfield(handles,'Spikes');
            end
            handles=CatStruct(handles,spikeData);
            %clear
            clear spikeData;
        end
        % check if we need to load recording info
        if ~isfield(handles,'rec_info')
            try
                if ~isempty(handles.spikeFile)
                    fileName=regexp(handles.spikeFile,'.+(?=_\w+.mat$)','match');
                    recInfo=load([fileName{:} '_info.mat']);
                elseif ~isempty(handles.offlineSort_SpikeFile)
                    exportDirListing=dir(cd);
                    if sum(~cellfun('isempty',cellfun(@(x) strfind(x,'_info.'),...
                {exportDirListing.name},'UniformOutput',false)))
                    fileName=exportDirListing(~cellfun('isempty',cellfun(@(x) strfind(x,'_info.'),...
                {exportDirListing.name},'UniformOutput',false))).name;
                    recInfo=load(fileName);
                    else %might be in folder above
                        exportDirListing=dir('..');
                        fileName=exportDirListing(~cellfun('isempty',cellfun(@(x) strfind(x,'_info.'),...
                            {exportDirListing.name},'UniformOutput',false))).name;
                        currDir=cd; cd('..');
                        recInfo=load(fileName); cd(currDir);
                    end 
                end 
            catch
                [fileInfoName,fileInfoDir,FilterIndex] = uigetfile({'*.mat;*.txt','File Info Formats';...
                    '*.*','All Files' },'Select file providing basic recording info, or cancel and enter info');
                if FilterIndex==0 % display prompt
                    prompt={'Sampling rate',...
                        'Number of channels in spike-sorted file',...
                        'uV/bit'};
                    name='Data file information';
                    numlines=1;
                    defaultanswer={'30000','32','0.25'};
                    options.Resize='on';
                    options.WindowStyle='normal';
                    info=inputdlg(prompt,name,numlines,defaultanswer,options);
                    info=cellfun(@(field) str2double(field),info,'UniformOutput',false);
                    info{2}=1:info{2};
                    recInfo.rec_info=cell2struct(info, {'samplingRate','exportedChan','bitResolution'}, 1);
                else
                    recInfo=load(fullfile(fileInfoDir,fileInfoName));
                end
            end
            % just in case the field is named differently
            field = fieldnames(recInfo);
            if ~strcmp(field{:},'rec_info')
                [recInfo.('rec_info')] = recInfo.(field{:});
                recInfo = rmfield(recInfo,field);
            end
            %concatenate
            if isfield(handles,'rec_info')
                handles = rmfield(handles,'rec_info');
            end
            handles=CatStruct(handles,recInfo);
            %clear
            clear recInfo;
        end
        if isfield(handles,'offlineSort_SpikeFile')
%             cd(handles.exportDir);
            if logical(regexp(handles.offlineSort_SpikeFile,'Ch\d+.')) % Spike2
                cd(handles.offlineSort_SpikeDir);
                handles.Spikes.Offline_Sorting=LoadSpikeData(handles.offlineSort_SpikeFile,...
                    handles.rec_info.numRecChan,handles.rec_info.samplingRate);
            elseif contains(handles.offlineSort_SpikeFile,'.hdf5') % Spyking-Circus
                % first load raw traces in memory
                fileName=regexp(handles.offlineSort_SpikeFile,'.+(?=\.\w+\.\w+$)','match');
                %                 tb_filename=[tb_filename{1} '.dat'];
                if ~exist([fileName{:} '.dat'],'file') && ~exist([fileName{:} '.mat'],'file')
                    % try one folder up
                    cd ..
                end
                if ~exist([fileName{:} '.dat'],'file') && ~exist([fileName{:} '.mat'],'file')
                    % then ask where it is
                    [handles.datFile,handles.datDir] = uigetfile({'*.dat;*.mat','Data Formats';...
                        '*.*','All Files' },'Select data file for spike waveform extraction');
                else
                    if exist([fileName{1} '.dat'],'file')
                        handles.datFile=[fileName{1} '.dat'];
                    elseif exist([fileName{1} '.mat'],'file')
                        handles.datFile=[fileName{1} '.mat'];
                    end
                    handles.datDir=cd;
                end
                if strfind(handles.datFile,'.mat') % .mat file contain rawData
                    load(handles.datFile);
                elseif strfind(handles.datFile,'.dat')
                    rawData = memmapfile(fullfile(handles.datDir,handles.datFile),'Format','int16');
                end
                %% then import spike Data
                cd(handles.offlineSort_SpikeDir);
                if ~isfield(handles.rec_info,'bitResolution') || isempty(handles.rec_info.bitResolution)
                    handles.rec_info.bitResolution=0.25; %default 0.25 bit per uV % 0.195 for Open Ephys
                end
                handles.Spikes.Offline_Sorting=LoadSpikeData_byElectrode(handles.offlineSort_SpikeFile,rawData,...
                    numel(handles.rec_info.exportedChan),handles.rec_info.samplingRate,handles.rec_info.bitResolution); %handles.tracesInfo.size(1)
                clear rawData;
            elseif contains(handles.offlineSort_SpikeFile,'rez.mat') % KiloSort
                fileName=exportDirListing(~cellfun('isempty',cellfun(@(x) strfind(x,'.dat'),...
                {exportDirListing.name},'UniformOutput',false))).name;
                if iscell(fileName) && length(fileName)>1
                    fileName=exportDirListing(~cellfun('isempty',cellfun(@(x) strfind(x,rec_info.exportname),...
                fileName,'UniformOutput',false))).name;
                end
                %                 tb_filename=[tb_filename{1} '.dat'];
                if ~exist(fileName,'file') % then ask where it is
                    [handles.datFile,handles.datDir] = uigetfile({'*.dat;*.mat','Data Formats';...
                        '*.*','All Files' },'Select data file for spike waveform extraction');
                else
                    handles.datFile=fileName;
                    handles.datDir=cd;
                end
                if strfind(handles.datFile,'.mat') % .mat file contain rawData
                    load(handles.datFile);
                elseif strfind(handles.datFile,'.dat')
                    rawData = memmapfile(fullfile(handles.datDir,handles.datFile),'Format','int16');
                end
                %then import
                cd(handles.offlineSort_SpikeDir);
                if ~isfield(handles.rec_info,'bitResolution') || isempty(handles.rec_info.bitResolution)
                    handles.rec_info.bitResolution=0.25; %default 0.25 uV/bit
                end
                handles.Spikes.Offline_Sorting=LoadSpikeData(handles.offlineSort_SpikeFile,rawData,...
                    numel(handles.rec_info.exportedChan),handles.rec_info.samplingRate,handles.rec_info.bitResolution); %handles.tracesInfo.size(1)
                clear rawData;
            elseif contains(handles.offlineSort_SpikeFile,'.csv') || ...
                    contains(handles.offlineSort_SpikeFile,'_jrc.mat') % assuming JRClust
                fileName=exportDirListing(~cellfun('isempty',cellfun(@(x) strfind(x,'.dat'),...
                {exportDirListing.name},'UniformOutput',false))).name;
                if iscell(fileName) && length(fileName)>1
                    fileName=exportDirListing(~cellfun('isempty',cellfun(@(x) strfind(x,rec_info.exportname),...
                fileName,'UniformOutput',false))).name;
                end
                %                 tb_filename=[tb_filename{1} '.dat'];
                if ~exist(fileName,'file') % then ask where it is
                    [handles.datFile,handles.datDir] = uigetfile({'*.dat;*.mat','Data Formats';...
                        '*.*','All Files' },'Select data file for spike waveform extraction');
                else
                    handles.datFile=fileName;
                    handles.datDir=cd;
                end
                if strfind(handles.datFile,'.mat') % .mat file contain rawData
                    load(handles.datFile);
                elseif strfind(handles.datFile,'.dat')
                    rawData = memmapfile(fullfile(handles.datDir,handles.datFile),'Format','int16');
                end
                %then import
                cd(handles.offlineSort_SpikeDir);
                if ~isfield(handles.rec_info,'bitResolution') || isempty(handles.rec_info.bitResolution)
                    handles.rec_info.bitResolution=0.25; %default 0.25 uV/bit
                end
                handles.Spikes.Offline_Sorting=LoadSpikeData(handles.offlineSort_SpikeFile,rawData); %,...
%                     ,numel(handles.rec_info.exportedChan),handles.rec_info.samplingRate,...
%                     handles.rec_info.bitResolution); %handles.tracesInfo.size(1)   
                clear rawData;
            end
            cd(handles.exportDir);
        end
        if strcmp(get(handles.Spikes_Th_Menu,'Checked'),'on')
            %extract waveforms if needed
            if ~isfield(handles.Spikes.Offline_Threshold,'Waveforms')
                % first load raw traces
                fileName=regexp(handles.spikeFile,'.+(?=_\w+.\w+$)','match');
                load([fileName{:} '_raw.mat']);
                %then import
                handles.Spikes=LoadSpikeData(handles.spikeFile,...
                    handles.Spikes.Offline_Threshold.electrode,...
                    handles.Spikes.Offline_Threshold.samplingRate(:,1),rawData);
                clear rawData
            end
        end
        
        if isfield(handles.Spikes,'HandSort')
            % all good
        elseif isfield(handles.Spikes,'Online_Sorting') || isfield(handles.Spikes,'Offline_Sorting') ...
                || isfield(handles.Spikes,'Offline_Threshold')
            if isfield(handles.Spikes,'Offline_Threshold') && strcmp(get(handles.Spikes_Th_Menu,'Checked'),'on')
                set(handles.Spikes_SortOff_Menu,'Checked','off');
                set(handles.Spikes_SortOn_Menu,'Checked','off');
                handles.Spikes.HandSort.Units=handles.Spikes.Offline_Threshold.Units;
                handles.Spikes.HandSort.SpikeTimes=handles.Spikes.Offline_Threshold.SpikeTimes;
                handles.Spikes.HandSort.Waveforms=handles.Spikes.Offline_Threshold.Waveforms;
                handles.Spikes.HandSort.samplingRate=handles.Spikes.Offline_Threshold.samplingRate;
            elseif isfield(handles.Spikes,'Offline_Sorting') || strcmp(get(handles.Spikes_SortOff_Menu,'Checked'),'on')
                set(handles.Spikes_SortOff_Menu,'Checked','on');
                set(handles.Spikes_SortOn_Menu,'Checked','off');
                try
                    handles.Spikes.HandSort.Units=handles.Spikes.Offline_Sorting.Units;
                catch
                    handles.Spikes.HandSort.Units=handles.Spikes.Offline_Sorting.unitID;
                end
                try
                    handles.Spikes.HandSort.SpikeTimes=handles.Spikes.Offline_Sorting.SpikeTimes;
                catch
                    handles.Spikes.HandSort.SpikeTimes=handles.Spikes.Offline_Sorting.times;
                end
                try
                    handles.Spikes.HandSort.Waveforms=handles.Spikes.Offline_Sorting.Waveforms;
                catch
                    handles.Spikes.HandSort.Waveforms=handles.Spikes.Offline_Sorting.waveforms;
                end
                handles.Spikes.HandSort.samplingRate=handles.Spikes.Offline_Sorting.samplingRate;
            elseif isfield(handles.Spikes,'Online_Sorting') || strcmp(get(handles.Spikes_SortOn_Menu,'Checked'),'on')
                set(handles.Spikes_SortOff_Menu,'Checked','off');
                set(handles.Spikes_SortOn_Menu,'Checked','on');
                handles.Spikes.HandSort.Units=handles.Spikes.Online_Sorting.Units;
                handles.Spikes.HandSort.SpikeTimes=handles.Spikes.Online_Sorting.SpikeTimes;
                handles.Spikes.HandSort.Waveforms=handles.Spikes.Online_Sorting.Waveforms;
                handles.Spikes.HandSort.samplingRate=handles.Spikes.Online_Sorting.samplingRate;
            end
            % set data classes
            handles.Spikes.HandSort.SpikeTimes=...

        cellfun(@(spktimes) uint32(spktimes), handles.Spikes.HandSort.SpikeTimes,'UniformOutput',false);  %safe for about 40 hours of recording
            handles.Spikes.HandSort.Units=...
                cellfun(@(units) int8(units), handles.Spikes.HandSort.Units,'UniformOutput',false); %256 units per channel max
            handles.Spikes.HandSort.Waveforms=...
                cellfun(@(wfs) int16(wfs), handles.Spikes.HandSort.Waveforms,'UniformOutput',false); %256 units per channel max
            handles.Spikes.HandSort.samplingRate=uint32(handles.Spikes.HandSort.samplingRate); %uint16 should be enough, but just to be on the safe side
            %squeeze in case of ND array
            squeezeIt=cellfun(@(x) numel(size(x))>2, handles.Spikes.HandSort.Waveforms,'UniformOutput',true);
            if sum(squeezeIt)
                handles.Spikes.HandSort.Waveforms=...
                    cellfun(@(x) squeeze(x), handles.Spikes.HandSort.Waveforms,'UniformOutput',false);
                handles.Spikes.HandSort.Waveforms(~squeezeIt)=...
                    cellfun(@(x) transpose(x), handles.Spikes.HandSort.Waveforms(~squeezeIt),'UniformOutput',false);
            end
            % re-order along time scale if needed
            timeOrder=cellfun(@(times) sum(diff(times)<0),...
                handles.Spikes.HandSort.SpikeTimes,'UniformOutput',true);
            if sum(timeOrder)
                [~,timeIndices]=cellfun(@(times) sort(times),...
                    handles.Spikes.HandSort.SpikeTimes,'UniformOutput',false);
                handles.Spikes.HandSort.Units=cellfun(@(units,times) units(times),...
                    handles.Spikes.HandSort.Units,timeIndices,'UniformOutput',false);
                handles.Spikes.HandSort.SpikeTimes=cellfun(@(spktimes,times) spktimes(times),...
                    handles.Spikes.HandSort.SpikeTimes,timeIndices,'UniformOutput',false);
                handles.Spikes.HandSort.Waveforms=cellfun(@(wf,times) wf(times,:),...
                    handles.Spikes.HandSort.Waveforms,timeIndices,'UniformOutput',false);
            end
            % check if empty cells
            nonEmpty=~cellfun('isempty',handles.Spikes.HandSort.Units);
            % if negative unit IDs (-1), shift all to get 0 as lowest ID
            negUnits=cellfun(@(units) sum(units<0), handles.Spikes.HandSort.Units(nonEmpty),...
                'UniformOutput',true);
            if sum(negUnits)
                handles.Spikes.HandSort.Units(nonEmpty)=...
                    cellfun(@(units) units+abs(min(units)),...
                    handles.Spikes.HandSort.Units(nonEmpty),'UniformOutput',false);
            end
            % permute dimensions if not the right orientation
            permuteIt=cellfun(@(x,y) find(size(x)==length(y))==1,...
                handles.Spikes.HandSort.Waveforms(nonEmpty),...
                handles.Spikes.HandSort.Units(nonEmpty),'UniformOutput',true);
            if sum(permuteIt)
                handles.Spikes.HandSort.Waveforms(nonEmpty)=...
                    cellfun(@(x) transpose(x), handles.Spikes.HandSort.Waveforms(nonEmpty),...
                    'UniformOutput',false);
            end
            permuteIt=cellfun(@(x) size(x,1)<numel(x), handles.Spikes.HandSort.Units(nonEmpty),...
                'UniformOutput',true);
            if sum(permuteIt)         % check units are not right as well
                handles.Spikes.HandSort.Units(nonEmpty)=...
                    cellfun(@(x) transpose(x), handles.Spikes.HandSort.Units(nonEmpty),'UniformOutput',false);
            end
            
            % downsample spike times to 1 millisecond bins
            %     spikeTimes=handles.Spikes.HandSort.SpikeTimes;
            for chNum=1:size(handles.Spikes.HandSort.samplingRate,1)
                handles.Spikes.HandSort.samplingRate(chNum,2)=1000;
                if ~isempty(handles.Spikes.HandSort.SpikeTimes{chNum,1})
                    handles.Spikes.HandSort.SpikeTimes{chNum,2}=...
                        uint32(handles.Spikes.HandSort.SpikeTimes{chNum,1})/...
                        (handles.Spikes.HandSort.samplingRate(chNum,1)/...
                        handles.Spikes.HandSort.samplingRate(chNum,2));
                else
                    handles.Spikes.HandSort.SpikeTimes{chNum,2}=[];
                end
                
                %         spikeTimeIdx=zeros(1,size(Spikes.Offline_Threshold.data{ChExN,1},2));
                %         spikeTimeIdx(Spikes.Offline_Threshold.data{ChExN,1})=1;
                %
                %         binSize=1;
                %         numBin=ceil(max(spikeTimes{chNum})/...
                %             (handles.Spikes.HandSort.samplingRate(chNum,1)/...
                %             handles.Spikes.HandSort.samplingRate(chNum,2))/binSize);
                %         % binspikeTime = histogram(double(spikeTimes), numBin); %plots directly histogram
                %         [data,binEdges] = histcounts(double(spikeTimes{chNum}),...
                %             linspace(0,max(double(spikeTimes{chNum})),numBin));
                %         data(data>1)=1; %no more than 1 spike per ms
            end
        end
    end
    %% Set number of electrodes and units,
    if strcmp(get(handles.SelectElectrode_LB,'string'),'none')
        % when opening GUI, selects cluster with most spikes and
        % channel with most spikes from that cluster
        allUnits=single(vertcat(handles.Spikes.HandSort.Units{:}));
        [unitOccurence,uniqueUnitIDs]=hist(allUnits,unique(allUnits)); unitOccurence=unitOccurence(uniqueUnitIDs>0);
        uniqueUnitIDs=uniqueUnitIDs(uniqueUnitIDs>0);
        selectedUnit=int8(uniqueUnitIDs(unitOccurence==max(unitOccurence),1));
        % find all occurences of this unit across channels
        channelIdx=cellfun(@(x) logical(sum(ismember(x,selectedUnit))), handles.Spikes.HandSort.Units);
        occurencePerChannel=cellfun(@(x) sum(ismember(x,selectedUnit)), handles.Spikes.HandSort.Units);
        bestChannel=find(occurencePerChannel==max(occurencePerChannel),1);   
        if diff(size(handles.rec_info.exportedChan))>0
            handles.rec_info.exportedChan=handles.rec_info.exportedChan';
        end
        set(handles.SelectElectrode_LB,'string',num2str(handles.rec_info.exportedChan(channelIdx)));
        channelNum=find(bestChannel==handles.rec_info.exportedChan(channelIdx));
        if isfield(handles.Spikes,'Online_Sorting') || isfield(handles.Spikes,'Offline_Sorting')
            set(handles.SelectElectrode_LB,'value',channelNum);
        else
            set(handles.SelectElectrode_LB,'value',1)
        end
    else
        channelNum=get(handles.SelectElectrode_LB,'value');
    end
    
    %% color electrode by shanks
    if isfield(handles.rec_info,'probeID')
        try
            handles.userinfo=UserDirInfo;
        catch
            handles.userinfo=[];
        end
        if ~isfield(handles.userinfo,'probemap')
            if ~exist('probemaps', 'dir')
                handles.userinfo.probemap = uigetdir(cd,'select probe maps location');
                path(handles.userinfo.probemap,path)
            else
                allPathDirs=strsplit(path,pathsep);
                handles.userinfo.probemap=allPathDirs{find(cellfun(@(pathdir) contains(pathdir,'probemaps'),allPathDirs),1)};
            end
            mapping=load([handles.userinfo.probemap filesep handles.rec_info.probeID '.mat']);
            handles.rec_info.differentShanks=logical(bwlabel(mod([mapping.(handles.rec_info.probeID).Shank]+1,2)));
        else
            shanksID=[handles.rec_info.probeLayout(~cellfun('isempty',...
                {handles.rec_info.probeLayout.Electrode})).Shank];
            handles.rec_info.differentShanks=logical(bwlabel(mod(shanksID+1,2)));    
        end
        channelNum=ReturnElectrodes(handles.SelectElectrode_LB);   
        colorChannels=cell(length(channelNum),1);
        for elNum=1:length(channelNum) % ASSUMING ALL EXPORTED CHANNELS ARE SORTED !
            if handles.rec_info.differentShanks(channelNum(elNum))>0
                colorChannels(elNum)=cellfun(@(thatChannel) sprintf(['<HTML><BODY bgcolor="%s">'...
                    '<FONT color="%s">%s</FONT></BODY></HTML>'],... %size="+1"
                    'black','white', thatChannel),{num2str(channelNum(elNum))},'UniformOutput',false);
            else
                colorChannels(elNum)=cellfun(@(thatChannel) sprintf(['<HTML><BODY bgcolor="%s">'...
                    '<FONT color="%s">%s</FONT></BODY></HTML>'],... %size="+1"
                    'white','black', thatChannel),{num2str(channelNum(elNum))},'UniformOutput',false);
            end
        end
        set(handles.SelectElectrode_LB, 'String', colorChannels);
        set(handles.SelectElectrode_LB ,'ListboxTop',max([1 numel(channelNum)-3]));
    end
    
    % namestr = cellstr(get(hObject, 'String'));
    % validx = get(hObject, 'Value');
    % newstr = regexprep(namestr{validx}, '"red"','"green"');
    % namestr{validx} = newstr;
    % set(hObject, 'String', namestr);
    
    %% initialize variables
    unitsIdx=vertcat(handles.Spikes.HandSort.Units{:}); %handles.Spikes.HandSort.Units{channelNum};
    waveForms=horzcat(handles.Spikes.HandSort.Waveforms{:}); %handles.Spikes.HandSort.Waveforms{channelNum};
    if ~isempty(unitsIdx)
        % how many units on that Channel?
        unitsID=unique(unitsIdx); %number of clustered units
        set(handles.SelectUnit_LB,'string',num2str(unitsID));
        handles=ClassificationColor(handles);
        set(handles.SelectUnit_LB,'value',find(unitsID>0));
        
        %% take out big ouliers
        %     if find(size(waveForms)==length(unitsIdx))==1
        %         waveForms=waveForms';
        %     end
        WFmeanZ=mean(abs(zscore(single(waveForms))),2);
        % if more than one, plot it and keep
        if sum(WFmeanZ>6)>1
            figure('name', 'Artifacts','position',[30   500   500   400]);
            plot(waveForms(:,WFmeanZ>6)','linewidth',2.5); hold on;
            plot(mean(waveForms,2),'linewidth',2.5);
            legend({'Artifact','Mean waveforms'});
            title('Potential artifacts removed, mean sigma > 6');
        else
            channelNum=ReturnElectrodes(handles.SelectElectrode_LB); 
            handles.Spikes.HandSort.Units{channelNum(get(handles.SelectElectrode_LB,'value'))}...
                (WFmeanZ>6)=-9;%artifacts
        end
        
        % here's how this can work:
        %     Spikes detected from threshold are the benchmark (unit code = 0).
        %     May want to extract waveforms, but most important is the time.
        %     Sorted units (whatever the source) "color" those units. One spike per ms max.
        handles=Plot_Unsorted_WF(handles);
        if get(handles.ShowWF_CB,'value')
            handles=Plot_Sorted_WF(handles);
        else
            cla(handles.SortedUnits_Axes);
        end
        Plot_Mean_WF(handles);
        Plot_ISI(handles);
        Plot_ACG(handles);
        Plot_XCG(handles);
    else
        %         unitsID=[];
        cla(handles.SortedUnits_Axes);
        cla(handles.MeanSortedUnits_Axes);
        cla(handles.ISI_Axes);
        cla(handles.ACG_Axes);
        cla(handles.XCorr_Axes);
        set(handles.SelectUnit_LB,'string',num2str(0));
        set(handles.SelectUnit_LB,'value',1);
    end
    handles.fileLoaded=1;
end

%% Load raw traces function
function handles=LoadRawData(handles)
try
    if ~(isfield(handles,'datFile') && ~isempty(strfind(handles.datFile,'.mat')))
        handles.datFile=regexp(handles.spikeFile,'.+(?=_\w+.\w+$)','match');
        handles.datFile=[handles.datFile{:} '_raw.mat'];
        handles.datDir=cd;
    end
    handles.tracesInfo=whos('-file',fullfile(handles.datDir,handles.datFile));
    handles.tracesInfo=rmfield(handles.tracesInfo,...
        {'bytes','global','sparse','complex','nesting','persistent'});
    if strfind(fullfile(handles.datDir,handles.datFile),'nopp')
        handles.tracesInfo.preproc=0;
    else
        % if "raw" data has been pre-processed already, do
        % not process later
        handles.tracesInfo.preproc=1;
    end
    handles.traces = matfile(fullfile(handles.datDir,handles.datFile));
catch % try to load from .dat file
    if ~(isfield(handles,'datFile') && contains(fullfile(handles.datDir,handles.datFile),'.dat'))
        if isfield(handles,'offlineSort_SpikeFile')
            handles.datFile=regexp(handles.offlineSort_SpikeFile,'.+?(?=\.)','match');
            handles.datFile=[handles.datFile{1} '.dat'];
            handles.datDir=cd;
        else
            % ask user for raw data file
        end
    end
    handles.traces = memmapfile(fullfile(handles.datDir,handles.datFile),'Format','int16');
    handles.tracesInfo= struct('name','rawData',...
        'size',size(handles.traces.Data),...
        'numChan',numel(handles.rec_info.exportedChan),...
        'source','dat');
    try 
    %check in params file if data was filtered, and find threshold value
    fid  = fopen([fullfile(handles.datDir,handles.datFile(1:end-4)) '.params'],'r');
    catch 
       fid=-1;
    end
    if fid~=-1
        params=fread(fid,'*char')';
        if contains(regexp(params,'(?<=filter_done\s+=\s)\w+(?=\s)','match','once'), 'True')
            handles.tracesInfo.preproc=1;
        else
            handles.tracesInfo.preproc=0;
        end
        handles.rec_info.threshold=regexp(params,'(?<=spike_thresh\s+=\s)\w+(?=\s)','match','once');
        fclose(fid);
    else
        handles.tracesInfo.preproc=0;
    end
end
handles.tracesInfo.excerptSize=handles.rec_info.samplingRate/2; %1 second as default (-:+ around loc)
if isa(handles.traces,'memmapfile')
    handles.tracesInfo.excerptLocation=round(max(handles.tracesInfo.size)/2/numel(handles.rec_info.exportedChan));  %mid-recording as default
    set(handles.TW_slider,'max',max(handles.tracesInfo.size)/numel(handles.rec_info.exportedChan));
else
    handles.tracesInfo.excerptLocation=round(max(handles.tracesInfo.size)/2);
    set(handles.TW_slider,'max',max(handles.tracesInfo.size));
end
% set(handles.TW_slider,'sliderstep',[0.01 max([0.01,...
%     handles.tracesInfo.excerptSize/max(handles.tracesInfo.size)])]);
% plot "raw" (filtered) trace
DisplayTraces(handles);
% plot spike rasters
DisplayMarkers(handles);

%% Plot unsorted spikes
function handles=Plot_Unsorted_WF(handles)
channelNum=get(handles.SelectElectrode_LB,'value');
waveForms=handles.Spikes.HandSort.Waveforms{channelNum};
unitsIdx=handles.Spikes.HandSort.Units{channelNum};
samplingRate=handles.Spikes.HandSort.samplingRate(channelNum,1);
%% Plot unsorted spikes
numWFtoPlot=str2double(get(handles.ShowHowManyUWF_ET,'string'));
if sum(unitsIdx==0)>numWFtoPlot %then only plot subset of waveforms
    subset=find(unitsIdx==0);
    handles.subset{1}=subset(1:ceil(sum(unitsIdx==0)/numWFtoPlot):end);
else
    handles.subset{1}=find(unitsIdx==0);
end
axes(handles.UnsortedUnits_Axes); hold on;
cla(handles.UnsortedUnits_Axes);
set(handles.UnsortedUnits_Axes,'Visible','on');
if find(size(waveForms)==length(unitsIdx))==1
    waveForms=waveForms';
end
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
set(gca,'XTick',linspace(0,size(waveForms(:,handles.subset{1}),1),5),...
    'XTickLabel',round(linspace(-round(size(waveForms(:,handles.subset{1}),1)/2),...
    round(size(waveForms(:,handles.subset{1}),1)/2),5)/(double(samplingRate)/1000),2),'TickDir','out');
% ordinateLabels=str2num(get(gca,'YTickLabel'))/4;
% set(gca,'YTickLabel',num2str(ordinateLabels));
% legend('Unclustered waveforms','location','southeast')
axis('tight');box off;
xlabel('Time (ms)');
ylabelh=ylabel('Voltage (\muV)');
set(ylabelh,'interpreter','tex');
set(gca,'Color','white','FontSize',10,'FontName','Calibri');

%% Plot clusters
function handles=Plot_Sorted_WF(handles)
% get units's data
[selectedUnits,selectedUnitsListIdx,unitsIDs,~,waveForms,samplingRate]=GetUnitData(handles);
% selected unit ids
axes(handles.SortedUnits_Axes); hold on;%colormap lines; cmap=colormap;
cla(handles.SortedUnits_Axes);
set(handles.SortedUnits_Axes,'Visible','on');

numWFtoPlot=str2double(get(handles.ShowHowManyWF_ET,'string'));
numUnits=nan(1,length(selectedUnits));
for unitP=1:length(selectedUnits)
    numUnits(unitP)=sum(unitsIDs==selectedUnits(unitP));
    %if there are too many waveforms to plot
    if sum(unitsIDs==selectedUnits(unitP))>numWFtoPlot %then only plot subset of waveforms
        subset=find(unitsIDs==selectedUnits(unitP));
        handles.subset{selectedUnitsListIdx(unitP)}=subset(1:round(numUnits(unitP)/numWFtoPlot):end);
    else
        handles.subset{selectedUnitsListIdx(unitP)}=find(unitsIDs==selectedUnits(unitP));
    end
    %make sure waveforms can be plotted
    if find(size(waveForms)==length(unitsIDs))==1
        waveForms=waveForms'; %figure; hold on; plot(waveForms(:,1:10))
    end
    plot(waveForms(:,handles.subset{selectedUnitsListIdx(unitP)}),...
        'linewidth',1,'Color',[handles.cmap(selectedUnits(unitP),:),0.4],...
        'Tag',num2str(selectedUnits(unitP)));
    %Change waveform tags from cluster number to unit ID
    lineH=flipud(findobj(gca,'Type', 'line'));
    skippedTags=0;
    if unitP==1 %define keepTags
        keepTags=cell(size(lineH,1),1);
    end
    for lineTag=1:size(lineH,1)
        % check if it's the right cluster
        if strcmp(lineH(lineTag).Tag,num2str(selectedUnits(unitP)))
            %record Tags
            keepTags{lineTag,1}=num2str(handles.subset{selectedUnitsListIdx(unitP)}(lineTag-skippedTags)); %Tag the unit ID
        else
            skippedTags=skippedTags+1;
        end
    end
    %     figure;hold on
    %     plot(lineH(lineTag).YData);
    %     plot(waveForms(:,handles.subset{selectedUnitsListIdx(unitP)}(lineTag-passedTags)))
end
%apply Tags
if exist('lineH','var')
    for lineTag=1:size(lineH,1)
        lineH(lineTag).Tag=keepTags{lineTag};
    end
    firstOfClus=[1, find(diff(cellfun(@(col) sum(col), {lineH.Color})))+1];
    legH=legend(lineH(firstOfClus),{num2str(numUnits')},'location','southeast');
    set(legH,'Box','Off','FontSize',7,'LineWidth',0.2,'Position',...
        [0.54 legH.Position(2)+legH.Position(3)-0.01 0 0]);
end
set(gca,'Ylim',[min(min(waveForms)) max(max(waveForms))],'XTick',linspace(0,...
    size(waveForms(:,handles.subset{selectedUnitsListIdx(unitP)}),1),5),...
    'XTickLabel',round(linspace(-round(size(waveForms,1)/2),...
    round(size(waveForms,1)/2),5)/(double(samplingRate)/1000),2),...
    'TickDir','out');
axis('tight');box off;
xlabel('Time (ms)');
ylabel('Voltage (\muV)');
set(gca,'Color','white','FontSize',10,'FontName','Calibri');
hold off

%% Plot mean waveforms
function Plot_Mean_WF(handles)
% get units's data
[selectedUnits,selectedUnitsListIdx,unitsIDs,~,waveForms,samplingRate]=GetUnitData(handles);
% channelNum=get(handles.SelectElectrode_LB,'value');
% waveForms=handles.Spikes.HandSort.Waveforms{channelNum};
% unitsIdx=handles.Spikes.HandSort.Units{channelNum};
% uniqueUnitIDs=unique(unitsIDs);
axes(handles.MeanSortedUnits_Axes); hold on;%colormap lines;
cla(handles.MeanSortedUnits_Axes);
set(handles.MeanSortedUnits_Axes,'Visible','on');

numWFtoPlot=str2double(get(handles.ShowHowManyWF_ET,'string'));
for unitP=1:length(selectedUnits)
    if selectedUnits(unitP)==0
        continue
    end
    %if subset is not defined
    if size(handles.subset,2)<selectedUnitsListIdx(unitP) || isempty(handles.subset{selectedUnitsListIdx(unitP)})
        if sum(unitsIDs==selectedUnits(unitP))>numWFtoPlot %then only plot subset of waveforms
            subset=find(unitsIDs==selectedUnits(unitP));
            handles.subset{selectedUnitsListIdx(unitP)}=subset(1:round(sum(unitsIDs==selectedUnits(unitP))/numWFtoPlot):end);
        else
            handles.subset{selectedUnitsListIdx(unitP)}=find(unitsIDs==selectedUnits(unitP));
        end
    end
    selectWF=single(waveForms(:,unitsIDs==selectedUnits(unitP))');
    if ~isnan(mean(selectWF))
        lineh(unitP)=plot(mean(selectWF),'linewidth',2,'Color',[handles.cmap(selectedUnits(unitP),:),0.7]);
        wfSEM=std(selectWF)/ sqrt(size(selectWF,2)); %standard error of the mean
        wfSEM = wfSEM * 1.96; % 95% of the data will fall within 1.96 standard deviations of a normal distribution
        patch([1:length(wfSEM),fliplr(1:length(wfSEM))],...
            [mean(selectWF)-wfSEM,fliplr(mean(selectWF)+wfSEM)],...
            handles.cmap(selectedUnits(unitP),:),'EdgeColor','none','FaceAlpha',0.2);
        %duplicate mean unit waveform over unsorted plot
        %             plot(handles.UnsortedUnits_Axes,mean(selectWF),'linewidth',2,'Color',handles.cmap(unitP,:));
        if unitP==1
            delete(findobj(handles.UnsortedUnits_Axes,'Type', 'patch'));
        end
        patch([1:length(wfSEM),fliplr(1:length(wfSEM))],...
            [mean(selectWF)-wfSEM,fliplr(mean(selectWF)+wfSEM)],...
            handles.cmap(selectedUnits(unitP),:),'EdgeColor','none','FaceAlpha',0.5,'Parent', handles.UnsortedUnits_Axes);
        set(handles.UnsortedUnits_Axes,'Color','white','FontSize',10,'FontName','Calibri');
    end
end
set(gca,'XTick',linspace(0,size(waveForms(:,handles.subset{selectedUnitsListIdx(unitP)}),1),5),...
    'XTickLabel',round(linspace(-round(size(waveForms,1)/2),...
    round(size(waveForms,1)/2),5)/(double(samplingRate)/1000),2),'TickDir','out');
% ordinateLabels=str2num(get(gca,'YTickLabel'))/4;
% set(gca,'YTickLabel',num2str(ordinateLabels));
if exist('selectWF','var')
    legH=legend(lineh,{num2str(selectedUnits)},'location','southeast');
    set(legH,'Box','Off');
end
axis('tight');box off;
xlabel('Time (ms)');
ylabel('Voltage (\muV)');
set(gca,'Color','white','FontSize',10,'FontName','Calibri');
hold off

%% Plot traces
function dataExcerpt=DisplayTraces(handles)
channelIndex=get(handles.SelectElectrode_LB,'value');
channelNum=ReturnElectrodes(handles.SelectElectrode_LB);
channelNum=channelNum(channelIndex);
if isa(handles.traces,'memmapfile')
    %     exwdw=[((30000*106)-15000)*28:((30000*106)+15000)*28-1];
    % % dataExcerpt=rawData.Data(exwdw);
    % dataExcerpt=handles.traces.Data(exwdw);
    % dataExcerpt=reshape(dataExcerpt,[28 30000]);
    % figure;
    % plot(dataExcerpt(10,:))
    % title('el10 106')
    winIdxStart=((handles.tracesInfo.excerptLocation-...
        double(handles.tracesInfo.excerptSize))*numel(handles.rec_info.exportedChan))+1;
    if mod(winIdxStart,numel(handles.rec_info.exportedChan))~=1 %set window index to correct point in data vector
        winIdxStart=winIdxStart-...
            mod(winIdxStart,numel(handles.rec_info.exportedChan))-...
            numel(handles.rec_info.exportedChan)+1;    % set index loc to first electrode
        %             (numel(handles.rec_info.exportedChan) - channelNum);     % set index loc to selected electrode
    end
    winIdxEnd=winIdxStart+...
        (2*handles.tracesInfo.excerptSize*numel(handles.rec_info.exportedChan));
    excerptWindow=winIdxStart:winIdxEnd-1;
    if size(excerptWindow,2)>(2*handles.tracesInfo.excerptSize*...
            numel(handles.rec_info.exportedChan)) %for some reason
        excerptWindow=excerptWindow(1:end-(size(excerptWindow,2)-...
            (2*handles.tracesInfo.excerptSize*numel(handles.rec_info.exportedChan))));
    end
    dataExcerpt=handles.traces.Data(excerptWindow);
    dataExcerpt=reshape(dataExcerpt,[numel(handles.rec_info.exportedChan)...
        handles.tracesInfo.excerptSize*2]);
    %         foo=handles.traces.Data;
    %         foo=reshape(foo,[numel(handles.rec_info.exportedChan)...
    %         size(foo,1)/numel(handles.rec_info.exportedChan)]);
else
    excerptWindow=handles.tracesInfo.excerptLocation-...
        handles.tracesInfo.excerptSize:handles.tracesInfo.excerptLocation+handles.tracesInfo.excerptSize-1;
    dataExcerpt=handles.traces.(handles.tracesInfo.name)(:,excerptWindow);
end
if handles.tracesInfo.preproc==0 % raw data is presumed bandpassed filtered at this point
    preprocOption={'CAR','all'};
    dataExcerpt=PreProcData(dataExcerpt,handles.rec_info.samplingRate,preprocOption);
end
dataExcerpt=int32(dataExcerpt(channelNum,:));
axes(handles.TimeRaster_Axes);
cla(handles.TimeRaster_Axes);
set(handles.TimeRaster_Axes,'Visible','on');
plot(handles.TimeRaster_Axes,dataExcerpt,'k','linewidth',0.1); hold on;
% threshold
if isfield(handles.rec_info,'threshold')
    threshold=str2double(handles.rec_info.threshold);
    plot(ones(1,size(dataExcerpt,2))*threshold*mad(single(dataExcerpt)),'--','Color',[0 0 0 0.3]);
    plot(ones(1,size(dataExcerpt,2))*-threshold*mad(single(dataExcerpt)),'--','Color',[0 0 0 0.3]);
else
    %     threshold=8;
end
hold off;
timeLabels=round(linspace(round(handles.tracesInfo.excerptLocation-...
    handles.tracesInfo.excerptSize)/(handles.rec_info.samplingRate/1000),...
    round(handles.tracesInfo.excerptLocation+handles.tracesInfo.excerptSize)/...
    (handles.rec_info.samplingRate/1000),4)./1000,3); % duration(X,'Format','h')
set(handles.TimeRaster_Axes,'xtick',round(linspace(0,max(get(handles.TimeRaster_Axes,'xtick')),4)),...
    'xticklabel',timeLabels); %,'TickDir','out');
set(handles.TimeRaster_Axes,'ytick',[],'yticklabel',[]); %'ylim'
axis('tight');box off;
set(handles.TimeRaster_Axes,'Color','white','FontSize',12,'FontName','calibri');
% if isa(handles.traces,'memmapfile')
%     set(handles.TW_slider,'value',handles.tracesInfo.excerptLocation/numel(handles.rec_info.exportedChan));
% else
set(handles.TW_slider,'value',handles.tracesInfo.excerptLocation);

%% Plot spike unit markers
function spkTimes=DisplayMarkers(handles)
% display color-coded markers for each identified spike
channelIndex=get(handles.SelectElectrode_LB,'value');
channelNum=ReturnElectrodes(handles.SelectElectrode_LB);
channelNum=channelNum(channelIndex);
if contains(get(handles.DisplayNtrodeMarkers_MenuItem,'Checked'),'on')
    % get channel values from same Ntrode
    shankNum=cumsum(logical([0 diff(handles.rec_info.differentShanks(1:handles.rec_info.numRecChan))]));
    channelNum=find(shankNum==shankNum(channelNum));
end
axes(handles.TimeRaster_Axes); hold on
for elNum=1:length(channelNum)
    % get which unit to plot
    if get(handles.ShowAllUnits_RB,'value') || length(channelNum)>1
        unitID=unique(handles.Spikes.HandSort.Units{channelNum(elNum), 1});%  unitID=ReturnUnits(handles.SelectUnit_LB);
        selectedUnitsListIdx=find(unitID>=0);
        selectedUnits=unitID(selectedUnitsListIdx);
    else
        unitID=ReturnUnits(handles.SelectUnit_LB);
        selectedUnitsListIdx=get(handles.SelectUnit_LB,'value');
        selectedUnits=unitID(selectedUnitsListIdx);
    end
    spkTimes=cell(size(selectedUnits,1),1);
    if isfield(handles.Spikes,'Online_Sorting') % plot above trace
        if ~isempty(handles.Spikes.Online_Sorting.SpikeTimes)
            for unitP=1:size(selectedUnits,1)
                spkTimes{unitP}=handles.Spikes.Online_Sorting.SpikeTimes{channelNum(elNum)}(...
                    (handles.Spikes.Online_Sorting.SpikeTimes{channelNum(elNum)}>=...
                    handles.tracesInfo.excerptLocation-handles.tracesInfo.excerptSize) &...
                    (handles.Spikes.Online_Sorting.SpikeTimes{channelNum(elNum)}<...
                    handles.tracesInfo.excerptLocation+handles.tracesInfo.excerptSize) &...
                    handles.Spikes.Online_Sorting.Units{channelNum(elNum)}==unitID(selectedUnitsListIdx(unitP)));
                if ~isempty(spkTimes{unitP})
                    rasterHeight=ones(1,size(spkTimes{unitP},2))*max(get(gca,'ylim'))/4*3;
                    wfWidthComp=round(size(handles.Spikes.Online_Sorting.Waveforms{channelNum(elNum)},1)); %will substract wf width to raster times
                    plot(spkTimes{unitP}-(handles.tracesInfo.excerptLocation-handles.tracesInfo.excerptSize)-wfWidthComp,...
                        rasterHeight,'Color',[handles.cmap(unitID(selectedUnitsListIdx(unitP)),:),0.4],...
                        'linestyle','none','Marker','v');
                end
            end
        end
    end
    if isfield(handles.Spikes,'Offline_Sorting') % plot below trace
        for unitP=1:size(selectedUnits,1)
            spkTimes{unitP}=handles.Spikes.Offline_Sorting.SpikeTimes{channelNum(elNum)}(...
                (handles.Spikes.Offline_Sorting.SpikeTimes{channelNum(elNum)}>=...
                handles.tracesInfo.excerptLocation-handles.tracesInfo.excerptSize) &...
                (handles.Spikes.Offline_Sorting.SpikeTimes{channelNum(elNum)}<...
                handles.tracesInfo.excerptLocation+handles.tracesInfo.excerptSize) &...
                handles.Spikes.Offline_Sorting.Units{channelNum(elNum)}==unitID(selectedUnitsListIdx(unitP)));
            if ~isempty(spkTimes{unitP})
                rasterHeight=ones(1,size(spkTimes{unitP},2))*(min(get(gca,'ylim'))/4*3);
                if unitID(selectedUnitsListIdx(unitP))==0 %"garbage spikes"
                    plot(single(spkTimes{unitP})-(handles.tracesInfo.excerptLocation-handles.tracesInfo.excerptSize),...
                        rasterHeight,'Color','k',...
                        'linestyle','none','Marker','*');
                else
                    plot(single(spkTimes{unitP})-(handles.tracesInfo.excerptLocation-handles.tracesInfo.excerptSize),...
                        rasterHeight,'Color',[handles.cmap(unitID(selectedUnitsListIdx(unitP)),:),0.4],...
                        'linestyle','none','Marker','^');
                end
            end
        end
    end
    if strcmp(get(handles.Spikes_CurrentVersion_Menu,'Checked'),'on') % also plot below trace, but circles
        for unitP=1:size(selectedUnits,1)
            spkTimes{unitP}=handles.Spikes.HandSort.SpikeTimes{channelNum(elNum)}(...
                (handles.Spikes.HandSort.SpikeTimes{channelNum(elNum)}>=...
                handles.tracesInfo.excerptLocation-handles.tracesInfo.excerptSize) &...
                (handles.Spikes.HandSort.SpikeTimes{channelNum(elNum)}<...
                handles.tracesInfo.excerptLocation+handles.tracesInfo.excerptSize) &...
                handles.Spikes.HandSort.Units{channelNum(elNum)}==unitID(selectedUnitsListIdx(unitP)));
            if ~isempty(spkTimes{unitP})
                rasterHeight=ones(1,size(spkTimes{unitP},2))*(min(get(gca,'ylim'))/4*3);
                if unitID(selectedUnitsListIdx(unitP))==0 %"garbage spikes"
                    plot(spkTimes{unitP}-(handles.tracesInfo.excerptLocation-handles.tracesInfo.excerptSize),...
                        rasterHeight,'Color','k',...
                        'linestyle','none','Marker','*');
                else
                    plot(spkTimes{unitP}-(handles.tracesInfo.excerptLocation-handles.tracesInfo.excerptSize),...
                        rasterHeight,'Color',[handles.cmap(unitID(selectedUnitsListIdx(unitP)),:),0.4],...
                        'linestyle','none','Marker','o');
                end
            end
        end
    end
end
hold off

%% --- Executes on slider movement.
function TW_slider_Callback(hObject, ~, handles)
handles.tracesInfo.excerptLocation=double(round(get(handles.TW_slider,'value')/...
    handles.rec_info.samplingRate)*handles.rec_info.samplingRate); %round to nearest second
if handles.tracesInfo.excerptLocation-handles.tracesInfo.excerptSize<1
    handles.tracesInfo.excerptLocation=handles.tracesInfo.excerptSize+1;
elseif handles.tracesInfo.excerptLocation+handles.tracesInfo.excerptSize>get(handles.TW_slider,'max')
    handles.tracesInfo.excerptLocation=get(handles.TW_slider,'max')-handles.tracesInfo.excerptSize;
end
% plot "raw" (filtered) trace
DisplayTraces(handles);
% plot spike rasters
DisplayMarkers(handles);
% update handles
guidata(hObject, handles);

%% --- Executes on button press in TWplus_PB.
function TWplus_PB_Callback(hObject, ~, handles)
handles.tracesInfo.excerptSize=handles.tracesInfo.excerptSize*2;
mem=memory;
if ((handles.tracesInfo.excerptSize*2+1)*2)/(0.5*mem.MemAvailableAllArrays)>1 %handles.tracesInfo.excerptSize>
    handles.tracesInfo.excerptSize=floor((0.5*mem.MemAvailableAllArrays)/2); % handles.tracesInfo.size(1);
end
if handles.tracesInfo.excerptSize*2+1>handles.tracesInfo.size(1)
    handles.tracesInfo.excerptSize=floor(handles.tracesInfo.size(1)/2);
end
% plot "raw" (filtered) trace
DisplayTraces(handles);
% plot spike rasters
DisplayMarkers(handles);
%  Update handles structure
guidata(hObject, handles);

%% --- Executes on button press in TWminus_PB.
function TWminus_PB_Callback(hObject, ~, handles)
handles.tracesInfo.excerptSize=round(handles.tracesInfo.excerptSize/2);
if handles.tracesInfo.excerptSize<handles.rec_info.samplingRate/100
    handles.tracesInfo.excerptSize=handles.rec_info.samplingRate/100;
end
% plot "raw" (filtered) trace
DisplayTraces(handles);
% plot spike rasters
DisplayMarkers(handles);
%  Update handles structure
guidata(hObject, handles);

%% --- Executes on button press in TWall_PB.
function TWall_PB_Callback(hObject, ~, handles)
handles.tracesInfo.excerptSize=(handles.rec_info.dur/2)-1;
mem=memory;
if ((handles.tracesInfo.excerptSize*2+1)*2)/(0.5*mem.MemAvailableAllArrays)>1 %handles.tracesInfo.excerptSize>
    handles.tracesInfo.excerptSize=floor((0.5*mem.MemAvailableAllArrays)/2); % handles.tracesInfo.size(1);
end
% plot "raw" (filtered) trace
DisplayTraces(handles);
% plot spike rasters
DisplayMarkers(handles);
%  Update handles structure
guidata(hObject, handles);

%% Plot ISI
function Plot_ISI(handles)
% get which unit to plot
[selectedUnits,~,unitsIDs,spikeTimes,~,samplingRate]=GetUnitData(handles);
if isempty(selectedUnits) || sum(selectedUnits)==0
    cla(handles.ISI_Axes);
    return
end
%keep the most numerous if more than one
if diff(size(selectedUnits))<0
    keepU=1;
    for uidx=1:size(selectedUnits,1)
        if sum(unitsIDs==selectedUnits(uidx))>sum(unitsIDs==selectedUnits(keepU))
            keepU=uidx;
        end
    end
    selectedUnits=selectedUnits(keepU);
end
%spike times for that unit
unitST=spikeTimes(unitsIDs==selectedUnits);
% compute interspike interval
if ~isempty(diff(unitST))
    ISI=diff(unitST)/(samplingRate/1000);
    axes(handles.ISI_Axes); hold on;
    cla(handles.ISI_Axes,'reset');
    set(handles.ISI_Axes,'Visible','on');
    ISIhist=histogram(double(ISI),logspace(0, 4, 50),'DisplayStyle','stairs','LineWidth',1.5);  %,'Normalization','probability'
%     ISIhist.FaceColor = handles.cmap(unitID(unitID==selectedUnits),:);
    ISIhist.EdgeColor = handles.cmap(selectedUnits,:); %'k';
    xlabel('Interspike Interval (ms)')
    axis('tight');box off; grid('on'); set(gca,'xscale','log','GridAlpha',0.25,'MinorGridAlpha',1);
    set(gca,'xlim',[0 10^4],... %'XTick',linspace(0,40,5),'XTickLabel',linspace(0,40,5),...
        'TickDir','out','Color','white','FontSize',10,'FontName','Calibri');
    hold off
end

%% Plot autocorrelogram
function Plot_ACG(handles)
% get which unit to plot
[selectedUnits,~,unitsIDs,spikeTimes,~,samplingRate]=GetUnitData(handles);
if isempty(selectedUnits) || sum(selectedUnits)==0
    cla(handles.ACG_Axes);
    return
end
%keep the most numerous if more than one
if diff(size(selectedUnits))<0
    keepU=1;
    for uidx=1:size(selectedUnits,1)
        if sum(unitsIDs==selectedUnits(uidx))>sum(unitsIDs==selectedUnits(keepU))
            keepU=uidx;
        end
    end
    selectedUnits=selectedUnits(keepU);
end
%get unit spike times
unitST=spikeTimes(unitsIDs==selectedUnits);
% change to ms timescale
unitST=unitST/(samplingRate/1000);
%get ISI
% ISI=diff(unitST)/(samplingRate/1000);
%bin
spikeTimeIdx=zeros(1,unitST(end));
spikeTimeIdx(unitST)=1;
binSize=1;
numBin=ceil(size(spikeTimeIdx,2)/binSize);
binUnits = histcounts(double(unitST), linspace(0,size(spikeTimeIdx,2),numBin));
binUnits(binUnits>1)=1; %no more than 1 spike per ms
% compute autocorrelogram
[ACG,lags]=xcorr(double(binUnits),200,'unbiased'); %'coeff'
ACG(lags==0)=0;
axes(handles.ACG_Axes); hold on;
cla(handles.ACG_Axes,'reset');
set(handles.ACG_Axes,'Visible','on');
ACGh=bar(lags,ACG);
ACGh.FaceColor = handles.cmap(selectedUnits,:);
ACGh.EdgeColor = 'none';
% axis('tight');
box off; grid('on'); %set(gca,'yscale','log','GridAlpha',0.25,'MinorGridAlpha',1);
xlabel('Autocorrelogram (1 ms bins)')
set(gca,'xlim',[-20 20],... %'ylim',[0 max([max(get(gca,'ylim')) 10^1])]
    'Color','white','FontSize',10,'FontName','Calibri','TickDir','out');
hold off

%% Plot cross-correlogram
function Plot_XCG(handles)
% get which unit to plot
[selectedUnits,~,unitsIDs,spikeTimes,~,samplingRate]=GetUnitData(handles);
if isempty(selectedUnits) || length(selectedUnits)~=2
    cla(handles.XCorr_Axes);
    return
end
%get units spike times
unitST{1}=spikeTimes(unitsIDs==selectedUnits(1));
unitST{2}=spikeTimes(unitsIDs==selectedUnits(2));
% change to ms timescale
unitST{1}=unitST{1}/(samplingRate/1000);
unitST{2}=unitST{2}/(samplingRate/1000);

%bin
spikeTimeIdx{1}=zeros(1,unitST{1}(end));
spikeTimeIdx{1}(unitST{1})=1;
spikeTimeIdx{2}=zeros(1,unitST{2}(end));
spikeTimeIdx{2}(unitST{2})=1;
binSize=5;
numBin=max(ceil(size(spikeTimeIdx{1},2)/binSize),ceil(size(spikeTimeIdx{2},2)/binSize));
binUnits{1} = histcounts(double(unitST{1}), linspace(0,size(spikeTimeIdx{1},2),numBin));
binUnits{1}(binUnits{1}>1)=1; %no more than 1 spike per ms
binUnits{2} = histcounts(double(unitST{2}), linspace(0,size(spikeTimeIdx{2},2),numBin));
binUnits{2}(binUnits{2}>1)=1; %no more than 1 spike per ms

% compute autocorrelogram
[XCG,lags]=xcorr(double(binUnits{1}),double(binUnits{2}),200,'unbiased'); %'coeff'
XCG(lags==0)=0;
axes(handles.XCorr_Axes); hold on;
cla(handles.XCorr_Axes);
set(handles.XCorr_Axes,'Visible','on');
XCGh=bar(lags,XCG);
XCGh.FaceColor = (handles.cmap(selectedUnits(1),:)+...
    handles.cmap(selectedUnits(2),:))/2;
XCGh.EdgeColor = 'none';
axis('tight');box off;
xlabel('CrossCorrelogram (5 ms bins)')
set(gca,'xlim',[-50 50],'Color','white','FontSize',10,'FontName','Calibri','TickDir','out');
hold off

% function  Plot_Raster_TW(handles)
% %% plot rasters
% channelNum=get(handles.SelectElectrode_LB,'value');
% spikeTimes=handles.Spikes.HandSort.SpikeTimes{channelNum,2};

% plot 10 sec or numWFtoPlot waveforms max

%% --- Executes on mouse press over unsorted units axes
function UnsortedUnits_Axes_ButtonDownFcn(hObject, ~, handles)
% left click to start selection line
% right click to end it

channelNum=get(handles.SelectElectrode_LB,'value');

%% initialize variables
unitsIdx=find(handles.Spikes.HandSort.Units{channelNum}==0);
% waveForms=handles.Spikes.HandSort.Waveforms{channelNum};
% spikeTimes=handles.Spikes.HandSort.SpikeTimes{channelNum,2};
% samplingRate=handles.Spikes.HandSort.samplingRate(channelNum,1);

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
% find similar waveforms that were not plotted
% if numel(find(handles.Spikes.HandSort.Units{channelNum}==0))>numel(clusterClasses) &...
%         sum(logical(clusterClasses))
%     meanSelectedWF=mean(handles.Spikes.HandSort.Waveforms{channelNum}(:,...
%         unitsIdx(logical(clusterClasses))),2);
%     allWF=handles.Spikes.HandSort.Waveforms{channelNum}(:,unitsIdx);
%     for wfNum=1:size(allWF,2)
%         ccVal(wfNum)=median(xcorr(double(allWF(:,wfNum)'),...
%             double(meanSelectedWF')),2);
%     end
%     figure;
%     plot(allWF(:,unitsIdx(logical(clusterClasses)))')
% end

handles.Spikes.HandSort.Units{channelNum}(unitsIdx(logical(clusterClasses)))=...
    clusterClasses(logical(clusterClasses));
unitsID=unique(handles.Spikes.HandSort.Units{channelNum});
%Check if unit selection still works
unitSelection=get(handles.SelectUnit_LB,'Value');
try
    unitsID(get(handles.SelectUnit_LB,'Value'));
catch
    unitSelection=unitSelection(ismember(unitSelection,(1:numel(unitsID))));
end
unitSelection(unitsID(unitSelection)==0)=unitSelection(unitsID(unitSelection)==0)+1;
set(handles.SelectUnit_LB,'Value',unitSelection);
set(handles.SelectUnit_LB,'String',num2str(unitsID(unitsID>=0)));
if max(get(handles.SelectUnit_LB,'value'))>length(ReturnUnits(handles.SelectUnit_LB))
    newSelection=get(handles.SelectUnit_LB,'value')-...
        (max(get(handles.SelectUnit_LB,'value'))-length(ReturnUnits(handles.SelectUnit_LB)));
    newSelection=newSelection(newSelection>0);
    set(handles.SelectUnit_LB,'value',newSelection);
end
if find(clusterClasses>0,1)
    if get(handles.ShowWF_CB,'value')
        handles=Plot_Sorted_WF(handles);
    else
        cla(handles.SortedUnits_Axes);
    end
    Plot_Mean_WF(handles);
    Plot_ISI(handles);
    Plot_ACG(handles);
    Plot_XCG(handles);
end
%  Update handles structure
guidata(hObject, handles);

%% --- Executes on mouse press over sorted units axes
function SortedUnits_Axes_ButtonDownFcn(hObject, ~, handles)
channelNum=get(handles.SelectElectrode_LB,'value');

%% initialize variables
unitID=ReturnUnits(handles.SelectUnit_LB);
% for uIdxNum=1:length(unitsID)
%     unitsIdx{uIdxNum}=find(handles.Spikes.HandSort.Units{channelNum}==unitsID(uIdxNum));
% end
unitsIdx=handles.Spikes.HandSort.Units{channelNum};

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
% foo=handles.Spikes.HandSort.Waveforms{channelNum}; foo=foo';
% figure;plot(foo(unitsIdx(logical(clusterClasses)),:)');hold on
% plot(lineH(flip(logical(clusterClasses))).YData)
%
% waveForms=handles.Spikes.HandSort.Waveforms{channelNum};
% figure; plot(waveForms(:,linesTags(lineSelecIdx)))
if lineSelecIdx==0
    return
else
    set(handles.Spikes_OriginalVersion_Menu,'Checked','off');
    set(handles.Spikes_CurrentVersion_Menu,'Checked','on');
end
handles.Spikes.HandSort.Units{channelNum}(linesTags(lineSelecIdx))=...
    clusterClasses(lineSelecIdx);
unitsID=unique(handles.Spikes.HandSort.Units{channelNum});
%Check if unit selection still works
unitSelection=get(handles.SelectUnit_LB,'Value');
try
    unitsID(get(handles.SelectUnit_LB,'Value'));
catch
    unitSelection=unitSelection(ismember(unitSelection,(1:numel(unitsID))));
end
unitSelection(unitsID(unitSelection)==0)=unitSelection(unitsID(unitSelection)==0)+1;
set(handles.SelectUnit_LB,'Value',unitSelection);
set(handles.SelectUnit_LB,'String',num2str(unitsID(unitsID>=0)))
if find(clusterClasses>0,1)
    if get(handles.ShowWF_CB,'value')
        handles=Plot_Sorted_WF(handles);
    else
        cla(handles.SortedUnits_Axes);
    end
end
if find(clusterClasses==0,1)
    handles=Plot_Unsorted_WF(handles);
end
Plot_Mean_WF(handles);
Plot_ISI(handles);
Plot_ACG(handles);
Plot_XCG(handles);
%  Update handles structure
guidata(hObject, handles);

%% --- Executes on mouse press over mean sorted units axes.
function MeanSortedUnits_Axes_ButtonDownFcn(hObject, ~, handles)
unitID=ReturnUnits(handles.SelectUnit_LB);
channelNum=get(handles.SelectElectrode_LB,'value');

%% initialize variables
unitsIdx=handles.Spikes.HandSort.Units{channelNum};

lineH=findobj(gca,'Type', 'line');
visibleLines=cellfun(@(x) strcmp(x,'on'), {lineH.Visible});

waveForms=fliplr(reshape([lineH(visibleLines).YData],...
    size([lineH(visibleLines).YData],2)/size(lineH(visibleLines),...
    1),size(lineH(visibleLines),1)));
waveForms=waveForms';%one waveform per row

[clusterClasses,viewClasses]=deal(flip(cellfun(@(x) str2double(x), {lineH(visibleLines).DisplayName})));
% figure;hold on
% plot(waveForms(1,:)','r')
% plot(waveForms(3,:)','b')

[newClasses,lineSelecIdx,groupReclass]=InteractiveClassification(waveForms,clusterClasses,viewClasses); % viewClasses=0
% foo=handles.Spikes.HandSort.Waveforms{channelNum}; foo=foo';
% figure;plot(foo(unitsIdx(logical(clusterClasses)),:)');hold on
% plot(lineH(flip(logical(clusterClasses))).YData)
%
% waveForms=handles.Spikes.HandSort.Waveforms{channelNum};
% figure; plot(waveForms(:,linesTags(lineSelecIdx)))

if lineSelecIdx==0
    return
else
    set(handles.Spikes_OriginalVersion_Menu,'Checked','off');
    set(handles.Spikes_CurrentVersion_Menu,'Checked','on');
end

changeUnits=clusterClasses(clusterClasses~=unique(newClasses(lineSelecIdx)) & lineSelecIdx');
for chgu=1:length(changeUnits)
    handles.Spikes.HandSort.Units{channelNum}(unitsIdx==changeUnits(chgu))=...
        unique(newClasses(lineSelecIdx));
end

if groupReclass==1 && sum(lineSelecIdx)==1
    
    selectedWF=double(waveForms(lineSelecIdx,:));
    WFpeak=find(abs(selectedWF)==max(abs(selectedWF)),1);
    wfRange=WFpeak-9:WFpeak+8;
    selectedWF=selectedWF(wfRange);
    allWf=handles.Spikes.HandSort.Waveforms{channelNum}; allWf=double(allWf(wfRange,:)');
    
    % figure; hold on
    % plot(selectedWF);
    allWf_wRef=[selectedWF;allWf];
    minDif=allWf_wRef-repmat(selectedWF,size(allWf_wRef,1),1);
    % figure; hist(median(abs(minDif')),20)
    % minResidual=find(median(abs(minDif'))==min(median(abs(minDif'))),1);
    maxPeaks=zscore(allWf_wRef')'.*repmat(selectedWF,size(allWf_wRef,1),1);
    % maxMultiplied=find(max(abs(maxPeaks'))==max(max(abs(maxPeaks'))),1);
    % plot(allWf_wRef(minResidual,:))
    % plot(allWf_wRef(maxMultiplied,:))
    
    roughCat=([median(abs(minDif'))' max(abs(maxPeaks'))']);
    [~,PrComps] = pca(roughCat);
    % figure
    % scatter(PrComps(:,1), PrComps(:,2), 'k.');
    clusterCat = gmdistribution(PrComps(1,:),...
        [abs(min(PrComps(:,1)))+abs(max(PrComps(:,1))),...
        -0; -0, abs(min(PrComps(:,2)))+abs(max(PrComps(:,2)))]);
    mahalDis = mahal(clusterCat,PrComps(2:end,:));
    % scatter(PrComps(:,1),PrComps(:,2),50,mahalDis(:,1),'.')
    % hb = colorbar;
    % ylabel(hb,'Mahalanobis Distance to Component 1')
    
    similIdx=find(mahalDis<max(mahalDis)/10);
    % figure
    % scatter(PrComps(similIdx,1), PrComps(similIdx,2), 'k.');
    
    subsetWF=allWf(similIdx,:);
    [cc, ccv]=deal(nan(size(subsetWF,1),1));
    for wf=1:size(subsetWF,1)
        cc(wf)=max(abs(xcorr((selectedWF),(subsetWF(wf,:)),2,'coeff'))); %'unbiased' 'coeff'
    end
    
    for wf=1:size(subsetWF,1)
        ccv(wf)=max(abs(xcov((selectedWF),(subsetWF(wf,:)),2,'biased'))); %'unbiased' 'coeff'
    end
    
    %     clustev = evalclusters([cc ccv], 'kmeans', 'silhouette', 'KList',
    %     2:3); this takes too long, do not use clustev.OptimalK
    [clusIdx ,Centroids]= kmeans([cc ccv],2,'Distance','cityblock',...
        'Replicates',5);
    
    % figure; hold on
    % plot(mean(subsetWF(clusIdx==2,:)))
    
    similWF=similIdx(clusIdx==find(Centroids(:,2)==max(Centroids(:,2))));
    % figure; plot(mean(allWf(similWF,:)))
    
    handles.Spikes.HandSort.Units{channelNum}(similWF)=...
        unique(newClasses(lineSelecIdx));
end

unitsID=unique(handles.Spikes.HandSort.Units{channelNum});
%Check if unit selection still works
unitSelection=get(handles.SelectUnit_LB,'Value');
if get(handles.ShowAllUnits_RB,'value')
    unitSelection=find(unitsID>0);
else
    unitSelection=unitSelection(ismember(unitSelection,(1:numel(unitsID))));
end
unitSelection(unitsID(unitSelection)==0)=unitSelection(unitsID(unitSelection)==0)+1;
set(handles.SelectUnit_LB,'Value',unitSelection);
set(handles.SelectUnit_LB,'String',num2str(unitsID(unitsID>=0)))
if find(clusterClasses>0,1)
    if get(handles.ShowWF_CB,'value')
        handles=Plot_Sorted_WF(handles);
    else
        cla(handles.SortedUnits_Axes);
    end
end
if find(clusterClasses==0,1)
    handles=Plot_Unsorted_WF(handles);
end
Plot_Mean_WF(handles);
Plot_ISI(handles);
Plot_ACG(handles);
Plot_XCG(handles);
%  Update handles structure
guidata(hObject, handles);

%% --- Executes on selection change in SelectElectrode_LB.
function SelectElectrode_LB_Callback(hObject, ~, handles)
if strcmp(get(gcf,'SelectionType'),'normal')
    handles=LoadSpikes(handles);
    % plot "raw" (filtered) trace
    DisplayTraces(handles);
    % plot spike rasters
    DisplayMarkers(handles);
    % Update handles structure
    guidata(hObject, handles);
end

%% Get channel and unit selection
%     channelMenu=get(handles.SelectElectrode_LB,'string');
%     channelSelected=get(handles.SelectElectrode_LB,'value');
%     channelSelected=channelMenu(channelSelected);
%
%     unitMenu=get(handles.SelectUnit_LB,'string');
%     unitsSelected=get(handles.SelectUnit_LB,'value');
%     unitsSelected=unitMenu(unitsSelected);

%% --- Executes on selection change in SelectUnit_LB.
function SelectUnit_LB_Callback(hObject, ~, handles)
set(handles.ShowAllUnits_RB,'value',0)
if strcmp(get(gcf,'SelectionType'),'normal')
    [~,~,~,~,~,~,~,channelIdx,bestChannel]=GetUnitData(handles);
    channelList=find(channelIdx);
    channelNum=find(bestChannel==channelList);
    set(handles.SelectElectrode_LB,'string',num2str(channelList),'value',channelNum);
    channelList=cellstr(get(handles.SelectElectrode_LB,'String'))';
    colorChannels=cell(length(channelList),1);
    for elNum=1:length(channelList)
        if handles.rec_info.differentShanks(str2double(channelList(elNum)))>0
            colorChannels(elNum)=cellfun(@(thatChannel) sprintf(['<HTML><BODY bgcolor="%s">'...
                '<FONT color="%s">%s</FONT></BODY></HTML>'],... %size="+1"
                'black','white', thatChannel),channelList(elNum),'UniformOutput',false);
        else
            colorChannels(elNum)=cellfun(@(thatChannel) sprintf(['<HTML><BODY bgcolor="%s">'...
                '<FONT color="%s">%s</FONT></BODY></HTML>'],... %size="+1"
                'white','black', thatChannel),channelList(elNum),'UniformOutput',false);
        end
    end
    set(handles.SelectElectrode_LB, 'String', colorChannels);
    set(handles.SelectElectrode_LB ,'ListboxTop',max([1 numel(channelNum)-3]));
    if get(handles.ShowWF_CB,'value')
        handles=Plot_Sorted_WF(handles);
    else
        cla(handles.SortedUnits_Axes);
    end
    Plot_Mean_WF(handles);
    Plot_ISI(handles);
    Plot_ACG(handles);
    Plot_XCG(handles);
    DisplayTraces(handles);
    DisplayMarkers(handles);
    
    % elseif strcmp(get(gcf,'SelectionType'),'open') % double click
    % else
    %  Update handles structure
    guidata(hObject, handles);
end

%% --- Executes on button press in ShowAllUnits_RB.
function ShowAllUnits_RB_Callback(hObject, ~, handles)
if get(hObject,'Value')
    if get(handles.ShowWF_CB,'value')
        handles=Plot_Sorted_WF(handles);
    else
        cla(handles.SortedUnits_Axes);
    end
    Plot_Mean_WF(handles);
    Plot_ISI(handles);
    Plot_ACG(handles);
    Plot_XCG(handles);
    guidata(hObject, handles);
end

%% --- Executes on button press in PreviewTh_PB.
function PreviewTh_PB_Callback(hObject, ~, handles)
if ~isfield(handles,'datFile')
    handles.datFile=[cell2mat(regexp(handles.spikeFile,'.+(?=\_\w+\_\w+\.)','match')) '_nopp'];
end
handles.PreviewTh_PB.ForegroundColor=[0.2 0.5 0.7];
RunSpykingCircus(cd,handles.datFile,{'previewspkc'});
handles.PreviewTh_PB.ForegroundColor=[0.3490 0.2000 0.3294];

%% --- Executes on button press in LauncherGUI_PB.
function LauncherGUI_PB_Callback(hObject, eventdata, handles)
if ~isfield(handles,'datFile')
    handles.datFile=[cell2mat(regexp(handles.spikeFile,'.+(?=\_\w+\_\w+\.)','match')) '_nopp'];
end
RunSpykingCircus(cd,handles.datFile,{'launcherGUI'});

%% --- Executes on button press in RefineSort_PB.
function RefineSort_PB_Callback(hObject, ~, handles)
if ~isfield(handles,'datFile')
    handles.datFile=[cell2mat(regexp(handles.spikeFile,'.+(?=\_\w+\_\w+\.)','match')) '_nopp'];
end
handles.RefineSort_PB.ForegroundColor=[0.2 0.5 0.7];
[status,cmdout]=RunSpykingCircus(cd,handles.datFile,{'startVisGUI'});
handles.RefineSort_PB.ForegroundColor=[0.3490 0.2000 0.3294];

%% --- Outputs from this function are returned to the command line.
function varargout = SpikeVisualizationGUI_OutputFcn(hObject, ~, handles)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in PeakReAlign_PB.
function PeakReAlign_PB_Callback(hObject, eventdata, handles)
% confirm re-alignment
confirm = questdlg('Realign Spikes Waveforms?', ...
    'Spike Alignment','Yes','No thank you','Yes');
switch confirm
    case 'Yes'
        % first load raw traces
        %         tb_filename=regexp(handles.spikeFile,'.+(?=_\w+.\w+$)','match');
        %         load([tb_filename{:} '_raw.mat']);
        selectedEl=get(handles.SelectElectrode_LB,'value');
        rawData=handles.traces.(handles.tracesInfo.name)(selectedEl,:);
        %then re-align (current waveforms, raw data, spike times, unit IDs)
        Spikes=ReAlignSpikes(handles.Spikes.HandSort.Waveforms{selectedEl,1},...
            rawData,...
            handles.Spikes.HandSort.SpikeTimes{selectedEl,1},...
            handles.Spikes.HandSort.Units{selectedEl,1});
        handles.Spikes.HandSort.Waveforms(selectedEl)=Spikes.Waveforms;
        handles.Spikes.HandSort.Waveforms(selectedEl)=...
            cellfun(@(wfs) int16(wfs), handles.Spikes.HandSort.Waveforms{selectedEl},'UniformOutput',false); %256 units per channel max
        handles.Spikes.HandSort.SpikeTimes(selectedEl,1)=Spikes.SpikeTimes;
        % downsample spike times to 1 millisecond bins
        for ElNum=1:length(selectedEl)
            if ~isempty(handles.Spikes.HandSort.SpikeTimes{ElNum,1})
                handles.Spikes.HandSort.SpikeTimes{ElNum,2}=...
                    uint32(handles.Spikes.HandSort.SpikeTimes{ElNum,1})/...
                    (handles.Spikes.HandSort.samplingRate(ElNum,1)/...
                    handles.Spikes.HandSort.samplingRate(ElNum,2));
            else
                handles.Spikes.HandSort.SpikeTimes{ElNum,2}=[];
            end
        end
    case 'No thank you'
        %OK
end
handles=Plot_Unsorted_WF(handles);
if get(handles.ShowWF_CB,'value')
    handles=Plot_Sorted_WF(handles);
else
    cla(handles.SortedUnits_Axes);
end
Plot_Mean_WF(handles);
% Update handles structure
guidata(hObject, handles);

%% --- Executes during object creation, after setting all properties.
function SelectElectrode_LB_CreateFcn(hObject, ~, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% --- Executes during object creation, after setting all properties.
function SelectUnit_LB_CreateFcn(hObject, ~, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% --- Executes during object creation, after setting all properties.
function TW_slider_CreateFcn(hObject, ~, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.839 .91 .851]);
end

%% --- Executes on button press in LoadFile_PB.
function LoadFile_PB_Callback(hObject, ~, handles)
if isfield(handles,'spikeFile')
    handles = rmfield(handles,'spikeFile');
end
try
    handles.userinfo=UserDirInfo;
catch
    handles.userinfo=[];
    handles.userinfo.user=getenv('username');
end
%% get most recently changed data folder (looks for a folder named "export")
if isfield(handles.userinfo,'directory')
    exportDir=regexprep(handles.userinfo.directory,'\\\w+$','\\export');
else
    [exportDir,handles.userinfo.directory]=deal(cd);
end
dataDirListing=dir(exportDir);
%removing dots
dataDirListing=dataDirListing(cellfun('isempty',cellfun(@(x) strfind(x,'.'),...
    {dataDirListing.name},'UniformOutput',false)));
%removing other folders
dataDirListing=dataDirListing(cellfun('isempty',cellfun(@(x)...
    regexp('list | all | unwanted | folders | here ',x),...
    {dataDirListing.name},'UniformOutput',false)));
[~,fDateIdx]=sort([dataDirListing.datenum],'descend');
recentDataFolder=[exportDir filesep dataDirListing(fDateIdx(1)).name filesep];

% open user input window
[handles.spikeFile,handles.exportDir] = uigetfile({'*.mat;*.hdf5;*.csv','Export Formats';...
    '*.dat','Raw data';'*.*','All Files' },'Most recent data',recentDataFolder);
if handles.spikeFile==0
    handles.spikeFile='';
    handles.exportDir=recentDataFolder;
else
    handles.fileLoaded=0;
    handles=LoadSpikes(handles);
end

% Update handles structure
guidata(hObject, handles);

%% --- Executes on button press in Reload_PB.
function Reload_PB_Callback(hObject, ~, handles)

%% --- Executes on button press in Save_PB.
function Save_PB_Callback(hObject, ~, handles)
if ~isfield(handles,'exportDir')
     handles.exportDir=handles.datDir;
end
if isfield(handles,'fname')
    outputName=[handles.exportDir handles.fname '_spikesResorted'];
elseif isfield(handles,'spikeFile') && ~isempty(handles.spikeFile) 
    outputName=[handles.exportDir cell2mat(regexp(handles.spikeFile,'.+(?=_spikes)','match'))...
            '_spikesResorted'];
elseif isfield(handles,'offlineSort_SpikeFile') && ~isempty(handles.offlineSort_SpikeFile) 
    exportDirListing=dir(handles.exportDir);
    infoFile=exportDirListing(~cellfun('isempty',cellfun(@(x) strfind(x,'_info.'),...
        {exportDirListing.name},'UniformOutput',false))).name;
    outputName=fullfile(handles.datDir,[cell2mat(regexp(infoFile,'.+(?=_info)','match'))...
        '_spikesResorted']);
elseif isfield(handles,'datFile')
    outputName=[handles.exportDir handles.datFile(1:end-4)  ...
        '_spikesResorted.mat'];
end

% find existing fields 
allFields=fieldnames(handles);
exportFields=allFields(cellfun(@(fldName) ...
        strcmp(fldName,'spikeFile') || strcmp(fldName,'exportDir') || ...
        strcmp(fldName,'datFile') || strcmp(fldName,'datDir') || ...
        strcmp(fldName,'classification') || strcmp(fldName,'Spikes') || ...
        strcmp(fldName,'rec_info') || strcmp(fldName,'subset') || ...
        strcmp(fldName,'userinfo') || strcmp(fldName,'rawData') || ...
        strcmp(fldName,'rawDataInfo') || strcmp(fldName,'fname'),...
        allFields,'UniformOutput',true));
exportFields=cellfun(@(x) [x '*'], exportFields,'UniformOutput',false); %placeholders for coma separation
exportFields=strrep([exportFields{:}],'*',[char(39) ',' char(39)]); %insert coma
% save file
eval(['save(outputName,''-struct'',''handles'',' [char(39) exportFields(1:end-2)] ',''-v7.3'')'])
        

% --- Executes on button press in ShowWF_CB.
function ShowWF_CB_Callback(hObject, eventdata, handles)

function ShowUWF_CB_Callback(hObject, eventdata, handles)

function ShowHowManyUWF_ET_Callback(hObject, eventdata, handles)

function ShowHowManyUWF_ET_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in ExportData_PB.
function ExportData_PB_Callback(hObject, eventdata, handles)
channelNum=get(handles.SelectElectrode_LB,'value');
channelList=ReturnElectrodes(handles.SelectElectrode_LB);
channelLabel=channelList(channelNum);

[selectedUnits,~,unitsIdx,spikeTimes,waveForms,samplingRate,~,channelIdx]=GetUnitData(handles);

allTraces=handles.traces; %allTraces.numChan=handles.tracesInfo.numChan;
traceInfo=handles.tracesInfo;
classification=handles.classification;
% get 'raw data' excerpt
traceExcerpt.data=DisplayTraces(handles);
traceExcerpt.xTicks=linspace(0,handles.rec_info.samplingRate*2,4);
traceExcerpt.xTicklabels=linspace(round(handles.tracesInfo.excerptLocation-...
    handles.tracesInfo.excerptSize)/handles.rec_info.samplingRate,...
    round(handles.tracesInfo.excerptLocation+handles.tracesInfo.excerptSize)...
    /handles.rec_info.samplingRate,4);
traceExcerpt.location=handles.tracesInfo.excerptLocation;
traceExcerpt.excerptSize=handles.tracesInfo.excerptSize;
traceExcerpt.spkTimes=DisplayMarkers(handles);

cd(handles.exportDir);
if isfield(handles,'Trials')
    TTLs=handles.Trials;
else
    TTLs=[]
end
channelList=channelList(ismember(channelList,find(channelIdx)));
if numel(selectedUnits)>1
    fName=[handles.datFile(1:end-4) '_SelectedData_Multi.mat'];
else
    fName=[handles.datFile(1:end-4) '_SelectedData_Ch' num2str(channelLabel) 'U' num2str(selectedUnits) '.mat'];
end
save(fName,'waveForms','spikeTimes','unitsIdx','samplingRate','channelList',...
    'selectedUnits','TTLs','traceExcerpt','allTraces','traceInfo','classification');

% --- Executes on button press in LoadTTL_PB.
function LoadTTL_PB_Callback(hObject, eventdata, handles)
[handles.TTLFile,handles.TTLFileDir] = uigetfile({'*.*','All Files' },...
    'Select file containing TTL data',handles.exportDir);% '*.mat;*.hdf5','Export Formats';'*.dat;*.nev','Raw data';
handles.Trials = LoadTTL([handles.TTLFileDir handles.TTLFile]);
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in ShowTTL_PB.
function ShowTTL_PB_Callback(hObject, eventdata, handles)
% hObject    handle to ShowTTL_PB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in ExportFigs_PB.
function ExportFigs_PB_Callback(hObject, eventdata, handles)
% hObject    handle to ExportFigs_PB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in ClassifySU_PB.
function ClassifySU_PB_Callback(hObject, eventdata, handles)
channelNum=get(handles.SelectElectrode_LB,'value');
channelList=ReturnElectrodes(handles.SelectElectrode_LB);
channelLabel=channelList(channelNum);

if get(handles.ShowAllUnits_RB,'value')
    helpdlg('All Clusters are selected. Please select a specific Cluster','SU classification');
    return;
else
    unitID=ReturnUnits(handles.SelectUnit_LB);
    if unitID==0
        return;
    end
    selectedUnitsListIdx=get(handles.SelectUnit_LB,'value');
    if isempty(selectedUnitsListIdx) || selectedUnitsListIdx(end)>length(unitID)
        selectedUnitsListIdx=length(unitID);
    end
    selectedUnits=unitID(selectedUnitsListIdx);
end

for unit=1:length(selectedUnits)
    %find previous records
    sortID=handles.classification.SortID(handles.classification.Channel==channelLabel &...
        handles.classification.UnitNumber==selectedUnits(unit));
    if isempty(sortID)
        sortID=size(handles.classification,1)+1;
    end
    handles.classification(sortID,:)={sortID,channelLabel,selectedUnits(unit),'SU',''};
end
handles=ClassificationColor(handles);
guidata(hObject, handles);

% --- Executes on button press in ClassifyMU_PB.
function ClassifyMU_PB_Callback(hObject, eventdata, handles)
channelNum=get(handles.SelectElectrode_LB,'value');
channelList=ReturnElectrodes(handles.SelectElectrode_LB);
channelLabel=channelList(channelNum);

if get(handles.ShowAllUnits_RB,'value')
    helpdlg('All Clusters are selected. Please select a specific Cluster','MU classification');
    return;
else
    unitID=ReturnUnits(handles.SelectUnit_LB);
    if unitID==0
        return;
    end
    selectedUnitsListIdx=get(handles.SelectUnit_LB,'value');
    if isempty(selectedUnitsListIdx) || selectedUnitsListIdx(end)>length(unitID)
        selectedUnitsListIdx=length(unitID);
    end
    selectedUnits=unitID(selectedUnitsListIdx);
end

for unit=1:length(selectedUnits)
    %find previous records
    sortID=handles.classification.SortID(handles.classification.Channel==channelLabel &...
        handles.classification.UnitNumber==selectedUnits(unit));
    if isempty(sortID)
        sortID=size(handles.classification,1)+1;
    end
    handles.classification(sortID,:)={sortID,channelLabel,selectedUnits(unit),'MU',''};
end
handles=ClassificationColor(handles);
guidata(hObject, handles);

% --- Executes on button press in ClassifyOther_PB.
function ClassifyOther_PB_Callback(hObject, eventdata, handles)
channelNum=get(handles.SelectElectrode_LB,'value');
channelList=ReturnElectrodes(handles.SelectElectrode_LB);
channelLabel=channelList(channelNum);

if get(handles.ShowAllUnits_RB,'value')
    helpdlg('All Clusters are selected. Please select a specific Cluster','Other classification');
    return;
else
    unitID=ReturnUnits(handles.SelectUnit_LB);
    if unitID==0
        return;
    end
    selectedUnitsListIdx=get(handles.SelectUnit_LB,'value');
    if isempty(selectedUnitsListIdx) || selectedUnitsListIdx(end)>length(unitID)
        selectedUnitsListIdx=length(unitID);
    end
    selectedUnits=unitID(selectedUnitsListIdx);
end

for unit=1:length(selectedUnits)
    %find previous records
    sortID=handles.classification.SortID(handles.classification.Channel==channelLabel &...
        handles.classification.UnitNumber==selectedUnits(unit));
    if isempty(sortID)
        sortID=size(handles.classification,1)+1;
    end
    handles.classification(sortID,:)={sortID,channelLabel,selectedUnits(unit),'Other',''};
end
handles=ClassificationColor(handles);
guidata(hObject, handles);

% --- Executes on button press in Note_PB.
function Note_PB_Callback(hObject, eventdata, handles)

prompt = {'Enter comment about special classification'};
dlg_title = '';
num_lines = 1;
defaultans = {''};
comment = inputdlg(prompt,dlg_title,num_lines,defaultans);

channelNum=get(handles.SelectElectrode_LB,'value');
channelList=ReturnElectrodes(handles.SelectElectrode_LB);
channelLabel=channelList(channelNum);

if get(handles.ShowAllUnits_RB,'value')
    helpdlg('All Clusters are selected. Please select a specific Cluster','Notes');
    return;
else
    unitID=ReturnUnits(handles.SelectUnit_LB);
    if unitID==0
        return;
    end
    selectedUnitsListIdx=get(handles.SelectUnit_LB,'value');
    if isempty(selectedUnitsListIdx) || selectedUnitsListIdx(end)>length(unitID)
        selectedUnitsListIdx=length(unitID);
    end
    selectedUnits=unitID(selectedUnitsListIdx);
end

for unit=1:length(selectedUnits)
    %find previous records
    sortID=handles.classification.SortID(handles.classification.Channel==channelLabel &...
        handles.classification.UnitNumber==selectedUnits(unit));
    if isempty(sortID)
        sortID=size(handles.classification,1)+1;
        currentClass=categorical();
    else
        currentClass=handles.classification(sortID,:).Classification;
    end
    handles.classification(sortID,:)={sortID,channelLabel,selectedUnits(unit),currentClass,comment};
end
guidata(hObject, handles);

% --- Executes on button press in PB_DeleteClass.
function PB_DeleteClass_Callback(hObject, eventdata, handles)
channelNum=get(handles.SelectElectrode_LB,'value');
channelList=ReturnElectrodes(handles.SelectElectrode_LB);
channelLabel=channelList(channelNum);

if get(handles.ShowAllUnits_RB,'value')
    helpdlg('All Clusters are selected. Please select a specific Cluster','Other classification');
    return;
else
    unitID=ReturnUnits(handles.SelectUnit_LB);
    if unitID==0
        return;
    end
    selectedUnitsListIdx=get(handles.SelectUnit_LB,'value');
    if isempty(selectedUnitsListIdx) || selectedUnitsListIdx(end)>length(unitID)
        selectedUnitsListIdx=length(unitID);
    end
    selectedUnits=unitID(selectedUnitsListIdx);
end

for unit=1:length(selectedUnits)
    %find previous records
    sortID=handles.classification.SortID(handles.classification.Channel==channelLabel &...
        handles.classification.UnitNumber==selectedUnits(unit));
    if isempty(sortID)
        sortID=size(handles.classification,1)+1;
    end
    handles.classification(sortID,:)={sortID,channelLabel,selectedUnits(unit),'',''};
end
handles=ClassificationColor(handles);
guidata(hObject, handles);

function handles=ClassificationColor(handles)
% channelNum=get(handles.SelectElectrode_LB,'value');
clusterList=ReturnUnits(handles.SelectUnit_LB);
colorClusters=cell(length(clusterList),1);
for clusNum=1:length(clusterList)
    if numel(unique(clusterList))==numel(clusterList) %global cluster IDs
        clusterNum=find(handles.classification.UnitNumber==clusterList(clusNum));
    else %each electrode has it's own cluster IDs
        clusterNum=find(handles.classification.Channel==channelNum & ...
            handles.classification.UnitNumber==clusterList(clusNum));
    end
    if isempty(clusterNum)
        colorClusters(clusNum)=cellfun(@(thatCluster) sprintf(['<HTML><BODY bgcolor="%s">'...
            '<FONT color="%s">%s</FONT></BODY></HTML>'],... %size="+1"
            'white','black', num2str(thatCluster)),{clusterList(clusNum)},'UniformOutput',false);
        continue;
    end
    if contains(char(handles.classification.Classification(clusterNum)),'SU')
        colorClusters(clusNum)=cellfun(@(thatCluster) sprintf(['<HTML><BODY bgcolor="%s">'...
            '<FONT color="%s">%s</FONT></BODY></HTML>'],... %size="+1"
            'white','red',  num2str(thatCluster)),{clusterList(clusNum)},'UniformOutput',false);
    elseif contains(char(handles.classification.Classification(clusterNum)),'MU')
        colorClusters(clusNum)=cellfun(@(thatCluster) sprintf(['<HTML><BODY bgcolor="%s">'...
            '<FONT color="%s">%s</FONT></BODY></HTML>'],... %size="+1"
            'white','blue',  num2str(thatCluster)),{clusterList(clusNum)},'UniformOutput',false);
    elseif contains(char(handles.classification.Classification(clusterNum)),'Other')
        colorClusters(clusNum)=cellfun(@(thatCluster) sprintf(['<HTML><BODY bgcolor="%s">'...
            '<FONT color="%s">%s</FONT></BODY></HTML>'],... %size="+1"
            'white','green',  num2str(thatCluster)),{clusterList(clusNum)},'UniformOutput',false);
    else
        colorClusters(clusNum)=cellfun(@(thatCluster) sprintf(['<HTML><BODY bgcolor="%s">'...
            '<FONT color="%s">%s</FONT></BODY></HTML>'],... %size="+1"
            'white','black', thatCluster),{clusterList(clusNum)},'UniformOutput',false);
    end
end
set(handles.SelectUnit_LB, 'String', colorClusters);


% --- Executes on button press in PB_SaveClass.
function PB_SaveClass_Callback(hObject, eventdata, handles)

if size(handles.classification,1)>0
    handles.classification = sortrows(handles.classification,2);
    if exist(fullfile(handles.datDir,[regexp(handles.datFile,'.+(?=\.)','match','once')...
            '_classification.xlsx']),'file')
        overwiteFile = questdlg('Classification file exits. Overwrite?', ...
            '','Yes','No','Yes');
        switch overwiteFile
            case 'Yes'
                delete([handles.datDir filesep handles.fname  '_classification.xlsx'])
            case 'No'
                disp('Classification file not saved')
                return
        end
    end
    writetable(handles.classification,fullfile(handles.datDir,...
        [regexp(handles.datFile,'.+(?=\.)','match','once')...
            '_classification.xlsx']));
end

% --------------------------------------------------------------------
function Help_Menu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function AddOptions_Menu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function DisplayNtrodeMarkers_MenuItem_Callback(hObject, eventdata, handles)
if contains(get(handles.DisplayNtrodeMarkers_MenuItem,'Checked'),'off')
    if isfield(handles.rec_info,'differentShanks')
        set(handles.DisplayNtrodeMarkers_MenuItem,'Checked','on')
    else
        disp('no information available about shanks');
        return;
    end
else
    set(handles.DisplayNtrodeMarkers_MenuItem,'Checked','off')
end

% --------------------------------------------------------------------
function SpikeSource_MenuItem_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function SortVersion_MenuItem_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function Spikes_OriginalVersion_Menu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function Spikes_CurrentVersion_Menu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function Spikes_SortOn_Menu_Callback(hObject, eventdata, handles)
set(handles.Spikes_SortOn_Menu,'Checked','on');
set(handles.Spikes_SortOff_Menu,'Checked','off');
set(handles.Spikes_Th_Menu,'Checked','off');
% [handles.spikeFile,handles.exportDir] = uigetfile({'*.mat;*.hdf5','All Data Formats';...
%     '*.*','All Files' },'Export folder',handles.exportDir);
handles.fileLoaded=0;
handles=LoadSpikes(handles);
% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function Spikes_SortOff_Menu_Callback(hObject, eventdata, handles)
set(handles.Spikes_SortOff_Menu,'Checked','on');
set(handles.Spikes_SortOn_Menu,'Checked','off');
set(handles.Spikes_Th_Menu,'Checked','off');

[handles.spikeFile,handles.exportDir] = uigetfile({'*.mat;*.hdf5','All Data Formats';...
    '*.*','All Files' },'Export folder',handles.exportDir);
handles.fileLoaded=0;
handles=LoadSpikes(handles);
% plot "raw" (filtered) trace
DisplayTraces(handles);
% plot spike rasters
DisplayMarkers(handles);
% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function Spikes_Th_Menu_Callback(hObject, eventdata, handles)
set(handles.Spikes_Th_Menu,'Checked','on');
set(handles.Spikes_SortOn_Menu,'Checked','off');
set(handles.Spikes_SortOff_Menu,'Checked','off');
if isempty(handles.spikeFile) || ~strcmp(handles.spikeFile(end-2:end),'mat')
    [handles.spikeFile,handles.exportDir] = uigetfile({'*.mat;*.hdf5','All Data Formats';...
        '*.*','All Files' },'Export folder',handles.exportDir);
end
handles.fileLoaded=0;
handles=LoadSpikes(handles);
% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function Spikes_PrevSorted_Menu_Callback(hObject, eventdata, handles)
set(handles.Spikes_SortOn_Menu,'Checked','off');
set(handles.Spikes_SortOff_Menu,'Checked','off');
set(handles.Spikes_Th_Menu,'Checked','off');
set(handles.Spikes_PrevSorted_Menu,'Checked','on');
if ~strcmp(handles.spikeFile(end-2:end),'mat')
    [handles.spikeFile,handles.exportDir] = uigetfile({'*.mat;*.hdf5','All Data Formats';...
        '*.*','All Files' },'Export folder',handles.exportDir);
end
handles.fileLoaded=0;
handles=LoadSpikes(handles);
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in PCA_PB.
function PCA_PB_Callback(hObject, eventdata, handles)
channelNum=get(handles.SelectElectrode_LB,'value');

if contains(get(handles.DisplayNtrodeMarkers_MenuItem,'Checked'),'on')
    % get channel values from same Ntrode
    shankNum=cumsum(logical([0 diff(handles.rec_info.differentShanks(1:handles.rec_info.numRecChan))]));
    channelNum=find(shankNum==shankNum(channelNum));
end

[spikeTimes,sortIdx]=sort(vertcat(handles.Spikes.HandSort.SpikeTimes{channelNum,1}));
clusterID={handles.Spikes.HandSort.Units{channelNum,1}};
clusterIncrement=cumsum(cellfun(@(x) length(unique(x)),clusterID)); clusterIncrement=clusterIncrement-min(clusterIncrement);
clusterID=cellfun(@(x,y) x+y, clusterID,num2cell(clusterIncrement,4),'UniformOutput',false);
clusterID=vertcat(clusterID{:}); clusterID=clusterID(sortIdx);

% waveForms=NaN(length(spikeTimes),length(channelNum),size(handles.Spikes.HandSort.Waveforms{1},1));
% for elNum=1:length(channelNum)
%     if isa(handles.traces,'memmapfile') % reading electrode data from .dat file
%         waveForms(:,elNum,:)=ExtractChunks(handles.traces.Data(channelNum(elNum):...
%             handles.rec_info.numRecChan:max(size(handles.traces.Data))),...
%             spikeTimes,50,'tshifted'); %'tzero' 'tmiddle' 'tshifted'
%     else
%         waveForms(:,elNum,:)=ExtractChunks(handles.traces(channelNum(elNum),:),...
%             spikeTimes,50,'tshifted'); %'tzero' 'tmiddle' 'tshifted'
%     end
% end

features = Create_FeatureSpace(waveForms);
figure;
plot3(features(:,1),features(:,2),features(:,8),'.')

[CluSep, m] = Cluster_Quality(Fet, clusterID);


[SpikeData,Y]=deal(cell(1,length(channelNum)));
for elNum=1:length(channelNum)
    SpikeData{elNum}=double(handles.Spikes.HandSort.Waveforms{channelNum(elNum),:});
    SpikeNum = size(SpikeData{elNum},2);
    SpikeData{elNum} = SpikeData{elNum} - repmat(mean(SpikeData{elNum},2),1,SpikeNum);
    [Y{elNum},Sig,~] = svd(SpikeData{elNum});
end
% sig = diag(Sig);
% figure; semilogy(sig(sig>1),'kx-') % plot the significant singular values
% xlabel('index','fontsize',14); ylabel('singular value','fontsize',14)
figure;
for elNum=1:length(channelNum)
    subplot(2,2,elNum)
plot3(SpikeData{elNum}'*Y{elNum}(:,1),SpikeData{elNum}'*Y{elNum}(:,2),SpikeData{elNum}'*Y{elNum}(:,3),'.')
end

function electrodes=ReturnElectrodes(unitsHandle)
elecString=get(unitsHandle,'string');
try
    electrodes=str2num(elecString);
catch
    electrodes=cellfun(@(x) str2double(regexp(x,'(?<=>)\s*\d+(?=</FONT>)','match')),elecString,'UniformOutput',true);
end

function units=ReturnUnits(unitsHandle)
unitString=get(unitsHandle,'string');
try
    units=str2num(unitString);
catch
    units=cellfun(@(x) str2double(regexp(x,'(?<=>)\s*\d+(?=</FONT>)','match')),unitString,'UniformOutput',true);
end

function [selectedUnits,selectedUnitsListIdx,unitsIDs,spikeTimes,waveForms,samplingRate,unitsIdx,channelIdx,bestChannel]=GetUnitData(handles)
% find which units are selected 
if get(handles.ShowAllUnits_RB,'value') %all units
    unitID=ReturnUnits(handles.SelectUnit_LB);
    selectedUnitsListIdx=find(unitID>0);
    selectedUnits=unitID(selectedUnitsListIdx);
else %or selected units
    unitID=ReturnUnits(handles.SelectUnit_LB);
    if unitID==0
        return;
    end
    selectedUnitsListIdx=get(handles.SelectUnit_LB,'value');
    if isempty(selectedUnitsListIdx) & ~isempty(unitID)
        set(handles.SelectUnit_LB,'value',1);
        selectedUnitsListIdx=1;
    end
%     if sum(~ismember(unique(unitsIdx(unitsIdx>=0)),unitID))>0 %Not sure when electrode's unit would not be part of listed units
%         unitID=unique(unitsIdx);
%     end
    if isempty(selectedUnitsListIdx) || selectedUnitsListIdx(end)>length(unitID)
        selectedUnitsListIdx=length(unitID);
    end
    selectedUnits=unitID(selectedUnitsListIdx);
    selectedUnitsListIdx=selectedUnitsListIdx(selectedUnits>0);
    selectedUnits=selectedUnits(selectedUnits>0);
end
% find all occurences of these units across channels
channelIdx=cellfun(@(x) logical(sum(ismember(x,selectedUnits))), handles.Spikes.HandSort.Units);
occurencePerChannel=cellfun(@(x) sum(ismember(x,selectedUnits)), handles.Spikes.HandSort.Units);
bestChannel=find(occurencePerChannel==max(occurencePerChannel),1);
allUnits=vertcat(handles.Spikes.HandSort.Units{:});
allSpikeTimes=vertcat(handles.Spikes.HandSort.SpikeTimes{:});
allWaveforms=horzcat(handles.Spikes.HandSort.Waveforms{:})';
unitsIdx=ismember(allUnits,selectedUnits);
unitsIDs=allUnits(unitsIdx,:);
spikeTimes=allSpikeTimes(unitsIdx,:);
waveForms=allWaveforms(unitsIdx,:)';
samplingRate=handles.Spikes.HandSort.samplingRate(bestChannel,1);
