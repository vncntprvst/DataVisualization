function varargout = SpikeVisualizationGUI(varargin)
% MATLAB code for SpikeVisualizationGUI.fig


% Last Modified by GUIDE v2.5 09-Jun-2016 17:06:38

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
    % open user input window
end

handles=LoadSpikes(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SpikeVisualizationGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

%% Load data function
function handles=LoadSpikes(handles)
% function declaration
axis_name= @(x) sprintf('Chan %.0f',x);
if strcmp(handles.fname,'')
    set(handles.FileName,'string','')
else
    cd(handles.exportdir);
    userinfo=UserDirInfo;
    exportDirListing=dir;
    handles.spikeFile={exportDirListing(~cellfun('isempty',cellfun(@(x) strfind(x,'_spikes'),...
        {exportDirListing.name},'UniformOutput',false))).name};
    if size(handles.spikeFile,2)>1
        nameComp=cellfun(@(x) sum(ismember(x,handles.fname)) ,handles.spikeFile);
        if abs(diff(nameComp))<2 %that can be tricky for some files
            %select the most recent
            fileDates=datenum({exportDirListing(~cellfun('isempty',cellfun(@(x) strfind(x,'_spikes'),...
                {exportDirListing.name},'UniformOutput',false))).date});
            handles.spikeFile=handles.spikeFile{fileDates==max(fileDates)};
        else
            handles.spikeFile=handles.spikeFile{nameComp==max(nameComp)};
        end
    else
        handles.spikeFile=handles.spikeFile{:};
    end
    set(handles.FileName,'string',[handles.exportdir userinfo.slash handles.spikeFile])
    
    %% Load spike data
    spikeData=load(handles.spikeFile);
    handles=catstruct(handles,spikeData);
    clear spikeData;
    
    %% Set number of electrodes and units, select electrode with most units and spikes
    set(handles.SelectChannel_LB,'string',num2str(handles.Spikes.Offline_Threshold.electrode'));
    if isfield(handles.Spikes,'Online_Sorting') || isfield(handles.Spikes,'Offline_Sorting')
        if isfield(handles.Spikes,'Offline_Sorting') %preferred
            set(handles.Spikes_SortOff_RB,'value',1);
            set(handles.Spikes_SortOn_RB,'value',0);
            handles.Units=handles.Spikes.Offline_Sorting.Units;
            handles.SpikeTimes=handles.Spikes.Offline_Sorting.SpikeTimes;
            handles.Waveforms=handles.Spikes.Offline_Sorting.Waveforms;
        else
            set(handles.Spikes_SortOff_RB,'value',0);
            set(handles.Spikes_SortOn_RB,'value',1);
            handles.Spikes.inGUI.Units=handles.Spikes.Online_Sorting.Units;
            handles.Spikes.inGUI.SpikeTimes=handles.Spikes.Online_Sorting.SpikeTimes;
            handles.Spikes.inGUI.Waveforms=handles.Spikes.Online_Sorting.Waveforms;
        end
        numUnits=cellfun(@(x) sum(length(x)*unique(x)), handles.Spikes.inGUI.Units);
        electrodeNum=find(numUnits==max(numUnits),1);
        set(handles.SelectChannel_LB,'value',electrodeNum);
    else
        set(handles.SelectChannel_LB,'value',1)
    end
    
    %% initialize variables
    unitsIdx=handles.Spikes.inGUI.Units{electrodeNum};
    waveForms=handles.Spikes.inGUI.Waveforms{electrodeNum};
    spikeTimes=handles.Spikes.inGUI.SpikeTimes{electrodeNum};
    samplingRate=handles.Spikes.Online_Sorting.samplingRate(electrodeNum);
    % how many units on that electrode?
    unitsID=unique(unitsIdx); %number of clustered units
    set(handles.SelectUnit_LB,'string',num2str(unitsID'));
    set(handles.SelectUnit_LB,'value',find(unitsID~=0));
    
    % here's how this can work:
    %     Spikes detected from threshold are the benchmark (unit code = 0).
    %     May want to extract waveforms, but most important is the time.
    %     Sorted units (whatever the source) "color" those units. One spike per ms max.
    
    %% Plot unsorted spikes
    if sum(unitsIdx==0)>2000 %then only plot subset of waveforms
        subset=find(unitsIdx==0);
        handles.subset{1}=subset(1:round(sum(unitsIdx==0)/2000):end);
    else
        handles.subset{1}=find(unitsIdx==0);
    end
    axes(handles.UnsortedUnits_Axes); hold on;colormap lines; cmap=colormap;
    cla(handles.UnsortedUnits_Axes);
    set(handles.UnsortedUnits_Axes,'Visible','on');
    plot(waveForms(:,handles.subset{1}),'linewidth',2,'Color','k');
    
    set(gca,'xtick',linspace(0,size(waveForms(:,handles.subset{1}),1),5),...
        'xticklabel',round(linspace(-round(size(waveForms(:,handles.subset{1}),1)/2),...
        round(size(waveForms(:,handles.subset{1}),1)/2),5)/(double(samplingRate)/1000),2),'TickDir','out');
    % legend('Unclustered waveforms','location','southeast')
    axis('tight');box off;
    xlabel('Time (ms)')
    ylabel('Voltage (?V)')
    set(gca,'Color','white','FontSize',10,'FontName','calibri');
    
    %% Plot clusters
    % selected unit ids
    axes(handles.SortedUnits_Axes); hold on;colormap lines; cmap=colormap;
    cla(handles.SortedUnits_Axes);
    set(handles.SortedUnits_Axes,'Visible','on');
    selectedUnits=get(handles.SelectUnit_LB,'value');
    for unitP=1:length(selectedUnits)
        if sum(unitsIdx==selectedUnits(unitP))>2000 %then only plot subset of waveforms
            subset=find(unitsIdx==selectedUnits(unitP));
            handles.subset{selectedUnits(unitP)}=subset(1:round(sum(unitsIdx==selectedUnits(unitP))/2000):end);
        else
            handles.subset{selectedUnits(unitP)}=find(unitsIdx==selectedUnits(unitP));
        end
        
        plot(waveForms(:,handles.subset{selectedUnits(unitP)}),'linewidth',1.2,'Color',cmap(unitP,:));
        
        set(gca,'xtick',linspace(0,size(waveForms(:,handles.subset{selectedUnits(unitP)}),1),5),...
            'xticklabel',round(linspace(-round(size(waveForms(:,handles.subset{selectedUnits(unitP)}),1)/2),...
            round(size(waveForms(:,handles.subset{selectedUnits(unitP)}),1)/2),5)/(double(samplingRate)/1000),2),'TickDir','out');
        % legend('Unclustered waveforms','location','southeast')
        axis('tight');box off;
        xlabel('Time (ms)')
        ylabel('Voltage (?V)')
        set(gca,'Color','white','FontSize',10,'FontName','calibri');
    end
    
    % mean waveforms
    axes(handles.MeanSortedUnits_Axes); hold on;colormap lines; cmap=colormap;
    cla(handles.MeanSortedUnits_Axes);
    set(handles.MeanSortedUnits_Axes,'Visible','on');
    selectedUnits=get(handles.SelectUnit_LB,'value');
    for unitP=1:length(selectedUnits)
        selectWF=single(waveForms(:,handles.subset{selectedUnits(unitP)})');
        plot(mean(selectWF),'linewidth',2,'Color',cmap(unitP,:));
            wfSEM=std(selectWF)/ sqrt(size(selectWF,2)); %standard error of the mean
            wfSEM = wfSEM * 1.96; % 95% of the data will fall within 1.96 standard deviations of a normal distribution
        patch([1:length(wfSEM),fliplr(1:length(wfSEM))],...
            [mean(selectWF)-wfSEM,fliplr(mean(selectWF)+wfSEM)],...
            cmap(unitP,:),'EdgeColor','none','FaceAlpha',0.2);
        %duplicate mean unit over unsorted plot
        plot(handles.UnsortedUnits_Axes,mean(selectWF),'linewidth',2,'Color',cmap(unitP,:));
        patch([1:length(wfSEM),fliplr(1:length(wfSEM))],...
            [mean(selectWF)-wfSEM,fliplr(mean(selectWF)+wfSEM)],...
            cmap(unitP,:),'EdgeColor','none','FaceAlpha',0.5,'Parent', handles.UnsortedUnits_Axes);
        set(gca,'xtick',linspace(0,size(waveForms(:,handles.subset{selectedUnits(unitP)}),1),5),...
            'xticklabel',round(linspace(-round(size(waveForms(:,handles.subset{selectedUnits(unitP)}),1)/2),...
            round(size(waveForms(:,handles.subset{selectedUnits(unitP)}),1)/2),5)/(double(samplingRate)/1000),2),'TickDir','out');
        % legend('Unclustered waveforms','location','southeast')
        axis('tight');box off;
        xlabel('Time (ms)')
        ylabel('Voltage (?V)')
        set(gca,'Color','white','FontSize',10,'FontName','calibri');
    end
    
    %% plot rasters
    % downsample to 1 millisecond bins
    %         Spikes.Offline_Threshold.samplingRate(ChExN,2)=1000;
    %         Spikes.Offline_Threshold.type{ChExN,2}='downSampled';
    %         spikeTimeIdx=zeros(1,size(Spikes.Offline_Threshold.data{ChExN,1},2));
    %         spikeTimeIdx(Spikes.Offline_Threshold.data{ChExN,1})=1;
    %         spikeTimes=find(Spikes.Offline_Threshold.data{ChExN,1});
    %         binSize=1;
    %         numBin=ceil(size(spikeTimeIdx,2)/(Spikes.Offline_Threshold.samplingRate(ChExN,1)/Spikes.Offline_Threshold.samplingRate(ChExN,2))/binSize);
    %         % binspikeTime = histogram(double(spikeTimes), numBin); %plots directly histogram
    %         [Spikes.Offline_Threshold.data{ChExN,2},Spikes.Offline_Threshold.binEdges{ChExN}] = histcounts(double(spikeTimes), linspace(0,size(spikeTimeIdx,2),numBin));
    %         Spikes.Offline_Threshold.data{ChExN,2}(Spikes.Offline_Threshold.data{ChExN,2}>1)=1; %no more than 1 spike per ms
    %
    % plot 10 sec or 2000 waveforms max
    
end

%% --- Outputs from this function are returned to the command line.
function varargout = SpikeVisualizationGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%% --- Executes on button press in Spikes_SortOff_RB.
function Spikes_SortOff_RB_Callback(hObject, eventdata, handles)
% hObject    handle to Spikes_SortOff_RB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Spikes_SortOff_RB


%% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


%% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


%% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%% --- Executes on selection change in SelectChannel_LB.
function SelectChannel_LB_Callback(hObject, eventdata, handles)
% hObject    handle to SelectChannel_LB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Get chanel and unit selection
%     channelMenu=get(handles.SelectChannel_LB,'string');
%     channelSelected=get(handles.SelectChannel_LB,'value');
%     channelSelected=channelMenu(channelSelected);
%
%     unitMenu=get(handles.SelectUnit_LB,'string');
%     unitsSelected=get(handles.SelectUnit_LB,'value');
%     unitsSelected=unitMenu(unitsSelected);


% Hints: contents = cellstr(get(hObject,'String')) returns SelectChannel_LB contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SelectChannel_LB


%% --- Executes during object creation, after setting all properties.
function SelectChannel_LB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SelectChannel_LB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


%% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


%% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%% --- Executes on button press in Spikes_Th_RB.
function Spikes_Th_RB_Callback(hObject, eventdata, handles)
% hObject    handle to Spikes_Th_RB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Spikes_Th_RB


%% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


%% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


%% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Executes on button press in Spikes_SortOn_RB.
function Spikes_SortOn_RB_Callback(hObject, eventdata, handles)
% hObject    handle to Spikes_SortOn_RB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Spikes_SortOn_RB


%% --- Executes on button press in PB_GetSortedSpikes.
function PB_GetSortedSpikes_Callback(hObject, eventdata, handles)
% hObject    handle to PB_GetSortedSpikes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in SelectUnit_LB.
function SelectUnit_LB_Callback(hObject, eventdata, handles)
% hObject    handle to SelectUnit_LB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SelectUnit_LB contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SelectUnit_LB


% --- Executes during object creation, after setting all properties.
function SelectUnit_LB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SelectUnit_LB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function TW_slider_Callback(hObject, eventdata, handles)
% hObject    handle to TW_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function TW_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TW_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in TWplus_PB.
function TWplus_PB_Callback(hObject, eventdata, handles)
% hObject    handle to TWplus_PB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in TWminus_PB.
function TWminus_PB_Callback(hObject, eventdata, handles)
% hObject    handle to TWminus_PB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in TWall_PB.
function TWall_PB_Callback(hObject, eventdata, handles)
% hObject    handle to TWall_PB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse press over axes background.
function UnsortedUnits_Axes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to UnsortedUnits_Axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

electrodeNum=get(handles.SelectChannel_LB,'value');

%% initialize variables
unitsIdx=handles.Spikes.inGUI.Units{electrodeNum};
waveForms=handles.Spikes.inGUI.Waveforms{electrodeNum};
% spikeTimes=handles.Spikes.inGUI.SpikeTimes{electrodeNum};
% samplingRate=handles.Spikes.Online_Sorting.samplingRate(electrodeNum);

waveForms=waveForms(:,unitsIdx==0);

%make cluster classes 0 (unsorted), and -1 (hidden)
clusterClasses=ones(size(waveForms,2),1)*-1;
clusterClasses(handles.subset{1})=0;%mark subset visible (might be all of them)
clusterClasses=clusterClasses(unitsIdx==0);%remove those already sorted
% find lines that were added

InteractiveClassification(waveForms,clusterClasses,0); % viewClasses=0


% unitsIdx
%  Update handles structure
guidata(hObject, handles);
