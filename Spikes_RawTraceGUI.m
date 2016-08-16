function varargout = Spikes_RawTraceGUI(varargin)
% MATLAB code for Spikes_RawTraceGUI.fig


% Last Modified by GUIDE v2.5 08-Jun-2016 17:55:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Spikes_RawTraceGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @Spikes_RawTraceGUI_OutputFcn, ...
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


%% --- Executes just before Spikes_RawTraceGUI is made visible.
function Spikes_RawTraceGUI_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for Spikes_RawTraceGUI
handles.output = hObject;

if ~isempty(varargin)
    handles=catstruct(handles,varargin{:}); % catstruct available here: 
    % http://www.mathworks.com/matlabcentral/fileexchange/7842-catstruct
    
    if isfield(handles,'fname')
        set(handles.FileName,'string',handles.fname(1:end-4));
    else
        set(handles.FileName,'string','');
    end
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Spikes_RawTraceGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


%% --- Outputs from this function are returned to the command line.
function varargout = Spikes_RawTraceGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%% --- Executes on button press in PB_GetSortedSpikes.
function PB_GetSortedSpikes_Callback(hObject, eventdata, handles)
% hObject    handle to PB_GetSortedSpikes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
