function varargout = BrainStateGUI(varargin)
% BRAINSTATEGUI MATLAB code for BrainStateGUI.fig
%      BRAINSTATEGUI, by itself, creates a new BRAINSTATEGUI or raises the existing
%      singleton*.
%
%      H = BRAINSTATEGUI returns the handle to a new BRAINSTATEGUI or the handle to
%      the existing singleton*.
%
%      BRAINSTATEGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BRAINSTATEGUI.M with the given input arguments.
%
%      BRAINSTATEGUI('Property','Value',...) creates a new BRAINSTATEGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BrainStateGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BrainStateGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BrainStateGUI

% Last Modified by GUIDE v2.5 24-Aug-2016 15:08:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BrainStateGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @BrainStateGUI_OutputFcn, ...
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


% --- Executes just before BrainStateGUI is made visible.
function BrainStateGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BrainStateGUI (see VARARGIN)

% Choose default command line output for BrainStateGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes BrainStateGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = BrainStateGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on selection change in SelectElectrode_LB.
function SelectElectrode_LB_Callback(hObject, eventdata, handles)
% hObject    handle to SelectElectrode_LB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SelectElectrode_LB contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SelectElectrode_LB


% --- Executes during object creation, after setting all properties.
function SelectElectrode_LB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SelectElectrode_LB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function TW_Slider_Callback(hObject, eventdata, handles)
% hObject    handle to TW_Slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function TW_Slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TW_Slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in loadPB.
function loadPB_Callback(hObject, eventdata, handles)
% hObject    handle to loadPB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in ReloadPB.
function ReloadPB_Callback(hObject, eventdata, handles)
% hObject    handle to ReloadPB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in SavePB.
function SavePB_Callback(hObject, eventdata, handles)
% hObject    handle to SavePB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in PlusTW.
function PlusTW_Callback(hObject, eventdata, handles)
% hObject    handle to PlusTW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in MinusTW.
function MinusTW_Callback(hObject, eventdata, handles)
% hObject    handle to MinusTW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in allTraceTW.
function allTraceTW_Callback(hObject, eventdata, handles)
% hObject    handle to allTraceTW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
