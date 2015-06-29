function varargout = sim_100G_single_laser(varargin)
% SIM_100G_SINGLE_LASER MATLAB code for sim_100G_single_laser.fig
%      SIM_100G_SINGLE_LASER, by itself, creates a new SIM_100G_SINGLE_LASER or raises the existing
%      singleton*.
%
%      H = SIM_100G_SINGLE_LASER returns the handle to a new SIM_100G_SINGLE_LASER or the handle to
%      the existing singleton*.
%
%      SIM_100G_SINGLE_LASER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SIM_100G_SINGLE_LASER.M with the given input arguments.
%
%      SIM_100G_SINGLE_LASER('Property','Value',...) creates a new SIM_100G_SINGLE_LASER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before sim_100G_single_laser_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to sim_100G_single_laser_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sim_100G_single_laser

% Last Modified by GUIDE v2.5 29-Jun-2015 13:24:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sim_100G_single_laser_OpeningFcn, ...
                   'gui_OutputFcn',  @sim_100G_single_laser_OutputFcn, ...
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


% --- Executes just before sim_100G_single_laser is made visible.
function sim_100G_single_laser_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sim_100G_single_laser (see VARARGIN)

% Choose default command line output for sim_100G_single_laser
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes sim_100G_single_laser wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = sim_100G_single_laser_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in mod_type_popup.
function mod_type_popup_Callback(hObject, eventdata, handles)
% hObject    handle to mod_type_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns mod_type_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from mod_type_popup


% --- Executes during object creation, after setting all properties.
function mod_type_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mod_type_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in system_popup.
function system_popup_Callback(hObject, eventdata, handles)
% hObject    handle to system_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns system_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from system_popup


% --- Executes during object creation, after setting all properties.
function system_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to system_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in results_popup.
function results_popup_Callback(hObject, eventdata, handles)
% hObject    handle to results_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns results_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from results_popup


% --- Executes during object creation, after setting all properties.
function results_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to results_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
