function varargout = GUIPlot(varargin)
% GUIPLOT MATLAB code for GUIPlot.fig
%      GUIPLOT, by itself, creates a new GUIPLOT or raises the existing
%      singleton*.
%
%      H = GUIPLOT returns the handle to a new GUIPLOT or the handle to
%      the existing singleton*.
%
%      GUIPLOT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIPLOT.M with the given input arguments.
%
%      GUIPLOT('Property','Value',...) creates a new GUIPLOT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUIPlot_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUIPlot_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUIPlot

% Last Modified by GUIDE v2.5 19-Dec-2017 09:29:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUIPlot_OpeningFcn, ...
                   'gui_OutputFcn',  @GUIPlot_OutputFcn, ...
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


% --- Executes just before GUIPlot is made visible.
function GUIPlot_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUIPlot (see VARARGIN)

% Choose default command line output for GUIPlot
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUIPlot wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUIPlot_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Plot_PPP.
function Plot_PPP_Callback(hObject, eventdata, handles)
% hObject    handle to Plot_PPP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global newpath;
global ylimppp;
ylimppp = str2double(get(handles.edit_ylimppp,'String'));
[filenames,pathname] = uigetfile('*.pos','open ppp files','MultiSelect','on',newpath);
if ~iscell(filenames)
    fullname = [pathname filenames];
    if judge_empytfile(fullname)==0
        return;
    end
    plot_neu(pathname,filenames,'ppp',ylimppp);
else
    for i=1:max(size(filenames))
        filename=filenames{i};
        fullname = [pathname filename];
        if judge_empytfile(fullname)==0
            continue;
        end
        plot_neu(pathname,filename,'ppp',ylimppp);
    end
end;
newpath = pathname;


function res = judge_empytfile(fullname)
sdir = dir(fullname);  
ss1 = sdir.bytes;
if ss1 == 0
    res = 0;
else
    res = 1;
end


function edit_ylimppp_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ylimppp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ylimppp as text
%        str2double(get(hObject,'String')) returns contents of edit_ylimppp as a double


% --- Executes during object creation, after setting all properties.
function edit_ylimppp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ylimppp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plot_nsat.
function plot_nsat_Callback(hObject, eventdata, handles)
% hObject    handle to plot_nsat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global newpath;
[filenames,pathname] = uigetfile('*.pdop','open nsat files','MultiSelect','on',newpath);
n= max(size(filenames));
if ~iscell(filenames)
    fullname = [pathname filenames];
    if judge_empytfile(fullname)==0
        return;
    end
    plot_nsat(pathname,filenames,'nsat+pdop');
else
    for i=1:n
        filename=filenames{i};
        fullname = [pathname filename];
        if judge_empytfile(fullname)==0
            continue;
        end
        plot_nsat(pathname,filename,'nsat+pdop');
    end
end;
newpath = pathname;


% --- Executes on button press in Plot_SPP.
function Plot_SPP_Callback(hObject, eventdata, handles)
% hObject    handle to Plot_SPP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global newpath;
global ylimspp;
ylimspp = str2double(get(handles.edit_ylimspp,'String'));
[filenames,pathname] = uigetfile('*.spp','open spp files','MultiSelect','on',newpath);
if ~iscell(filenames)
    fullname = [pathname filenames];
    if judge_empytfile(fullname)==0
        return;
    end
    plot_neu(pathname,filenames,'spp',ylimspp);
else
    for i=1:max(size(filenames))
        filename=filenames{i};
        fullname = [pathname filename];
        if judge_empytfile(fullname)==0
            continue;
        end
        plot_neu(pathname,filename,'spp',ylimspp);
    end
end;
newpath = pathname;


function edit_ylimspp_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ylimspp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ylimspp as text
%        str2double(get(hObject,'String')) returns contents of edit_ylimspp as a double


% --- Executes during object creation, after setting all properties.
function edit_ylimspp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ylimspp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Plot_resx.
function Plot_resx_Callback(hObject, eventdata, handles)
% hObject    handle to Plot_resx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global newpath;
[filenames,pathname] = uigetfile('*.res*','open resx files','MultiSelect','on',newpath);
n= max(size(filenames));
if ~iscell(filenames)
    fullname = [pathname filenames];
    if judge_empytfile(fullname)==0
        return;
    end
    if filenames(n)=='c'
        plot_resx(pathname,filenames,'code');
    else
        plot_resx(pathname,filenames,'phase');
    end
else
    for i=1:n
        filename=filenames{i};
        fullname = [pathname filename];
        if judge_empytfile(fullname)==0
            continue;
        end
        n1 = max(size(filename));
        if filename(n1-1)=='c'
            plot_resx(pathname,filename,'code');
        else
            plot_resx(pathname,filename,'phase');
        end
    end
end;
newpath = pathname;


% --- Executes on button press in Plot_ztd.
function Plot_ztd_Callback(hObject, eventdata, handles)
% hObject    handle to Plot_ztd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global newpath;
global ylimztd;
ylimztd = str2double(get(handles.edit_ylimztd,'String'));
[filenames,pathname] = uigetfile('*.dtrp','open ztd files','MultiSelect','on',newpath);
if ~iscell(filenames)
        fullname = [pathname filenames];
        if judge_empytfile(fullname)==0
            return;
        end
    plot_ztd(pathname,filenames,'ztd',ylimztd);
else
    for i=1:max(size(filenames))
        filename=filenames{i};
        fullname = [pathname filename];
        if judge_empytfile(fullname)==0
            continue;
        end
        plot_ztd(pathname,filename,'ztd',ylimztd);
    end
end;
newpath = pathname;


function edit_ylimztd_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ylimztd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ylimztd as text
%        str2double(get(hObject,'String')) returns contents of edit_ylimztd as a double


% --- Executes during object creation, after setting all properties.
function edit_ylimztd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ylimztd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Plot_wlamb.
function Plot_wlamb_Callback(hObject, eventdata, handles)
% hObject    handle to Plot_wlamb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global newpath;
[filenames,pathname] = uigetfile('*.wlamb*','open resx files','MultiSelect','on',newpath);
n= max(size(filenames));
if ~iscell(filenames)
    fullname = [pathname filenames];
    if judge_empytfile(fullname)==0
        return;
    end
    if filenames(n)=='o'
        plot_wlamb(pathname,filenames,'no');
    else
        plot_wlamb(pathname,filenames,'yes');
    end
else
    for i=1:n
        filename=filenames{i};
        fullname = [pathname filename];
        if judge_empytfile(fullname)==0
            continue;
        end
        n1 = max(size(filename));
        if filename(n1)=='o'
            plot_resx(pathname,filename,'no');
        else
            plot_resx(pathname,filename,'yes');
        end
    end
end;
newpath = pathname;


% --- Executes on button press in Plot_ion.
function Plot_ion_Callback(hObject, eventdata, handles)
% hObject    handle to Plot_ion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global newpath;
[filenames,pathname] = uigetfile('*.stec','open sTEC files','MultiSelect','on',newpath);
if ~iscell(filenames)
        fullname = [pathname filenames];
        if judge_empytfile(fullname)==0
            return;
        end
    plot_ion(pathname,filenames,'stec');
else
    for i=1:max(size(filenames))
        filename=filenames{i};
        fullname = [pathname filename];
        if judge_empytfile(fullname)==0
            continue;
        end
        plot_ion(pathname,filename,'stec');
    end
end;
newpath = pathname;


% --- Executes on button press in Analysis.
function Analysis_Callback(hObject, eventdata, handles)
% hObject    handle to Analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global newpath;
[filenames,pathname] = uigetfile('*.pos','open ppp files','MultiSelect','on',newpath);
anafile = [pathname , 'analysis.ana'];
if exist(anafile,'file') > 0
    delete(anafile);
end;
anafile1 = [pathname , 'analysis.ana1'];
if exist(anafile1,'file') > 0
    delete(anafile1);
end;
n= max(size(filenames));
if ~iscell(filenames)
    fullname = [pathname filenames];
    if judge_empytfile(fullname)==0
        return;
    end
    [cvg_n,cvg_e,cvg_u,rms_n,rms_e,rms_u] = analysis(pathname,filenames,'ppp');
    disp(filenames(1:4));
else
    for i=1:n
        filename=filenames{i};
        fullname = [pathname filename];
        if judge_empytfile(fullname)==0
            continue;
        end
        [cvg_n,cvg_e,cvg_u,rms_n,rms_e,rms_u] = analysis(pathname,filename,'ppp');
        disp(filename(1:4));
    end
end;
disp('***end***');
% plot_bar(pathname,'analysis.ana1');
newpath = pathname;


% --- Executes on button press in Plot_any.
function Plot_any_Callback(hObject, eventdata, handles)
% hObject    handle to Plot_any (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global newpath;
global ylimppp;
global ylimspp;
global ylimztd;
ylimppp = str2double(get(handles.edit_ylimppp,'String'));
ylimspp = str2double(get(handles.edit_ylimspp,'String'));
ylimztd = str2double(get(handles.edit_ylimztd,'String'));
%clear;
[filenames,pathname] = uigetfile({'*.spp;*.pos;*.dtrp;*.stec;*.resc*;*.resp*;*.nspp;*.pdop;*.wlamb_no;*.wlamb_yes', ...
    'All result Files'; '*.*','All Files'},'open ALL files','MultiSelect','on',newpath);
judgefilenames(pathname,filenames)
newpath = pathname;

function judgefilenames(pathname,filenames)
global ylimppp;
global ylimspp;
global ylimztd;
n= max(size(filenames));
if ~iscell(filenames)
    %empty file?
    fullname = [pathname filenames];
    if judge_empytfile(fullname)==0
        return;
    end
    if strcmp(filenames(n-3:n),'.pos')
        plot_neu(pathname,filenames,'ppp',ylimppp);
    elseif strcmp(filenames(n-3:n),'.spp')
        plot_neu(pathname,filenames,'spp',ylimspp);
    elseif strcmp(filenames(n-4:n),'.dtrp')
        plot_ztd(pathname,filenames,'ztd',ylimztd);
    elseif strcmp(filenames(n-4:n),'.stec')
        plot_ion(pathname,filenames,'stec');
    elseif strcmp(filenames(n-5:n-1),'.resc')
        plot_resx(pathname,filenames,'code');
    elseif strcmp(filenames(n-5:n-1),'.resp')
        plot_resx(pathname,filenames,'phase');
    elseif strcmp(filenames(n-4:n),'.nspp')
        plot_nsat(pathname,filenames,'nspp');
    elseif strcmp(filenames(n-4:n),'.pdop')
        plot_nsat(pathname,filenames,'nppp');
    elseif strcmp(filenames(n-8:n),'.wlamb_no')
        plot_wlamb(pathname,filenames,'no');
    elseif strcmp(filenames(n-9:n),'.wlamb_yes')
        plot_wlamb(pathname,filenames,'yes');
    end
else
    for i=1:n
        filename=filenames{i}; 
        %empty file?
        fullname = [pathname filename];
        if judge_empytfile(fullname)==0
            continue;
        end
        n1 = max(size(filename));
        if strcmp(filename(n1-3:n1),'.pos')
            plot_neu(pathname,filename,'ppp',ylimppp);
        elseif strcmp(filename(n1-3:n1),'.spp')
            plot_neu(pathname,filename,'spp',ylimspp);
        elseif strcmp(filename(n1-4:n1),'.dtrp')
            plot_ztd(pathname,filename,'ztd',ylimztd);
        elseif strcmp(filename(n1-4:n1),'.stec')
            plot_ion(pathname,filename,'stec');
        elseif strcmp(filename(n1-5:n1-1),'.resc')
            plot_resx(pathname,filename,'code');
        elseif strcmp(filename(n1-5:n1-1),'.resp')
            plot_resx(pathname,filename,'phase');
        elseif strcmp(filename(n1-4:n1),'.nspp')
            plot_nsat(pathname,filename,'nspp');
        elseif strcmp(filename(n1-4:n1),'.pdop')
            plot_nsat(pathname,filename,'nppp');
        elseif strcmp(filename(n1-8:n1),'.wlamb_no')
            plot_wlamb(pathname,filename,'no');
        elseif strcmp(filename(n1-9:n1),'.wlamb_yes')
            plot_wlamb(pathname,filename,'yes');
        end
    end
end
