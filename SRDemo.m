function varargout = SRDemo(varargin)
% SRDEMO M-file for SRDemo.fig
%      SRDEMO, by itself, creates a new SRDEMO or raises the existing
%      singleton*.
%
%      H = SRDEMO returns the handle to a new SRDEMO or the handle to
%      the existing singleton*.
%
%      SRDEMO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SRDEMO.M with the given input arguments.
%
%      SRDEMO('Property','Value',...) creates a new SRDEMO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SRDemo_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SRDemo_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SRDemo

% Last Modified by GUIDE v2.5 04-Aug-2007 16:12:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
  'gui_Singleton',  gui_Singleton, ...
  'gui_OpeningFcn', @SRDemo_OpeningFcn, ...
  'gui_OutputFcn',  @SRDemo_OutputFcn, ...
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


% --- Executes just before SRDemo is made visible.
function SRDemo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SRDemo (see VARARGIN)

addpath([pwd '\LKOFlow']);

% Choose default command line output for SRDemo
handles.output = hObject;
axis(handles.axesLR,'off');
axis(handles.axesHR,'off');

handles.prevHR = [];
handles.HR = [];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SRDemo wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SRDemo_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in cmdLoad.
function cmdLoad_Callback(hObject, eventdata, handles)
% hObject    handle to cmdLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Open dialog to retreive the filename
[FileName,PathName] = uigetfile('*.mat','Select the movie MAT file');

if FileName ~= 0

  handles.LR=LoadVideoMat([PathName FileName]);
  handles.LRDisplayI = 1;

  set(handles.cmdRegister, 'enable', 'on');
  set(handles.cmdSR, 'enable', 'off');
  set(handles.cmdClear, 'enable', 'off');
  set(handles.cmdSave, 'enable', 'off');
  set(handles.cmdSaveLR, 'enable', 'on');
  set(handles.cbShowReg, 'enable', 'off');

  % Update LR display
  UpdateLRDisplay(hObject, handles)

end

function LRScanButtonEnable(hObject, handles)

if handles.LRDisplayI == 1
  set(handles.cmdPrev, 'enable', 'off');
else
  set(handles.cmdPrev, 'enable', 'on');
end
if handles.LRDisplayI == size(handles.LR, 3)
  set(handles.cmdNext, 'enable', 'off');
else
  set(handles.cmdNext, 'enable', 'on');
end

% Update handles structure
guidata(hObject, handles);

function UpdateLRDisplay(hObject, handles)

axes(handles.axesLR);
if get(handles.cbShowReg, 'Value')
  imagesc(handles.LR(:,:,handles.LRDisplayI));colormap('gray')
else
  tform = maketform('affine',[1 0 0; 0 1 0; handles.D(handles.LRDisplayI, :) 1]);
  sz = size(handles.LR(:,:,handles.LRDisplayI));
  img = imtransform(squeeze(handles.LR(:,:,handles.LRDisplayI)),tform,'bilinear','XData',[1 sz(2)], 'YData',[1 sz(1)]);
  imagesc(img);colormap('gray')
end

axis(handles.axesLR,'off');

set(handles.lblLRImg, 'String', sprintf('Low Resolution Image: %u out of %u', handles.LRDisplayI, size(handles.LR,3)));

% Update scan buttons
LRScanButtonEnable(hObject, handles)

% --- Executes on button press in cmdNext.
function cmdNext_Callback(hObject, eventdata, handles)
% hObject    handle to cmdNext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.LRDisplayI = handles.LRDisplayI + 1;

% Update LR display
UpdateLRDisplay(hObject, handles)


% --- Executes on button press in cmdPrev.
function cmdPrev_Callback(hObject, eventdata, handles)
% hObject    handle to cmdPrev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.LRDisplayI = handles.LRDisplayI - 1;

% Update LR display
UpdateLRDisplay(hObject, handles)


% --- Executes on button press in cmdSave.
function cmdSave_Callback(hObject, eventdata, handles)
% hObject    handle to cmdSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uiputfile('*.jpg','Save image file');

if FileName ~= 0

  if get(handles.cbSelectHR, 'Value')
    imwrite(uint8(handles.HR), [PathName FileName]);
  else
    imwrite(uint8(handles.prevHR), [PathName FileName]);
  end
end

% --- Executes on button press in cmdClear.
function cmdClear_Callback(hObject, eventdata, handles)
% hObject    handle to cmdClear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla(handles.axesHR);

set(handles.cmdClear, 'enable', 'off');
set(handles.cmdSave, 'enable', 'off');
set(handles.rbRobust, 'enable', 'off');

if get(handles.gbSRType, 'SelectedObject') == handles.rbRobust
  set(handles.gbSRType, 'SelectedObject', handles.rbFast);
end

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in cmdRegister.
function cmdRegister_Callback(hObject, eventdata, handles)
% hObject    handle to cmdRegister (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Register the image sequence to base image (frame 1)
handles.D=RegisterImageSeq(handles.LR);

set(handles.cmdSR, 'enable', 'on');
set(handles.cbShowReg, 'enable', 'on');

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in cmdSR.
function cmdSR_Callback(hObject, eventdata, handles)
% hObject    handle to cmdSR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[props, resFactor, D, LR, Hpsf] = CollectParms(hObject, handles);

handles.prevHR = handles.HR;

% Check the selected algorithm 
switch get(handles.gbSRType, 'SelectedObject')
  
  case handles.rbSpline
   
    handles.HR=SplineSRInterp(LR, resFactor, Hpsf, props);
    
  case handles.rbRobust

    handles.HR=RobustSR(LR(3:end-2,3:end-2,:), D, handles.HR, resFactor, Hpsf, props);

  case handles.rbFast

    handles.HR=FastRobustSR(LR(3:end-2,3:end-2,:), D, resFactor, Hpsf, props);
    
end

DisplayHRImage(hObject, handles);

set(handles.cmdClear, 'enable', 'on');
set(handles.cmdSave, 'enable', 'on');
set(handles.rbRobust, 'enable', 'on');

% Update handles structure
guidata(hObject, handles);


function [props, resFactor, D, LR, Hpsf] = CollectParms(hObject, handles)

try

  resFactor = str2double(get(handles.txtResFactor, 'String'));

  psfSize = str2double(get(handles.txtPsfSize, 'String'));
  psfSig =  str2double(get(handles.txtPsfSig, 'String'));
  Hpsf = fspecial('gaussian', [psfSize psfSize], psfSig);


  props.alpha = str2double(get(handles.txtAlpha, 'String'));
  props.beta = str2double(get(handles.txtBeta, 'String'));
  props.lambda = str2double(get(handles.txtLambda, 'String'));
  props.P = str2double(get(handles.txtP, 'String'));
  props.maxIter = str2double(get(handles.txtIterNum, 'String'));

  % Round translation to nearest neighbor
  D=round(handles.D.*resFactor);

  % Shift all images so D is bounded from 0-resFactor
  Dr=floor(D/resFactor);
  D=mod(D,resFactor)+resFactor;

  LR = handles.LR;
  [X,Y]=meshgrid(1:size(LR, 2), 1:size(LR, 1));

  for i=1:size(LR, 3)

    LR(:,:,i)=interp2(X+Dr(i,1), Y+Dr(i,2), LR(:,:,i), X, Y, '*nearest');

  end


catch

  err = lasterror;
  errordlg(err.message,'Parsing error');
end




function txtResFactor_Callback(hObject, eventdata, handles)
% hObject    handle to txtResFactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtResFactor as text
%        str2double(get(hObject,'String')) returns contents of txtResFactor as a double


% --- Executes during object creation, after setting all properties.
function txtResFactor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtResFactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end



function txtPsfSize_Callback(hObject, eventdata, handles)
% hObject    handle to txtPsfSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtPsfSize as text
%        str2double(get(hObject,'String')) returns contents of txtPsfSize as a double


% --- Executes during object creation, after setting all properties.
function txtPsfSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtPsfSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end



function txtAlpha_Callback(hObject, eventdata, handles)
% hObject    handle to txtAlpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtAlpha as text
%        str2double(get(hObject,'String')) returns contents of txtAlpha as a double


% --- Executes during object creation, after setting all properties.
function txtAlpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtAlpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end



function txtBeta_Callback(hObject, eventdata, handles)
% hObject    handle to txtBeta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtBeta as text
%        str2double(get(hObject,'String')) returns contents of txtBeta as a double


% --- Executes during object creation, after setting all properties.
function txtBeta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtBeta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end



function txtLambda_Callback(hObject, eventdata, handles)
% hObject    handle to txtLambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtLambda as text
%        str2double(get(hObject,'String')) returns contents of txtLambda as a double


% --- Executes during object creation, after setting all properties.
function txtLambda_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtLambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end



function txtP_Callback(hObject, eventdata, handles)
% hObject    handle to txtP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtP as text
%        str2double(get(hObject,'String')) returns contents of txtP as a double


% --- Executes during object creation, after setting all properties.
function txtP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end



function txtIterNum_Callback(hObject, eventdata, handles)
% hObject    handle to txtIterNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtIterNum as text
%        str2double(get(hObject,'String')) returns contents of txtIterNum as a double


% --- Executes during object creation, after setting all properties.
function txtIterNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtIterNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end





function txtPsfSig_Callback(hObject, eventdata, handles)
% hObject    handle to txtPsfSig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtPsfSig as text
%        str2double(get(hObject,'String')) returns contents of txtPsfSig as a double


% --- Executes during object creation, after setting all properties.
function txtPsfSig_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtPsfSig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end

function DisplayHRImage(hObject, handles)
axes(handles.axesHR);

if get(handles.cbSelectHR, 'Value')
  imagesc(handles.HR);
else
  imagesc(handles.prevHR);
end

colormap('gray')
axis(handles.axesHR,'off');


% --- Executes on button press in cbSelectHR.
function cbSelectHR_Callback(hObject, eventdata, handles)
% hObject    handle to cbSelectHR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbSelectHR
DisplayHRImage(hObject, handles)

if get(handles.cbSelectHR, 'Value')
  set(handles.cbSelectHR, 'String', 'Displaying Current Image');
else
  set(handles.cbSelectHR, 'String', 'Displaying Previous Image');
end


% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in cmdSaveLR.
function cmdSaveLR_Callback(hObject, eventdata, handles)
% hObject    handle to cmdSaveLR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[FileName,PathName] = uiputfile('*.jpg','Save image file');

if FileName ~= 0

  imwrite(uint8(handles.LR(:,:,handles.LRDisplayI)), [PathName FileName]);
end



% --- Executes on button press in cbShowReg.
function cbShowReg_Callback(hObject, eventdata, handles)
% hObject    handle to cbShowReg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbShowReg
UpdateLRDisplay(hObject, handles)

if get(handles.cbShowReg, 'Value')
  set(handles.cbShowReg, 'String', 'Displaying Unregistered Images');
else
  set(handles.cbShowReg, 'String', 'Displaying Registered Images');
end


% Update handles structure
guidata(hObject, handles);


