function varargout = kakuroGUI(varargin)
% KAKUROGUI M-file for kakuroGUI.fig
%      KAKUROGUI, by itself, creates a new KAKUROGUI or raises the existing
%      singleton*.
%
%      H = KAKUROGUI returns the handle to a new KAKUROGUI or the handle to
%      the existing singleton*.
%
%      KAKUROGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in KAKUROGUI.M with the given input arguments.
%
%      KAKUROGUI('Property','Value',...) creates a new KAKUROGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before kakuroGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to kakuroGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help kakuroGUI

% Last Modified by GUIDE v2.5 10-Mar-2011 17:09:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @kakuroGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @kakuroGUI_OutputFcn, ...
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


% --- Executes just before kakuroGUI is made visible.
function kakuroGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to kakuroGUI (see VARARGIN)

% Choose default command line output for kakuroGUI
handles.output = hObject;
handles.rowSize = 3;   % Kakuro row size
handles.colSize = 3;   % Kakuro column size

global pMat;
global inputErr;
pMat = [];  % Kakuro problem matrix
inputErr = 0;   % input error count

% create Edit fields
global gap;
global width;
global height;

temp1 = get(handles.rowMenu, 'String');
temp2 = get(handles.colMenu, 'String');
maxRow = str2double(temp1{end});
maxCol = str2double(temp2{end});
handles.editFields = zeros(maxRow, maxCol);

gap = 5;    % gap between edit fields
height = 25;
width = 40; % width of edit fields

yField = gap + 30;
for ri = 1:maxRow,
    xField = gap;
    for ci = 1:maxCol,
        handles.editFields(ri, ci) = uicontrol(...%'Parent', pane, ... 
                        'Style', 'edit', 'Position', [xField, yField, width, height], ...
                        'CallBack', @validateInput, 'BackgroundColor', 'white', ...
                        'FontSize', 10, 'Visible', 'off', 'ButtonDownFcn', @toggleEnable, ...
                        'Tag', sprintf('%.0f,%.0f', ri, ci));  % use Tag to refer to row and col indices
        xField = xField + width;
    end;
    
    yField = yField + width;
end;

handles = clearAndShowEditFields(handles, [0, 0]);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes kakuroGUI wait for user response (see UIRESUME)
% uiwait(handles.mainFig);


% --- Outputs from this function are returned to the command line.
function varargout = kakuroGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function handles = clearAndShowEditFields(handles, prevSize)
global inputErr;
global pMat;

global height;
global width;
global gap;

clearPos = get(handles.clearButton, 'Position');
solvePos = get(handles.solveButton, 'Position');

% minimum main window width taken to be the width of the controls
if (~isfield(handles, 'minMainWidth') || isempty(handles.minMainWidth)),
    tempPos = get(handles.rowMenu, 'Position');
    handles.minMainWidth = tempPos(3);
    tempPos = get(handles.colMenu, 'Position');
    handles.minMainWidth = handles.minMainWidth + tempPos(3);
    tempPos = get(handles.label, 'Position');
    handles.minMainWidth = handles.minMainWidth + tempPos(3);

    handles.minMainWidth = handles.minMainWidth + clearPos(3);

    handles.minMainWidth = handles.minMainWidth + solvePos(3) + 6*gap;
end;

gridWidth = handles.colSize*width + 2*gap;
extraWidth = handles.minMainWidth - gridWidth;

if (extraWidth > 0),
    mainWidth = handles.minMainWidth;
    gridXOffset = extraWidth/2;
else
    mainWidth = gridWidth;
    gridXOffset = 0;
end;

% reposition clear button
set(handles.clearButton, 'Position', [(mainWidth-solvePos(3)-clearPos(3)-2*gap), gap, clearPos(3), clearPos(4)]);

% reposition solve button
set(handles.solveButton, 'Position', [(mainWidth-solvePos(3)-gap), gap, solvePos(3),solvePos(4)]);

% resize main window
fig = handles.mainFig;
pos = get(fig, 'Position');
% pos(1) = pos(1) + mainWidth - pos(3);   % reposition the origin upon window resize

mainHeight = handles.rowSize*height+gap+30;
pos(2) = pos(2) - mainHeight + pos(4);  % reposition the y-component of origin upon window resize

set(fig, 'Position', [pos(1), pos(2), mainWidth, mainHeight]);

yField = mainHeight - gap - 25;
for ri = 1:handles.rowSize,
    xField = gap + gridXOffset;
    for ci = 1:handles.colSize,
        set(handles.editFields(ri, ci), 'Visible', 'on', 'String', '', 'FontWeight', 'normal', 'BackgroundColor', 'white', ...
            'Position', [xField, yField, width, height], 'Enable', 'on'); 
        
        xField = xField + width;
    end;
    
    yField = yField - height;
end;

if (prevSize(1) > handles.rowSize), % if previous rowSize > current rowSize, make the edit field invisible
    for ri = (handles.rowSize+1):prevSize(1),
        for ci = 1:handles.colSize,
            set(handles.editFields(ri, ci), 'Visible', 'off');
        end;
    end;
end;

if (prevSize(2) > handles.colSize), % if previous colSize > current colSize, make the edit field invisible
    for ri = 1:handles.rowSize,
        for ci = (handles.colSize+1):prevSize(2),
            set(handles.editFields(ri, ci), 'Visible', 'off');
        end;
    end;
end;

inputErr = 0;
pMat = zeros(handles.rowSize, handles.colSize);

% mouse right click to disable/enable edit fields
function toggleEnable(hObject, eventdata)

global pMat;

tag = get(hObject, 'Tag');
ind = sscanf(tag, '%d,%d');
    
status = get(hObject, 'Enable');
if (strcmpi(status, 'on')),
    set(hObject, 'Enable', 'off', 'String', '', 'BackgroundColor', 'white', 'FontWeight', 'normal');
    pMat(ind(1), ind(2)) = -1;
else
    set(hObject, 'Enable', 'on');
    pMat(ind(1), ind(2)) = 0;
end;



% function used for edit fields callback for validating inputs
function validateInput(hObject, eventdata)
global inputErr;
global pMat;

val = get(hObject, 'String');

color = get(hObject, 'BackgroundColor');
if (~isempty(val)),
    val = regexp(val,'^[0-9 \\]*$', 'match');   % only these characters are allowed
    if (~isempty(val)),
        val = val{1};
        a = regexp(val,'\s*(\d*)\s*\\?\s*(\d*)\s*','tokens');
        num = str2double(a{1});
        nan = isnan(num);
    end;
    
    if (isempty(val) || sum(nan) == 2), % if both tokens are invalid
        if (sum(color ~= [1,0,0]) >= 1),    % if bg color is not red
            inputErr = inputErr + 1;
            set(hObject, 'BackgroundColor', 'red', 'FontWeight', 'bold');
        end;
    else
        % correct error
        if (sum(color == [1,0,0]) == 3),    % if bg color is red
            inputErr = inputErr - 1;
            set(hObject, 'BackgroundColor', 'white', 'FontWeight', 'normal');
        end;
        
        if (sum(num == 0) >= 1),
            set(hObject, 'BackgroundColor', 'white', 'FontWeight', 'normal');
        else
            set(hObject, 'BackgroundColor', [0.9, 0.99, 0.99], 'FontWeight', 'bold');
        end;

        aVal = 0;
        if (~nan(1)),   % vertical sum available
            aVal = aVal + num(1);
        end;

        if (~nan(2)),   % horizontal suym avaiable
            aVal = aVal + 1i*num(2);
        end;

        % update input matrix
        tag = get(hObject, 'Tag');
        ind = sscanf(tag, '%d,%d');
        pMat(ind(1), ind(2)) = aVal;
    end;
else
    set(hObject, 'BackgroundColor', 'white', 'FontWeight', 'normal');
end;



% --- Executes on button press in clearButton.
function clearButton_Callback(hObject, eventdata, handles)
% hObject    handle to clearButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% rowSize = str2double();
% colSize = str2double(get(handles.colMenu, 'String'));

handles = clearAndShowEditFields(handles, [handles.rowSize, handles.colSize]);

guidata(hObject, handles);

% --- Executes on button press in solveButton.
function solveButton_Callback(hObject, eventdata, handles)
% hObject    handle to solveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global pMat;
global inputErr;

if (inputErr == 0),
    sol = solveKakuro(pMat);
    [rI, cI] = find(sol > 0 & imag(sol) == 0);

    for ii = 1:length(rI),
        set(handles.editFields(rI(ii), cI(ii)), 'String', num2str(sol(rI(ii), cI(ii))));
    end;
else
    str = sprintf('Invalid input(s) in specified Kakuro problem. Please check.');
    msgbox(str, 'Error in Kakuro problem!', 'error');
end;

% --- Executes on selection change in rowMenu.
function rowMenu_Callback(hObject, eventdata, handles)
% hObject    handle to rowMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = get(hObject,'String');
prevRow = handles.rowSize;
handles.rowSize = str2double(contents{get(hObject,'Value')});

handles = clearAndShowEditFields(handles, [prevRow, handles.colSize]);

guidata(hObject, handles);


% Hints: contents = get(hObject,'String') returns rowMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from rowMenu


% --- Executes during object creation, after setting all properties.
function rowMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rowMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in colMenu.
function colMenu_Callback(hObject, eventdata, handles)
% hObject    handle to colMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns colMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from colMenu
contents = get(hObject,'String');
prevCol = handles.colSize;
handles.colSize = str2double(contents{get(hObject,'Value')});

handles = clearAndShowEditFields(handles, [handles.rowSize, prevCol]);

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function colMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to colMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function file_Callback(hObject, eventdata, handles)
% hObject    handle to file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function about_Callback(hObject, eventdata, handles)
% hObject    handle to about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str = sprintf('Written by: LCC Corps \nData: Mar 10, 2011 \nDedicated to: ...');
msgbox(str, 'About', 'help');

% --------------------------------------------------------------------
function exit_Callback(hObject, eventdata, handles)
% hObject    handle to exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcf);
delete(gcf);

