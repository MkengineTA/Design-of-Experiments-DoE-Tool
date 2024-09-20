function FinalDesign = createExperimentDesign()
% CREATEEXPERIMENTDESIGN adds a UI query for chemical experiments
% to the CALCULATEOPTIMALREDUCEDDESIGN function.
%
% The chemical variables can be any of the following (case-insensitive) strings:
%     - 'Temperature'
%     - 'Time'
%     - 'Amount'
%     - 'Species'
%
% There is no limit for the number of factors or factor-levels.
%
%   See also CALCULATEOPTIMALREDUCEDDESIGN, ROWEXCH, CANDEXCH, UITABLE
%
%   Copyright © 2020 by Tobias Averbeck, tobiasav92@gmail.com
% ====================================================================

clc
clearvars

% Check whether the buffer file exists and delete it if so.
if isfile('Cache.fig')
    delete Cache.fig
end

% Start UI.
z = InputQuery();
uiwait(z);
fig = openfig('Cache.fig','invisible');

% Get data from the buffer file.
dataObjs = findobj(fig);
FactorTable = dataObjs(end);
Data = FactorTable.Data;
[~,NumberOfFactors] = size(Data);
levels = zeros(1,NumberOfFactors);

Positions = zeros(NumberOfFactors,1);
for k = 1:NumberOfFactors
    CheckBox = dataObjs(k+2);
    Positions(k) = CheckBox.Position(1);
end
% Sometimes the positions in the buffer file are mixed up. These lines fix that.
if Positions(1) > Positions(end)
    mirror = 1;
else
    mirror = 0;
end

retain = zeros(NumberOfFactors,1);
Cache_Matrix = cell(1,NumberOfFactors);
for k = 1:NumberOfFactors
    CheckBox = dataObjs(k+2);
    if mirror == 1
        retain(end-k+1,1) = CheckBox.Value;
    else
        retain(k,1) = CheckBox.Value;
    end
    levels(1,k) = length(find(~cellfun(@isempty,Data(:,k))));
    index = ~cellfun(@isempty,Data(:,k));
    Cache_Matrix{1,k} = Data(index,k);
end

ExperimentMatrix = cell(max(levels),NumberOfFactors);
for m = 1:NumberOfFactors
    ExperimentMatrix(1:levels(1,m),m) = Cache_Matrix{1,m};
end

levels_cell = num2cell(levels);

% Clear the factors to be retained from the DOE input.
for n = 1:NumberOfFactors
    if retain(n) == 1
        levels_cell(n) = [];
    end
end

levels_reduced = cell2mat(levels_cell);

% Start the DOE algorithm.
[ChosenDesign,~] = calculateOptimalReducedDesign(levels_reduced);

% Add the factors to be retained back to the design matrix.
repdesign = ChosenDesign;
indexRetain = find(retain);
for q = 1:sum(retain)
    repdesign = repmat(repdesign,levels(indexRetain(q)),1);
    newColumn = sort(repmat((1:levels(indexRetain(q)))',size(repdesign,1)/levels(indexRetain(q)),1));
    repdesign(:,size(repdesign,2)+1) = newColumn;
end

sortedMat = sortrows(repdesign,NumberOfFactors:-1:1);
FinalDesign = sortedMat;
for a = 1:size(levels,2)
    for b = levels(a):-1:1
        idx = FinalDesign(:,a) == b;
        FinalDesign(idx,a) = ExperimentMatrix{b,a};
    end
end

% Delete the buffer file.
if isfile('Cache.fig')
    delete Cache.fig
end

end

function z = InputQuery()
% User Interface

% Initial query
prompt = {'Enter the number of levels of the factor with the most levels:'};
dlgtitle = 'Maximum Level Number';
dims = [1 30];
maxLevelNumber = inputdlg(prompt,dlgtitle,dims);
NumberOfLevels = str2double(maxLevelNumber{1});

prompt = {'Enter the number of factors:'};
dlgtitle = 'Input Factor Number';
dims = [1 30];
answer = inputdlg(prompt,dlgtitle,dims);
NumberOfFactors = str2double(answer{1});
prompt = cell(1,NumberOfFactors);
for k = 1:NumberOfFactors
    prompt{k} = sprintf('Enter the %d. variable type',k);
end
dlgtitle = 'Input Factor Names';
dims = [1 35];
answerOK = zeros(NumberOfFactors,5);

ScreenSize = get(0,'ScreenSize');
width = 400;
height = 161;
startX = ScreenSize(3)/2 - width/2;
startY = ScreenSize(4)/2 - height/2;

answer2 = inputdlg(prompt,dlgtitle,dims);
% Check if the input is correct
% while true
%     answer2 = inputdlg(prompt,dlgtitle,dims);  
%     for k = 1:NumberOfFactors
%         answerOK(k,1) = count(answer2{k},'amount','IgnoreCase',true);
%         answerOK(k,2) = count(answer2{k},'species','IgnoreCase',true);
%         answerOK(k,3) = count(answer2{k},'temperature','IgnoreCase',true);
%         answerOK(k,4) = count(answer2{k},'time','IgnoreCase',true);
%         answerOK(k,5) = sum(answerOK(k,1:4));
%     end
%     if sum(answerOK(1:NumberOfFactors,5)) ~= NumberOfFactors
%         f = uifigure('Position',[startX   startY   width   height]);
%         uialert(f,['Please name the factors correctly (each name must'...
%             ' contain exactly one of the following words: ''Temperature'''...
%             ',''Time'',''Amount'',''Species'''],...
%             'Input Error','CloseFcn','uiresume(f)');
%         uiwait(f);
%         try
%             close(f);
%         catch
%         end
%     else
%         break
%     end
% end

ScreenSize = get(0,'ScreenSize');

rowNames = cell(1,NumberOfLevels);
for a = 1:NumberOfLevels
    rowNames{a} = sprintf('Level %d',a);
end

f = figure('visible','off');
g = uicontrol(f,'Style', 'text');
set(g,'units','pixel');
ExtentMatrix_col = zeros(NumberOfFactors,1);
ExtentMatrix_row = zeros(NumberOfLevels,1);
for k = 1:NumberOfFactors
    set(g,'String',answer2{k});
    Extent = get(g,'Extent');
    ExtentMatrix_col(k) = Extent(3);
end
for m = 1:NumberOfLevels
    set(g,'String',rowNames{m});
    Extent = get(g,'Extent');
    ExtentMatrix_row(m) = Extent(3);
end
Extent_row = max(ExtentMatrix_row);
close

ExtentMatrix_adj_col = num2cell(1.1962*ExtentMatrix_col + 16.325)';
sumExtent_col = sum([ExtentMatrix_adj_col{:}]);

BorderSize = 20;
ButtonHeight = 25;

if NumberOfLevels*18+4*BorderSize+ButtonHeight < ScreenSize(4)/2
    extraSpace = 5/6*Extent_row+7.5;
else
    extraSpace = 5/6*Extent_row+24.5;
end

widthFigure = extraSpace + Extent_row + sumExtent_col + 2*BorderSize;
heightFigure = min(NumberOfLevels*18+4*BorderSize+ButtonHeight,ScreenSize(4)/2);
startXFigure = ScreenSize(3)/2 - widthFigure/2;
startYFigure = ScreenSize(4)/2 - heightFigure/2;

heightTable = heightFigure - 83;
widthTable = (extraSpace + Extent_row + sumExtent_col);
startXTable = BorderSize;
startYTable = heightFigure-heightTable-BorderSize;

close
z = figure('Position',[startXFigure startYFigure widthFigure heightFigure],'Name', 'Factor Levels');
columnnames = answer2';

dat = cell(NumberOfLevels,NumberOfFactors);

t = uitable('Parent',z,'Data',dat,'ColumnName',columnnames,... 
            'RowName',rowNames,'Position',[startXTable startYTable...
            widthTable heightTable],...
            'ColumnWidth', ExtentMatrix_adj_col);
        
set(t,'ColumnEditable',true(1,NumberOfFactors));
t.ColumnFormat = repmat({'numeric'},1,NumberOfFactors);
set(z, 'MenuBar', 'none');
ButtonWidth = 60;
startXButton = widthFigure/2 - ButtonWidth/2;
startYButton = 20;

uicontrol('Parent',z,'Style','text','String',{'Retain:'},'Position',...
    [BorderSize+20 heightFigure-BorderSize+4 45 15],'FontSize',10);

sizeBox = 14;
yCheckBox = heightFigure-BorderSize+4;
xCheckBox = zeros(NumberOfFactors,1);
xExtra = 0;
for p = 1:NumberOfFactors
    xCheckBox(p,1) = BorderSize + extraSpace + Extent_row + ExtentMatrix_adj_col{p}/2 - sizeBox/2 + xExtra;
    xExtra = xExtra + ExtentMatrix_adj_col{p};
    uicontrol('Parent',z,'Style','checkbox',...
    'Position',[xCheckBox(p,1) yCheckBox sizeBox sizeBox],'Visible','on');
end

uicontrol('Parent',z,'Style','pushbutton','String','Save',...
    'Position',[startXButton startYButton ButtonWidth ButtonHeight],'Visible',...
    'on','Callback',@SaveButtonPushed);

function SaveButtonPushed(src,event)
    saveas(z,'Cache.fig')
    close
end

end