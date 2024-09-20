function [ChosenDesign,CandidateAssessment] = calculateOptimalReducedDesign(levels,varargin)
%CALCULATEOPTIMALREDUCEDDESIGN constructs model-variance-optimized experimental designs.
%   CHOSENDESIGN, CANDIDATEASSESSMENT = CALCULATEOPTIMALREDUCEDDESIGN(LEVELS)
%   uses up to four different variance criteria in combination to construct
%   optimal experimental designs. The following criteria are supported:
%
%   (1) A-optimality
%   (2) D-optimality
%   (3) E-optimality
%   (4) G-optimality
%
%   LEVELS is a 1-by-N vector containing the factors and levels. Each
%   number stands for a different factor, while the magnitude of the number
%   represents the number of levels the particular factor. 
%   Example: LEVELS = [2 2 3 7 9].
%   
%   CHOSENDESIGN, CANDIDATEASSESSMENT = CALCULATEOPTIMALREDUCEDDESIGN(LEVELS,...,'PARAM1',VALUE1,'PARAM2',VALUE2,...)
%   provides more control over the model generation through a set of
%   parameter/value pairs (case insensitive).
%   Valid parameters are the following:
%
%      Parameter     Value
%
%      'OptiCond'    Determines which optimality criteria are used. Must be
%                    passed as a string, seperated by '+'. Example:
%                    'A+D+G' or 'E' (default = 'A+D+E+G').
%
%      'startValue'  Specifies the starting trial number for the algorithm.
%                    (default = sum(levels))
%
%      'excludefun'  Handle to a function that excludes undesirable runs. 
%                    If the function is f, it must support the syntax b = f(S),
%                    where S is a matrix of treatments with nfactors columns 
%                    and b is a vector of Boolean values with the same number 
%                    of rows as S. b(i) is true if the ith row S should be 
%                    excluded (default = []).
%
%      'Options'     A structure that contains options specifying whether to
%                    compute multiple tries in parallel, and specifying how
%                    to use random numbers when generating the starting points
%                    for the tries. This argument can be created by a call to 
%                    STATSET. CALCULATEOPTIMALREDUCEDDESIGN uses the following 
%                    fields:
%                        'UseParallel'
%                        'UseSubstreams'
%                        'Streams'
%                    For information on these fields see PARALLELSTATS.
%                    (default = FALSE).
%                    NOTE: If 'UseParallel' is TRUE and 'UseSubstreams' is 
%                    FALSE, then the length of 'Streams' must equal the 
%                    number of workers used by CALCULATEOPTIMALREDUCEDDESIGN.
%                    If a parallel pool is already open, this will be the 
%                    size of the parallel pool. If a parallel pool is not
%                    already open, then MATLAB may try to open a pool for 
%                    you (depending on your installation and preferences).
%                    To ensure more predictable results, it is best 
%                    to use the PARPOOL command and explicitly create a 
%                    parallel pool prior to invoking ROWEXCH with 
%                    'UseParallel' set to TRUE.
%
%   Example:
%      [ChosenDesign,CandidateAssessment] = calculateOptimalReducedDesign(...
%          [2 2 3 7 9],'OptiCond','A+D+G','startValue',24,'Options',...
%          statset('UseParallel',true))
%
%   NOTE: All local maxima are found at divisors or integer fractions of
%   the full-factorial number of trials, thus the algorithm only
%   investigates only these trial numbers. The Design matrices are
%   generated using standardized orthogonal contrast coding.
%
%   See also ROWEXCH, CANDEXCH
%
%   Copyright © 2020 by Tobias Averbeck, tobiasav92@gmail.com
% ====================================================================

% set defaults
[pnames, design_options] = set_default(levels);

if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
end

% process optional input
idxOpt = cell(length(pnames),2); 
idxOpt(1:length(pnames),1) = pnames';
for a = 1:length(pnames)
    if ~isempty(find(cellfun(@(x) strcmp(x,pnames{a}), varargin), 1))
        idxOpt{a,2} = varargin{1,find(cellfun(@(x) strcmp(x,pnames{a}), varargin), 1)+1};
    end
end

design_options(~cellfun(@isempty,idxOpt(:,2)),2)...
    = idxOpt(~cellfun(@isempty,idxOpt(:,2)),2);
design_options(~cellfun(@isempty,idxOpt(:,2)),3) = {1};

NumberOfFactors = size(levels,2);

% Generate the Design- and Fisher-Information matrix, 
% using orthogonal contrast coding.
[D_full,~,X_orth_full,~] = calcX(levels);

% Disable warning for bad starting design.
warning('off', 'stats:candexch:BadStartingDesign')
if ~isempty(idxOpt{1,2})
    if ~isequal(idxOpt{1,2},'on') && ~isequal(idxOpt{1,2},'off')
        m = message('stats:doptargchk:BadDisplay');
        throwAsCaller(MException(m.Identifier,'%s',getString(m)));
    end
end

% Check exclusion function, if any.
if ~isempty(idxOpt{6,2})
   if ~isa(idxOpt{6,2},'function_handle') && ...
      ~(ischar(idxOpt{6,2}) && exist(idxOpt{6,2},'file')>0)
      m = message('stats:doptargchk:BadExcludeFun');
      throwAsCaller(MException(m.Identifier,'%s',getString(m)));
   end
end

% Check levels array.
if ~isempty(levels)
   if ~isscalar(levels) && ~isvector(levels)
      m = message('stats:doptargchk:BadLevelsSize');
      throwAsCaller(MException(m.Identifier,'%s',getString(m)));
   elseif ~isscalar(levels) && length(levels)~=NumberOfFactors
      m = message('stats:doptargchk:BadLevelsLength');
      throwAsCaller(MException(m.Identifier,'%s',getString(m)));
   elseif any(levels<2 | levels~=round(levels))
      m = message('stats:doptargchk:BadLevelsValues');
      throwAsCaller(MException(m.Identifier,'%s',getString(m)));
   end
end

% Delete buffer file if existing.
if isfile('selection.txt')
    delete selection.txt
end

% Check 'OptiCond' argument.
idxStr = cell(4,2);
idxStr(:,2) = {0};
idxStr(:,1) = {'A' 'D' 'E' 'G'};
OptiC = idxStr(:,1);
for b = 1:length(OptiC)
    if ~isempty(idxOpt{3,2})
        idxStr{b,2} = logical(strfind(idxOpt{3,2},OptiC{b}));
    else
        idxStr{b,2} = logical(strfind(design_options{3,2},OptiC{b}));
    end
end
if sum(find(cellfun(@isempty,idxStr(:,2))))
    idxEmptyConds = find(cellfun(@isempty,idxStr(:,2)));
    for k = 1:length(idxEmptyConds)
        idxStr{idxEmptyConds(k),2} = 0;
    end
end

CreateStruct.Interpreter = 'tex';
CreateStruct.WindowStyle = 'modal';

% Initial information window.
createMessage_1(CreateStruct)

% Query for the algorithm parameter.
prompt = {'Maximum computation time (in hours):',...
    'Upper bound for number of trials:',...
    'Minimum number of calculations per trial:'};
        dlgtitle = 'Calculation parameters';
dims = [1 30];
defInput = {'16',num2str(min(1000,prod(levels))),'10000'};
Parameter = inputdlg(prompt,dlgtitle,dims,defInput);

MaxCompTime = str2double(Parameter{1});
UpperTrialBound = str2double(Parameter{2});
MinNumberOfCalcsPerTrial = str2double(Parameter{3});

% Define size for the progress window.
ScreenSize = get(0,'ScreenSize');
widthFigure = ScreenSize(3)*0.6;
heightFigure = ScreenSize(4)*0.8;
startXFigure = ScreenSize(3)/2 - widthFigure/2;
startYFigure = ScreenSize(4)/2 - heightFigure/2;

% Calculate divisors and integer fractions of the full-factorial size.
TrialNumberOfFullDesign = size(D_full,1);
NumberOfDivisors = 1:ceil(sqrt(TrialNumberOfFullDesign));
PartialDivisors = NumberOfDivisors(rem(TrialNumberOfFullDesign,NumberOfDivisors)==0);
AllDivisors = [PartialDivisors sort(TrialNumberOfFullDesign./PartialDivisors)]';

ExtraDivisor = zeros(max(levels),size(levels,2));
for a = 1:size(levels,2)
    for b = 1:levels(a)
        ExtraDivisor(b,a) = size(D_full,1)*(b/levels(a));
    end
end
ExtraDivisor = unique(ExtraDivisor);
ExtraDivisor = ExtraDivisor(ExtraDivisor~=0);
AllDivisors = unique([AllDivisors;ExtraDivisor]);

% Use only divisor greater than the starting trial number.
RelevantDivisorsIdx = AllDivisors >= (design_options{4,2});
RelevantDivisors = AllDivisors(RelevantDivisorsIdx);
RelevantDivisors = RelevantDivisors(1:size(RelevantDivisors,1)-1,1);

unique_levels = unique(levels);

% If the number of unique factor sizes (= number of levels) is too small,
% add additional points to look for optima.
if length(unique_levels) < 4
    LinearCombinations = cell(length(unique_levels),length(unique_levels)+1);
    for a = 1:length(unique_levels)
        LinearCombinations{a,1} = unique_levels(a):unique_levels(a):UpperTrialBound;
    end

    for b = 1:length(unique_levels)
        for c = 1:length(unique_levels)
            LinearCombinations{b,c+1} = log(kron(exp(LinearCombinations{b,1}),exp(LinearCombinations{c,1})));
        end
    end

    UniqueLinearCombinations = unique([LinearCombinations{:}])';

    RelevantDivisors = [RelevantDivisors; UniqueLinearCombinations];
    RelevantDivisors = sort(RelevantDivisors);

    RelevantDivisorsIdx2 = RelevantDivisors >= (design_options{4,2});
    RelevantDivisors = RelevantDivisors(RelevantDivisorsIdx2);
    RelevantDivisors = unique(RelevantDivisors);
end

if isempty(RelevantDivisors)
    error('The full factorial matrix is too small, either increase the number of factors or the number of levels!');
end

Fitness = zeros(length(RelevantDivisors),MinNumberOfCalcsPerTrial);
maxFitness = zeros(length(RelevantDivisors),1);
Runs2 = length(RelevantDivisors);

% Progess window.
f = figure('units','pixels','Position',[startXFigure...
    startYFigure widthFigure heightFigure],'name',...
    'Visual presentation of the calculation progress',...
    'menubar','none','numbertitle','off','resize','off','visible','off');

axBar1 = axes('Units','pix','Position',...
    [widthFigure*0.175 heightFigure*0.18 widthFigure*0.65 30]);
set(axBar1,'Xtick',[],'Ytick',[],'Xlim',[0 1000]);
box on;

axBar2=axes('Units','pix','Position',...
    [widthFigure*0.175 heightFigure*0.0231 widthFigure*0.65 30]);
set(axBar2,'Xtick',[],'Ytick',[],'Xlim',[0 1000]);
box on;

ax = axes(f);
ax.Units = 'pixels';
ax.Position = [widthFigure*0.1 heightFigure*0.15+0.2894*heightFigure...
    widthFigure*0.8 heightFigure*0.5];

ylim([0 100]);
xlim([RelevantDivisors(1,1)-5 RelevantDivisors(size(RelevantDivisors,1),1)+5])
title('Design efficiency of calculated reduced designs')
xlabel('Number of Runs [-]') 
ylabel('Design efficiency [%]')
allFontSizes = findall(f,'-property','FontSize');
set(allFontSizes(1:3),'FontSize',20);

uicontrol('Parent',f,'Style','text','String',...
    {'Calculation progress for current design (Number of trials: 0)'},...
    'Position',[widthFigure/2-400 heightFigure*0.18+50 800 32],...
    'FontSize',20);

uicontrol('Parent',f,'Style','text','String',...
    {'Overall calculation progress'},...
    'Position',[widthFigure/2-400 heightFigure*0.0231+50 800 32],...
    'FontSize',20);

% Calculation of the required computation time.
% MinimumNumberOfTriesForTimeCalculation = 5;
% EstimatedTimeForTimeCalculationForOneTry = -1.621706 + 0.01411077*AllDivisors(end-1)...
%     - 0.000005801121*AllDivisors(end-1)^2 + 1.444878e-9*AllDivisors(end-1)^3; %empirical formula
% EstimatedTimeForTimeCalculationForOneTry = max(0,EstimatedTimeForTimeCalculationForOneTry);
% MaximumTimeForTimeCalculation = 600; %seconds
% NumberOfTriesAtMax = MaximumTimeForTimeCalculation/...
%     EstimatedTimeForTimeCalculationForOneTry;
% if isinf(NumberOfTriesAtMax)
%     NumberOfTriesAtMax = 1000;
% end
% NumberOfTriesForTimeCalculation = round(max(NumberOfTriesAtMax,...
%     MinimumNumberOfTriesForTimeCalculation));
% 
% elapsedTime = calculateRequiredTime(Runs2,idxStr,RelevantDivisors,...
%     NumberOfTriesForTimeCalculation,NumberOfFactors,maxFitness,...
%     D_full,X_orth_full,design_options,levels);
% 
% load elapsedTime
% 
% NeededTime = zeros(size(elapsedTime,1),3);
% NeededTime(:,1) = RelevantDivisors;
% NeededTime(:,2) = elapsedTime;
% 
% NeededTime(:,3) = cumsum(elapsedTime.*MinNumberOfCalcsPerTrial)/3600;
% RoundToDecimalPlace = 3;
% 
% if MaxCompTime > 12
%     BufferTime = 1;
% else
%     BufferTime = 0;
% end

% Calculate if the entered values are feasible.
% newUpperTrialBound = 0;
% newUpperTrialBound2 = 0;
% while 1
%     NewNeededTime = NeededTime(:,1:3);
%     while 1
%         x = 1;
%         if newUpperTrialBound
%             NewNeededTime(:,3) = cumsum(elapsedTime(size(NewNeededTime(:,3),1),1).*x)/3600;
%             NewNeededTime = NewNeededTime(1:find(NewNeededTime(:,1) >= UpperTrialBound,1),:);
%             newUpperTrialBound = 0;
%         end
%         while 1
%             disp(MinNumberOfCalcsPerTrial);
%             NewNeededTime(:,3) = cumsum(NewNeededTime(:,2).*x)/3600;
%             if NewNeededTime(end,3) < MaxCompTime-BufferTime
%                 x = x+1;
%             else
%                 break
%             end
%         end
%         if x < MinNumberOfCalcsPerTrial 
%             NewNeededTime = NewNeededTime(1:end-1,:);
%             if isempty(NewNeededTime)
%                 NewNeededTime = NeededTime;
%                 if newUpperTrialBound2
%                     NewNeededTime = NewNeededTime(1:find(NewNeededTime >= UpperTrialBound,1),1:3);
%                 end
%                 
%                 Sum_sqr = (NewNeededTime(:,1) - UpperTrialBound).^2;
%                 [~,idx_min] = min(Sum_sqr);
%                 
%                 createWarn_1(x,RoundToDecimalPlace,NewNeededTime,idx_min,CreateStruct);
%                 
%                 prompt = {'Maximum computation time (in hours):',...
%                     'Upper bound for number of trials:',...
%                     'Minimum number of calculations per trial:'};
%                         dlgtitle = 'Calculation parameters';
%                 dims = [1 30];
%                 defInput = {sprintf(['%.' sprintf('%d',RoundToDecimalPlace) 'f'],ceil((1.1*NewNeededTime(idx_min,3))*...
%                     10^RoundToDecimalPlace)/10^RoundToDecimalPlace),num2str(UpperTrialBound),'10000'};
%                 Parameter2 = inputdlg(prompt,dlgtitle,dims,defInput);
% 
%                 MaxCompTime = str2double(Parameter2{1});
%                 UpperTrialBound = str2double(Parameter2{2});
%                 if UpperTrialBound > prod(levels)
%                     UpperTrialBound = prod(levels);
%                 end
%                 MinNumberOfCalcsPerTrial = str2double(Parameter2{3});
%                 if str2double(Parameter2{2}) ~= str2double(Parameter{2})
%                     newUpperTrialBound = 1;
%                     newUpperTrialBound2 = 1;
%                 end
%             end
%         else
%             Sum_sqr = (NewNeededTime(:,1) - UpperTrialBound).^2;
%                 [~,idx_min] = min(Sum_sqr);
%             if NewNeededTime(idx_min,3) < NewNeededTime(end,3)
%                 MinNumberOfCalcsPerTrial = MinNumberOfCalcsPerTrial + 100;
%             else
%                 break
%             end
%         end
%     end
%     if NewNeededTime(end,1) < UpperTrialBound
%         TestNeededTime = NeededTime(:,1:2);
%         test_x = 1;
%         while 1
%             TestNeededTime(:,3) = cumsum(TestNeededTime(:,2).*test_x)/3600;
%             if TestNeededTime(end,3) < MaxCompTime-BufferTime
%                 test_x = test_x+1;
%             else
%                 break
%             end
%         end
%         
%         if x >= MinNumberOfCalcsPerTrial
%             if ~isempty(NewNeededTime(find(NewNeededTime(:,1) >= UpperTrialBound,1),3)+1 ...
%                 <= MaxCompTime+0.05) 
%                 if NewNeededTime(find(NewNeededTime(:,1) >= UpperTrialBound,1),3)+1 ...
%                     <= MaxCompTime+0.05 || NewNeededTime(end,3) <= MaxCompTime+0.05
%                     break
%                 end
%             elseif NeededTime(end,3) <= MaxCompTime+0.05
%                 NewNeededTime(:,3) = cumsum(elapsedTime(size(NewNeededTime(:,3),1),1).*x)/3600;
%                 break
%             end
%         end
%         
%         time = NeededTime(find(NeededTime(:,1) >= UpperTrialBound,1),3)+1;
%         if NeededTime(end,1) < UpperTrialBound && NeededTime(end,3) > MaxCompTime
%             time = NeededTime(end,3)+BufferTime;
%         end      
%         
%         mh = createWarn_2(MaxCompTime,NewNeededTime,MinNumberOfCalcsPerTrial,...
%             test_x,UpperTrialBound,time,CreateStruct);
%         
%         th = findall(mh, 'Type', 'Text');
%         deltaWidth = sum(th.Extent([1,3]))-mh.Position(3) + th.Extent(1);
%         deltaHeight = sum(th.Extent([2,4]))-mh.Position(4) + 10;
%         mh.Position([3,4]) = mh.Position([3,4]) + [deltaWidth, deltaHeight];        
%         mh.Resize = 'on';  
%         
%         uiwait(mh);
%         
%         prompt = {'Minimum number of calculations per trial:','Upper bound for number of trials:'...
%             'Maximum computation time (in hours):'};
%         dlgtitle = 'Calculation parameters';
%         dims = [1 30];
%         defInput = {sprintf('%d',MinNumberOfCalcsPerTrial),sprintf('%d',UpperTrialBound)...
%             sprintf('%.2f',time)};
%         newParameter = inputdlg(prompt,dlgtitle,dims,defInput);
%         MinNumberOfCalcsPerTrial = str2double(newParameter{1});
%         UpperTrialBound = str2double(newParameter{2});
%         MaxCompTime = str2double(newParameter{3});
%         
%         NeededTime(:,3) = cumsum(elapsedTime.*MinNumberOfCalcsPerTrial)/3600;
%     else
%         break
%     end
% end
% 
% if abs(x - MinNumberOfCalcsPerTrial) <= 1
%     x = MinNumberOfCalcsPerTrial;
% end
% 
% message2 = ['\rmNumber of calculations per Trial: \bf' sprintf('%d',x)]; 
% if NewNeededTime(find(NewNeededTime(:,1) >= UpperTrialBound,1),3) <= MaxCompTime+0.05
%     time = NewNeededTime(find(NewNeededTime(:,1) >= UpperTrialBound,1),3)+BufferTime;
%     trials = NewNeededTime(find(NewNeededTime(:,1) >= UpperTrialBound,1),1);
% elseif NewNeededTime(end,3) <= MaxCompTime+0.05
%     time = NewNeededTime(end,3)+BufferTime;
%     trials = NewNeededTime(end,1);
% end  
% message1 = ['\fontsize{12}Upper bound for number of trials: \bf'...
%     sprintf('%d',trials)];
% message3 = ['\rmComputation time: \bf' sprintf(['%.' sprintf('%d',...
%     RoundToDecimalPlace) 'f'],round(time,3)) ' hours\rm'];
% j = msgbox({message1;message2;message3},'Calculation Information','help',CreateStruct);
% uiwait(j);

% numberOfTries = x;
% trials = RelevantDivisors;

% RelevantDivisors = RelevantDivisors(1:find(RelevantDivisors == trials));

numberOfTries = MinNumberOfCalcsPerTrial;
[~,idx]=min(abs(RelevantDivisors-UpperTrialBound));
RelevantDivisors = RelevantDivisors(1:idx);
reducedDesigns = cell(length(RelevantDivisors(:,1)),numberOfTries);

NumberOfDivisorsTries = 100;
TryDivisors = zeros(NumberOfDivisorsTries,1);
for b = 1:NumberOfDivisorsTries
    TryDivisors(b) = round(b/NumberOfDivisorsTries*numberOfTries);
end
Runs2 = zeros(1,1);
Fitness = zeros(length(RelevantDivisors),MinNumberOfCalcsPerTrial);
maxFitness = zeros(length(RelevantDivisors),1);

%--------------------------------------------
%---------------- Loop Start ----------------
%--------------------------------------------
% The optimality criteria are calculated for each new design.
for k = 1:length(RelevantDivisors)
    
    uicontrol('Parent',f,'Style','text','String',...
    {sprintf('Calculation progress for current design (Number of trials: %d)',...
    RelevantDivisors(k))},'Position',[widthFigure/2-400 heightFigure*0.18+50 ...
    800 32],'FontSize',20);

    for m = 1:numberOfTries
        if ismember(m,TryDivisors) || m == 1
            axes(axBar2)
            cla
            rectangle('Position',[0,0,(round(1000*(numberOfTries*(k-1)+m)...
                /(length(RelevantDivisors)*numberOfTries)))+1,20],'FaceColor','g'); 
            text(480,10,[num2str(round(100*(numberOfTries*(k-1)+m)...
                /(length(RelevantDivisors)*numberOfTries))),'%']);

            axes(axBar1)
            cla
            rectangle('Position',[0,0,(round(1000*m/numberOfTries)),20],'FaceColor','g');
            text(482,10,[num2str(round(100*m/numberOfTries)),'%']);
        end
        
        fprintf('Divisor: %d - Try: %d \n',RelevantDivisors(k),m);
        
        try
            [reducedDesigns{k,m},~,X_orth_red_m] = rowexchVarInput(RelevantDivisors(k),NumberOfFactors,...
                D_full,X_orth_full,levels,design_options);  
        catch
            X_orth_red_m = 0;
        end
        % A-efficiency
        if idxStr{1,2} 
            try
                Runs2(k,1) = real(100/(size(X_orth_red_m,1)*...
                    trace((X_orth_red_m'*X_orth_red_m)\...
                    eye(size(X_orth_red_m,2)))/size(X_orth_red_m,2)));
                if isnan(Runs2(k,1))
                    Runs2(k,1) = 0;
                end
            catch
                Runs2(k,1) = 0;
            end
        end

        % D-efficiency
        if idxStr{2,2}
            try
                Runs2(k,2) = real(100/(size(X_orth_red_m,1)*...
                    det((X_orth_red_m'*X_orth_red_m)\eye(size(X_orth_red_m,2)))^...
                    (1/size(X_orth_red_m,2))));
                if isnan(Runs2(k,2))
                    Runs2(k,2) = 0;
                end
            catch
                Runs2(k,2) = 0;
            end
        end

        % E-efficiency
        if idxStr{3,2}
            try
                Runs2(k,3) = real(100/(size(X_orth_red_m,1)*...
                    eigs((X_orth_red_m'*X_orth_red_m)\eye(size(X_orth_red_m,2)),1)));
                if isnan(Runs2(k,3))
                    Runs2(k,3) = 0;
                end
            catch
                Runs2(k,3) = 0;
            end
        end

        % G-efficiency    
        if idxStr{4,2}
            try
                Runs2(k,4) = real(100*sqrt(size(X_orth_red_m,2)/(size(X_orth_red_m,1)*...
                    max(diag(X_orth_red_m*((X_orth_red_m'*X_orth_red_m)\...
                    eye(size(X_orth_red_m,2)))*X_orth_red_m')))));
                if isnan(Runs2(k,4))
                    Runs2(k,4) = 0;
                end
            catch
                Runs2(k,4) = 0;
            end
        end

        % Calculate Fitness as Average of all criteria.
        Fitness(k,m) = (idxStr{1,2}*Runs2(k,1) + idxStr{2,2}*Runs2(k,2) + idxStr{3,2}...
            *Runs2(k,3) + idxStr{4,2}*Runs2(k,4))/sum([idxStr{:,2}]);

    end
% Store the best fitness.
maxFitness(k) = max(Fitness(k,:));

% Update the plot.
axes(ax)
plot(RelevantDivisors(1:k),maxFitness(1:k),'Color','b','linestyle',...
    'none','marker','x','MarkerSize',20,'LineWidth',2);
title('Design efficiency of calculated reduced designs')
ylim([maxFitness(1,1)-5 100]);
xlim([RelevantDivisors(1,1)-5 RelevantDivisors(size(RelevantDivisors,1),1)+5])
% xlim([0 540])
grid on
xlabel('Number of Runs [-]') 
ylabel('Design efficiency [%]') 
allFontSizes = findall(f,'-property','FontSize');
set(allFontSizes(1:end-2),'FontSize',20);
drawnow
end
close all

% Assess all candidates.
CandidateAssessment = zeros(size(RelevantDivisors,1),3);
CandidateAssessment(:,1:2) = [RelevantDivisors maxFitness];
minFitnessFirstSurvey = min(maxFitness);
CandidateAssessment(:,3) = ((1-RelevantDivisors/size(D_full,1))*100 +...
    (maxFitness - minFitnessFirstSurvey)/(100 - minFitnessFirstSurvey)*100)/2;

% Show a table with the results.
z = newInput(CandidateAssessment);
uiwait(z);

% Get data from buffer file.
fileID = fopen('selection.txt','r');
DesignNumber = fscanf(fileID,'%s');
fclose(fileID);
delete selection.txt

% Store chosen design for output.
index = find(~isletter(DesignNumber));
DesignNumberNumerical = str2double(DesignNumber(index(1):index(end)));

[~,q] = find(Fitness(DesignNumberNumerical,:) == max(Fitness(DesignNumberNumerical,:)));

ChosenDesign = reducedDesigns{DesignNumberNumerical,q};

end

function createMessage_1(CreateStruct)

msg = ['\fontsize{12}\bfFor the following calculations, 3 parameters '...
    'have to be defined, which are:\rm' newline newline '\bfMaximum '...
    'computation time\rm:' newline 'The upper time limit until all '...
    'calculations are finished. \bf16 hours\rm would be '...
    'recommended, so the algorithm can run overnight and is finished '...
    'in the morning.' newline newline '\bfUpper bound for the trial '...
    'number\rm:' newline 'This is the maximum number of trials which '...
    'are considered for the calculations. This number depends on personal '...
    'preference and depends on the duration of each trial. If one trial '...
    'of the whole experiments takes you 2 hours to finish, you should '...
    'set a smaller number than if it would take just 2 minutes. Also '...
    'keep in mind, that the computation time increases exponentially '...
    'with increasing trial numbers. With regard to the computation time '...
    'a value around \bf1000\rm or the product of the all factor levels '...
    'would be recommended, whatever is lower. Example factor levels: '...
    '[2 3 4 5]\rightarrow 2\cdot3\cdot4\cdot5 = \bf120\rm. ' newline newline...
    '\bfNumber of Calculations '...
    'per trial number\rm:' newline 'Due to the extremely high number of '...
    'possible combinations the possibility to find THE best solution '...
    'in a timely manner is next to zero. Because of this, the algorithm '...
    'iteratively searches for near-optimal solutions and by increasing '...
    'this number, the probability to find near-optimal solutions '...
    'increases linearly, as does the computation time. '...
    'Recommended would be a value around \bf10,000\rm.' newline...
    newline '\bfPlease define these parameters in the following '...
    'window.\rm' newline newline];

b = msgbox(msg,'Calculation Information','none',CreateStruct);
uiwait(b);

end

function createWarn_1(x,RoundToDecimalPlace,NewNeededTime,idx_min,CreateStruct)

TooSmallWarning = warndlg(['\fontsize{12}The \itminimum calculation number per trial\rm is '...
    'never reached while staying below the \itmaximum computation time\rm. '...
    'Please pick a lower value for the \itminimum calculation number per '...
    'trial\rm or a higher value for the \itmaximum computation time\rm!' newline...
    newline 'Maximum value for the \itminimum calculation number per '...
    'trial\rm while retaining the same \itmaximum computation time\rm: \bf'...
    sprintf('%d',x) '.' newline 'Using this value is not recommended, as it '...
    'reduces the number of design candidates to 1. Increasing the '...
    '\itmaximum computation time\rm \bfis recommended, as the number of '...
    'design candidates are retained:' newline newline '\rmMinimum '...
    'value for the \itmaximum computation time\rm to retain both the same '...
    '\itminimum calculation number per trial\rm and \itnumber of design candidates\rm: \bf'...
    sprintf(['%.' sprintf('%d',RoundToDecimalPlace) 'f'],ceil((1.1*NewNeededTime(idx_min,3))*...
    10^RoundToDecimalPlace)/10^RoundToDecimalPlace) '.\rm'],...
    'Parameter Warning',CreateStruct);
uiwait(TooSmallWarning);
                
end

function mh = createWarn_2(MaxCompTime,NewNeededTime,MinNumberOfCalcsPerTrial,test_x,UpperTrialBound,time,CreateStruct)

mh = warndlg( ['\fontsize{12}\fontname{Calibri}The set value for the '...
    '\itmaximum computation time\rm \bf('  sprintf('%.2f',MaxCompTime) ...
    ' hours)\rm, admits only an \itupper bound for the trial number\rm'...
    ' of \bf' sprintf('%d',NewNeededTime(end,1)) '\rm with a \itminimum'...
    ' calculation number per trial\rm of \bf' sprintf('%d',MinNumberOfCalcsPerTrial)...
    '\rm OR a \itcalculation number per trial\rm of \bf' sprintf('%d',test_x)...
    '\rm with a \itminimum upper bound for the trial number\rm of \bf'...
    sprintf('%d',UpperTrialBound) '\rm.' newline newline 'With a set \itminimum'...
    ' calculation number per trial\rm of \bf' sprintf('%d',MinNumberOfCalcsPerTrial)...
    '\rm and a \itset upper bound for the trial number\rm of \bf'...
    sprintf('%d',UpperTrialBound) '\rm, the \itrequired computation time\rm'...
    ' would be \bf' sprintf('%.2f',time) ' hours\rm. If that is too long, please adjust'...
    ' one or more parameters in the following window.'],'Parameter Warning',...
    CreateStruct);

end

% function [pnames, design_options] = set_default(levels)
% 
% pnames = {'display' 'Options' 'OptiCond' 'weightOptiCond' 'startValue' 'excludefun' 'model' 'safetyDistance' 'weightJ2' 'Plot'};
% changed = num2cell(zeros(length(pnames),1));
% design_options = cell(length(pnames),3);
% design_options(1:length(pnames),1) = pnames';
% design_options(1:length(pnames),3) = changed;
% design_options{1,2} = 'off';
% design_options{2,2} = [];
% % design_options{2,2} = statset('UseParallel',true);
% design_options{3,2} = 'A+D+E+G';
% design_options{4,2} = 0.5;
% design_options{5,2} = sum(levels);
% design_options{6,2} = [];
% design_options{7,2} = 'linear';
% design_options{8,2} = 0;
% design_options{9,2} = [];
% design_options{10,2} = 'combined';
% 
% end

function [pnames, design_options] = set_default(levels)
% set default values.

pnames = {'display' 'Options' 'OptiCond' 'startValue' 'excludefun' 'model'};
changed = num2cell(zeros(length(pnames),1));
design_options = cell(length(pnames),3);
design_options(1:length(pnames),1) = pnames';
design_options(1:length(pnames),3) = changed;
design_options{1,2} = 'off';
design_options{2,2} = [];
design_options{3,2} = 'A+D+E+G';
design_options{4,2} = sum(levels);
design_options{5,2} = [];
design_options{6,2} = 'linear';

end

function elapsedTime = calculateRequiredTime(Runs2,idxStr,RelevantDivisors,numberOfTries,NumberOfFactors,maxFitness,D_full,X_orth_full,design_options,levels)
% Perform a test run to determine required computation time.
Runs2 = zeros(length(RelevantDivisors),NumberOfFactors);
ScreenSize = get(0,'ScreenSize');
widthFigure = 400;
heightFigure = 100;
startXFigure = ScreenSize(3)/2 - widthFigure/2;
startYFigure = ScreenSize(4)/2 - heightFigure/2;

f = uifigure('Position',[startXFigure startYFigure widthFigure heightFigure]);
d = uiprogressdlg(f,'Title','Calculate computation time',...
    'Message','This usually takes a maximum of 10 minutes.');

reducedDesigns = cell(length(RelevantDivisors(:,1)),numberOfTries);
elapsedTime = zeros(size(RelevantDivisors,1),1);
timerVal = zeros(size(RelevantDivisors,1),1,'uint64');

for k = 1:length(RelevantDivisors)
    timerVal(k) = tic;
    
    for m = 1:numberOfTries
        d.Value = (numberOfTries*(k-1)+m)/(length(RelevantDivisors)*numberOfTries); 
        drawnow
        try
            [reducedDesigns{k,m},~,X_orth_red_m] = rowexchVarInput(RelevantDivisors(k),NumberOfFactors,...
                D_full,X_orth_full,levels,design_options);  
        catch
            X_orth_red_m = 0;
        end
        if idxStr{1,2}
            try
                Runs2(k,1) = real(100/(size(X_orth_red_m,1)*...
                    trace((X_orth_red_m'*X_orth_red_m)\...
                    eye(size(X_orth_red_m,2)))/size(X_orth_red_m,2)));
                if isnan(Runs2(k,1))
                    Runs2(k,1) = 0;
                end
            catch
                Runs2(k,1) = 0;
            end
        end


        if idxStr{2,2}
            try
                Runs2(k,2) = real(100/(size(X_orth_red_m,1)*...
                    det((X_orth_red_m'*X_orth_red_m)\eye(size(X_orth_red_m,2)))^...
                    (1/size(X_orth_red_m,2))));
                if isnan(Runs2(k,2))
                    Runs2(k,2) = 0;
                end
            catch
                Runs2(k,2) = 0;
            end
        end


        if idxStr{3,2}
            try
                Runs2(k,3) = real(100/(size(X_orth_red_m,1)*...
                    eigs((X_orth_red_m'*X_orth_red_m)\eye(size(X_orth_red_m,2)),1)));
                if isnan(Runs2(k,3))
                    Runs2(k,3) = 0;
                end
            catch
                Runs2(k,3) = 0;
            end
        end


        if idxStr{4,2}
            try
                Runs2(k,4) = real(100*sqrt(size(X_orth_red_m,2)/(size(X_orth_red_m,1)*...
                    max(diag(X_orth_red_m*((X_orth_red_m'*X_orth_red_m)\...
                    eye(size(X_orth_red_m,2)))*X_orth_red_m')))));
                if isnan(Runs2(k,4))
                    Runs2(k,4) = 0;
                end
            catch
                Runs2(k,4) = 0;
            end
        end

        Fitness(k,m) = (idxStr{1,2}*Runs2(k,1) + idxStr{2,2}*Runs2(k,2) + idxStr{3,2}...
            *Runs2(k,3) + idxStr{4,2}*Runs2(k,4))/sum([idxStr{:,2}]);

    end
    
maxFitness(k) = max(Fitness(k,:));
elapsedTime(k) = toc(timerVal(k))/numberOfTries;
end

close(f)
end

function [D_red,X_orth_red,X_orth_red_m] = rowexchVarInput(k,NumberOfFactors,D_full,X_orth_full,levels,design_option)
% Calculate fisher information matrix using the ROWEXCH function.

indexRand = randperm(prod(levels));
indexRand = indexRand(1:k);

if ~isempty(design_option{2,2}) && isempty(design_option{5,2})
    [D_red,~] = rowexch(NumberOfFactors,k,design_option{6,2},'categorical',...
        1:NumberOfFactors,'levels',levels,'display',design_option{1,2},...
        'Options',design_option{2,2},'init',D_full(indexRand,:));
elseif isempty(design_option{2,2}) && ~isempty(design_option{5,2})
    [D_red,~] = rowexch(NumberOfFactors,k,design_option{6,2},'categorical',...
        1:NumberOfFactors,'levels',levels,'display',design_option{1,2},...
        'excludefun',design_option{5,2},'init',D_full(indexRand,:));
elseif ~isempty(design_option{2,2}) && ~isempty(design_option{5,2})
    [D_red,~] = rowexch(NumberOfFactors,k,design_option{6,2},'categorical',...
        1:NumberOfFactors,'levels',levels,'display',design_option{1,2},...
        'excludefun',design_option{5,2},'Options',design_option{2,2},...
        'init',D_full(indexRand,:));
else
    [D_red,~] = rowexch(NumberOfFactors,k,design_option{6,2},'categorical',...
        1:NumberOfFactors,'levels',levels,'display',design_option{1,2},...
        'init',D_full(indexRand,:));
end

X_orth_red = X_orth_full(ismember(D_full,D_red,'rows'),:);

[~,I,~] = unique(D_red, 'rows', 'first');
ixDupRows = setdiff(1:size(D_red,1), I)';

X_orth_red_m = [X_orth_red; X_orth_full(ixDupRows,:)];

end

function [X,X_eff,X_orth,NumberOfTerms] = calcX(levels)
% Calculate full factorial fisher information matrix using 
% standardized orthogonal contrast coding.

NumberOfFactors = size(levels,2);

X = fullfact(levels);
X(:,NumberOfFactors+1) = randperm(size(X,1))';
X_tbl = array2table(X);

formula = cell(NumberOfFactors+2,1);
formula{1,1} = sprintf('X%d ~ ',NumberOfFactors+1);
formula{2,1} = ' + ';
formula(3:NumberOfFactors+2,1) = {''};
formulaString = formula{1,1};
for k = 1:NumberOfFactors
    formula{k+2,1} = sprintf('X%d',k);
    formulaString = [formulaString formula{k+2,1} formula{2,1}];
end
formulaString = formulaString(1:length(formulaString)-2);

varnames = X_tbl.Properties.VariableNames;
X_tbl = [X_tbl(:,~ismember(X_tbl.Properties.VariableNames,varnames(1:NumberOfFactors)))...
   varfun(@categorical,X_tbl,'inputvariables',varnames(1:NumberOfFactors))];
X_tbl = [X_tbl(:,2:NumberOfFactors+1), X_tbl(:,1)]; 
X_tbl.Properties.VariableNames = varnames;

% Calculate effects matrix
lme = fitlme(X_tbl,formulaString,...
    'FitMethod','ML','DummyVarCoding','effects');
X_eff = designMatrix(lme,'Fixed');

NumberOfTerms = size(X_eff,2);

X = X(:,1:NumberOfFactors);

% Standardized orthogonal contrast coding.
[X_orth,~] = gson(X_eff);
% X_orth = X_orth * fzero(@(x) sumsqr(X_orth*x) - size(X_orth,1)*size(X_orth,2),[0, 10^100]);
X_orth = X_orth * sqrt((size(X_orth,1)*size(X_orth,2))/sumsqr(X_orth));

end

function [Q, R] = gson(X)
% Gram-Schmidt orthonormalization

[d,n] = size(X);
m = min(d,n);
R = zeros(m,n);
Q = zeros(d,m);

for i = 1:m
    R(1:i-1,i) = Q(:,1:i-1)'*X(:,i);
    v = X(:,i)-Q(:,1:i-1)*R(1:i-1,i);
    R(i,i) = norm(v);
    Q(:,i) = v/R(i,i);
end

R(:,m+1:n) = Q'*X(:,m+1:n);

end

function z = newInput(CandidateAssessment)
% Final graphical representation using UITABLE.

ScreenSize = get(0,'ScreenSize');
NumberOfRows = size(CandidateAssessment,1);

columnNames = cell(3,1);
columnNames{1} = 'Number of Trials';
columnNames{2} = 'Calculated Fitness';
columnNames{3} = 'Objective Function';

f = figure('visible','off');
g = uicontrol(f,'Style', 'text');
set(g,'units','pixel');
ExtentMatrix = zeros(size(columnNames,1),1);
for k = 1:size(columnNames,1)
    set(g,'String',columnNames{k});
    Extent = get(g,'Extent');
    ExtentMatrix(k) = Extent(3);
end

close

ExtentMatrix_adj = num2cell(1.1962*ExtentMatrix + 16.325)';
sumExtent = sum([ExtentMatrix_adj{:}]);

widthFigure = 130 + sumExtent + double(logical(max(0,size(CandidateAssessment,1)-9)))*11;
heightFigure = 62+NumberOfRows*18+100;
startXFigure = ScreenSize(3)/2 - widthFigure/2;
startYFigure = ScreenSize(4)/2 - heightFigure/2;

z = figure('Position',[startXFigure startYFigure widthFigure heightFigure],'Name', 'Optimum Trial Number');

rowNames = cell(1,NumberOfRows);
for a = 1:NumberOfRows
    rowNames{a} = sprintf('Design %d',a);
end

cmp = [248,105,107;248,109,107;248,113,108;248,117,109;248,...
    121,110;249,125,110;249,129,111;249,133,112;249,...
    138,113;250,142,114;250,146,114;250,150,115;250,...
    154,116;250,158,117;251,162,118;251,166,118;251,...
    171,119;251,175,120;252,179,121;252,183,122;252,...
    187,122;252,191,123;252,195,124;253,199,125;253,...
    204,126;253,208,126;253,212,127;254,216,128;254,...
    220,129;254,224,130;254,228,130;254,232,131;253,...
    235,132;248,233,132;243,232,132;238,230,131;233,...
    229,131;228,228,131;223,226,131;218,225,130;213,...
    223,130;208,222,130;204,221,130;199,219,129;194,...
    218,129;189,216,129;184,215,128;179,213,128;174,...
    212,128;169,210,127;164,209,127;159,208,127;154,...
    206,127;149,205,126;144,203,126;139,202,126;134,...
    201,126;129,199,125;124,198,125;119,196,125;114,...
    195,124;109,193,124;104,192,124;99,190,123];

ColorsForTableO = zeros(size(CandidateAssessment,1),5);
ColorsForTableO(:,1) = CandidateAssessment(:,1);
ColorsForTableO(:,2) = round(((CandidateAssessment(:,3) - min(CandidateAssessment(:,3)))...
    /(max(CandidateAssessment(:,3))- min(CandidateAssessment(:,3))))*(size(cmp,1)-1)+1);

if size(ColorsForTableO,1) > 1
    ColorsForTableO(:,3:5) = cmp(ColorsForTableO(:,2),:);
else
    ColorsForTableO(:,3:5) = cmp(64,:);
end


ColorsForTableQ = zeros(size(CandidateAssessment,1),5);
ColorsForTableQ(:,1) = CandidateAssessment(:,1);
ColorsForTableQ(:,2) = round(((CandidateAssessment(:,2) - min(CandidateAssessment(:,2)))...
    /(max(CandidateAssessment(:,2))- min(CandidateAssessment(:,2))))*(size(cmp,1)-1)+1);

if size(ColorsForTableQ,1) > 1
    ColorsForTableQ(:,3:5) = cmp(ColorsForTableQ(:,2),:);
else
    ColorsForTableQ(:,3:5) = cmp(64,:);
end 

ColorsForTableT = zeros(size(CandidateAssessment,1),5);
ColorsForTableT(:,1) = CandidateAssessment(:,1);
ColorsForTableT(:,2) = (size(cmp,1)+1) - round(((CandidateAssessment(:,1) - min(CandidateAssessment(:,1)))...
    /(max(CandidateAssessment(:,1))- min(CandidateAssessment(:,1))))*(size(cmp,1)-1)+1);

if size(ColorsForTableT,1) > 1
    ColorsForTableT(:,3:5) = cmp(ColorsForTableT(:,2),:);
else
    ColorsForTableT(:,3:5) = cmp(64,:);
end 

ColorsForTable = cell(3,1);
ColorsForTable(3,1) = {ColorsForTableO};
ColorsForTable(2,1) = {ColorsForTableQ};
ColorsForTable(1,1) = {ColorsForTableT};

dat = num2cell(CandidateAssessment);
dat(:,1) = cellfun(@(x)sprintf('%d',x),dat(:,1), 'UniformOutput',false);
dat(:,2:end) = cellfun(@(x)sprintf('%.2f',x),dat(:,2:end), 'UniformOutput',false);
dat4 = cell(size(dat,1),size(dat,2));
for b = 1:3
    for d = 1:NumberOfRows
        hex = rgb2hex(ColorsForTable{b,1}(d,3:5));
        dat4(d,b) = strcat(sprintf('<html><table bgcolor=%s><tr align=center><td width=%d>',hex,...
            1.1962*ExtentMatrix(b) + 16.325), dat(d,b));
    end
end

t = uitable('Parent',z,'Data',dat4,'ColumnName',columnNames,... 
            'RowName',rowNames,'Position',[20 heightFigure-(22+18*NumberOfRows+18)...
            (widthFigure-40) 22+18*NumberOfRows],'ColumnWidth', ExtentMatrix_adj);
        
t.ColumnFormat = repmat({'numeric'},1,size(columnNames,1));
set(z, 'MenuBar', 'none');

ButtonWidth = 60;
ButtonHeight = 25;
startXButton = widthFigure/2 - ButtonWidth/2;
startYButton = 20;

uicontrol('Parent',z,'Style','pushbutton','String','Save',...
    'Position',[startXButton startYButton ButtonWidth ButtonHeight],'Visible',...
    'on','Callback',@SaveButtonPushed);

popMenu = cell(1,NumberOfRows);
for n = 1:NumberOfRows
    popMenu{1,n} = sprintf('Design %d',n);
end

c = uicontrol('Parent',z,'Style','popupmenu','String',popMenu,'Position',[startXButton-8 75 75 20]);

uicontrol('Parent',z,'Style','text','String',{'Please choose a design:'},'Position',[66 71 122 20]);

function SaveButtonPushed(src,event)
    val = c.Value;
    str = c.String;
    str{val};
    fileID = fopen('selection.txt','w');
    fprintf(fileID,'%s',str{val});
    fclose(fileID);
    close
end

end

function hex = rgb2hex(rgb)
% Convert RGB to HEX.

% Check inputs: 
assert(nargin==1,'This function requires an RGB input.') 
assert(isnumeric(rgb)==1,'Function input must be numeric.') 

sizergb = size(rgb); 
assert(sizergb(2)==3,'rgb value must have three components in the form [r g b].')
assert(max(rgb(:))<=255& min(rgb(:))>=0,'rgb values must be on a scale of 0 to 1 or 0 to 255')

% If no value in RGB exceeds unity, scale from 0 to 255: 
if max(rgb(:))<=1
    rgb = round(rgb*255); 
else
    rgb = round(rgb); 
end

hex(:,2:7) = reshape(sprintf('%02X',rgb.'),6,[]).'; 
hex(:,1) = '#';

end