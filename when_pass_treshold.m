
% This script gives you times after which the OD data from the plate reader
% reach a certain value. 
%
% It sorts these time per wellName.
% 
% Output is structure like 
% all_my_thresholds = {name, times_pass, idxs_pass, nameidx};
% Where 
% - name: name of well group
% - times_pass: times when threshold was passed
% - idxs_pass: indices of when time was passed in sortedData.times() var
% - nameidx: index of name in wellNames (convenient to retreieve idx of
% specific well, use membersOfGroups for this).

% -------------------------------------------------------------------------
% CONFIG SETTINGS
% -------------------------------------------------------------------------
% Set the threshold value for which the first time this value is 
% encountered should be determined
myThreshold = 7*10^-2;
USESMOOTH = 1;

% specify folder, date, etc (also done in
% ExtractFitPlateReaderData_General.m).

myRootDir='U:\EXPERIMENTAL_DATA\platereader\';
myScriptDir='platereader_scripts\'; % leave empty if scripts are in root
myDateDir='2014_04_20\';

myFullDir=[myRootDir myDateDir];
% -------------------------------------------------------------------------

% Some additional wells to be ignored aside from those marked with 
% realData=0.
toIgnore = {'karlblank','H2O','15A','15B','15C','16A','16B','16C','17A','17B','17C'};

% For this script to work, the following objects should be present:
if (exist('sortedData') ~= 1 || ~exist('membersOfGroups') ~= 1)
    
    error('Need objects sortedData and membersOfGroups for this function to work. Run ExtractFitPlateReaderData_General.m first or load data.');
    
end

all_my_thresholds = {};

% Loop over wellNames
for nameidx = 1:length(wellNames)
          
    % Get current group name
    name = char(wellNames(nameidx));

    % Skip the loop if this wellName should be ignored
    if ismember(name,toIgnore)
        continue;
    end
    
    % Collect indices of wells that belong to this group
    currentDataIdxs = cell2mat(membersOfGroup(nameidx));
    
    % Times associated with crossing the threshold 
    idxs_pass = [];
    times_pass = [];
    
    for j = currentDataIdxs
        
        if USESMOOTH
            % Find index where threshold is passed for this well (from
            % smoothed line)
            current_idx_pass = find(sortedData(j).OD_subtr_smooth>myThreshold,1);
        else
            % Find index where threshold is passed for this well
            current_idx_pass = find(sortedData(j).OD_subtr>myThreshold,1);
        end
        
        % Determine corresponding time
        current_time_pass = sortedData(j).time(current_idx_pass);        
        
        if isempty(current_idx_pass)
            current_idx_pass = -1;
            current_time_pass = -1;
        end        
        
        % Add these to lists
        idxs_pass = [idxs_pass current_idx_pass];
        times_pass = [times_pass current_time_pass]; 
    end

    % store thresholds for current wells
    all_my_thresholds{end+1} = {name, times_pass, idxs_pass, nameidx};
    
end

% Tell user where to find data.
disp('Done. See all_my_thresholds for data.');

% Store data
filename = ['thresholds_' num2str(myThreshold)];
save([myFullDir filename '.mat'],'all_my_thresholds');

disp(['Writing to ' filename ' (.mat and .xls).'])

% Exporting to Excel ------------------------------------------------------

% Select the threshold pass times only
thresholds_pass_times = cellfun(@(x) x(2), all_my_thresholds);
% Select accompanying cell names
names = cellfun(@(x) x(1), all_my_thresholds);

% Determine averages and 
average_thresholds_pass_times = mean(cell2mat(thresholds_pass_times')');
std_thresholds_pass_times = std(cell2mat(thresholds_pass_times')');

% Export to Excel sheet
filename = [myFullDir 'thresholds_' num2str(myThreshold) '.xlsx'];
%filename = ['d:\' 'thresholds_' num2str(myThreshold) '.xlsx'];
myThresholdTable=cell([names', num2cell(average_thresholds_pass_times'),num2cell(std_thresholds_pass_times')]);
%myThresholdTable=cell([wellNames,num2cell(average_thresholds_pass_times')])
xlswrite(filename,myThresholdTable,'Thresholdvalues','B2');

% -------------------------------------------------------------------------

%{
% Plot data
h = figure();
barh(average_thresholds_pass_times);
set(gca, 'YTick', [1:length(wellNames)]);
set(gca, 'YTickLabel', wellNames);

% save with (moving) averages on linear scale
figFullName=[myJustPlotDir 'plateauvalues' ];
saveas(h,[figFullName '.fig'], 'fig');
saveas(h,[figFullName '.png'], 'png');
%}






