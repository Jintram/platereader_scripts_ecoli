
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
%% CONFIG SETTINGS
% -------------------------------------------------------------------------
% Set the threshold value for which the first time this value is 
% encountered should be determined
myThreshold = 0.01
USESMOOTH = 1;


myFullDir=[myRootDir myDateDir];

% Some additional wells to be ignored aside from those marked with 
% realData=0.
toIgnore = {'karlblank','H2O','15A','15B','15C','16A','16B','16C','17A','17B','17C'};

% -------------------------------------------------------------------------
% Checks whether other script has run..
if ~mainscriptsettingran==1
    disp('ERROR: Run first section ExtractFitPlateReaderData_General.m first!'); % error('')?
end
% Already done in ExtractFitPlateReaderData_General.m; and this needs to be
% run anyways.
%{
myRootDir='U:\EXPERIMENTAL_DATA\platereader\';
myScriptDir='platereader_scripts\'; % leave empty if scripts are in root
myDateDir='2014_06_14\';
%}

% For this script to work, the following objects should be present:
if (exist('sortedData') ~= 1 || exist('membersOfGroup') ~= 1) % || exist('sortedData.OD_subtr_smooth') ~= 1)
    
    disp('ERROR: Need objects sortedData, membersOfGroups and sortedData.OD_subtr_smooth for this function to work. Run ExtractFitPlateReaderData_General.m first or load data.');
    
end

%% Find treshold times

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
        
        % Store if found
        if ~isempty(current_idx_pass)            
            % Add these to lists
            idxs_pass = [idxs_pass current_idx_pass];
            times_pass = [times_pass current_time_pass]; 
        end        
        
        
    end

    % To flag groups where no treshold was passed in all duplicates, put
    % value of -1 instead of []
    if isempty(times_pass)
        times_pass = -1;
        idxs_pass = -1;
    end
    
    % store thresholds for current wells
    all_my_thresholds{end+1} = {name, times_pass, idxs_pass, nameidx};
    
end

% Tell user where to find data.
disp('Done. See all_my_thresholds for data.');

% Store data
filename = [currentdate 'thresholds_' num2str(myThreshold)];
save([myFullDir filename '.mat'],'all_my_thresholds');

disp(['Writing to ' filename ' (.mat and .xls).'])

% Exporting to Excel ------------------------------------------------------

% Select the threshold pass times only
thresholds_pass_times = cellfun(@(x) x(2), all_my_thresholds);
% Select accompanying cell names
names = cellfun(@(x) x(1), all_my_thresholds);

% Determine averages, std and number of duplicates for each group
average_thresholds_pass_times = [];
std_thresholds_pass_times = [];
nr_duplicates = [];
for i = [1:length(all_my_thresholds)]
    average_thresholds_pass_times = ...
        [average_thresholds_pass_times mean(cell2mat(thresholds_pass_times(i)))];
    std_thresholds_pass_times = ...
        [std_thresholds_pass_times std(cell2mat(thresholds_pass_times(i)))];
    nr_duplicates = ...
        [nr_duplicates length(cell2mat(thresholds_pass_times(i)))];
end

%% Export data

% Export to Excel sheet
filename = [myFullDir currentdate 'thresholds_' num2str(myThreshold) '.xlsx'];
%filename = ['d:\' 'thresholds_' num2str(myThreshold) '.xlsx'];
myThresholdTable=cell([names', num2cell(average_thresholds_pass_times'),num2cell(std_thresholds_pass_times'),num2cell(nr_duplicates')]);
%myThresholdTable=cell([wellNames,num2cell(average_thresholds_pass_times')])
xlswrite(filename,myThresholdTable,'Thresholdvalues','B2');

% -------------------------------------------------------------------------


% Plot data
h = figure();
barh(average_thresholds_pass_times);
set(gca, 'YTick', [1:length(wellNames)]);
set(gca, 'YTickLabel', wellNames);

% save with (moving) averages on linear scale
figFullName=[myJustPlotDir 'plateauvalues' ];
saveas(h,[figFullName '.fig'], 'fig');
saveas(h,[figFullName '.png'], 'png');







