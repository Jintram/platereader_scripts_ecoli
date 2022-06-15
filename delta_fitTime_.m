
% File is dependent on ExtractFitPlateReaderData_General.m.
%
% DEPENDENCIES
% This file needs first section of ExtractFitPlateReaderData_General.m
% to have been run.
% Also needs sortedData and membersOfGroups to exist (can be generated
% by ExtractFitPlateReaderData_General or loaded from .mat).
%
% USE
% Get delta time between manual fit range. 
% Convenient to estimate length of period of exponential growth.

% Some additional wells to be ignored aside from those marked with 
% realData=0.
toIgnore = {'karlblank','H2O','15A','15B','15C','16A','16B','16C','17A','17B','17C'};

% For this script to work, the following objects should be present:
if (exist('sortedData') ~= 1 || ~exist('membersOfGroups') ~= 1)    
    error('Need objects sortedData and membersOfGroups for this function to work. Run ExtractFitPlateReaderData_General.m first or load data.');    
end

all_my_deltas = {};

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
    
    % Vector with delta time values
    deltas = [];
    
    for j = currentDataIdxs        
               
        % Store if found
        if ~isempty(sortedData(j).fitTimeManual)            
            % Determine delta time
            current_delta = sortedData(j).fitTimeManual(2)-sortedData(j).fitTimeManual(1);        
            
            % Add these to lists
            deltas = [deltas current_delta]; 
        end                
        
    end

    % To flag groups where no treshold was passed in all duplicates, put
    % value of -1 instead of []
    if isempty(deltas)
        deltas = -1;
    end
    
    % store thresholds for current wells
    all_my_deltas{end+1} = {name, deltas, nameidx};
    
end

% Select the deltas only
deltas_only = cellfun(@(x) x(2), all_my_deltas);
% Select accompanying cell names
names = cellfun(@(x) x(1), all_my_deltas);

% Determine averages, std and number of duplicates for each group
average_deltas = [];
std_deltas = [];
nr_duplicates = [];
for i = [1:length(all_my_deltas)]
    average_deltas = ...
        [average_deltas mean(cell2mat(deltas_only(i)))];
    std_deltas = ...
        [std_deltas std(cell2mat(deltas_only(i)))];
    nr_duplicates = ...
        [nr_duplicates length(cell2mat(deltas_only(i)))];
end

% Store data --------------------------------------------------------------
filename = ['deltas'];
save([myFullDir filename '.mat'],'all_my_deltas');

disp(['Writing to ' filename ' (.mat and .xls).'])

% Tell user where to find data.
disp('Done. See all_my_deltas for data.');

% Exporting to Excel ------------------------------------------------------

% Export to Excel sheet
filename = [myFullDir filename '.xlsx'];
myDeltaTable=cell([names', num2cell(average_deltas'),num2cell(std_deltas'),num2cell(nr_duplicates')]);
%myThresholdTable=cell([wellNames,num2cell(average_thresholds_pass_times')])
xlswrite(filename,myDeltaTable,'Deltas','B2');




