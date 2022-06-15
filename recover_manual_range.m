

% Recover manual range into sortedData from saved 
% fitTimeManual_ranges[DATE].mat file 
myfilename = 'fitTimeManual_ranges2014-4-11_20-7.mat';

% -------------------------------------------------------------------------
% specify folder, date, etc (also done in
% ExtractFitPlateReaderData_General.m).

myRootDir='U:\PROJECTS\Temperature_Mutants\platereader\';
myScriptDir='platereader_scripts\'; % leave empty if scripts are in root
myDateDir='2014_03_29\';

myFullDir=[myRootDir myDateDir];
% -------------------------------------------------------------------------

% Load appropiate file
load([myFullDir myfilename],'myfitTimeManual','myfitRangeManual');

% put in sortedData
for i = 1:length(sortedData)
    % get current value fitTimeManual
    if (myfitTimeManual(i,:)==[0, 0])
        sortedData(i).fitTimeManual = [];   
    else
        sortedData(i).fitTimeManual = myfitTimeManual(i,:);
    end
    if (cell2mat(myfitRangeManual(i)) == 0)
        sortedData(i).fitRangeManual = [];
    else
        sortedData(i).fitRangeManual = cell2mat(myfitRangeManual(i));
    end
end

display('Done.')