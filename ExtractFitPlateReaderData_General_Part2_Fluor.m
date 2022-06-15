
%% INSTRUCTIONS
% 
% See ExtractFitPlateReaderData_General_Part1
%
% For the used fluor measurement script, the parameters below should be OK.
%
% When finishing, continue with 
% ExtractFitPlateReaderData_General_Part3

%% USER PARAMETERS

if USERSETTINGS.ODorFluor==1
    timeField    = 'timeOD';
    yField       = 'OD';
    yField_subtr = 'OD_subtr';    
elseif USERSETTINGS.ODorFluor==2
    timeField    = 'timeFluor';
    yField       = 'fluor';
    yField_subtr = 'fluor_subtr';
else
    error('Something went wrong. USERSETTINGS.ODorFluor should be 1 or 2.');
end

if ~isfield(USERSETTINGS,'platereader')
    warning('Platereader not set, assuming it is the new one..');
    USERSETTINGS.platereader = 'new';
end

%%

% THIS PARAMETER MIGHT NEED TO BE SET DEPENDING ON THE PROTOCOL USED
% define function to extrapolate which ODfields belong to a certain fluor 
% time.
if strcmp(lower(USERSETTINGS.platereader),'new')
    disp('Using preset new reader protocol OD/fluor index conversion');
    fnODINDEXTOFLUORINDEX = @(i) (i-1)*3+1+[-1 1];
elseif strcmp(lower(USERSETTINGS.platereader),'old')
    disp('Using preset old reader protocol OD/fluor index conversion');
    fnODINDEXTOFLUORINDEX = @(i) (i)*3;
else
    error('No platereader defined (USERSETTINGS.platereader ~= ''old''/''new'')..');
end


%% Now divide by OD.

if isfield('fitTimeManual',sortedData);
    TIMERANGEFIELD = 'fitTimeManual';
else
    TIMERANGEFIELD = 'fitTime';
end

for i = 1:numel(sortedData)
    
    % MW REMOVE TODO
    %i=38 
    
    %myTimesOD    = sortedData(i).timeOD;
    %myTimesFluor = sortedData(i).timeFluor;

    if sortedData(i).realData == 1
    
        % Create fitRange for fluo
        startTime = sortedData(i).(TIMERANGEFIELD)(1);
        endTime   = sortedData(i).(TIMERANGEFIELD)(2);

        sortedData(i).fitRangeFluor = ...
            find(sortedData(i).(timeField)>=startTime & sortedData(i).(timeField)<= endTime);
        

        % Calculate 
        sortedData(i).fluorDivOD = [];
        for j = 2:numel(sortedData(i).timeFluor)
            % matching OD indices
            ODindicesThisPoint = fnODINDEXTOFLUORINDEX(j);
            % Check whether indices indeed correct (commented out as only valid speficic protocol)
            %{
            if ~( (sortedData(i).timeOD(ODindicesThisPoint(1)) < sortedData(i).timeFluor(j)) && ...
                  (sortedData(i).timeOD(ODindicesThisPoint(2)) > sortedData(i).timeFluor(j)) ) 
                  error('End points do not straddle. Index selection went wrong!');
                  % The times from which the ODs have been taken, should
                  % straddle the timepoint from which the fluor value comes.
            end
            %}
            sortedData(i).fluorDivOD(j) = ...
               sortedData(i).fluor_subtr(j) / mean(sortedData(i).OD_subtr(ODindicesThisPoint));
        end
    
    end
    
end

%% Save this data
if ~exist('DONTSAVE','var')
    myFilePath = [myFullDir currentdate 'CompleteAnalyzedData' USERSETTINGS.customSuffix '.mat'];
    save(myFilePath);
    disp(['Saved data in: ' 10 myFilePath]);
end

%% Automatically call the plotting script
% (This is also useful for backwards compatiblity);

%{
if ~exist('DONTPLOT','var')
    ExtractFitPlateReaderData_General_Part3_Plotting
end
%}









