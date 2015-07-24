
%% INSTRUCTIONS
% 
% See ExtractFitPlateReaderData_General_Part1
%
% For the used fluor measurement script, the parameters below should be OK.

%% USER PARAMETERS

% define function to extrapolate which ODfields belong to a certain fluor 
% time.
fnODINDEXTOFLUORINDEX = @(i) (i-1)*3+1+[-1 1];


%% Now divide by OD.

for i = 1:numel(sortedData)
    
    % MW REMOVE TODO
    %i=38 
    
    %myTimesOD    = sortedData(i).timeOD;
    %myTimesFluor = sortedData(i).timeFluor;

    if sortedData(i).realData == 1
    
        % Create fitRange for fluo
        startTime = sortedData(i).fitTime(1);
        endTime   = sortedData(i).fitTime(2);

        sortedData(i).fitRangeFluor = ...
            find(sortedData(i).(timeField)>=startTime & sortedData(i).(timeField)<= endTime);

        % Calculate 
        sortedData(i).fluorDivOD = [];
        for j = 2:numel(sortedData(i).timeFluor)
            % matching OD indices
            ODindicesThisPoint = fnODINDEXTOFLUORINDEX(j);
            % Check whether indices indeed correct
            if ~( (sortedData(i).timeOD(ODindicesThisPoint(1)) < sortedData(i).timeFluor(j)) && ...
                  (sortedData(i).timeOD(ODindicesThisPoint(2)) > sortedData(i).timeFluor(j)) ) 
                  error('End points do not straddle. Index selection went wrong!');
                  % The times from which the ODs have been taken, should
                  % straddle the timepoint from which the fluor value comes.
            end
            sortedData(i).fluorDivOD(j) = ...
               sortedData(i).fluor_subtr(j) / mean(sortedData(i).OD_subtr(ODindicesThisPoint));
        end
    
    end
    
end

%% Now plot

if isfield(USERSETTINGS, 'wellNamesToPlot')

    output = struct; 
    for i = 1:numel(USERSETTINGS.wellNamesToPlot)

        % 
        wellNameToPlot = USERSETTINGS.wellNamesToPlot{i};

        plotGroupIdx = find(ismember(wellNames,wellNameToPlot));

        toPlot = cell2mat(membersOfGroup(plotGroupIdx));

        % All data
        figure(i), clf, hold on
        for j = toPlot
            myTimes   = sortedData(j).timeFluor;
            myYvalues = sortedData(j).fluorDivOD;
            plot(myTimes, myYvalues,'LineWidth',2,'Color',[.6 .6 .6])
        end
        % Selected window data
        for j = toPlot
            if isfield(sortedData,'fitRangeManualFluor')
                fitRangeFluor = sortedData(j).fitRangeManualFluor;
            else
                fitRangeFluor = sortedData(j).fitRangeFluor;
            end
            mySelectedTimes   = sortedData(j).timeFluor(fitRangeFluor);
            mySelectedYvalues = sortedData(j).fluorDivOD(fitRangeFluor);
            plot(mySelectedTimes, mySelectedYvalues,'r','LineWidth',4)
        end
        
        % plot markup
        title(USERSETTINGS.wellNamesToPlot{i});    
        xlabel('time (hr)');
        ylabel('fluor / OD');
        MW_makeplotlookbetter(14);

        % create summary var
        output(i).groupName = USERSETTINGS.wellNamesToPlot{i};
        output(i).meanFluor = mean(mySelectedYvalues);
        output(i).stdFluor = std(mySelectedYvalues);
    end

    % Go back per figure and set y max
    for i = 1:numel(USERSETTINGS.wellNamesToPlot)
        figure(i), ylim([0, max([output.meanFluor])*1.1]);
    end
    
    % make errorbar plot
    % ===
    figure, clf, 
    [h,hErrorbar]=barwitherr([output.stdFluor],[output.meanFluor]);
    set(hErrorbar, 'LineWidth', 3) %MW added
    set(h, 'FaceColor', [.6 .6 .6]);
    MW_makeplotlookbetter(14);
    title('mean Fluor values');
    set(gca, 'XTickLabel',USERSETTINGS.wellNamesToPlot, 'XTick',1:numel(USERSETTINGS.wellNamesToPlot))
    ylabel('fluor / OD (a.u.)')
else
    error('wellNamesToPlot was not set!');
end

%{
sortedData.fluorNormalized = ...
    sortedData.fluor_subtr ./ sortedData.fluor_subtr
%}

%% Save this data
myFilePath = [myFullDir currentdate 'CompleteAnalyzedData' USERSETTINGS.customSuffix '.mat'];
save(myFilePath,'sortedData','muAvStdev','membersOfGroup','wellNames');
disp(['Saved data in: ' 10 myFilePath]);


