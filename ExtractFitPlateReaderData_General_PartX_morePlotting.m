
%% INSTRUCTIONS
% 
% See ExtractFitPlateReaderData_General_Part1
%
% For the used fluor measurement script, the parameters below should be OK.


%% Plot some more: mean fluor values

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
            myMu   = sortedData(j).mu;
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


