
%% INSTRUCTIONS
% 
% See ExtractFitPlateReaderData_General_Part1
%
% For the used fluor measurement script, the parameters below should be OK.
%

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

%%

% define function to extrapolate which ODfields belong to a certain fluor 
% time.
fnODINDEXTOFLUORINDEX = @(i) (i-1)*3+1+[-1 1];


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
% Things that can be adjusted change before calling script
% TIMEFIELD
% YFIELD
% RANGEFIELD
if ~exist ('TIMEFIELD','var'), TIMEFIELD = 'timeFluor'; end
if ~exist ('YFIELD','var'), YFIELD = 'fluorDivOD'; end
if ~exist ('RANGEFIELD','var'), RANGEFIELD = 'fitRangeFluor'; end

if isfield(USERSETTINGS, 'wellNamesToPlot')

    output = struct; 
    for i = 1:numel(USERSETTINGS.wellNamesToPlot)

        % current well to plot
        wellNameToPlot = USERSETTINGS.wellNamesToPlot{i};
        % find out which group number is associated w. this label
        plotGroupIdx = find(ismember(wellNames,wellNameToPlot));
        % now find out which datasets contain info on this label
        toPlot = cell2mat(membersOfGroup(plotGroupIdx));

        % Simply plot all separate datasets into one figure
        figure(i), clf, hold on
        for j = toPlot
            myTimes   = sortedData(j).(TIMEFIELD);
            myYvalues = sortedData(j).(YFIELD);
            plot(myTimes, myYvalues,'LineWidth',2,'Color',[.6 .6 .6])
        end
        
        % Now plot (in same plot; overlay) idem, but only within timewindow
        mySelectedTimes = {};
        mySelectedYvalues = {};
        traceYmeans = [];
        loopcount = 1;
        for j = toPlot
            % get range
            theRange = sortedData(j).(RANGEFIELD);
                        
            % store data for later usage
            mySelectedTimes{loopcount}   = sortedData(j).timeFluor(theRange);
            mySelectedYvalues{loopcount} = sortedData(j).fluorDivOD(theRange);
            traceYmeans(loopcount)       = mean(sortedData(j).fluorDivOD(theRange));
            
            % plot data (in thick red)
            plot(mySelectedTimes{loopcount}, mySelectedYvalues{loopcount},'r','LineWidth',4)
            
            % loop counter
            loopcount = loopcount + 1;
        end
        
        % plot markup
        title(USERSETTINGS.wellNamesToPlot{i});    
        xlabel('time (hr)');
        ylabel('fluor / OD');
        MW_makeplotlookbetter(14);

        % create summary var
        output(i).groupName = USERSETTINGS.wellNamesToPlot{i};
        % actual datapoints
        output(i).rawY = mySelectedYvalues;
        % mean values of traces:
        output(i).traceYmeans = traceYmeans;
        % mean and std of above means
        output(i).meanFluor = mean(traceYmeans);
        output(i).stdFluor = std(traceYmeans);
    end

    % Go back per figure and set y max
    for i = 1:numel(USERSETTINGS.wellNamesToPlot)
        figure(i), ylim([0, max([output.meanFluor])*1.1]);
    end
    
    % make errorbar plot
    % ===
    figure, clf, hold on
    [h,hErrorbar]=barwitherr([output.stdFluor],[output.meanFluor]);
    set(hErrorbar, 'LineWidth', 3) %MW added
    set(h, 'FaceColor', [.6 .6 .6]);
    MW_makeplotlookbetter(14);
    title('mean Fluor values');
    set(gca, 'XTickLabel',USERSETTINGS.wellNamesToPlot, 'XTick',1:numel(USERSETTINGS.wellNamesToPlot))
    ylabel('fluor / OD (a.u.)')
    
    % Add values per trace
    figure, clf, hold on
    for i = 1:numel(USERSETTINGS.wellNamesToPlot) % loop over labels        
        % raw
        plot(ones(1,numel(output(i).traceYmeans))*i,output(i).traceYmeans,'x','MarkerSize',15,'LineWidth',2,'Color',[.6 .6 .6])
        % mean
        currentMean = mean(output(i).traceYmeans);
        plot(i,currentMean,'ok','LineWidth',3)        
        % std        
        currentStd = std(output(i).traceYmeans);
        plot(ones(1,2)*i,[currentMean-currentStd , currentMean+currentStd],'-','LineWidth',3,'Color','k')        
    end
    ylim([0, max([output.traceYmeans])*1.1]);
    xlim([0,numel(USERSETTINGS.wellNamesToPlot)+1])
    set(gca, 'XTickLabel',USERSETTINGS.wellNamesToPlot, 'XTick',1:numel(USERSETTINGS.wellNamesToPlot))
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

%% Create a plot with determined growth rates
if isfield(USERSETTINGS, 'wellNamesToPlot')

    % Set flag for whether manual determination of timewindow was done
    if isfield(sortedData, 'muManual')
        manualDetermined = 1;
    else
        manualDetermined = 0;
    end
    
    for i = 1:numel(USERSETTINGS.wellNamesToPlot)

        % current well to plot
        wellNameToPlot = USERSETTINGS.wellNamesToPlot{i};
        % find out which group number is associated w. this label
        plotGroupIdx = find(ismember(wellNames,wellNameToPlot));
        % now find out which datasets contain info on this label
        toPlot = cell2mat(membersOfGroup(plotGroupIdx));                        
        
        % Now plot (in same plot; overlay) idem, but only within timewindow
        muValues = [];
        manualMuValues = [];
        loopcount = 1;
        for j = toPlot                        
            % store data for later usage
            muValues(end+1)         = sortedData(j).mu;
            if manualDetermined
                manualMuValues(end+1)   = sortedData(j).muManual;
            end
        end                
        
        % create summary var
        output(i).muValues = muValues;
        if manualDetermined
            output(i).manualMuValues = manualMuValues;            
        end
        
    end


% Plot values per label
if manualDetermined
    muFieldNameToPlot = 'manualMuValues';
else
    muFieldNameToPlot = 'muValues';
end
figure, clf, hold on

theXlim = [0,numel(USERSETTINGS.wellNamesToPlot)+1]
plot(theXlim, [0 0],'k-');
for i = 1:numel(USERSETTINGS.wellNamesToPlot) % loop over labels        
   
    % raw
    plot(ones(1,numel(output(i).(muFieldNameToPlot)))*i,output(i).(muFieldNameToPlot),'x','MarkerSize',15,'LineWidth',2,'Color',[.6 .6 .6])
    % mean
    currentMean = mean(output(i).(muFieldNameToPlot));
    plot(i,currentMean,'ok','LineWidth',3)        
    % std        
    currentStd = std(output(i).(muFieldNameToPlot));
    plot(ones(1,2)*i,[currentMean-currentStd , currentMean+currentStd],'-','LineWidth',3,'Color','k')        
    
end

%ylim([0, max([output.(muFieldNameToPlot)])*1.1]);
xlim(theXlim)
set(gca, 'XTickLabel',USERSETTINGS.wellNamesToPlot, 'XTick',1:numel(USERSETTINGS.wellNamesToPlot))
MW_makeplotlookbetter(14);
title('growth rate values');
set(gca, 'XTickLabel',USERSETTINGS.wellNamesToPlot, 'XTick',1:numel(USERSETTINGS.wellNamesToPlot))
ylabel('fitted growth rate (dbl/hr)')

end
%% Scatter plot fluor signal against growth rate
if isfield(USERSETTINGS, 'wellNamesToPlot')
    
figure, clf, hold on
for i = 1:numel(USERSETTINGS.wellNamesToPlot) % loop over labels        
   
    
    plot( output(i).traceYmeans, ...
          output(i).(muFieldNameToPlot), ...
          'x','MarkerSize',15,'LineWidth',2,'Color',[.6 .6 .6])
        
end

%ylim([0, max([output.(muFieldNameToPlot)])*1.1]);
%xlim([0,numel(USERSETTINGS.wellNamesToPlot)+1])

MW_makeplotlookbetter(14);
title(['For ' USERSETTINGS.wellNamesToPlot{1} ' etc.']);
xlabel('fluor intensity (a.u./OD)')
ylabel('fitted growth rate (dbl/hr)')
    
end

%% Save this data
if ~exist('DONTSAVE','var')
    myFilePath = [myFullDir currentdate 'CompleteAnalyzedData' USERSETTINGS.customSuffix '.mat'];
    save(myFilePath,'sortedData','muAvStdev','membersOfGroup','wellNames');
    disp(['Saved data in: ' 10 myFilePath]);
end

