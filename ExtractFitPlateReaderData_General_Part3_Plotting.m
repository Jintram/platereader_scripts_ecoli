
%% README
%
% See ExtractFitPlateReaderData_General_Part1 for more instructions.
%
% First execute 
% ===
% ExtractFitPlateReaderData_General_Part1
% Then
% ExtractFitPlateReaderData_General_Part2_OD 
% OR
% ExtractFitPlateReaderData_General_Part2_GFP
% 
% Data is then saved. And can be loaded also, in which case one can
% immediately call the plotting script if correct USERSETTINGS have been
% set.
%
% Then continue with this
% ExtractFitPlateReaderData_General_Part3_Plotting


%% For plotting when data has been loaded from file------------------------

if ~exist('myFullDir','var')
    % Directory with datafiles
    myFullDir=[USERSETTINGS.myRootDir USERSETTINGS.myDateDir];
end
if ~exist('myPlotsSaveDir','var')
    % Output directory for plots
    myPlotsSaveDir=[myFullDir 'Plots' USERSETTINGS.customSuffix '\'];
end

if ~exist('myJustPlotDir','var')
    % in case script is called directly
    myJustPlotDir=[myPlotsSaveDir MYCATEGORIEPLOTDIRNAME];
    % make if doesn't exist
    if exist(myJustPlotDir)~=7
      [status,msg,id] = mymkdir([myJustPlotDir]);
      if status == 0
        disp(['Warning: unable to mkdir ' myJustPlotDir ' : ' msg]);
        return;
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

saveas(gcf,[myJustPlotDir 'groupedfluor' USERSETTINGS.wellNamesToPlot{1} 'etc.png'],'png')
saveas(gcf,[myJustPlotDir 'groupedfluor' USERSETTINGS.wellNamesToPlot{1} 'etc.eps'],'epsc')


%{
sortedData.fluorNormalized = ...
    sortedData.fluor_subtr ./ sortedData.fluor_subtr
%}

%% Create a plot with determined growth rates
if isfield(USERSETTINGS, 'wellNamesToPlot')

    manualDetermined = USERSETTINGS.fitManual;    
    
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

saveas(gcf,[myJustPlotDir 'groupedgrowth' USERSETTINGS.wellNamesToPlot{1} 'etc.png'],'png')
saveas(gcf,[myJustPlotDir 'groupedgrowth' USERSETTINGS.wellNamesToPlot{1} 'etc.eps'],'epsc')


% Plot of distribution of growth rate values
% ===
figure, clf, hold on;
[n,centers]=hist([output(:).(muFieldNameToPlot)]);
l=bar(centers,n,'FaceColor',[.6 .6 .6],'LineWidth',1);
%set(l,'FaceColor',[.6 .6 .6],'LineWidth',1);
xlabel('Growth rate (dbl/hr)'); ylabel('Count'); title(['For ' USERSETTINGS.wellNamesToPlot{1} ' etc.']); MW_makeplotlookbetter(20);
if exist('MYGROWTHRATEXMAX','var')
    xlim([0,MYGROWTHRATEXMAX]);
else, xlim([0,1]);
end

saveas(gcf,[myJustPlotDir 'groupedgrowthdistr' USERSETTINGS.wellNamesToPlot{1} 'etc.png'],'png')
saveas(gcf,[myJustPlotDir 'groupedgrowthdistr' USERSETTINGS.wellNamesToPlot{1} 'etc.eps'],'epsc')


end
%% Scatter plot fluor signal against growth rate
% Optional parameters:
% - MYXMAX  X axis max of fluor PDF

if isfield(USERSETTINGS, 'wellNamesToPlot')
    
figure, clf, hold on
for i = 1:numel(USERSETTINGS.wellNamesToPlot) % loop over labels        
       
    plot( output(i).traceYmeans, ...
          output(i).(muFieldNameToPlot), ...
          'x','MarkerSize',15,'LineWidth',2,'Color','k')
        
end

ylim([0, max([output.(muFieldNameToPlot)])*1.1]);
%xlim([0,numel(USERSETTINGS.wellNamesToPlot)+1])

MW_makeplotlookbetter(14);
title(['For ' USERSETTINGS.wellNamesToPlot{1} ' etc.']);
xlabel('fluor intensity (a.u./OD)')
ylabel('fitted growth rate (dbl/hr)')
    
saveas(gcf,[myJustPlotDir 'groupedscatter' USERSETTINGS.wellNamesToPlot{1} 'etc.png'],'png')
saveas(gcf,[myJustPlotDir 'groupedscatter' USERSETTINGS.wellNamesToPlot{1} 'etc.eps'],'epsc')

%set(gca, 'Xscale', 'log')

% Plot of distribution of fluor values
% ===
figure, clf, hold on;
[n,centers]=hist([output(:).traceYmeans]);
l=bar(centers,n,'FaceColor',[.6 .6 .6],'LineWidth',1);
%set(l,'FaceColor',[.6 .6 .6],'LineWidth',1);
xlabel('Fluor signal (a.u.)'); ylabel('Count'); title(['For ' USERSETTINGS.wellNamesToPlot{1} ' etc.']); MW_makeplotlookbetter(20);
if exist('MYXMAX','var')
    xlim([0,MYXMAX]);
else, xlim([0,300000]);
end

saveas(gcf,[myJustPlotDir 'groupedfluordistr' USERSETTINGS.wellNamesToPlot{1} 'etc.png'],'png')
saveas(gcf,[myJustPlotDir 'groupedfluordistr' USERSETTINGS.wellNamesToPlot{1} 'etc.eps'],'epsc')

end