%% (4) choose fitTime according to OD thresholds 
% Set USERSETTINGS.ODmin and ODmax to define tresholds

%reset all actual data to 'real data' -> also "bad wells"are considered for
% fitting as real data. only background and blank are not considered.
% bad data can be excluded at point (8)
for i=1:length(sortedData)
    realData=(strcmp(sortedData(i).DescriptionPos,'blank')==0 & strcmp(sortedData(i).DescriptionPos,'x')==0);
    sortedData(i).realData=realData;
end

%create subSaveDirectory according to fit OD
myPlotsSaveDirODsub=[myPlotsSaveDir 'ODmin'  num2str(USERSETTINGS.ODmin) 'ODmax' num2str(USERSETTINGS.ODmax) '\'];
if exist(myPlotsSaveDirODsub)~=7
  [status,msg,id] = mkdir([myPlotsSaveDirODsub]);
  if status == 0
    disp(['Warning: unable to mkdir ' myPlotssSaveDirODsub ' : ' msg]);
    return;
  end
end



for i=1:length(sortedData)
    if sortedData(i).realData==1 % ignore blanks and empty wells
        
        % obtain all values above and below specified values ODmin ODmax
        idxODminOrHigher=find(sortedData(i).movingAverage>USERSETTINGS.ODmin);
        idxODmaxOrLower=find(sortedData(i).movingAverage<USERSETTINGS.ODmax);
        
        % get min and max of those
        idxMin=min(idxODminOrHigher);
        idxMax=max(idxODmaxOrLower);
        
        % moving avg has different timerange, so get that one
        myTimes = sortedData(i).(timeField)(sortedData(i).rangeMovingAverage);
        
        % use these as base for time window
        startTime = myTimes(idxMin);
        endTime   = myTimes(idxMax);
        if isempty(startTime)
            startTime=min(sortedData(i).(timeField));
            disp('Warning: Couldn''t find start fit time (taking min).')
        end
        if isempty(endTime)
            endTime=max(sortedData(i).(timeField));
            disp('Warning: Couldn''t find end fit time (taking max).')
        end
        sortedData(i).fitTime=[startTime, endTime];
        sortedData(i).fitRange = ...
            find(sortedData(i).(timeField)>=startTime & sortedData(i).(timeField)<= endTime);
        
    else sortedData(i).fitTime=[]; sortedData(i).fitRange=[];
    end
end
clear i idxODmaxOrLower idxODminOrHigher idxMin idxMax status msg id

%% (5) Select manual fitranges
% ************************************************
% Select manual fitranges
% ************************************************
if USERSETTINGS.fitManual

    % These are important for plotting later
    ignoreList(end+1)={'blank'};
    ignoreList(end+1)={'x'};

    % Set figure
    m=figure(1);

    for i=1:length(sortedData) 

            % unless not to be set manually
            if ~ismember(sortedData(i).DescriptionPos,ignoreList) & sortedData(i).realData

                % Plot the current well
                clf;
                hold on
                set(gca, 'Yscale', 'log');
                xlabel('time [h]')
                ylabel('OD')            
                title([sortedData(i).DescriptionPos ' - PLEASE DETERMINE {X1,X2} FOR FITRANGE'])

                % Plot linear scale            
                %plot(sortedData(i).(timeField),sortedData(i).OD_subtr','x','Color',myColor(colorcounter,:)','Linewidth',2);
                % Plot log scale
                semilogy(sortedData(i).(timeField),sortedData(i).OD_subtr','x','Color',[.5 .5 .5],'Linewidth',2);
                semilogy(sortedData(i).(timeField)(sortedData(i).rangeMovingAverage),sortedData(i).movingAverage','-','Color','r','Linewidth',2);

                % set manual range by using ginput
                myxy = ginput();

                % set fitTime
                fitTimeManual = myxy(:,1)'
                sortedData(i).fitTimeManual = fitTimeManual;

                % also update fitrange
                sortedData(i).fitRangeManual = ...
                    find(sortedData(i).(timeField)>=fitTimeManual(1) & sortedData(i).(timeField)<= fitTimeManual(2));

            end
    end

    close(m);

    % Save this manual selection
    myfitTimeManual = [];
    myfitRangeManual = {};
    for i = 1:length(sortedData)
        % get current value fitTimeManual
        currentmyfitTimeManual = sortedData(i).fitTimeManual
        % make zero to signal it doesnt exist
        if isempty(currentmyfitTimeManual) currentmyfitTimeManual = [0, 0]; end

        % get current value fitRangeManual
        currentmyfitRangeManual = sortedData(i).fitRangeManual;    
        % make zero to signal it doesnt exist
        if isempty(currentmyfitRangeManual) currentmyfitRangeManual = 0; end

        % Add to array
        myfitTimeManual = [myfitTimeManual; currentmyfitTimeManual(1:2)];
        myfitRangeManual{end+1} = currentmyfitRangeManual;
    end

    % save this data 
    thetime=clock();
    mytime=[num2str(thetime(1)) '-' num2str(thetime(2)) '-' num2str(thetime(3)) '_' num2str(thetime(4)) '-' num2str(thetime(5))];
    save([myFullDir currentdate 'fitTimeManual_ranges' mytime '_' USERSETTINGS.customSuffix '.mat'],'myfitTimeManual','myfitRangeManual');

    %{
    %% To restore saved data use:
    for i = 1:length(sortedData)
              
        if ~(myfitTimeManual(i,:)==[0,0])
            sortedData(i).fitTimeManual  =myfitTimeManual(i,:);
            sortedData(i).fitRangeManual =myfitRangeManual{i};
        end
        
    end        
    %}
    
end

%% (6)
% ************************************************
% fit growth rate to each graph
% ************************************************

% If smoothing desired, overwrite original data (TODO MW - change this?!)
% MW BUG! overwriting original data is not allowed, as it might result in 
% smoothing the data multiple times!!
if USERSETTINGS.useSmooth
    % Copy data structure of OD_subtr to OD_subtr_smooth
    for i = [1:length(sortedData)]        
            sortedData(i).OD_subtr_smooth = sortedData(i).OD_subtr        
    end
        
    % Smooth the newly created dataset
    for i = [1:length(sortedData)]
        for j = [1:length(sortedData(i).rangeMovingAverage)]            
            sortedData(i).OD_subtr_smooth(sortedData(i).rangeMovingAverage(j)) = sortedData(i).movingAverage(j);
        end
    end
end


for i=1:length(sortedData)
        
    % fit growth rate according to fitTimeManual
    if (sortedData(i).realData==1 & ~isempty(sortedData(i).fitTimeManual)) %bad data has '0' as entry in fitTimeManual
        
        if USERSETTINGS.useSmooth
            [muManual,x0Manual]=NW_ExponentialFit_fitTime(sortedData(i).(timeField),sortedData(i).OD_subtr_smooth,sortedData(i).fitTimeManual);            
        else
            [muManual,x0Manual]=NW_ExponentialFit_fitTime(sortedData(i).(timeField),sortedData(i).OD_subtr,sortedData(i).fitTimeManual);
        end
        
        sortedData(i).muManual=muManual;
        sortedData(i).x0Manual=x0Manual;
    else
        sortedData(i).muManual=[];
        sortedData(i).x0Manual=[];
    end
    
    % fit growth rate according to fitTime (determined by OD thresholds)
    if sortedData(i).realData==1 & length(sortedData(i).fitTime)==2 %exclude failed wells where OD threshold is not reached
        % test if fitTime range contains enough data points (2) to perform
        % fitting
        idx1=find(sortedData(i).(timeField)==sortedData(i).fitTime(1));
        idx2=find(sortedData(i).(timeField)==sortedData(i).fitTime(2));
        if idx2>=idx1+1
            if USERSETTINGS.useSmooth
                [mu,x0]=NW_ExponentialFit_fitTime(sortedData(i).(timeField),sortedData(i).OD_subtr_smooth,sortedData(i).fitTime);
            else
                [mu,x0]=NW_ExponentialFit_fitTime(sortedData(i).(timeField),sortedData(i).OD_subtr,sortedData(i).fitTime);
            end
            sortedData(i).mu=mu;
            sortedData(i).x0=x0;
        else
            sortedData(i).mu=[];
            sortedData(i).x0=[];  
        end
    else
        sortedData(i).mu=[];
        sortedData(i).x0=[];  
    end
end
clear i muManual x0Manual mu x0 idx1 idx2

%% (7) (use (8) to change with data should be used)
% ************************************************
% plot growth curves for each group of wells (that contains same
% strain/same medium etc -> has same name)
% plot also fitted growth curve
% calculate average growth rate and its stddev
% save plots and .txt files of growth rates and save sortedData.mat
% ************************************************

% names of wells, excluding 'x' and 'blank'
currentTime=unique(DescriptionPlateCoordinates);
wellNames={};
for i=1:length(currentTime);
    if strcmp(currentTime(i),'x')==0 & strcmp(currentTime(i),'blank')==0  % HERE possibility to exclude 
                                                              % more names (e.g. contaminated wells)
        wellNames{end+1,1}=char(currentTime(i));
    end
end

muAvStdev=zeros(length(wellNames),6); 
% mu -- stdev(mu) --#repetitions --  muManual -- stdev(muManual) -- #repetitions(Manual)

% NW's colors
% myColor=[1 0 0 ; 0 1 0 ; 0 0 1; 1 0.6 0.2;  0 1 1; 0 0.5 0.5; 0 0.6 0; 0.6 0 0.4; 0.8 0.5 0; 0.7 0 0; 0.4 0.2 0.6; 0.6 0.2 0; 1 0 0 ; 0 1 0 ; 0 0 1; 1 0.6 0.2;  0 1 1; 0 0.5 0.5; 0 0.6 0; 0.6 0 0.4; 0.8 0.5 0; 0.7 0 0; 0.4 0.2 0.6; 0.6 0.2 0];
% MW's colors loaded in section (1)

% -----------------------------------------------
%loop over different groups in well (wellNames)
% -----------------------------------------------
for nameidx=1:length(wellNames)    %blubb
    name=char(wellNames(nameidx));
    muAccum=[];
    muManualAccum=[];
    xlimfit=[10000 0]; %[min(sortedData(i).(timeField)) max(sortedData(i).(timeField))]; %start with extreme values
    ylimfit=[10000 0]; %[min(sortedData(i).OD_subtr) max(sortedData(i).OD_subtr)]; %start with extreme values
    colorcounter=0;
  %  usedColors=[]; % needed for legend in correct color
  %  dataidx=[]; % array with indices of sortedData that contain 'name'
    mylegendText=[];
    
    % initiate PLOTS
    if SHOW_FIG_FITMANUAL
        g=figure;
        clf
        title([name '  muManual of fitTimeManual'])
        hold on
        xlabel('time [h]')
        ylabel('OD')
        
    end
    if USERSETTINGS.showBigFit
        h=figure('Position',[100 100 900 700]);
        clf
        title([name '.  mu fitted between OD=' num2str(USERSETTINGS.ODmin) ' and ' num2str(USERSETTINGS.ODmax)])
        hold on
        xlabel('time [h]')
        ylabel('OD')
        %plot OD ranges as horizontal lines                
        ODmaxline=plot(sortedData(1).(timeField),(zeros(size(sortedData(1).(timeField)))+USERSETTINGS.ODmax),'-k');
        ODminline=plot(sortedData(1).(timeField),(zeros(size(sortedData(1).(timeField)))+USERSETTINGS.ODmin),'-k');
        set(get(get(ODmaxline,'Annotation'),'LegendInformation'),...
                            'IconDisplayStyle','off'); % Exclude line from legend
        set(get(get(ODminline,'Annotation'),'LegendInformation'),...
                            'IconDisplayStyle','off'); % Exclude line from legend
    end   
    
    ylimMaxDescriptionPos=0; %get max OD (y axis) of repetitive measurements to adjust full axis range
    
    % -----------------------------------------------
    %  loop over all wells and search for matching description (same
    %  'name')
    % MW TODO: use membersOfGroup for this.
    % -----------------------------------------------
    for i=1:length(sortedData) 
        if strcmp(sortedData(i).DescriptionPos,name)==1
            colorcounter=colorcounter+1;
            % add data for average growth rate if realData is set to =1 )and
            % if manualFitTime exists)
            if sortedData(i).realData==1
                if length(sortedData(i).fitTimeManual==2)
                    muManualAccum=[muManualAccum, sortedData(i).muManual];
                end
                muAccum=[muAccum, sortedData(i).mu];
            end
            
            % PLOT 
            if SHOW_FIG_FITMANUAL % not kept uptodate
                figure(g)
                plot(sortedData(i).(timeField),sortedData(i).OD_subtr,'Color',myColor(colorcounter,:),'Linewidth',2);
                if length(sortedData(i).fitTimeManual)==2
                    fitTimeManualext=[sortedData(i).fitTimeManual(1)-1:0.01:sortedData(i).fitTimeManual(2)+1];  %in [h]
                    ODcalcManual=sortedData(i).x0Manual*2.^(sortedData(i).muManual*fitTimeManualext);
                    plot(fitTimeManualext,ODcalcManual,'-.','Color',myColor(colorcounter,:))
                end
            end
            
            if USERSETTINGS.showBigFit
                figure(h)
                % plot ignored data in gray and used data in color
                if sortedData(i).realData==1
                    currentColor=myColor(colorcounter,:);
                    plot(sortedData(i).(timeField),sortedData(i).OD_subtr,'Color',currentColor,'Linewidth',2);
                else
                    currentColor=[1 1 1]*colorcounter/(0.5+colorcounter);
                    plot(sortedData(i).(timeField),sortedData(i).OD_subtr,'--','Color',currentColor,'Linewidth',1);
                end
                
                
                    
                % plot fitted growth rate
                if length(sortedData(i).fitTime)==2 & ~isempty(sortedData(i).mu) 
                    if length(sortedData(i).fitTime)==2 %exclude failed wells where OD threshold is not reached
                        fitTimeext=[sortedData(i).fitTime(1)-1:0.01:sortedData(i).fitTime(2)+1];  %in [h]
                        ODcalc=sortedData(i).x0*2.^(sortedData(i).mu*fitTimeext);
                        fitline=plot(fitTimeext,ODcalc,'-.','Color',currentColor);                        
                        set(get(get(fitline,'Annotation'),'LegendInformation'),...
                            'IconDisplayStyle','off'); % Exclude line from legend
                    end
                end
                
                % get correct xlim and ylim. Also takes ignored data into
                % account!! can be changed with a if ... .realdata
                % condition
                ylimMaxDescriptionPos=max(ylimMaxDescriptionPos,max(sortedData(i).OD_subtr)*1.05);
                if length(sortedData(i).fitTime)==2 %exclude failed wells where OD threshold is not reached
                    xlimfit(1)=min(xlimfit(1),fitTimeext(1)*0.8); xlimfit(2)=max(xlimfit(2),fitTimeext(end)*1.05);
                    ylimfit(1)=min(ylimfit(1),ODcalc(1))*0.8; ylimfit(2)=max(ylimfit(2),ODcalc(end)*1.05);
                else 
                    xlimfit=[sortedData(i).(timeField)(1), sortedData(i).(timeField)(end)];
                    ylimfit(1)=0;
                    ylimfit(2)=max(sortedData(i).OD_subtr*1.05);
                    
                end
                if (xlimfit(1)>xlimfit(2)) | ylimfit(1)>ylimfit(2)
                    xlimfit=[min(sortedData(i).(timeField)) max(sortedData(i).(timeField))];
                    ylimfit=[min(sortedData(i).OD_subtr) max(sortedData(i).OD_subtr)]; 
                end
                % Set limits
                xlim(xlimfit);
                ylim(ylimfit);
                
                %collect data for legend
                if isempty(mylegendText)
                    mylegendText=['''idx=', num2str(i), ', mu=' num2str(sortedData(i).mu) '''' ];
                else
                    mylegendText=[mylegendText, ', ''idx=', num2str(i), ', mu=' num2str(sortedData(i).mu) '''' ];
                end
               % dataidx=[dataidx,i];
               % usedColors=[usedColors;currentColor];
                
            end
            
            
        end
    end
    % end loop over all data and search for repetitions with same 'name'
    % -----------------------------------------------
    muAvStdev(nameidx,:)=[mean(muAccum),std(muAccum),length(muAccum), mean(muManualAccum),...
        std(muManualAccum),length(muManualAccum)];
    
    %create legend and save image
    if USERSETTINGS.showBigFit
        eval(['legend(', mylegendText, ',''Location'',''NW'')']);
        figFullName=[myPlotsSaveDirODsub currentdate 'GrowthCurves_' name '_automaticFitTime'];
        saveas(h,[figFullName '.fig'], 'fig');
        saveas(h,[figFullName '.png'], 'png');
        %save also image with full axis range
        xlim([sortedData(1).(timeField)(1) sortedData(1).(timeField)(end)]);
        ylim([0 ylimMaxDescriptionPos]);
        figFullName=[myPlotsSaveDirODsub currentdate 'Full_GrowthCurves_' name '_automaticFitTime'];
        saveas(h,[figFullName '.fig'], 'fig');
        saveas(h,[figFullName '.png'], 'png');
        close(h)
    end
    
      
end
% and loop over all names (different exp's)
% -----------------------------------------------

% give some output for growthrates and save them in .txt file. 
%-------------------------------------------------
fid = fopen([myFullDir 'FittedGrowthRateData_ODrange' num2str(USERSETTINGS.ODmin) '_' num2str(USERSETTINGS.ODmax) '.txt'],'wt');
% function from Daan
disp(['  ']); disp(['   '])
dispAndWrite(fid, ['--------------------------------------------------------------']);
dispAndWrite(fid, ['Averages Of Fitted Growth Rates']);
dispAndWrite(fid, ['--------------------------------------------------------------']);
dispAndWrite(fid,['                 mu      stdev(mu)  #repetitions      muManual  stdev(muManual)  #repetitions(Manual)']);
for i=1:size(muAvStdev,1)
    str = sprintf('%10s     %1.4f     %1.4f        %1.0f              %1.4f        %1.4f         %1.0f', ...
        char(wellNames(i)), muAvStdev(i,1), muAvStdev(i,2), muAvStdev(i,3), muAvStdev(i,4), muAvStdev(i,5), muAvStdev(i,6));
    dispAndWrite(fid,str);
end
fclose(fid);

% Output above table also to Excel file
if USERSETTINGS.useSmooth smoothyesnow='SMOOTHED_'; else smoothyesnow=''; end
if USERSETTINGS.fitManual manualyesno='inclManual_'; else manualyesno=''; end
filename = [myFullDir currentdate 'FittedGrowthRateData_' smoothyesnow manualyesno 'ODrange' num2str(USERSETTINGS.ODmin) '_' num2str(USERSETTINGS.ODmax) '.xlsx'];
%disp(['Saving ' filename]);
%myPlateauTable=table(wellNames,myPlateauValues'; 
% Sheet w. averages and stds mu
TitleLine={'name','mu','stdev(mu)','#repetitions','muManual','stdev(muManual)','#repetitions(Manual)'}
myDataTable=cell([TitleLine;wellNames,num2cell(muAvStdev)])
xlswrite(filename,myDataTable,'FittedGrowthRateData','B2');

clear dummy nameidx name muAccum muManualAccum
clear xlimfit ylimfit colorcounter mylegendText g h fitTimeManualext ODcalcManual  ODcalc
clear fitline figFullName ans currentColor fid i str SHOW_FIG_FIT ODmaxline ODminline

%% (MW)
% Save data to matlab file for later use

%save 'sortedData' and growthrate data  'muAvStdev'
save([myFullDir currentdate 'CompleteAnalyzedData' USERSETTINGS.customSuffix '.mat']);

%% (7b) - ONLY IF MANUALLY FITTED!
% -------------------------
% Plot half-logarithmic plots for manual fit ranges
% -------------------------

if USERSETTINGS.fitManual

    % -----------------------------------------------
    %loop over different groups in well (wellNames)
    % -----------------------------------------------
    for nameidx=1:length(wellNames)           

        name=char(wellNames(nameidx));    

        % exclude well groups for which manual range is not set + above
        if ismember(wellNames(nameidx),ignoreList)
            disp(['Skipping group ' name]);
            continue
        end

        muAccum=[];
        muManualAccum=[];
        xlimfit=[10000 0]; %[min(sortedData(i).(timeField)) max(sortedData(i).(timeField))]; %start with extreme values
        ylimfit=[10000 0]; %[min(sortedData(i).OD_subtr) max(sortedData(i).OD_subtr)]; %start with extreme values
        colorcounter=0;
      %  usedColors=[]; % needed for legend in correct color
      %  dataidx=[]; % array with indices of sortedData that contain 'name'
        mylegendText=[];

        % initiate PLOTS
        if USERSETTINGS.showBigFit
            h=figure('Position',[100 100 900 700]);
            clf
            set(gca, 'Yscale', 'log');
            hold on
            title(['Log-Plot. ' name '.  mu fitted manually, USESMOOTH=' num2str(USERSETTINGS.useSmooth)]);
            %hold on  % "hold on" only after first semilogy-plot
            xlabel('time [h]');
            ylabel('log_10(OD)','Interpreter','None');
        end   

        ylimMaxDescriptionPos=0; %get max OD (y axis) of repetitive measurements to adjust full axis range  

        % obtain indices of data with wellname(i)
        currentDataIdx = cell2mat(membersOfGroup(nameidx));
        % see if a manual fit exists here
        %sortedData(currentDataIdx).fitTimeManual

        for i = currentDataIdx

            if ~isempty(sortedData(i).fitTimeManual)

                colorcounter=colorcounter+1;
                % add data for average growth rate if realData is set to =1 )and
                % if manualFitTime exists)
                if sortedData(i).realData==1
                    muAccum=[muAccum, sortedData(i).muManual];
                end

                % PLOT 
                if USERSETTINGS.showBigFit

                    figure(h)
                    % plot those graphs for which manual fitTime is chosen
                    currentColor=myColor(colorcounter,:);
                    %plot(sortedData(i).(timeField),log(sortedData(i).OD_subtr)/log(2),'o','Color',currentColor,'Markersize',3);
                    plot(sortedData(i).(timeField),sortedData(i).OD_subtr,'o','Color',currentColor,'MarkerSize',3);
                    % Highlight datapoint used for fit
                    fitRangeManual=sortedData(i).fitRangeManual;
                    plot(sortedData(i).(timeField)(fitRangeManual),sortedData(i).OD_subtr(fitRangeManual),'x','Color',currentColor,'MarkerSize',6,'LineWidth',3);

                    % plot fitted growth rate
                    if ~isempty(sortedData(i).muManual) 

                        fitTimeext=[sortedData(i).fitTimeManual(1)-1:0.01:sortedData(i).fitTimeManual(2)+1];  %in [h]
                        ODcalc=sortedData(i).x0Manual*2.^(sortedData(i).muManual*fitTimeext);
                        %fitline=plot(fitTimeext,log(ODcalc)/log(2),'-','Color',currentColor);
                        fitline=plot(fitTimeext,ODcalc,'-','Color',currentColor,'LineWidth',2);
                            set(get(get(fitline,'Annotation'),'LegendInformation'),...
                            'IconDisplayStyle','off'); % Exclude line from legend                            

                        % get correct xlim and ylim. Also takes ignored data into
                        % account!! can be changed with a if ... .realdata
                        % condition
                        ylimMaxDescriptionPos=max(ylimMaxDescriptionPos,max(sortedData(i).OD_subtr)*1.05);
                        if length(sortedData(i).fitTime)==2 %exclude failed wells where OD threshold is not reached
                            xlimfit(1)=min(xlimfit(1),fitTimeext(1)*0.8); xlimfit(2)=max(xlimfit(2),fitTimeext(end)*1.05);
                            ylimfit(1)=min(ylimfit(1),ODcalc(1))*0.8; ylimfit(2)=max(ylimfit(2),ODcalc(end)*1.05);
                        else 
                            xlimfit=[sortedData(i).(timeField)(1), sortedData(i).(timeField)(end)];
                            ylimfit(1)=0;
                            ylimfit(2)=max(sortedData(i).OD_subtr*1.05);

                        end
                        if (xlimfit(1)>xlimfit(2)) & ylimfit(1)>ylimfit(2)
                            xlimfit=[min(sortedData(i).(timeField)) max(sortedData(i).(timeField))];
                            ylimfit=[min(sortedData(i).OD_subtr) max(sortedData(i).OD_subtr)]; 
                        end
                        xlim([xlimfit(1) xlimfit(2)+1]);  %blubb
                        %ylim([log(ODcalc(1))/log(2)-1    log(ODcalc(end))/log(2)+1]);
                        ylim1 = ylimfit(1);
                        ylim2 = ylimfit(2)*2;
                        if ylim2>ylim1
                            ylim([ylim1 ylim2]);
                        end

                       % dataidx=[dataidx,i];
                       % usedColors=[usedColors;currentColor];
                    end

                    %collect data for legend
                    if isempty(mylegendText)
                        mylegendText=['''idx=', num2str(i), ', mu=' num2str(sortedData(i).muManual) , ', datapoints = ',num2str(length(sortedData(i).fitRangeManual)),'' '''' ];
                    else
                        mylegendText=[mylegendText, ', ''idx=', num2str(i), ', mu=' num2str(sortedData(i).muManual), ', datapoints = ',num2str(length(sortedData(i).fitRangeManual)),'' '''' ];
                    end

                end

            end

        end

        %    end
        %end
        % end loop over all data and search for repetitions with same 'name'
        % -----------------------------------------------
        muAvStdev(nameidx,:)=[mean(muAccum),std(muAccum),length(muAccum), mean(muManualAccum),...
            std(muManualAccum),length(muManualAccum)];

        %create legend and save image
        if USERSETTINGS.showBigFit
            %create save Directory for log images
            myPlotsSaveDirLogODsub=[myPlotsSaveDir 'LogPlot_manualRange\'];        
            if exist(myPlotsSaveDirLogODsub)~=7
            [status,msg,id] = mkdir([myPlotsSaveDirLogODsub]);
                if status == 0
                    disp(['Warning: unable to mkdir ' myPlotssSaveDirLogODSub ' : ' msg]);
                    return;
                end
            end

            eval(['legend(', mylegendText, ',''Location'',''Best'')']);
            figFullName=[myPlotsSaveDirLogODsub currentdate 'GrowthCurves_' name '_automaticFitTime'];
            saveas(h,[figFullName '.fig'], 'fig');
            saveas(h,[figFullName '.png'], 'png');
            %save also image with full axis range
            xlim([sortedData(1).(timeField)(1) sortedData(1).(timeField)(end)]);
            %ylim([log(0.01)/log(2) log(ylimMaxDescriptionPos)/log(2)]);
            ylim([0 ylimMaxDescriptionPos]);
            figFullName=[myPlotsSaveDirLogODsub currentdate 'Full_GrowthCurves_' name '_automaticFitTime'];
            saveas(h,[figFullName '.fig'], 'fig');
            saveas(h,[figFullName '.png'], 'png');
            close(h)
        end


    end
    % and loop over all names (different exp's)
end

%% (7c)----------------------------
% -------------------------
% Plot half-logarithmic plots to check range of exponential phase
% -------------------------
% Manual Fit range not implemented in log-plots!

% NW's colors:
%myColor=[1 0 0 ; 0 1 0 ; 0 0 1; 1 0.6 0.2;  0 1 1; 0 0.5 0.5; 0 0.6 0; 0.6 0 0.4; 0.8 0.5 0; 0.7 0 0; 0.4 0.2 0.6; 0.6 0.2 0; 1 0 0 ; 0 1 0 ; 0 0 1; 1 0.6 0.2;  0 1 1; 0 0.5 0.5; 0 0.6 0; 0.6 0 0.4; 0.8 0.5 0; 0.7 0 0; 0.4 0.2 0.6; 0.6 0.2 0];
% MW's colors loaded at section (1)

% -----------------------------------------------
%loop over different groups in well (wellNames)
% -----------------------------------------------
for nameidx=1:length(wellNames)    %blubb
    name=char(wellNames(nameidx));
    muAccum=[];
    muManualAccum=[];
    xlimfit=[10000 0]; %[min(sortedData(i).(timeField)) max(sortedData(i).(timeField))]; %start with extreme values
    ylimfit=[10000 0]; %[min(sortedData(i).OD_subtr) max(sortedData(i).OD_subtr)]; %start with extreme values
    colorcounter=0;
  %  usedColors=[]; % needed for legend in correct color
  %  dataidx=[]; % array with indices of sortedData that contain 'name'
    mylegendText=[];
    
    % initiate PLOTS
    if USERSETTINGS.showBigFit
        h=figure('Position',[100 100 900 700]);
        clf
        %plot OD ranges as horizontal lines                
        ODmaxline=semilogy(sortedData(1).(timeField),(zeros(size(sortedData(1).(timeField)))+USERSETTINGS.ODmax),'-k');
        hold on
        ODminline=plot(sortedData(1).(timeField),(zeros(size(sortedData(1).(timeField)))+USERSETTINGS.ODmin),'-k');
        set(get(get(ODmaxline,'Annotation'),'LegendInformation'),...
                            'IconDisplayStyle','off'); % Exclude line from legend
        set(get(get(ODminline,'Annotation'),'LegendInformation'),...
                            'IconDisplayStyle','off'); % Exclude line from legend
        title(['Log-Plot. ' name '.  mu fitted between OD=' num2str(USERSETTINGS.ODmin) ' and ' num2str(USERSETTINGS.ODmax)])
        %hold on  % "hold on" only after first semilogy-plot
        xlabel('time [h]')
        ylabel('log_10(OD)','Interpreter','None')
        
    end   
    
    ylimMaxDescriptionPos=0; %get max OD (y axis) of repetitive measurements to adjust full axis range
    
    % -----------------------------------------------
    %  loop over all wells and search for matching description (same
    %  'name')
    % MW TODO: use membersOfGroup for this.
    % -----------------------------------------------
    for i=1:length(sortedData) 
        if strcmp(sortedData(i).DescriptionPos,name)==1
            colorcounter=colorcounter+1;
            % add data for average growth rate if realData is set to =1 )and
            % if manualFitTime exists)
            if sortedData(i).realData==1
                muAccum=[muAccum, sortedData(i).mu];
            end
            
            % PLOT 
            if USERSETTINGS.showBigFit
                figure(h)
                % plot ignored data in gray and used data in color
                if sortedData(i).realData==1
                    currentColor=myColor(colorcounter,:);
                    %plot(sortedData(i).(timeField),log(sortedData(i).OD_subtr)/log(2),'o','Color',currentColor,'Markersize',3);
                    plot(sortedData(i).(timeField),sortedData(i).OD_subtr,'o','Color',currentColor,'MarkerSize',3);
                    % Highlight datapoint used for fit
                    fitRange=sortedData(i).fitRange;
                    plot(sortedData(i).(timeField)(fitRange),sortedData(i).OD_subtr(fitRange),'x','Color',currentColor,'MarkerSize',6,'LineWidth',3);
                else
                    currentColor=[1 1 1]*colorcounter/(0.5+colorcounter);
                    %plot(sortedData(i).(timeField),log(sortedData(i).OD_subtr)/log(2),'o','Color',currentColor,'MarkerSize',2);
                    plot(sortedData(i).(timeField),sortedData(i).OD_subtr,'o','Color',currentColor,'MarkerSize',2);
                end
                
                    
                % plot fitted growth rate
                if length(sortedData(i).fitTime)==2 & ~isempty(sortedData(i).mu) 
                    if length(sortedData(i).fitTime)==2 %exclude failed wells where OD threshold is not reached
                        fitTimeext=[sortedData(i).fitTime(1)-1:0.01:sortedData(i).fitTime(2)+1];  %in [h]
                        ODcalc=sortedData(i).x0*2.^(sortedData(i).mu*fitTimeext);
                        %fitline=plot(fitTimeext,log(ODcalc)/log(2),'-','Color',currentColor);
                        fitline=plot(fitTimeext,ODcalc,'-','Color',currentColor);
                            set(get(get(fitline,'Annotation'),'LegendInformation'),...
                            'IconDisplayStyle','off'); % Exclude line from legend                            
                    end
                end
                
                % get correct xlim and ylim. Also takes ignored data into
                % account!! can be changed with a if ... .realdata
                % condition
                ylimMaxDescriptionPos=max(ylimMaxDescriptionPos,max(sortedData(i).OD_subtr)*1.05);
                if length(sortedData(i).fitTime)==2 %exclude failed wells where OD threshold is not reached
                    xlimfit(1)=min(xlimfit(1),fitTimeext(1)*0.8); xlimfit(2)=max(xlimfit(2),fitTimeext(end)*1.05);
                    ylimfit(1)=min(ylimfit(1),ODcalc(1))*0.8; ylimfit(2)=max(ylimfit(2),ODcalc(end)*1.05);
                else 
                    xlimfit=[sortedData(i).(timeField)(1), sortedData(i).(timeField)(end)];
                    ylimfit(1)=0;
                    ylimfit(2)=max(sortedData(i).OD_subtr*1.05);
                    
                end
                if (xlimfit(1)>xlimfit(2)) & ylimfit(1)>ylimfit(2)
                    xlimfit=[min(sortedData(i).(timeField)) max(sortedData(i).(timeField))];
                    ylimfit=[min(sortedData(i).OD_subtr) max(sortedData(i).OD_subtr)]; 
                end
                xlim([xlimfit(1) xlimfit(2)+1]);  %blubb
                %ylim([log(ODcalc(1))/log(2)-1    log(ODcalc(end))/log(2)+1]);
                ylim1 = ylimfit(1);
                ylim2 = ylimfit(2)*2;
                if ylim2>ylim1
                    ylim([ylim1 ylim2]);
                end
                
                %collect data for legend
                if isempty(mylegendText)
                    mylegendText=['''idx=', num2str(i), ', mu=' num2str(sortedData(i).mu) , ', datapoints = ',num2str(length(sortedData(i).fitRange)),'' '''' ];
                else
                    mylegendText=[mylegendText, ', ''idx=', num2str(i), ', mu=' num2str(sortedData(i).mu), ', datapoints = ',num2str(length(sortedData(i).fitRange)),'' '''' ];
                end
               % dataidx=[dataidx,i];
               % usedColors=[usedColors;currentColor];
                
            end
            
            
        end
    end
    % end loop over all data and search for repetitions with same 'name'
    % -----------------------------------------------
    muAvStdev(nameidx,:)=[mean(muAccum),std(muAccum),length(muAccum), mean(muManualAccum),...
        std(muManualAccum),length(muManualAccum)];
    
    %create legend and save image
    if USERSETTINGS.showBigFit
        %create save Directory for log images
        myPlotsSaveDirLogODsub=[myPlotsSaveDir 'LogPlot_ODmin'  num2str(USERSETTINGS.ODmin) 'ODmax' num2str(USERSETTINGS.ODmax) '\'];        
        if exist(myPlotsSaveDirLogODsub)~=7
        [status,msg,id] = mkdir([myPlotsSaveDirLogODsub]);
            if status == 0
                disp(['Warning: unable to mkdir ' myPlotssSaveDirLogODSub ' : ' msg]);
                return;
            end
        end
               
        eval(['legend(', mylegendText, ',''Location'',''Best'')']);
        figFullName=[myPlotsSaveDirLogODsub currentdate 'GrowthCurves_' name '_automaticFitTime'];
        saveas(h,[figFullName '.fig'], 'fig');
        saveas(h,[figFullName '.png'], 'png');
        %save also image with full axis range
        xlim([sortedData(1).(timeField)(1) sortedData(1).(timeField)(end)]);
        %ylim([log(0.01)/log(2) log(ylimMaxDescriptionPos)/log(2)]);
        ylim([0 ylimMaxDescriptionPos]);
        figFullName=[myPlotsSaveDirLogODsub currentdate 'Full_GrowthCurves_' name '_automaticFitTime'];
        saveas(h,[figFullName '.fig'], 'fig');
        saveas(h,[figFullName '.png'], 'png');
        close(h)
    end
    
      
end
% and loop over all names (different exp's)



% MW - want to keep data.
%{
clear dummy SHOW_FIG_FIT SHOW_FIG_FITMANUAL myColor nameidx name muAccum muManualAccum
clear xlimfit ylimfit colorcounter mylegendText g h fitTimeManualext ODcalcManual fitTimeext ODcalc
clear fitline figFullName ans currentColor fid i str ODmaxline ODminline
%}



%% (8)
% ************************************************
% set which data should be used and which should be ignored,
% i.e. set sortedData(i).realData to =0 or =1
% to determine the idx of the data, run (7) and consult the legend
% ************************************************
%{
disp('If data should not be used, sortedData(i).realData is set to =0,')
disp('If data should be used, it is set to =1. Blank and not used (''x'')')
disp('data can never be set to =1. Only enter numbers where you wish')
disp('to change the status.')
disp('If e.g. the data of well 4 and 11 should not be used, enter:')
disp('[4 11]')
DontUseMe=input('NOT to be used: ');
UseMe=input('TO be used: ');

for i=1:length(DontUseMe)
    sortedData(DontUseMe(i)).realData=0;
end
for i=1:length(UseMe)
    % data can be set to=1 at all
    if strcmp(sortedData(UseMe(i)).DescriptionPos,'x')==0 & strcmp(sortedData(UseMe(i)).DescriptionPos,'blank')==0
        sortedData(UseMe(i)).realData=1;
    end
end

clear i DontUseMe UseMe
%}