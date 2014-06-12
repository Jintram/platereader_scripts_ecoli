% THIS IS A SCRIPT TO EXTRACT DATA FROM A EXCEL FILE THAT CONTAINS GROWTH
% MEASUREMENTS FROM THE PLATE READER (VICTOR)
% DATA IS FITTED WITH EXPONENTIAL
% copy piece by piece into the command window or run single cells
% Note: file can be used for measurements >24h without alterations (.xls
% file contains intrinsically information that time is >24h)
%
% ************************************************
% REQUIRED INPUT:
% - Excel file form Plate reader (e.g. '120620.xls')
% - Excel file with description of each well (e.g. '#553' or 'glucose' in 
%                                  'plateCoordinates.xlsx')
%          ************************************************************************
%          ********** ADJUST plateCoordinates.xlsx AND SAVE IN DATEFOLDER *********
%          ************************************************************************
%
% OUTPUT:
% - .txt file with fitted growth rates
% - plot with fitted growth rates for each group of repeated experiments
% - sortedData.mat: matlab file containing all info about the wells and
%     growth rates
% ************************************************
%
%
% Some additional notes (MW):
% - "Smoothing" and "moving average" are used as synonyms here.

%% (1) Setting directories and such
% ************************************************
% specify folder, date, etc
% ************************************************
myRootDir='U:\EXPERIMENTAL_DATA\platereader\'; % should also contain folder with scripts
myDateDir='2014_06_09\';
myFullDir=[myRootDir myDateDir];

datafile='2014_06_09_directM9_failH2Orefill';
myPlotsSaveDir=[myFullDir 'Plots\'];

% Location of scripts
myScriptDir='platereader_scripts\'; % leave empty if scripts are in root

% Depends on data but needed for general functioning script
load([myRootDir myScriptDir 'myColor.mat'],'myColor'); % load MW colors

% *************************************************
% Some general setting
% *************************************************

% Some additional wells to be ignored aside from those marked with 
% realData=0.
ignoreList = {'H2O','karlblank'};
USEMARKER = 0; % Whether to use markers when plotting - only applies (3b)
yLimMinLog = 10^-3; % minimum for log scale

% Some not-to-be-changed config
availableMarkerSpecifiers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};

% ************************************************
% import data. names: 'data' 'textdata'. create saveDirectory
% ************************************************

[data textdata]=xlsread([myFullDir, datafile, '.xls']);
[~, DescriptionPlateCoordinates]=xlsread([myFullDir 'plateCoordinates.xlsx'],'C4:N11'); %what is tested on which position
load([myRootDir myScriptDir 'PositionNames.mat']); % cell array with 'A1' 'B1' etc
%
% nb: the timefield is now in IEEE format and can be converted into minutes
% by DJK_getMinutesfromTimestamp(timeIeee). Timefield will be converted to
% hours further below

%% (1b) Loading data and setting up structures
% DO NOT RUN if you have already processed data!

% OD range to fit is set at comment "% Change OD range here".

% Data structure of Excel sheet: (X=not important)
% data:     X -- X -- Nan   -- Nan -- time1 -- OD1 -- time2 -- OD2 -- time3 -- OD3 -- time4 -- OD4
% textdata: X -- X -- A01   --  X --    X --  xxxx
%                    (well nr)

% if header line was imported into 'textdata', delete this line
if size(textdata,1)==size(data,1)+1
    textdata=textdata(2:end,:);
    disp('length was corrected')
elseif size(textdata,1)==size(data,1)
    disp('length already correct')
else
    disp('''data'' and ''textdata'' length do not fit')
end
        
% ************************************************
% sort data by well position
% create one matrix for each well position and save position and matrix in struct 'sortedData'
% ************************************************

% --------------------------------------------------------------------
% ----- CELLS OF 'SORTEDDATA': -----
% wellCoordinate: A1, A2, .., B1 etc
% time: time points of complete measurement (typically 24h)
% OD: corresponding OD measuremnt
% OD_subtr: background subtracted OD (at each time point the mean of the
%           blanks from this time point is subtracted)
% fitTimeManual: manually chosen fit Time according to plots. Chosen by eye which
%          time range seems to contain most steady exponential growth
% fitTime: fitTimeManual [h] corresponding to OD range 0.03 till 0.08
%           (recommendation by Marjon. range can be adjusted further below)
% realData: =1 if real bacterial measurement. =0 if not used (x) or blank.
%           Can be set to =0 if measured data shall be ignored
% DescriptionPos: what has been measured in this well, e.g. '#553' or
%           'glucose'. Serves to group repeated measurement.
% muManual: fitted growthrate, determined with fitTimeManual range
% x0Manual: fitted initial OD, determined with fitTimeManual range
% mu: fitted growthrate, determined with fitTime range
% x0: fitted initial OD, determined with fitTime range
% --------------------------------------------------------------------


% get names of well positions
wellCoordinates=unique(textdata(:,3));

% prepare for sorting well plate data
clear sortedData;
sortedData=struct;
for i=1:length(wellCoordinates)
    sortedData(i).wellCoordinate=char(wellCoordinates(i));
    sortedData(i).time=[];
    sortedData(i).OD=[];
    sortedData(i).OD_subtr=[];
    sortedData(i).fitTime=[];    
    sortedData(i).fitTimeManual=[];
    sortedData(i).fitRange=[];
    sortedData(i).fitRangeManual=[];    
    [idx_row,idx_col]=find(strcmp(PositionNames,wellCoordinates(i))==1);
    sortedData(i).DescriptionPos=char(DescriptionPlateCoordinates(idx_row,idx_col));
    realData=(strcmp(sortedData(i).DescriptionPos,'blank')==0 & strcmp(sortedData(i).DescriptionPos,'x')==0);
    sortedData(i).realData=realData;
    sortedData(i).muManual=[];
    sortedData(i).x0Manual=[];
    sortedData(i).mu=[];
    sortedData(i).x0=[];
    sortedData(i).movingAverage=[];
    sortedData(i).rangeMovingAverage=[];
end

% sort data
for i=1:length(sortedData)
    idx=find(strcmp(textdata(:,3),sortedData(i).wellCoordinate)==1);
    clear dummy
    dummy(:,1)=[data(idx,5); data(idx,7); data(idx,9); data(idx,11)]; %time
    dummy(:,2)=[data(idx,6); data(idx,8); data(idx,10); data(idx,12)]; %OD
    % bring into right order (strange original excel format)
    dummy=sortrows(dummy,1);
    % convert time to hours
    dummy(:,1)=DJK_getMinutesFromTimestamp(dummy(:,1))/60;    
    sortedData(i).time=dummy(:,1);
    sortedData(i).OD=dummy(:,2);
    % delete NaN values in 'time' and 'OD' (can happen if experiment is
    % aborted)
    idxtime=~isnan(sortedData(i).time);
    sortedData(i).time=sortedData(i).time(idxtime);
    sortedData(i).OD=sortedData(i).OD(idxtime);
    idxod=~isnan(sortedData(i).OD);
    sortedData(i).time=sortedData(i).time(idxod);
    sortedData(i).OD=sortedData(i).OD(idxod);
    
end

% Create list w. names of wells, excluding 'x' and 'blank'
dummy=unique(DescriptionPlateCoordinates);
wellNames={};
for i=1:length(dummy);
    if strcmp(dummy(i),'x')==0 & strcmp(dummy(i),'blank')==0  % HERE possibility to exclude 
                                                              % more names (e.g. contaminated wells)
        wellNames{end+1,1}=char(dummy(i));
    end
end

% Group wells together by wellNames.
% Each cell of membersOfGroups contains a list of indices that point to
% sortedData entries that belong to the same group (i.e. the same
% wellName). The index of membersOfGroups itself corresponds to the index
% of wellNames.
membersOfGroup={};
for i = [1:length(wellNames)]
    currentGroupMembers=[];
    for j = [1:length(sortedData)]
        if strcmp(wellNames(i),sortedData(j).DescriptionPos)
            currentGroupMembers = [currentGroupMembers j];
        end
    end
    membersOfGroup(end+1) = {currentGroupMembers};
end

%----------------------------------------------------------
%create save directory
if exist(myPlotsSaveDir)~=7
  [status,msg,id] = mymkdir([myPlotsSaveDir]);
  if status == 0
    disp(['Warning: unable to mkdir ' myPlotssSaveDir ' : ' msg]);
    return;
  end
end

clear i dummy realData idx idx_col idx_row data textdata status msg id idxtime idxod

%% (2)
% ************************************************
% subtract background (blank)
% ************************************************
% x = empty tube  %for the moment, leave empty values just be
%[x_row,x_col]=find(strcmp(DescriptionPlateCoordinates,'x')==1);

% 'blank'= medium without bacteria -> subtract from all other entries
% all blank data is averaged (over all positions, but individually for each
% time point)
SHOW_FIG_BLANK=1;
if SHOW_FIG_BLANK
    figure
    clf
    title('blanks')
    hold on
    xlabel('time [h]')
    ylabel('OD (600)')
end
[blank_row,blank_col]=find(strcmp(DescriptionPlateCoordinates,'blank')==1); %find all blank positions
numBlanks=length(blank_row); % number of blanks
avBlankPerTime=zeros(size(sortedData,1),1);
avPerBlank = []; stdPerBlank = []; % MW
for i=1:numBlanks
    wellName=PositionNames(blank_row(i),blank_col(i)); % e.g. 'A1'
    idx=find(strcmp(wellCoordinates,wellName)==1); % which entry in sortedData corresponds to wellName
    avBlankPerTime=avBlankPerTime+sortedData(idx).OD;
    if SHOW_FIG_BLANK
        plot(sortedData(idx).time,sortedData(idx).OD,'x','Color', 0.8*i/numBlanks*[1 1 1])
    end
    % determine avg and std per blank
    avPerBlank = [avPerBlank mean(sortedData(idx).OD)]; % MW
    stdPerBlank = [stdPerBlank std(sortedData(idx).OD)]; % MW
end
avBlankPerTime=avBlankPerTime/numBlanks;
if SHOW_FIG_BLANK
        % plot averages per timepoint
        plot(sortedData(1).time,avBlankPerTime,'xr')
        % plot averages per blank
        figure
        errorbar(avPerBlank,stdPerBlank) % MW
        axis([-.5 numBlanks+.5 0 max(avPerBlank)*1.1])
end
totalAvBlank=mean(avBlankPerTime); % maybe use complete average instead of av per time point. Todo

%add fields with blank subtracted to 'sortedData' (only useful for fields with actual
%bacteria in it, but never bothers)
for i=1:length(sortedData)
    sortedData(i).OD_subtr=sortedData(i).OD-avBlankPerTime;
end
clear SHOW_FIG_BLANK blank_row blank_col numBlanks i idx wellName

%% (3)
% ************************************************ 
% Plot all graphs grouped by category, without fitting
% ************************************************

HIDE_GRAPHS=1; % TODO MW - doesn't work
windowSize=21; % needs to be odd number! - windowsize for moving average

% For determining plateau values of plots
myPlateauValues     = [];
myPlateauValues_std = [];
PLATEAUSTART = 0.95; % fraction of data after which averaging is performed 
                  % to estimate plateau value.

%create subSaveDirectory for these plots
myJustPlotDir=[myPlotsSaveDir 'Rawplots\'];
if exist(myJustPlotDir)~=7
  [status,msg,id] = mymkdir([myJustPlotDir]);
  if status == 0
    disp(['Warning: unable to mkdir ' myJustPlotDir ' : ' msg]);
    return;
  end
end

% -----------------------------------------------
%loop over different groups in well (wellNames)
% -----------------------------------------------
for nameidx=1:length(wellNames)    
       
    name=char(wellNames(nameidx));

    colorcounter=0;
  %  usedColors=[]; % needed for legend in correct color
  %  dataidx=[]; % array with indices of sortedData that contain 'name'
    mylegendText=[];
       
    
    % initiate PLOTS
    % normal scale plots
    h=figure(1);    
    clf
    title([name ' OD values over time'])
    hold on
    xlabel('time [h]')
    ylabel('OD')      
    if HIDE_GRAPHS
        set(h,'Visible','off'); % TODO MW - doesn't work
    end         
    % logplots
    hlog=figure(2); 
    clf    
    set(gca, 'Yscale', 'log')
    title(['log ' name ' OD values over time'])
    hold on
    xlabel('time [h]')
    ylabel('OD')        
    % hide plots if desired  
    if HIDE_GRAPHS
        set(hlog,'Visible','off'); % TODO MW - doesn't work
    end     
    
    % -----------------------------------------------
    %  loop over all wells and search for matching description (same
    %  'name')
    % MW TODO: use membersOfGroup for this.
    % -----------------------------------------------
    myCurrentDataTime = []; myCurrentDataOD_substr = []; myCurrentColors = [];
    for i=1:length(sortedData) 
        
        % if name matches current group
        if strcmp(sortedData(i).DescriptionPos,name)==1
            
            colorcounter=colorcounter+1;           
            
            % Plot the current well
            % Plot linear scale
            figure(h);
            plot(sortedData(i).time,sortedData(i).OD_subtr','x','Color',myColor(colorcounter,:)','Linewidth',2);
            % Plot log scale
            figure(hlog); % also plot on logarithmic scale
            semilogy(sortedData(i).time,sortedData(i).OD_subtr','x','Color',myColor(colorcounter,:)','Linewidth',2);
           
            % Determine moving averages for this group
            %[sortedData(i).movingAverage,sortedData(i).rangeMovingAverage]=movingaverage(sortedData(i).OD_subtr,windowSize);
            [movingAverage,rangeMovingAverage]=movingaverage(sortedData(i).OD_subtr,windowSize);
            sortedData(i).movingAverage = movingAverage;
            sortedData(i).rangeMovingAverage=rangeMovingAverage;
                       
            myColorNow=floor(myColor(colorcounter,:)*1.2); % quick n dirty edit color
            figure(h); plot(sortedData(i).time(rangeMovingAverage),movingAverage','-','Color',myColorNow,'Linewidth',2);
            figure(hlog); semilogy(sortedData(i).time(rangeMovingAverage),movingAverage','-','Color',myColorNow,'Linewidth',2);

            % Collect all data of this group in time and OD vector
            %myCurrentDataTime = [myCurrentDataTime sortedData(i).time];
            %myCurrentDataOD_substr = [myCurrentDataOD_substr sortedData(i).OD_subtr];         
            
            % Collect moving averages of this group to be able to determine
            % plateau value later.
            myCurrentDataOD_substr_movavg = [myCurrentDataOD_substr sortedData(i).OD_subtr];
            
        end
    end
    % end loop over all data and search for repetitions with same 'name'
    % -----------------------------------------------
        
    % save with (moving) averages on linear scale
    figFullName=[myJustPlotDir 'GrowthCurves_' name];
    saveas(h,[figFullName '.fig'], 'fig');
    saveas(h,[figFullName '.png'], 'png');
    % save with (moving) averages on log scale
    figFullName=[myJustPlotDir 'log_GrowthCurves_' name];
    saveas(hlog,[figFullName '.fig'], 'fig');
    saveas(hlog,[figFullName '.png'], 'png');

    % close figures
    close(h); close(hlog);
    
    % Determine plateau values for this group
    % (determined here as average of last 
    % [fraction PLATEAUSTART:totallength] values.)
    
    % First average different moving averages, if there are more lines
    % available.
    if size(myCurrentDataOD_substr_movavg,2) > 1
        myMeanCurrentDataOD_substr_movavg = mean(myCurrentDataOD_substr_movavg);
    else
        % simply take the 1 line
        myMeanCurrentDataOD_substr_movavg = myCurrentDataOD_substr_movavg;
    end
    
    % Determine plateauvalue    
    range = [ceil(size(myMeanCurrentDataOD_substr_movavg,1)*PLATEAUSTART) size(myMeanCurrentDataOD_substr_movavg,1)];
    if range(1) == 0, range(1)=1; end % this might be done prettier.. MW TODO
    currentPlateauValues     = mean(myMeanCurrentDataOD_substr_movavg(range));
    currentPlateauValues_std = std(myMeanCurrentDataOD_substr_movavg(range));

    myPlateauValues     = [myPlateauValues      currentPlateauValues];
    myPlateauValues_std = [myPlateauValues_std  currentPlateauValues_std];
      
end
% and loop over all names (different exp's)
% -----------------------------------------------

% Save estimates of plateau values
h = figure();
%barh(myPlateauValues);
barwitherr(myPlateauValues_std,myPlateauValues,'FaceColor',[0.8,0.8,0.8]);
xlabel('Strain/medium');
ylabel('OD value');
title(['Plateau values determined from ' num2str(PLATEAUSTART) '-1.00 interval'])
set(gca, 'XTick', [1:length(wellNames)]);
set(gca, 'XTickLabel', wellNames);

% save with (moving) averages on linear scale
figFullName=[myJustPlotDir 'plateauvalues' ];
saveas(h,[figFullName '.fig'], 'fig');
saveas(h,[figFullName '.png'], 'png');

% Output plateau values to Excel file
filename = [myJustPlotDir 'plateauvalues.xlsx'];
%myPlateauTable=table(wellNames,myPlateauValues'; 
myPlateauTable=cell([wellNames,num2cell(myPlateauValues'),num2cell(myPlateauValues_std')])
xlswrite(filename,myPlateauTable,'Plateauvalues','B2');

clear dummy nameidx name muAccum muManualAccum
clear xlimfit ylimfit colorcounter mylegendText g h fitTimeManualext ODcalcManual  ODcalc
clear fitline figFullName ans currentColor fid i str SHOW_FIG_FIT ODmaxline ODminline

%% (3b)
% ************************************************ 
% Plot all graphs in one plot 
% TODO MW: and also in subplot figure
% ************************************************

HIDE_GRAPHS = 0;

%create subSaveDirectory for these plots
myJustPlotDir=[myPlotsSaveDir 'Summaryplots\'];
if exist(myJustPlotDir)~=7
  [status,msg,id] = mymkdir([myJustPlotDir]);
  if status == 0
    disp(['Warning: unable to mkdir ' myJustPlotDir ' : ' msg]);
    return;
  end
end

% Determine subplots necessary
nrSubplots = ceil(sqrt(length(wellNames)));

% Initiate plots
% ===
% normal scale plots
h=figure(1);    
clf
title([' OD values over time'])
hold on
xlabel('time [h]')
ylabel('OD')      
if HIDE_GRAPHS
    set(h,'Visible','off'); % TODO MW - doesn't work
end         
% logplots
hlog=figure(2); 
clf    
set(gca, 'Yscale', 'log')
title(['log OD values over time'])
hold on
xlabel('time [h]')
ylabel('OD')        
% hide plots if desired  
if HIDE_GRAPHS
    set(hlog,'Visible','off'); % TODO MW - doesn't work
end   

% -----------------------------------------------
%loop over different groups in well (wellNames)
% -----------------------------------------------
colorcounter=0;
mylegendText=[];
myCurrentMarker = '';
maxima = [];
for nameidx=1:length(wellNames)    
       
    % administration
    name=char(wellNames(nameidx));
    colorcounter=colorcounter+1;
    if USEMARKER 
        myCurrentMarker = availableMarkerSpecifiers{nameidx}; % TODO MW: will crash if large number groups -> use mod()
    end
    
    % Skip this group if in ignorelist
    if ismember(name,ignoreList)
        disp(['FYI: ' name ' group was ignored.']);
        continue;
    end
    
  %  usedColors=[]; % needed for legend in correct color
  %  dataidx=[]; % array with indices of sortedData that contain 'name'
             
    % -----------------------------------------------
    %  loop over members of each group
    % -----------------------------------------------
    myCurrentDataTime = []; myCurrentDataOD_substr = []; myCurrentColors = [];
    countInGroup=0;
    for i=cell2mat(membersOfGroup(nameidx))
        
        % Counter for legend (see below)
        countInGroup = countInGroup+1;
        
        % Plot the current well
        rangeMA=sortedData(i).rangeMovingAverage; % range moving average
        % Plot linear scale
        figure(h);
        lineh = plot(sortedData(i).time(rangeMA),sortedData(i).movingAverage',['-' myCurrentMarker],'Color',myColor(colorcounter,:)','Linewidth',1);
        % Plot log scale
        figure(hlog); % also plot on logarithmic scale
        linehlog = semilogy(sortedData(i).time(rangeMA),sortedData(i).movingAverage',['-' myCurrentMarker],'Color',myColor(colorcounter,:)','Linewidth',1);

        % Determine maxima
        current_max = max(sortedData(i).movingAverage);
        maxima = [maxima current_max];
        
        %collect data for legend
        if (countInGroup == 1)
            if isempty(mylegendText)
                mylegendText=['''name=' name '''' ];
            else
                mylegendText=[mylegendText, ', ''name=' name '''' ];
            end        
        else
            set(get(get(lineh,'Annotation'),'LegendInformation'),...
                        'IconDisplayStyle','off'); % Exclude line from legend
            set(get(get(linehlog,'Annotation'),'LegendInformation'),...
                'IconDisplayStyle','off'); % Exclude line from legend
        end
        
    end
    % end loop over all data and search for repetitions with same 'name'
    % -----------------------------------------------
      
end
% end loop over all names (different exp's)
% -----------------------------------------------

% Set legends
figure(h)
eval(['legend(', mylegendText, ',''Location'',''NorthEastOutside'')']); % Best
figure(hlog)
eval(['legend(', mylegendText, ',''Location'',''NorthEastOutside'')']); % NorthEastOutside
yLimMaxLog = max(maxima); % hail to the queen!
ylim([yLimMinLog yLimMaxLog])

% save with (moving) averages on linear scale
figFullName=[myJustPlotDir 'All_GrowthCurves'];
saveas(h,[figFullName '.fig'], 'fig');
saveas(h,[figFullName '.png'], 'png');
% save with (moving) averages on log scale
figFullName=[myJustPlotDir 'All_Log_GrowthCurves'];
saveas(hlog,[figFullName '.fig'], 'fig');
saveas(hlog,[figFullName '.png'], 'png');

% close figures
close(h); close(hlog);

%{
clear dummy nameidx name muAccum muManualAccum
clear xlimfit ylimfit colorcounter mylegendText g h fitTimeManualext ODcalcManual  ODcalc
clear fitline figFullName ans currentColor fid i str SHOW_FIG_FIT ODmaxline ODminline
%}

%% (4)
% ************************************************
%choose fitTime according to OD thresholds 
% ************************************************
% XXXXXXXXXXXXXXXX
% Change OD range here (standard = [0.03, 0.08]):
% XXXXXXXXXXXXXXXX
ODmin=4*10^-3; ODmax=13*10^-3; % does not take into account sudden random umps over threshold (e.g. avoid by averaging)
ODmin=0.03; ODmax=0.07;

%reset all actual data to 'real data' -> also "bad wells"are considered for
% fitting as real data. only background and blank are not considered.
% bad data can be excluded at point (8)
for i=1:length(sortedData)
    realData=(strcmp(sortedData(i).DescriptionPos,'blank')==0 & strcmp(sortedData(i).DescriptionPos,'x')==0);
    sortedData(i).realData=realData;
end

%create subSaveDirectory according to fit OD
myPlotsSaveDirODsub=[myPlotsSaveDir 'ODmin'  num2str(ODmin) 'ODmax' num2str(ODmax) '\'];
if exist(myPlotsSaveDirODsub)~=7
  [status,msg,id] = mymkdir([myPlotsSaveDirODsub]);
  if status == 0
    disp(['Warning: unable to mkdir ' myPlotssSaveDirODsub ' : ' msg]);
    return;
  end
end



for i=1:length(sortedData)
    if sortedData(i).realData==1 % ignore blanks and empty wells
        
        % obtain all values above and below specified values ODmin ODmax
        idxODminOrHigher=find(sortedData(i).movingAverage>ODmin);
        idxODmaxOrLower=find(sortedData(i).movingAverage<ODmax);
        
        % get min and max of those
        idxMin=min(idxODminOrHigher);
        idxMax=max(idxODmaxOrLower);
        
        % moving avg has different timerange, so get that one
        myTimes = sortedData(i).time(sortedData(i).rangeMovingAverage);
        
        % use these as base for time window
        startTime = myTimes(idxMin);
        endTime   = myTimes(idxMax);
        if isempty(startTime)
            startTime=min(sortedData(i).time);
            disp('Warning: Couldn''t find start fit time (taking min).')
        end
        if isempty(endTime)
            endTime=max(sortedData(i).time);
            disp('Warning: Couldn''t find end fit time (taking max).')
        end
        sortedData(i).fitTime=[startTime, endTime];
        sortedData(i).fitRange = ...
            find(sortedData(i).time>=startTime & sortedData(i).time<= endTime);
        
    else sortedData(i).fitTime=[]; sortedData(i).fitRange=[];
    end
end
clear i idxODmaxOrLower idxODminOrHigher idxMin idxMax status msg id

%% (5)
% ************************************************
% Select manual fitranges
% ************************************************

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
            %plot(sortedData(i).time,sortedData(i).OD_subtr','x','Color',myColor(colorcounter,:)','Linewidth',2);
            % Plot log scale
            semilogy(sortedData(i).time,sortedData(i).OD_subtr','x','Color',[.5 .5 .5],'Linewidth',2);
            semilogy(sortedData(i).time(sortedData(i).rangeMovingAverage),sortedData(i).movingAverage','-','Color','r','Linewidth',2);
                                  
            % set manual range by using ginput
            myxy = ginput();
            
            % set fitTime
            fitTimeManual = myxy(:,1)'
            sortedData(i).fitTimeManual = fitTimeManual;
            
            % also update fitrange
            sortedData(i).fitRangeManual = ...
                find(sortedData(i).time>=fitTimeManual(1) & sortedData(i).time<= fitTimeManual(2));
            
        end
end

close(m);

% Save this manual selection
myfitTimeManual = [];
myfitRangeManual = {};
for i = 1:length(sortedData)
    % get current value fitTimeManual
    currentmyfitTimeManual = sortedData(i).fitTimeManual;
    % make zero to signal it doesnt exist
    if isempty(currentmyfitTimeManual) currentmyfitTimeManual = [0, 0]; end
    
    % get current value fitRangeManual
    currentmyfitRangeManual = sortedData(i).fitRangeManual;    
    % make zero to signal it doesnt exist
    if isempty(currentmyfitRangeManual) currentmyfitRangeManual = 0; end
    
    % Add to array
    myfitTimeManual = [myfitTimeManual; currentmyfitTimeManual];
    myfitRangeManual{end+1} = currentmyfitRangeManual;
end

% save this data 
thetime=clock();
mytime=[num2str(thetime(1)) '-' num2str(thetime(2)) '-' num2str(thetime(3)) '_' num2str(thetime(4)) '-' num2str(thetime(5))];
save([myFullDir 'fitTimeManual_ranges' mytime '.mat'],'myfitTimeManual','myfitRangeManual');

%% (6)
% ************************************************
% fit growth rate to each graph
% ************************************************

% Set whether data should be smoothed first (i.e. whether moving avg should
% be used as input).
USESMOOTH = 1;

% If smoothing desired, overwrite original data (TODO MW - change this?!)
% MW BUG! overwriting original data is not allowed, as it might result in 
% smoothing the data multiple times!!
if USESMOOTH
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
        
        if USESMOOTH
            [muManual,x0Manual]=NW_ExponentialFit_fitTime(sortedData(i).time,sortedData(i).OD_subtr_smooth,sortedData(i).fitTimeManual);            
        else
            [muManual,x0Manual]=NW_ExponentialFit_fitTime(sortedData(i).time,sortedData(i).OD_subtr,sortedData(i).fitTimeManual);
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
        idx1=find(sortedData(i).time==sortedData(i).fitTime(1));
        idx2=find(sortedData(i).time==sortedData(i).fitTime(2));
        if idx2>=idx1+1
            if USESMOOTH
                [mu,x0]=NW_ExponentialFit_fitTime(sortedData(i).time,sortedData(i).OD_subtr_smooth,sortedData(i).fitTime);
            else
                [mu,x0]=NW_ExponentialFit_fitTime(sortedData(i).time,sortedData(i).OD_subtr,sortedData(i).fitTime);
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
dummy=unique(DescriptionPlateCoordinates);
wellNames={};
for i=1:length(dummy);
    if strcmp(dummy(i),'x')==0 & strcmp(dummy(i),'blank')==0  % HERE possibility to exclude 
                                                              % more names (e.g. contaminated wells)
        wellNames{end+1,1}=char(dummy(i));
    end
end

muAvStdev=zeros(length(wellNames),6); 
% mu -- stdev(mu) --#repetitions --  muManual -- stdev(muManual) -- #repetitions(Manual)


SHOW_FIG_FIT=0; % Default 1 - MW
SHOW_FIG_FITMANUAL=0; % no options to save! not kept uptodate!

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
    xlimfit=[10000 0]; %[min(sortedData(i).time) max(sortedData(i).time)]; %start with extreme values
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
    if SHOW_FIG_FIT
        h=figure('Position',[100 100 900 700]);
        clf
        title([name '.  mu fitted between OD=' num2str(ODmin) ' and ' num2str(ODmax)])
        hold on
        xlabel('time [h]')
        ylabel('OD')
        %plot OD ranges as horizontal lines                
        ODmaxline=plot(sortedData(1).time,(zeros(size(sortedData(1).time))+ODmax),'-k');
        ODminline=plot(sortedData(1).time,(zeros(size(sortedData(1).time))+ODmin),'-k');
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
                plot(sortedData(i).time,sortedData(i).OD_subtr,'Color',myColor(colorcounter,:),'Linewidth',2);
                if length(sortedData(i).fitTimeManual)==2
                    fitTimeManualext=[sortedData(i).fitTimeManual(1)-1:0.01:sortedData(i).fitTimeManual(2)+1];  %in [h]
                    ODcalcManual=sortedData(i).x0Manual*2.^(sortedData(i).muManual*fitTimeManualext);
                    plot(fitTimeManualext,ODcalcManual,'-.','Color',myColor(colorcounter,:))
                end
            end
            
            if SHOW_FIG_FIT
                figure(h)
                % plot ignored data in gray and used data in color
                if sortedData(i).realData==1
                    currentColor=myColor(colorcounter,:);
                    plot(sortedData(i).time,sortedData(i).OD_subtr,'Color',currentColor,'Linewidth',2);
                else
                    currentColor=[1 1 1]*colorcounter/(0.5+colorcounter);
                    plot(sortedData(i).time,sortedData(i).OD_subtr,'--','Color',currentColor,'Linewidth',1);
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
                    xlimfit=[sortedData(i).time(1), sortedData(i).time(end)];
                    ylimfit(1)=0;
                    ylimfit(2)=max(sortedData(i).OD_subtr*1.05);
                    
                end
                if (xlimfit(1)>xlimfit(2)) & ylimfit(1)>ylimfit(2)
                    xlimfit=[min(sortedData(i).time) max(sortedData(i).time)];
                    ylimfit=[min(sortedData(i).OD_subtr) max(sortedData(i).OD_subtr)]; 
                end
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
    if SHOW_FIG_FIT
        eval(['legend(', mylegendText, ',''Location'',''NW'')']);
        figFullName=[myPlotsSaveDirODsub 'GrowthCurves_' name '_automaticFitTime'];
        saveas(h,[figFullName '.fig'], 'fig');
        saveas(h,[figFullName '.png'], 'png');
        %save also image with full axis range
        xlim([sortedData(1).time(1) sortedData(1).time(end)]);
        ylim([0 ylimMaxDescriptionPos]);
        figFullName=[myPlotsSaveDirODsub 'Full_GrowthCurves_' name '_automaticFitTime'];
        saveas(h,[figFullName '.fig'], 'fig');
        saveas(h,[figFullName '.png'], 'png');
        close(h)
    end
    
      
end
% and loop over all names (different exp's)
% -----------------------------------------------

% give some output for growthrates and save them in .txt file. 
%-------------------------------------------------
fid = fopen([myFullDir 'FittedGrowthRateData_ODrange' num2str(ODmin) '_' num2str(ODmax) '.txt'],'wt');
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
if USESMOOTH smoothyesnow='SMOOTHED'; else smoothyesnow=''; end
filename = [myFullDir 'FittedGrowthRateData_' smoothyesnow 'ODrange' num2str(ODmin) '_' num2str(ODmax) '.xlsx'];
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
save([myFullDir 'CompleteAnalyzedData.mat'],'sortedData','muAvStdev','membersOfGroup','wellNames');

%% (7b)----------------------------
% -------------------------
% Plot half-logarithmic plots for manual fit ranges
% -------------------------
% Manual Fit range not implemented in log-plots!
SHOW_FIG_FIT=1;


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
    xlimfit=[10000 0]; %[min(sortedData(i).time) max(sortedData(i).time)]; %start with extreme values
    ylimfit=[10000 0]; %[min(sortedData(i).OD_subtr) max(sortedData(i).OD_subtr)]; %start with extreme values
    colorcounter=0;
  %  usedColors=[]; % needed for legend in correct color
  %  dataidx=[]; % array with indices of sortedData that contain 'name'
    mylegendText=[];
   
    % initiate PLOTS
    if SHOW_FIG_FIT
        h=figure('Position',[100 100 900 700]);
        clf
        set(gca, 'Yscale', 'log');
        hold on
        title(['Log-Plot. ' name '.  mu fitted manually, USESMOOTH=' num2str(USESMOOTH)]);
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
            if SHOW_FIG_FIT
                
                figure(h)
                % plot those graphs for which manual fitTime is chosen
                currentColor=myColor(colorcounter,:);
                %plot(sortedData(i).time,log(sortedData(i).OD_subtr)/log(2),'o','Color',currentColor,'Markersize',3);
                plot(sortedData(i).time,sortedData(i).OD_subtr,'o','Color',currentColor,'MarkerSize',3);
                % Highlight datapoint used for fit
                fitRangeManual=sortedData(i).fitRangeManual;
                plot(sortedData(i).time(fitRangeManual),sortedData(i).OD_subtr(fitRangeManual),'x','Color',currentColor,'MarkerSize',6,'LineWidth',3);

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
                        xlimfit=[sortedData(i).time(1), sortedData(i).time(end)];
                        ylimfit(1)=0;
                        ylimfit(2)=max(sortedData(i).OD_subtr*1.05);

                    end
                    if (xlimfit(1)>xlimfit(2)) & ylimfit(1)>ylimfit(2)
                        xlimfit=[min(sortedData(i).time) max(sortedData(i).time)];
                        ylimfit=[min(sortedData(i).OD_subtr) max(sortedData(i).OD_subtr)]; 
                    end
                    xlim([xlimfit(1) xlimfit(2)+1]);  %blubb
                    %ylim([log(ODcalc(1))/log(2)-1    log(ODcalc(end))/log(2)+1]);
                    ylim([ylimfit(1) ylimfit(2)*2]);

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
    if SHOW_FIG_FIT
        %create save Directory for log images
        myPlotsSaveDirLogODsub=[myPlotsSaveDir 'LogPlot_manualRange\'];        
        if exist(myPlotsSaveDirLogODsub)~=7
        [status,msg,id] = mymkdir([myPlotsSaveDirLogODsub]);
            if status == 0
                disp(['Warning: unable to mkdir ' myPlotssSaveDirLogODSub ' : ' msg]);
                return;
            end
        end
               
        eval(['legend(', mylegendText, ',''Location'',''Best'')']);
        figFullName=[myPlotsSaveDirLogODsub 'GrowthCurves_' name '_automaticFitTime'];
        saveas(h,[figFullName '.fig'], 'fig');
        saveas(h,[figFullName '.png'], 'png');
        %save also image with full axis range
        xlim([sortedData(1).time(1) sortedData(1).time(end)]);
        %ylim([log(0.01)/log(2) log(ylimMaxDescriptionPos)/log(2)]);
        ylim([0 ylimMaxDescriptionPos]);
        figFullName=[myPlotsSaveDirLogODsub 'Full_GrowthCurves_' name '_automaticFitTime'];
        saveas(h,[figFullName '.fig'], 'fig');
        saveas(h,[figFullName '.png'], 'png');
        close(h)
    end
    
      
end
% and loop over all names (different exp's)

%% (7c)----------------------------
% -------------------------
% Plot half-logarithmic plots to check range of exponential phase
% -------------------------
% Manual Fit range not implemented in log-plots!
SHOW_FIG_FIT=1;

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
    xlimfit=[10000 0]; %[min(sortedData(i).time) max(sortedData(i).time)]; %start with extreme values
    ylimfit=[10000 0]; %[min(sortedData(i).OD_subtr) max(sortedData(i).OD_subtr)]; %start with extreme values
    colorcounter=0;
  %  usedColors=[]; % needed for legend in correct color
  %  dataidx=[]; % array with indices of sortedData that contain 'name'
    mylegendText=[];
    
    % initiate PLOTS
    if SHOW_FIG_FIT
        h=figure('Position',[100 100 900 700]);
        clf
        %plot OD ranges as horizontal lines                
        ODmaxline=semilogy(sortedData(1).time,(zeros(size(sortedData(1).time))+ODmax),'-k');
        hold on
        ODminline=plot(sortedData(1).time,(zeros(size(sortedData(1).time))+ODmin),'-k');
        set(get(get(ODmaxline,'Annotation'),'LegendInformation'),...
                            'IconDisplayStyle','off'); % Exclude line from legend
        set(get(get(ODminline,'Annotation'),'LegendInformation'),...
                            'IconDisplayStyle','off'); % Exclude line from legend
        title(['Log-Plot. ' name '.  mu fitted between OD=' num2str(ODmin) ' and ' num2str(ODmax)])
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
            if SHOW_FIG_FIT
                figure(h)
                % plot ignored data in gray and used data in color
                if sortedData(i).realData==1
                    currentColor=myColor(colorcounter,:);
                    %plot(sortedData(i).time,log(sortedData(i).OD_subtr)/log(2),'o','Color',currentColor,'Markersize',3);
                    plot(sortedData(i).time,sortedData(i).OD_subtr,'o','Color',currentColor,'MarkerSize',3);
                    % Highlight datapoint used for fit
                    fitRange=sortedData(i).fitRange;
                    plot(sortedData(i).time(fitRange),sortedData(i).OD_subtr(fitRange),'x','Color',currentColor,'MarkerSize',6,'LineWidth',3);
                else
                    currentColor=[1 1 1]*colorcounter/(0.5+colorcounter);
                    %plot(sortedData(i).time,log(sortedData(i).OD_subtr)/log(2),'o','Color',currentColor,'MarkerSize',2);
                    plot(sortedData(i).time,sortedData(i).OD_subtr,'o','Color',currentColor,'MarkerSize',2);
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
                    xlimfit=[sortedData(i).time(1), sortedData(i).time(end)];
                    ylimfit(1)=0;
                    ylimfit(2)=max(sortedData(i).OD_subtr*1.05);
                    
                end
                if (xlimfit(1)>xlimfit(2)) & ylimfit(1)>ylimfit(2)
                    xlimfit=[min(sortedData(i).time) max(sortedData(i).time)];
                    ylimfit=[min(sortedData(i).OD_subtr) max(sortedData(i).OD_subtr)]; 
                end
                xlim([xlimfit(1) xlimfit(2)+1]);  %blubb
                %ylim([log(ODcalc(1))/log(2)-1    log(ODcalc(end))/log(2)+1]);
                ylim([ylimfit(1) ylimfit(2)*2]);
                
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
    if SHOW_FIG_FIT
        %create save Directory for log images
        myPlotsSaveDirLogODsub=[myPlotsSaveDir 'LogPlot_ODmin'  num2str(ODmin) 'ODmax' num2str(ODmax) '\'];        
        if exist(myPlotsSaveDirLogODsub)~=7
        [status,msg,id] = mymkdir([myPlotsSaveDirLogODsub]);
            if status == 0
                disp(['Warning: unable to mkdir ' myPlotssSaveDirLogODSub ' : ' msg]);
                return;
            end
        end
               
        eval(['legend(', mylegendText, ',''Location'',''Best'')']);
        figFullName=[myPlotsSaveDirLogODsub 'GrowthCurves_' name '_automaticFitTime'];
        saveas(h,[figFullName '.fig'], 'fig');
        saveas(h,[figFullName '.png'], 'png');
        %save also image with full axis range
        xlim([sortedData(1).time(1) sortedData(1).time(end)]);
        %ylim([log(0.01)/log(2) log(ylimMaxDescriptionPos)/log(2)]);
        ylim([0 ylimMaxDescriptionPos]);
        figFullName=[myPlotsSaveDirLogODsub 'Full_GrowthCurves_' name '_automaticFitTime'];
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
