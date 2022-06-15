% THIS IS A SCRIPT TO EXTRACT DATA FROM A EXCEL FILE THAT CONTAINS GROWTH
% MEASUREMENTS FROM THE PLATE READER (VICTOR)
% DATA IS FITTED WITH EXPONENTIAL
% ===
% Written by Noreen Walker and Martijn Wehrens.
%
% Note: file can be used for measurements >24h without alterations (.xls
% file contains intrinsically information that time is >24h)
%
% If the parameter USERSETTINGS is set correctly (see below), you can just
% run the whole script at once. It will ask for the required input.
% The script is cut up in two parts, to accomodate for both OD and
% fluorescent signal measurements. Continue with the appropriate part II to
% finish the analysis (either 
% ExtractFitPlateReaderData_General_Part2_OD.m or 
% ExtractFitPlateReaderData_General_Part2_Fluorescence.m).
%
% See README_example_execution.m for an example how to run part 1 and 2.
% Note that part 1 needs to be executed both for OD and Fluor.
% For fluor, plotting can be done from 
% ExtractFitPlateReaderData_General_Part3_Plotting (this is automatically
% called from part 2 fluor).
% 
%
% ************************************************
% REQUIRED INPUT (with examples):
% - Adjust platreader coordinates in plateCoordinates.xlsx. 
% - USERSETTINGS.myRootDir = 'T:\TRAVELING_DATA\00_PLATEREADER\'; 
%       directory where all platereader datasets are stored
% - USERSETTINGS.myDateDir = '2015-06-19\';
%       directory with platereader dataset to be analyzed (usually this
%       directory is named with the data of when the data was taken).
% - USERSETTINGS.datafile = '2015_06_19_CRPcAMP_plasmids_repeat.xls';
%       name of the excel datafile created by platereader software
%       (created with export function of that software).
% - TIMEINDEXES
%       optional parameter for special case;
%       sets column indices which are regarded as time measurements. 
%       parameter needs to be altered when e.g. one of measurements is GFP.
% - YINDEXES
%       optional parameter for special case;
%       sets column indices which are regarded as OD measurements. 
%       parameter needs to be altered when e.g. one of measurements is GFP.
% - There are also some hard-coded parameters, see first section below.
% - USERSETTINGS.fitManual
%       important parameter, as when set to 1 it activates a part of the
%       code that allows for manual selection of growth rate fit ranges.
%       can be 0 or 1.
% REQUIRED FOR FLUO ANALYSIS ONLY:
% - USERSETTINGS.wellNamesToPlot  
%       cell with the names of which wells to plot
% REQUIRED FILES:
% - Excel file form Plate reader (e.g. '120620.xls')
% - Excel file with description of each well (e.g. '#553' or 'glucose' in 
%                                  'plateCoordinates.xlsx')
%
% OUTPUT:
% - .txt file with fitted growth rates
% - plot with fitted growth rates for each group of repeated experiments
% - sortedData.mat: matlab file containing all info about the wells and
%     growth rates
% (- wellNames and membersOfGroup)
% ************************************************
%
%
% Some additional notes:
% - "Smoothing" and "moving average" are used as synonyms here.
%
%
% Data structure of Excel sheet: (X=not important)
% data:     X -- X -- Nan   -- Nan -- time1 -- OD1 -- time2 -- OD2 -- time3 -- OD3 -- time4 -- OD4
% textdata: X -- X -- A01   --  X --    X --  xxxx
%                    (well nr)

%% User settings

% User USERSETTINGS (see comment above for explanation)
if ~exist('USERSETTINGS')
    disp('WARNING: No usersettings found! Please give USERSETTINGS. Now using USERSETTINGS stored in script.')
    disp('See comments for more information.');
    
    USERSETTINGS.myRootDir='X:\Plate Reader\Clone\';
    USERSETTINGS.myDateDir='20161118\';
    USERSETTINGS.datafile='20161118.xls';
end
if ~isfield(USERSETTINGS, 'showBlankFig')
    USERSETTINGS.showBlankFig=1;
    disp('USERSETTINGS.showBlankFig set to default value, 1');
end
if ~isfield(USERSETTINGS, 'hideGraphs')
    USERSETTINGS.hideGraphs = 0;
    disp('USERSETTINGS.hideGraphs set to default value, 0'); 
end
if ~isfield(USERSETTINGS, 'ODmin')
    USERSETTINGS.ODmin=0.05; 
    USERSETTINGS.ODmax=0.10;
    disp('USERSETTINGS.ODmin and ODmax set to defaults');
end
if ~isfield(USERSETTINGS, 'useSmooth')
    % Set whether data should be smoothed first (i.e. whether moving avg should
    % be used as input).
    USERSETTINGS.useSmooth = 0;
end
if ~isfield(USERSETTINGS, 'showBigFit')
    USERSETTINGS.showBigFit=1; % Default 1 - MW
end
if ~isfield(USERSETTINGS, 'fitManual');
    USERSETTINGS.fitManual = 1;
    disp('USERSETTINGS.fitManual set to default, 0');
end
if ~isfield(USERSETTINGS, 'customSuffix')
    USERSETTINGS.customSuffix='_OD';
end
if ~isfield(USERSETTINGS, 'ODorFluor')
    USERSETTINGS.ODorFluor=1; % 1 = OD, 2 = Fluor
end

% Some hard-coded configuration options:
windowSize=21; % needs to be odd number! - windowsize for moving average

% For determining plateau values of plots
PLATEAUSTART = 0.95; % fraction of data after which averaging is performed 
                  % to estimate plateau value.

SUMMARYPLOTDIRNAME     = ['Summaryplots\'];
MYCATEGORIEPLOTDIRNAME = ['categoriePlots\'];

SHOW_FIG_FITMANUAL = 0; % This is rather unused
   
FONTSIZE=15;

% Some parameters for special cases
% If default measurement was used, default values can be used for 
% TIMEINDEXES and YINDEXES, but when e.g. platereader also measured GFP
% signal, then the values might be in different fields.
% E.g., code for 1st measurement being fluor:
%{ TIMEINDEXES=[7,9,11], YINDEXES   = [8, 10, 12] %}
if ~exist('TIMEINDEXES'), TIMEINDEXES = [5, 7, 9]; end
if ~exist('YINDEXES'), YINDEXES   = [6, 8, 10]; end

% Creating parameters, directories 
% -------------------------------------------------------------------------

% Create correct field names depending OD or fluor measurement
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

% Directory with datafiles
myFullDir=[USERSETTINGS.myRootDir USERSETTINGS.myDateDir];
% Output directory for plots
myPlotsSaveDir=[myFullDir 'Plots' USERSETTINGS.customSuffix '\'];

% Get current date to label output files
currentdate=date();

% Depends on data but needed for general functioning script
load(['myColor.mat'],'myColor'); % load MW colors

% Create subsave dir if doesn't exist.
%name for subSaveDirectory for categorie plots
myJustPlotDir=[myPlotsSaveDir MYCATEGORIEPLOTDIRNAME];
if exist(myJustPlotDir)~=7
  [status,msg,id] = mkdir([myJustPlotDir]);
  if status == 0
    disp(['Warning: unable to mkdir ' myJustPlotDir ' : ' msg]);
    return;
  end
end

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

%% Load the data itself
%[data textdata]=xlsread([myFullDir, datafile]);
% More extended version to allow for long measurement files 
% (These are spread over multiple sheets by the platereader software. A
% data sheet can be recognized by the word "List" in the name.)
[myinfo, sheets]=xlsfinfo([myFullDir, USERSETTINGS.datafile]);
data=[]; textdata=[];
sheet_idx=0;
for name=sheets
    sheet_idx=sheet_idx+1;
    if ~isempty(strfind(name{1},'List'))
        disp(['Loading sheet ' name{1}])
        [tempData tempTextdata]=xlsread([myFullDir, USERSETTINGS.datafile],sheet_idx);
        
        % if header line was imported into 'tempTextdata', delete this line
        if size(tempTextdata,1)==size(tempData,1)+1
            tempTextdata=tempTextdata(2:end,:);
            disp('length was corrected')
        elseif size(tempTextdata,1)==size(tempData,1)
            disp('length already correct')
        else
            disp('''data'' and ''textdata'' length do not fit')
        end

        data =      [data;      tempData];
        textdata =  [textdata;  tempTextdata];        
        
        %size(tempData)
        %size(tempTextdata)
        %size(data)
        %size(textdata)
    end
end

[~, DescriptionPlateCoordinates]=xlsread([myFullDir 'plateCoordinates.xlsx'],'C4:N11'); %what is tested on which position
load(['PositionNames.mat']); % cell array with 'A1' 'B1' etc
%
% nb: the timefield is now in IEEE format and can be converted into minutes
% by getMinutesFromTimestamp(timeIeee). Timefield will be converted to
% hours further below

mainscriptsettingran=1; % Flag for other scripts

%% (1b) Loading data and setting up structures
% DO NOT RUN if you have already processed data!

% OD range to fit is set at comment "% Change OD range here".
        
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
if USERSETTINGS.ODorFluor == 1
    clear sortedData;
    sortedData=struct;
end
for i=1:length(wellCoordinates)

    % Field that y-data is written to, either "OD" and "OD_substr", or 
    % "fluor" and "fluor_subtr".
    
    sortedData(i).(yField)=[];
    sortedData(i).(yField_subtr)=[];    
    
    % Only required for OD measurement
    if USERSETTINGS.ODorFluor==1
        sortedData(i).wellCoordinate=char(wellCoordinates(i));
        sortedData(i).(timeField)=[];    

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
end

% sort data (now stored in data)
for i=1:length(sortedData)
    % each line of sortedData corresponds to a well, find indices of well
    % corresponding to current line of sortedData.
    idx=find(strcmp(textdata(:,3),sortedData(i).wellCoordinate)==1);

    % Collect data for this well for processing
    currentTimes=[]; currentYvalues=[];
    for j = 1:numel(TIMEINDEXES)
        currentTimes   = [currentTimes; data(idx,TIMEINDEXES(j))]; 
        currentYvalues = [currentYvalues; data(idx,YINDEXES(j))];
        
    end
    currentTimesAndODs=[currentTimes(:), currentYvalues(:)]; %time
    % Process dummy var   
    % bring into right order (strange original excel format)
    currentTimesAndODs=sortrows(currentTimesAndODs,1);
    % convert time to hours
    currentTimesAndODs(:,1)=getMinutesFromTimestamp(currentTimesAndODs(:,1))/60;    
    
    % Put dummy var into sortedData structure
    sortedData(i).(timeField)=currentTimesAndODs(:,1);
    sortedData(i).(yField)=currentTimesAndODs(:,2);
    
    % delete NaN values in 'time' and 'OD' (can happen if experiment is
    % aborted)
    idxtime=~isnan(sortedData(i).(timeField));
    sortedData(i).(timeField)=sortedData(i).(timeField)(idxtime);
    sortedData(i).(yField)=sortedData(i).(yField)(idxtime);
    idxod=~isnan(sortedData(i).(yField));
    sortedData(i).(timeField)=sortedData(i).(timeField)(idxod);
    sortedData(i).(yField)=sortedData(i).(yField)(idxod);
    
end

% Create list w. names of wells, excluding 'x' and 'blank'
currentTime=unique(DescriptionPlateCoordinates);
wellNames={};
for i=1:length(currentTime);
    if strcmp(currentTime(i),'x')==0 & strcmp(currentTime(i),'blank')==0  % HERE possibility to exclude 
                                                              % more names (e.g. contaminated wells)
        wellNames{end+1,1}=char(currentTime(i));
    end
end

% Group wells together by wellNames.
% Each cell of membersOfGroups contains a list of indices that point to
% sortedData entries that belong to the same group (i.e. the same
% wellName). The index of membersOfGroups itself corresponds to the index
% of wellNames.
if USERSETTINGS.ODorFluor == 1 % only when OD is processed
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
end

%----------------------------------------------------------
%create save directory
if exist(myPlotsSaveDir)~=7
  [status,msg,id] = mkdir([myPlotsSaveDir]);
  if status == 0
    disp(['Warning: unable to mkdir ' myPlotssSaveDir ' : ' msg]);
    return;
  end
end


%% (2)
% ************************************************
% subtract background (blank)
% ************************************************
% x = empty tube  %for the moment, leave empty values just be
%[x_row,x_col]=find(strcmp(DescriptionPlateCoordinates,'x')==1);

close all;

% 'blank'= medium without bacteria -> subtract from all other entries
% all blank data is averaged (over all positions, but individually for each
% time point)
if USERSETTINGS.showBlankFig
    figure(1); clf; hold on;
    title('time trace blanks')
    xlabel('time [h]')
    ylabel('OD 600nm')
end
[blank_row,blank_col]=find(strcmp(DescriptionPlateCoordinates,'blank')==1); %find all blank positions
numBlanks=length(blank_row); % number of blanks
avBlankPerTime=zeros(size(sortedData,1),1);
avPerBlank = []; stdPerBlank = []; plateauPerBlank = []; plateauStdPerBlank = [];% MW
for i=1:numBlanks
    wellName=PositionNames(blank_row(i),blank_col(i)); % e.g. 'A1'
    idx=find(strcmp(wellCoordinates,wellName)==1); % which entry in sortedData corresponds to wellName
    avBlankPerTime=avBlankPerTime+sortedData(idx).(yField);
    if USERSETTINGS.showBlankFig
        % dataseries of this blank sample
        plot(sortedData(idx).(timeField),sortedData(idx).(yField),'x','Color', 0.8*i/numBlanks*[1 1 1],'LineWidth',2)
        blankmaxvalues(i)=max(sortedData(idx).(yField));
        blanknames{i}=char(wellName);
    end
    % determine avg and std per blank
    avPerBlank = [avPerBlank mean(sortedData(idx).(yField))]; % MW
    stdPerBlank = [stdPerBlank std(sortedData(idx).(yField))]; % MW
    plateauPerBlank = [plateauPerBlank mean(sortedData(idx).(yField)(end-10:end))];
    plateauStdPerBlank = [plateauStdPerBlank std(sortedData(idx).(yField)(end-10:end))];
end
blankmaxvalues, blanknames
avBlankPerTime=avBlankPerTime/numBlanks;
if USERSETTINGS.showBlankFig
        % plot averages per timepoint
        plot(sortedData(1).(timeField),avBlankPerTime,'or','LineWidth',2)
                
        set(findall(gcf,'type','text'),'FontSize',FONTSIZE,'fontWeight','normal');
        set(gca,'FontSize',FONTSIZE);
        
        saveas(gcf,[myJustPlotDir currentdate '_blanksraw.png'], 'png');
        saveas(gcf,[myJustPlotDir currentdate '_blanksraw.eps'], 'epsc');
    
        % plot averages per blank
        figure(2)
        [h,hErrorbar]= barwitherr(stdPerBlank,avPerBlank,'LineWidth',1,'FaceColor',[.7 .7 .7]); % MW
        set(hErrorbar, 'LineWidth', 3) %MW added
        %errorbar(avPerBlank,stdPerBlank,'LineWidth',2) % MW
        axis([.5 numBlanks+.5 min(0-stdPerBlank) max(avPerBlank+stdPerBlank)*1.1])
        title('averages per blank'); xlabel('blank #'); ylabel('OD (600)')
        set(gca, 'XTickLabel',blanknames, 'XTick',1:numel(avPerBlank));
        
        set(findall(gcf,'type','text'),'FontSize',FONTSIZE,'fontWeight','normal');
        set(gca,'FontSize',FONTSIZE);
        
        saveas(gcf,[myJustPlotDir currentdate '_blanksmean.png'], 'png');
        saveas(gcf,[myJustPlotDir currentdate '_blanksmean.eps'], 'epsc');
        
        % plot plateaus per blank
        figure(3)
        [h,hErrorbar]= barwitherr(plateauStdPerBlank,plateauPerBlank,'LineWidth',1,'FaceColor',[.7 .7 .7]); % MW
        set(hErrorbar, 'LineWidth', 3) %MW added
        %errorbar(avPerBlank,stdPerBlank,'LineWidth',2) % MW
        axis([.5 numBlanks+.5 min(0-plateauStdPerBlank) max(plateauPerBlank+plateauStdPerBlank)*1.1])
        title('plateau values per blank (avg last 10)'); xlabel('blank #'); ylabel('OD (600)')
        set(gca, 'XTickLabel',blanknames, 'XTick',1:numel(avPerBlank));
        
        set(findall(gcf,'type','text'),'FontSize',FONTSIZE,'fontWeight','normal')
        set(gca,'FontSize',FONTSIZE)
        
        saveas(gcf,[myJustPlotDir currentdate '_blanksplateau.png'], 'png');
        saveas(gcf,[myJustPlotDir currentdate '_blanksplateau.eps'], 'epsc');
        plateauPerBlank
end
totalAvBlank=mean(avBlankPerTime); % maybe use complete average instead of av per time point. Todo

%add fields with blank subtracted to 'sortedData' (only useful for fields with actual
%bacteria in it, but never bothers)
for i=1:length(sortedData)
    sortedData(i).(yField_subtr) = sortedData(i).(yField)-avBlankPerTime;
end

%clear blank_row blank_col numBlanks i idx wellName

%% (3)
% ************************************************ 
% Plot all graphs grouped by category, without fitting
% ************************************************

currentPlateauFieldName = [yField 'Plateaus'];
currentPlateauFieldName_std = [yField 'Plateaus_std'];

% -----------------------------------------------
%loop over different groups in well (wellNames)
% -----------------------------------------------
myPlateauValues = []; myPlateauValues_std = [];
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
    title(name)
    hold on
    xlabel('time [hrs]')
    ylabel('OD600 (n.u.)')      
    if USERSETTINGS.hideGraphs
        set(h,'Visible','off'); % TODO MW - doesn't work
    end         
    % logplots
    hlog=figure(2); 
    clf    
    set(gca, 'Yscale', 'log')
    title(name)
    hold on
    xlabel('time [hrs]')
    ylabel('log_{2}(OD 600) (n.u.)')        
    % hide plots if desired  
    if USERSETTINGS.hideGraphs
        set(hlog,'Visible','off'); % TODO MW - doesn't work
    end     
    
    % -----------------------------------------------
    %  loop over all wells and search for matching description (same
    %  'name')
    % MW TODO: use membersOfGroup for this.
    % -----------------------------------------------
    myCurrentDataTime = []; myCurrentDataOD_substr = []; myCurrentColors = [];
    groupIdxs=[];
    for i=1:length(sortedData) 
        
        % if name matches current group
        if strcmp(sortedData(i).DescriptionPos,name)==1
            
            groupIdxs=[groupIdxs i];
            
            colorcounter=colorcounter+1;           
            
            % Plot the current well
            % Plot linear scale
            figure(h);
            plot(sortedData(i).(timeField),sortedData(i).(yField_subtr)','x','Color',myColor(colorcounter,:)','Linewidth',2);
            % Plot log scale
            figure(hlog); % also plot on logarithmic scale
            semilogy(sortedData(i).(timeField),sortedData(i).(yField_subtr)','x','Color',myColor(colorcounter,:)','Linewidth',2);
           
            % Determine moving averages for this group
            %[sortedData(i).movingAverage,sortedData(i).rangeMovingAverage]=movingaverage(sortedData(i).(yField_subtr),windowSize);
            [movingAverage,rangeMovingAverage]=movingaverage(sortedData(i).(yField_subtr),windowSize);
            sortedData(i).movingAverage = movingAverage;
            sortedData(i).rangeMovingAverage=rangeMovingAverage;
                       
            myColorNow=floor(myColor(colorcounter,:)*1.2); % quick n dirty edit color
            figure(h); plot(sortedData(i).(timeField)(rangeMovingAverage),movingAverage','-','Color',myColorNow,'Linewidth',2);
            figure(hlog); semilogy(sortedData(i).(timeField)(rangeMovingAverage),movingAverage','-','Color',myColorNow,'Linewidth',2);

            % Collect all data of this group in time and OD vector
            %myCurrentDataTime = [myCurrentDataTime sortedData(i).(timeField)];
            %myCurrentDataOD_substr = [myCurrentDataOD_substr sortedData(i).(yField_subtr)];         
            
            % Collect moving averages of this group to be able to determine
            % plateau value later.
            myCurrentDataOD_substr_movavg = [myCurrentDataOD_substr sortedData(i).(yField_subtr)];
            
            
            % Determine plateau value for this well 
            range = [ceil(size(sortedData(i).(yField_subtr),1)*PLATEAUSTART) size(sortedData(i).(yField_subtr),1)];
            if range(1) == 0, range(1)=1; end % this might be done prettier.. MW TODO
            sortedData(i).(currentPlateauFieldName)     = mean(sortedData(i).(yField_subtr)(range));
            sortedData(i).(currentPlateauFieldName_std) = std(sortedData(i).(yField_subtr)(range));
        end
    end
    % end loop over all data and search for repetitions with same 'name'
    % -----------------------------------------------
        
    % save with (moving) averages on linear scale
    figFullName=[myJustPlotDir currentdate 'GrowthCurves_' name];
    saveas(h,[figFullName '.fig'], 'fig');
    saveas(h,[figFullName '.png'], 'png');
    saveas(h,[figFullName '.eps'], 'epsc');
    % save with (moving) averages on log scale
    figFullName=[myJustPlotDir currentdate 'log_GrowthCurves_' name];
    saveas(hlog,[figFullName '.fig'], 'fig');
    saveas(hlog,[figFullName '.png'], 'png');
    saveas(hlog,[figFullName '.eps'], 'epsc');

    % close figures
    close(h); close(hlog);
    
    % Now determine the averages per group of the plateau values    
    myPlateauValues(nameidx) = mean([sortedData(groupIdxs).(currentPlateauFieldName)]);
    myPlateauValues_std(nameidx) = std([sortedData(groupIdxs).(currentPlateauFieldName)]); % std over means!
     
end
% and loop over all names (different exp's)
% -----------------------------------------------

% Save estimates of plateau values
h = figure();
%barh((currentPlateauFieldName));
barwitherr(myPlateauValues_std,myPlateauValues,'FaceColor',[0.8,0.8,0.8]);
xlabel('Strain/medium');
ylabel([yField ' value']);
title(['Plateau values determined from ' num2str(PLATEAUSTART) '-1.00 interval'])
set(gca, 'XTick', [1:length(wellNames)]);
set(gca, 'XTickLabel', wellNames);

set(findall(gcf,'type','text'),'FontSize',FONTSIZE,'fontWeight','normal');
set(gca,'FontSize',FONTSIZE);

% save with (moving) averages on linear scale
figFullName=[myJustPlotDir currentdate 'plateauvalues' ];
saveas(h,[figFullName '.fig'], 'fig');
saveas(h,[figFullName '.png'], 'png');

% Output plateau values to Excel file
filename = [myFullDir currentdate 'plateauvalues.xlsx'];
%myPlateauTable=table(wellNames,(currentPlateauFieldName)'; 
myPlateauTable=cell([wellNames,num2cell(myPlateauValues'),num2cell(myPlateauValues_std')])
xlswrite(filename,myPlateauTable,'Plateauvalues','B2');

clear dummy nameidx name muAccum muManualAccum
clear xlimfit ylimfit colorcounter mylegendText g h fitTimeManualext ODcalcManual  ODcalc
clear fitline figFullName ans currentColor fid i str SHOW_FIG_FIT ODmaxline ODminline

%% horizontal bar plot
% Save estimates of plateau values
h = figure();
barh(myPlateauValues);
ylabel('Strain/medium');
xlabel([yField ' value']);
title(['Plateau values determined from ' num2str(PLATEAUSTART) '-1.00 interval'])
set(gca, 'YTick', [1:length(wellNames)]);
set(gca, 'YTickLabel', wellNames);

set(findall(gcf,'type','text'),'FontSize',FONTSIZE,'fontWeight','normal');
set(gca,'FontSize',FONTSIZE);

%% (3b)
% ************************************************ 
% Plot all graphs in one plot 
% TODO MW: and also in subplot figure
% ************************************************

%create subSaveDirectory for these plots
myJustPlotDir=[myPlotsSaveDir SUMMARYPLOTDIRNAME];
if exist(myJustPlotDir)~=7
  [status,msg,id] = mkdir([myJustPlotDir]);
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
xlabel('time [hrs]')
ylabel(yField)      
if USERSETTINGS.hideGraphs
    set(h,'Visible','off'); % TODO MW - doesn't work
end         
% logplots
hlog=figure(2); 
clf    
set(gca, 'Yscale', 'log')
title(['log OD values over time'])
hold on
xlabel('time [hrs]')
ylabel(yField)        
% hide plots if desired  
if USERSETTINGS.hideGraphs
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
        lineh = plot(sortedData(i).(timeField)(rangeMA),sortedData(i).movingAverage',['-' myCurrentMarker],'Color',myColor(colorcounter,:)','Linewidth',1);
        % Plot log scale
        figure(hlog); % also plot on logarithmic scale
        linehlog = semilogy(sortedData(i).(timeField)(rangeMA),sortedData(i).movingAverage',['-' myCurrentMarker],'Color',myColor(colorcounter,:)','Linewidth',1);

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
figFullName=[myJustPlotDir currentdate 'All_GrowthCurves'];
saveas(h,[figFullName '.fig'], 'fig');
saveas(h,[figFullName '.png'], 'png');
% save with (moving) averages on log scale
figFullName=[myJustPlotDir currentdate 'All_Log_GrowthCurves'];
saveas(hlog,[figFullName '.fig'], 'fig');
saveas(hlog,[figFullName '.png'], 'png');

% close figures
%close(h); close(hlog);

%{
clear dummy nameidx name muAccum muManualAccum
clear xlimfit ylimfit colorcounter mylegendText g h fitTimeManualext ODcalcManual  ODcalc
clear fitline figFullName ans currentColor fid i str SHOW_FIG_FIT ODmaxline ODminline
%}


%% Tell user what to do

instructionString = ...
 ['------------------------------------------------------------------' 10 ...
  'Continue with the appropriate part II to ' 10 ...
  'finish the analysis: ' 10 ...
  '>> ExtractFitPlateReaderData_General_Part2_OD.m' 10 ...
  '>> ExtractFitPlateReaderData_General_Part2_Fluorescence.m' 10 ...
  '------------------------------------------------------------------' ];
disp(instructionString);











