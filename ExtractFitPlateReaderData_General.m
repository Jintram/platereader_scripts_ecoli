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


%% (1)
% ************************************************
% specify folder, date, etc
% ************************************************
myRootDir='U:\PROJECTS\Temperature_Mutants\platereader\';
myDateDir='2014_03_29\';
datafile='2014_03_29_results_temperature_mutants';

% ************************************************
% import data. names: 'data' 'textdata'. create saveDirectory
% ************************************************
myFullDir=[myRootDir myDateDir];
[data textdata]=xlsread([myFullDir, datafile, '.xls']);
[~, DescriptionPlateCoordinates]=xlsread([myFullDir 'plateCoordinates.xlsx'],'C4:N11'); %what is tested on which position
load([myRootDir 'PositionNames.mat']); % cell array with 'A1' 'B1' etc
%
% nb: the timefield is now in IEEE format and can be converted into minutes
% by DJK_getMinutesfromTimestamp(timeIeee). Timefield will be converted to
% hours further below


% Data structure: (X=not important)
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
    [idx_row,idx_col]=find(strcmp(PositionNames,wellCoordinates(i))==1);
    sortedData(i).DescriptionPos=char(DescriptionPlateCoordinates(idx_row,idx_col));
    realData=(strcmp(sortedData(i).DescriptionPos,'blank')==0 & strcmp(sortedData(i).DescriptionPos,'x')==0);
    sortedData(i).realData=realData;
    sortedData(i).muManual=[];
    sortedData(i).x0Manual=[];
    sortedData(i).mu=[];
    sortedData(i).x0=[];
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

%----------------------------------------------------------
%create save directory
myPlotsSaveDir=[myFullDir 'Plots\'];
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
avPerBlank = [], stdPerBlank = []; % MW
for i=1:numBlanks
    wellName=PositionNames(blank_row(i),blank_col(i)); % e.g. 'A1'
    idx=find(strcmp(wellCoordinates,wellName)==1); % which entry in sortedData corresponds to wellName
    avBlankPerTime=avBlankPerTime+sortedData(idx).OD;
    if SHOW_FIG_BLANK
        plot(sortedData(idx).time,sortedData(idx).OD,'x','Color', 0.8*i/numBlanks*[1 1 1])
    end
    avPerBlank = [avPerBlank mean(sortedData(idx).OD)]; % addition MW
    stdPerBlank = [stdPerBlank std(sortedData(idx).OD)]; % addition MW
end
avBlankPerTime=avBlankPerTime/numBlanks;
if SHOW_FIG_BLANK
        % plot averages per timepoint
        plot(sortedData(1).time,avBlankPerTime,'xr')
        % plot averages per blank
        figure(2)
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
%choose fitTime according to OD thresholds
% ************************************************
ODmin=0.03; ODmax=0.08; % does not take into account sudden random umps over threshold (e.g. avoid by averaging)

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
        idxODminOrHigher=find(sortedData(i).OD_subtr>ODmin);
        idxODmaxOrLower=find(sortedData(i).OD_subtr<ODmax);
        idxMin=min(idxODminOrHigher);
        idxMax=max(idxODmaxOrLower);
        sortedData(i).fitTime=[sortedData(i).time(idxMin), sortedData(i).time(idxMax)];
    else sortedData(i).fitTime=[];
    end
end
clear i idxODmaxOrLower idxODminOrHigher idxMin idxMax status msg id

%% (4) optional
% ************************************************
% loop over all wells and choose fitTimeManual
% if 'x','blank' or bad data, fitTimeManual is set to =0
% ************************************************
% this cell is older and because of that a little akwardly written. it is
% correct but not written in the most optimal way
%
stopIt=0;

startindex=input('Startindex=');
disp('enter 2 numbers seperated by space for each fitTimeManual. Use 0 (1 single number) if data is')
disp('bad and shall be omitted. Pressing space (empty fitTimeManual) is interpreted as not reviewed.')
disp('To stop, press ''q''. If you want to keep old fitTimeManual, press ''f''.') %f=forward 
for row=1:size(DescriptionPlateCoordinates,1)
    for col=1:size(DescriptionPlateCoordinates,2)
        idx=find(strcmp(wellCoordinates,PositionNames(row,col))==1);
        if idx>=startindex & stopIt==0
            if strcmp(DescriptionPlateCoordinates(row,col),'x')==1 | strcmp(DescriptionPlateCoordinates(row,col),'blank')==1
                sortedData(idx).fitTimeManual=0;
            else
                fitTimeManualstr='review';
                while strcmp(fitTimeManualstr,'f')==0  & stopIt==0
                    figure(2)
                    clf
                    semilogy(sortedData(idx).time*60,sortedData(idx).OD_subtr)
                    hold on
                    title(DescriptionPlateCoordinates(row,col))
                    xlabel('time [min]')
                    ylabel('OD')
                    grid on
                    %semilogy(sortedData(idx).time*60,sortedData(idx).OD,'r')
                    %if fitTimeManual already exists, plot it
                    y_lim=get(gca,'ylim');
                    if length(sortedData(idx).fitTimeManual)==2 % not foolproof                       
                        plot([sortedData(idx).fitTimeManual(1)*60 sortedData(idx).fitTimeManual(1)*60],y_lim,'-k','LineWidth',2);
                        plot([sortedData(idx).fitTimeManual(2)*60 sortedData(idx).fitTimeManual(2)*60],y_lim,'-k','LineWidth',2);
                    end
                    % plot alternative fitTimeManual according to OD thresholds
                    plot([sortedData(idx).fitTime(1)*60 sortedData(idx).fitTime(1)*60],y_lim,'-m');
                    plot([sortedData(idx).fitTime(2)*60 sortedData(idx).fitTime(2)*60],y_lim,'-m');
                    

                    fitTimeManualstr=input(['sortedData('  num2str(idx)  ...
                        '). fitTimeManual='],'s');
                    if strcmp(fitTimeManualstr,'q')==1
                        stopIt=1;
                        continue
                    end
                    if strcmp(fitTimeManualstr,'f')~=1
                        fitT=str2num(fitTimeManualstr);
                        fitT=fitT/60; %convert to hours
                        % invert if fitTimeManual(1)>fitTimeManual(2). check if fitTimeManual
                        % consists of 2 entries
                        if (fitT~=0 & ~isempty(fitT)) & (length(fitT)~=2 | fitT(1)==fitT(2))
                            disp(['redo fitTimeManual'])
                            fitT=[];
                        elseif (fitT~=0 & ~isempty(fitT)) & fitT(1)>fitT(2)
                            disp(['inverted start & end of fitTimeManual'])
                            dummy=fitT(2); fitT(2)=fitT(1); fitT(1)=dummy;
                        end
                        sortedData(idx).fitTimeManual=fitT;

                    end
                end
            end
        end
    end
end
clear fitTimeManualstr fitT col fitT fitTimeManualstr idx row startindex y_lim stopIt

%% (5) optional (use if (4) was used, otherwise ignore)
% ************************************************
% check if some fitTimeManuals are accidently empty
% ************************************************
for i=1:length(sortedData)
    if sortedData(i).realData==1 & isempty(sortedData(i).fitTimeManual)
        disp(['fitTimeManual of sortedData(' num2str(i) ') is missing']);
    end
end         
clear i

%% (6)
% ************************************************
% fit growth rate to each graph
% ************************************************
for i=1:length(sortedData)
    % fit growth rate according to fitTimeManual
    if (sortedData(i).realData==1 & length(sortedData(i).fitTimeManual)==2) %bad data has '0' as entry in fitTimeManual
        [muManual,x0Manual]=NW_ExponentialFit_fitTime(sortedData(i).time,sortedData(i).OD_subtr,sortedData(i).fitTimeManual);
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
            [mu,x0]=NW_ExponentialFit_fitTime(sortedData(i).time,sortedData(i).OD_subtr,sortedData(i).fitTime);
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


SHOW_FIG_FIT=1;
SHOW_FIG_FITMANUAL=0; % no options to save! not kept uptodate!

myColor=[1 0 0 ; 0 1 0 ; 0 0 1; 1 0.6 0.2;  0 1 1; 0 0.5 0.5; 0 0.6 0; 0.6 0 0.4; 0.8 0.5 0; 0.7 0 0; 0.4 0.2 0.6; 0.6 0.2 0; 1 0 0 ; 0 1 0 ; 0 0 1; 1 0.6 0.2;  0 1 1; 0 0.5 0.5; 0 0.6 0; 0.6 0 0.4; 0.8 0.5 0; 0.7 0 0; 0.4 0.2 0.6; 0.6 0.2 0];
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

%save 'sortedData' and growthrate data  'muAvStdev'
save([myFullDir 'CompleteAnalyzedData.mat'],'sortedData','muAvStdev');

clear dummy nameidx name muAccum muManualAccum
clear xlimfit ylimfit colorcounter mylegendText g h fitTimeManualext ODcalcManual  ODcalc
clear fitline figFullName ans currentColor fid i str SHOW_FIG_FIT ODmaxline ODminline

%% (7b)----------------------------
% -------------------------
% Plot half-logarithmic plots to check range of exponential phase
% -------------------------
% Manual Fit range not implemented in log-plots!
SHOW_FIG_FIT=1;
myColor=[1 0 0 ; 0 1 0 ; 0 0 1; 1 0.6 0.2;  0 1 1; 0 0.5 0.5; 0 0.6 0; 0.6 0 0.4; 0.8 0.5 0; 0.7 0 0; 0.4 0.2 0.6; 0.6 0.2 0; 1 0 0 ; 0 1 0 ; 0 0 1; 1 0.6 0.2;  0 1 1; 0 0.5 0.5; 0 0.6 0; 0.6 0 0.4; 0.8 0.5 0; 0.7 0 0; 0.4 0.2 0.6; 0.6 0.2 0];
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
                    plot(sortedData(i).time,sortedData(i).OD_subtr,'o','Color',currentColor,'Markersize',3);
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
        %create save Directory for log images
        myPlotsSaveDirLogODsub=[myPlotsSaveDir 'LogPlot_ODmin'  num2str(ODmin) 'ODmax' num2str(ODmax) '\'];
        if exist(myPlotsSaveDirLogODsub)~=7
        [status,msg,id] = mymkdir([myPlotsSaveDirLogODsub]);
            if status == 0
                disp(['Warning: unable to mkdir ' myPlotssSaveDirLogODSub ' : ' msg]);
                return;
            end
        end
        
        eval(['legend(', mylegendText, ',''Location'',''NW'')']);
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

clear dummy SHOW_FIG_FIT SHOW_FIG_FITMANUAL myColor nameidx name muAccum muManualAccum
clear xlimfit ylimfit colorcounter mylegendText g h fitTimeManualext ODcalcManual fitTimeext ODcalc
clear fitline figFullName ans currentColor fid i str ODmaxline ODminline



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
