% Example how to execute script.

% This example applies to a combined OD and fluor measurement, therefor
% TIMEINDEXES and YINDEXES are set to non-default values, which correspond
% to the OD fields.
%

% Last update: MW, 2015/11

% =========================================================================

% Example of how script can be executed:

%% Part 1: OD values

USERSETTINGS.myRootDir='U:\ZZ_EXPERIMENTAL_DATA\C_Platereader\';
USERSETTINGS.myDateDir='2015_11_23\';
USERSETTINGS.datafile= '2015_11_23_CRP_oldreader.xls';
USERSETTINGS.customSuffix = '_OD';
USERSETTINGS.ODorFluor = 1;

USERSETTINGS.ODmin=0.05;
USERSETTINGS.ODmax=0.10;

USERSETTINGS.fitManual = 0;

% From which column to get values, ATTENTION! This changes for old/new
% platereader.
% All columns with data, either OD or fluor:   
%    TIMEINDEXES=[5,7,9,11], YINDEXES   = [6,8,10,12]
% One of them is fluor usually, three others OD.
%    TIMEINDEXES=[7,9,11], YINDEXES   = [8,10,12] % new reader
%    TIMEINDEXES=[5,7,9], YINDEXES   = [6,8, 10] % old
% reader
% Different for old platereader protocol!
TIMEINDEXES=[7,9,11], YINDEXES   = [8,10,12] % values for old reader (open excel file w. data to check)

ExtractFitPlateReaderData_General_Part1
ExtractFitPlateReaderData_General_Part2_OD
  
%% Part 2: Repeat analysis to add GFP values.

USERSETTINGS.customSuffix = '_GFP';
USERSETTINGS.ODorFluor = 2;
USERSETTINGS.platereader = 'OLD';

% for fluor: only if you want to REDO the manual range based on fluor
USERSETTINGS.fitManual = 0; 

% Different for old platereader protocol!
%TIMEINDEXES=[5], YINDEXES   = [6] % new platereader
TIMEINDEXES=[11], YINDEXES   = [12] % old platereader

ExtractFitPlateReaderData_General_Part1
ExtractFitPlateReaderData_General_Part2_Fluor
% ExtractFitPlateReaderData_General_Part3_Plotting.m
%========================================================================

%% Part 3: Plot specific wells if desired
DONTSAVE=1;
USERSETTINGS.wellNamesToPlot = ...
{ 'aTc00000s926r1','aTc40000s926r1','aTc11599s926r1','aTc03363s926r1','aTc00975s926r1','aTc00282s926r1','aTc00082s926r1','aTc00024s926r1','aTc00007s926r1'};
ExtractFitPlateReaderData_General_Part2_Fluor


