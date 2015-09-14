% Example how to execute script.

% This example applies to a combined OD and fluor measurement, therefor
% TIMEINDEXES and ODINDEXES are set to non-default values, which correspond
% to the OD fields.
%
% This script processes the OD part of the data.

%% Dataset 2015-07-02, OD values

USERSETTINGS.myRootDir='U:\ZZ_EXPERIMENTAL_DATA\C_Platereader\';
USERSETTINGS.myDateDir='2015-07-02_FULL\';
USERSETTINGS.datafile='2015_07_02_CRPcAMP_plasmids_repeat3_FULL.xls';
USERSETTINGS.customSuffix = '_OD';
USERSETTINGS.ODorFluor = 1;

USERSETTINGS.ODmin=0.05;
USERSETTINGS.ODmax=0.10;

USERSETTINGS.fitManual = 0;

TIMEINDEXES=[7,9,11], ODINDEXES   = [8, 10, 12]

ExtractFitPlateReaderData_General_Part1
ExtractFitPlateReaderData_General_Part2_OD

%% Dataset 2015-07-02, GFP values

USERSETTINGS.myRootDir='U:\ZZ_EXPERIMENTAL_DATA\C_Platereader\';
USERSETTINGS.myDateDir='2015-07-02_FULL\';
USERSETTINGS.datafile='2015_07_02_CRPcAMP_plasmids_repeat3_FULL.xls';
USERSETTINGS.customSuffix = '_GFP';
USERSETTINGS.ODorFluor = 2;

USERSETTINGS.fitManual = 0; % for fluor: only if you want to REDO the manual range based on fluor

TIMEINDEXES=[5], ODINDEXES   = [6]

USERSETTINGS.wellNamesToPlot = ...
    {'lac838', 'lac852', 'lac853', 'lac854', 'lac855'};

ExtractFitPlateReaderData_General_Part1
ExtractFitPlateReaderData_General_Part2_Fluor

% see T:\TRAVELING_DATA\00_PLATEREADER\2015-07-02\Plots\LogPlot_manualRange for output.

