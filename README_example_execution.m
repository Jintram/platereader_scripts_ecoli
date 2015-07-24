% Example how to execute script.

% This example applies to a combined OD and fluor measurement, therefor
% TIMEINDEXES and ODINDEXES are set to non-default values, which correspond
% to the OD fields.
%
% This script processes the OD part of the data.

% Dataset 2015-07-02

USERSETTINGS.myRootDir='U:\ZZ_EXPERIMENTAL_DATA\C_Platereader\';
USERSETTINGS.myDateDir='2015-07-02_FULL\';
USERSETTINGS.datafile='2015_07_02_CRPcAMP_plasmids_repeat3_FULL.xls';
USERSETTINGS.customSuffix = 'OD';

USERSETTINGS.ODmin=0.05;
USERSETTINGS.ODmax=0.22;

USERSETTINGS.fitManual = 1;

TIMEINDEXES=[7,9,11], ODINDEXES   = [8, 10, 12]

ExtractFitPlateReaderData_General_Part1
ExtractFitPlateReaderData_General_Part2_OD
