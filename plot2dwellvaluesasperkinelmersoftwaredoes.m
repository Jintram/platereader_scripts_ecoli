
% Visual representation of the 96 well plate for a certain point
% (TIMEPOINT) in time.
%
% Note that you need to run
% ExtractFitPlateReaderData_General_Part1 
% etc first.

% Time of which to show wells
TIMEPOINT=100;

figure; hold on;

% Max value to normalize data with
maxValue = max(max([sortedData(:).OD]));

% loop over 
for i=1:12 % columns of plate
    for j = 0:7 % rows of plate
               
        % Get value of it 
        currentWellValue = sortedData(i+j*12).OD(TIMEPOINT);
        % plot it
        plot(i,j,'o','Color',[ones(1,3).*(1-currentWellValue./maxValue)])
                
    end
end

