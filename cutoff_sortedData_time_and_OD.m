

% Cuts out part of data
% Convenient when measurement was aborted during measurement run.

% To be executed inbetween steps (1) and (2) [i.e. before substr.
% background].

%cutoff=664
cutoff=664

for i = 1:length(sortedData)
    
    sortedData(i).time = sortedData(i).time(1:cutoff);
    sortedData(i).OD   =sortedData(i).OD(1:cutoff);
    
end
    
