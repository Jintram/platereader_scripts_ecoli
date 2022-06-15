function [movingaverage,xrange] = movingaverage(y,windowSize)
% Determine moving average for vector x, using a moving window of 
% windowsize length.
% Edited from: http://stackoverflow.com/questions/3453663/computing-running-averages-in-matlab
%
% MW 2014

% Create vector with weigths to determine moving average
F = ones(1,windowSize)/windowSize;

% Determine moving average using convolution
movingaverage = conv(y,F,'valid');

% Also return range to use for x
halfSize = floor(windowSize/2);
xrange = [halfSize+1:length(y)-halfSize]';

% Throw away exess points - only required if not used 'valid' option in conv()
%movingaverage = movingaverage(halfSize+1:end-halfSize);

end