function [power x0] = DJK_ExponentialFit( t, y,myfitTime);
% Performs an exponential fit by linearization of log values with base 2.
% only takes values of t within myfitTime into account
%
%
idx=find(t>=myfitTime(1) & t<= myfitTime(2));
t_sub=t(idx);
y_sub=y(idx);

log2_y_sub = log2(y_sub);

% apparently sometimes data gets processed to be imaginary - 
% prevent this:
if ~isreal(log2_y_sub)
    disp('Note: Converting imaginary data to real.');    
    log2_y_sub = real(log2_y_sub);
end

p = polyfit(t_sub,log2_y_sub,1);
power = p(1);
x0 = 2^p(2);
