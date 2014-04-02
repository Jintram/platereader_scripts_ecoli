function [power x0] = DJK_ExponentialFit( t, y,myfitTime);
% Performs an exponential fit by linearization of log values with base 2.
% only takes values of t within myfitTime into account
%
%
idx=find(t>=myfitTime(1) & t<= myfitTime(2));
t_sub=t(idx);
y_sub=y(idx);

p = polyfit(t_sub,log2(y_sub),1);
power = p(1);
x0 = 2^p(2);
