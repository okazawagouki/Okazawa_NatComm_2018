
function label = makeTickLabel(YTicks, Interval)

% 
% function label = makeTickLabel(YTicks, Interval)
% 
% creates tick labels by making a cell array of strings of size YTicks. the values of
% cell array are empty unless they correspond to the specified Interval
% Example:
%   makeTickLabel(100:100:1000, 200) returns
%   {'100', '', '300', '', '500', '', '700', '', '900', ''}
% 

label = cell(1, length(YTicks));
    %mod does not function well, if it was we could just use the following line
    % L = find(mod(YTicks-YTicks(1),Interval) < 1e-6);
    %instead we have to rewrite the code for mod
d = (YTicks-YTicks(1))/Interval;
L = find(abs(d-round(d))<1e-6);
YTicks(abs(YTicks)<1e-6) = 0;
for i = L
    label{i} = num2str(YTicks(i));
end




