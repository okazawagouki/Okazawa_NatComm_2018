
function [version, release, date]  = getMatlabVersion
% 
% function [version, release, date]  = getMatlabVersion
% 
% get Matlab version
% Outputs:
%   version         matlab version, for example 8.5 for Matlab2015a
%   release         matlab release, for example '2015a'
%   date            release date in datenum format, use datestr to convert 
% 

% RK, 5/19/2015

version = [];
release = [];
date = [];

v = ver;
for i = 1 : length(v)
    if strcmpi(v(i).Name,'MATLAB')
        version = str2double(v(i).Version);
        release = strtok(v(i).Release,'()');
        date = datenum(v(i).Date);
    end
end



