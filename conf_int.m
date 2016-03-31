function CI = conf_int(x,perc,add_mean,varargin)
% calculate confidence interval
%
% Input:
% x             data
% perc          percentage for CI as 0 < perc < 1 (default = 0.95)
% add_mean      add the mean to the output? 0 = no, 1 = yes (default = 1)
%
% Output:
% CI        confidence interval
%
% Author: 
% Frantisek Vasa (fv247@cam.ac.uk) - February 2016

% default settings for inputs
if nargin < 3; add_mean = 1; end    % CI relative to mean
if nargin < 2; perc = 0.95; end     % 95% confidence interval

SEM = std(x)/sqrt(length(x));                               % Standard Error
ts = tinv([(1-perc)/2  perc+(1-perc)/2],length(x)-1);       % T-Score
if add_mean == 1
    CI = mean(x) + ts*SEM;
elseif add_mean == 0
    CI = ts*SEM;
end

end