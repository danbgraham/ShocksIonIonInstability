function [P,x] = hist_dg(varargin)
%HIST_DG Compute 1D histogram
%
%   [P,x] = hist_dg(inp,'range',range,'nbins',nbins)
%
% Inputs:
%     inp - Array of data to make histogram
%
% Options:
%     range - two element array of min/max value [min max]
%     nbins - number of bins
% 
% Written by D. B. Graham

if nargin == 0
    help hist_dg;
    return;
end

data=varargin{1};
args=varargin(2:end);
if numel(args)>0
    flag_have_options=1;
else
    flag_have_options=0;
end

% set default values
nbins = 20;
range = [min(data) max(data)];

while flag_have_options
    l = 2;
    switch(lower(args{1}))
        case 'nbins'
            if numel(args)>1 && isnumeric(args{2})
                nbins = args{2};
            end
        case 'range'
            if numel(args)>1 && isnumeric(args{2})
                range = args{2};
            end
        otherwise
            irf_log('fcal',['Unknown flag: ' args{1}])
            l=1;
            break
    end
    args = args(l+1:end);
    if isempty(args), flag_have_options=0; end
end



xtemp = linspace(range(1), range(2), nbins+1);
dx = median(diff(xtemp))/2;
x = xtemp(1:end-1)+dx;
P = zeros(1,length(x));

for ii = 1:length(x)
    P(ii) = numel(find(data >= xtemp(ii) & data < xtemp(ii+1)));
end

end