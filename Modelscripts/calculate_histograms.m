function histstruct = calculate_histograms(varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%   Input: 
%     xrange -  range of positions normal to the shock
%     dx -      spacing between positions in normal direction
%     vxrange - range of velocities considered in the direction normal to
%               the shock
%     vyrange - range of velocities considered in the direction aligned
%               with B component
%     vzrange - range of velocities considered in the direction
%               perpendicular to B and the normal
%     dv -      velocity spacing of calculated distributions
%     xposarr - arrays of particle positions
%     vxarr -   arrays of 
%
%

if (nargin < 10)
  help calculate_histograms;
  histstruct = NaN;
  return;
end

xrange = varargin{1};
dx = varargin{2};
vxrange = varargin{3};
vyrange = varargin{4};
vzrange = varargin{5};
dv = varargin{6};
xposarr = varargin{7};
vxarr = varargin{8};
vyarr = varargin{9};
vzarr = varargin{10};

if (nargin == 11)
  weights = varargin{11};
else
  weights = ones(size(vxarr));
end

weightssize = size(weights);
if weightssize(2) == 1
  numtimes = length(vxarr(1,:));
  timesones = ones(1,numtimes);
  weights = weights*timesones;
end

% Define domain
xpositions = xrange(1)+dx/2:dx:xrange(2);
numpositions = length(xpositions);
  
% Define ranges
vxpositions = vxrange(1)+dv/2:dv:vxrange(2);
vypositions = vyrange(1)+dv/2:dv:vyrange(2);
vzpositions = vzrange(1)+dv/2:dv:vzrange(2);
numvxs = length(vxpositions);
numvys = length(vypositions);
numvzs = length(vzpositions);
  
n1Dvx = zeros(numpositions,numvxs);
n1Dvy = zeros(numpositions,numvys);
n1Dvz = zeros(numpositions,numvzs);
  
n2Dvxvy = zeros(numpositions,numvxs,numvys);
n2Dvxvz = zeros(numpositions,numvxs,numvzs);
n2Dvyvz = zeros(numpositions,numvys,numvzs);
  
  for ii = 1:numpositions
    % 1D reduced histograms
    idxx = xposarr > xpositions(ii)-dx/2 & xposarr < xpositions(ii)+dx/2; 
    vxarrdx = vxarr(idxx);
    vyarrdx = vyarr(idxx);
    vzarrdx = vzarr(idxx);
    weightsx = weights(idxx);
    
    for jj=1:numvxs
      % 1D vx
      idxvx = vxarrdx > vxpositions(jj)-dv/2 & vxarrdx < vxpositions(jj)+dv/2;
      weightsx1Dvx = weightsx(idxvx);
      n1Dvx(ii,jj) = sum(weightsx1Dvx);
      vyarrdx2 = vyarrdx(idxvx);
      vzarrdx2 = vzarrdx(idxvx);
      % 2D vxvy
      for mm=1:numvys
        idxvxvy = vyarrdx2 > vypositions(mm)-dv/2 & vyarrdx2 < vypositions(mm)+dv/2;
        weightsx2Dvxvy = weightsx1Dvx(idxvxvy);
        n2Dvxvy(ii,jj,mm) = sum(weightsx2Dvxvy);
      end
      % 2D vxvz
      for nn=1:numvzs
        idxvxvz = vzarrdx2 > vzpositions(nn)-dv/2 & vzarrdx2 < vzpositions(nn)+dv/2;
        weightsx2Dvxvz = weightsx1Dvx(idxvxvz);
        n2Dvxvz(ii,jj,nn) = sum(weightsx2Dvxvz);
      end
    end
    % 1D vy
    for kk=1:numvys
      idxvy = vyarrdx > vypositions(kk)-dv/2 & vyarrdx < vypositions(kk)+dv/2;
      weightsx1Dvy = weightsx(idxvy);
      n1Dvy(ii,kk) = sum(weightsx1Dvy);
      vzarrdx2 = vzarrdx(idxvy);
      % 2D vyvz
      for oo=1:numvzs
        idxvyvz = vzarrdx2 > vzpositions(oo)-dv/2 & vzarrdx2 < vzpositions(oo)+dv/2;
        weightsx2Dvyvz = weightsx1Dvy(idxvyvz);
        n2Dvyvz(ii,kk,oo) = sum(weightsx2Dvyvz);
      end
    end
    % 1D vz
    for ll=1:numvzs
      idxvz = vzarrdx > vzpositions(ll)-dv/2 & vzarrdx < vzpositions(ll)+dv/2;
      weightsx1Dvz = weightsx(idxvz);
      n1Dvz(ii,ll) = sum(weightsx1Dvz);
    end
  end
  

  histstruct = struct('n1Dvx',n1Dvx,'n1Dvy',n1Dvy,'n1Dvz',n1Dvz,'n2Dvxvy',n2Dvxvy,'n2Dvxvz',n2Dvxvz,'n2Dvyvz',n2Dvyvz,... 
    'xpositions',xpositions,'vxpositions',vxpositions,'vypositions',vypositions,'vzpositions',vzpositions);

end

