function histstruct = calculate_hist2D(xrange,dx,vxrange,dv,xposarr,vxarr)
% calculate_hist2D Routine used for calculating 2D histrograms for x-v_x
% distribution plots
%
%
%
% Written by D. B. Graham

  % Define domain
  xpositions = xrange(1)+dx/2:dx:xrange(2);
  numpositions = length(xpositions);
  
  % Define ranges
  vxpositions = vxrange(1)+dv/2:dv:vxrange(2);
  numvxs = length(vxpositions);
  
  n1Dvx = zeros(numpositions,numvxs);
  
  for ii = 1:numpositions
    % 1D reduced histograms
    idxx = xposarr > xpositions(ii)-dx/2 & xposarr < xpositions(ii)+dx/2; 
    vxarrdx = vxarr(idxx);
    
    for jj=1:numvxs
      % 1D vx
      idxvx = vxarrdx > vxpositions(jj)-dv/2 & vxarrdx < vxpositions(jj)+dv/2;
      n1Dvx(ii,jj) = sum(idxvx);
    end
  end
  

  histstruct = struct('n1Dvx',n1Dvx','xpositions',xpositions,'vxpositions',vxpositions);

end

