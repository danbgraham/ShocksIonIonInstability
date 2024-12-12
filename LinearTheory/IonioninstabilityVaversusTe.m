% Routine to calculate peak growth rate and associated mode properties of
% the electrostatic proton-alpha streaming instability as a function of
% alpha speed and electron temperature. The routine checks all angles, so 
% given the number of Te and V_\alpha considered it is very slow. 
% Written by D. B. Graham.
%% Define arrays

Tearr = [0:5:200];
Vabarr = [0:5:250]*1e3;

Valength = length(Vabarr);
Telength = length(Tearr);

wmaxarr = zeros(Telength,Valength);
gmaxarr = zeros(Telength,Valength);
thetaarr = zeros(Telength,Valength);
kmaxarr = zeros(Telength,Valength);
vphmaxarr = zeros(Telength,Valength);
csarr = zeros(Telength,Valength);
ldarr = zeros(Telength,Valength);

%% Run main loop
wmaxarr(1,:) = NaN;
gmaxarr(1,:) = NaN;
thetaarr(1,:) = NaN;
vphmaxarr(1,:) = NaN;
kmaxarr(1,:) = NaN;
csarr(1,:) = NaN;
ldarr(1,:) = NaN;

tic;

for ii=2:1:Telength
  for jj = 1:1:Valength
    [wmaxwpp,gmaxwpp,thetakamax,kldmax,vphmax,wppval,csval,ldval] = getmaxparameters(Tearr(ii),Vabarr(jj)); 
    wmaxarr(ii,jj) = wmaxwpp;
    gmaxarr(ii,jj) = gmaxwpp;
    thetaarr(ii,jj) = thetakamax;
    kmaxarr(ii,jj) = kldmax;
    vphmaxarr(ii,jj) = vphmax;
    csarr(ii,jj) = csval;
    ldarr(ii,jj) = ldval;
    jj
  end
  ii
end

toc;

wpp = wppval;

%%

disprelTeVa = struct('wmaxarr',wmaxarr,'gmaxarr',gmaxarr,'thetaarr',thetaarr,'kmaxarr',kmaxarr,'vphmaxarr',vphmaxarr,...
  'csarr',csarr,'ldarr',ldarr,'wpp',wpp,'Tearr',Tearr,'Vabarr',Vabarr);

save('disprelTeVa.mat','disprelTeVa');