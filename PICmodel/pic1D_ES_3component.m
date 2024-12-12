% Simple 1D ES PIC code used to model the proton-alpha streaming
% instability. 
%
% Written by D. B. Graham
%% Define initial conditions
Units = irf_units;
qe = Units.e;
me = Units.me;
mp = Units.mp;
ma = mp*4;
eps0 = Units.eps0;

qeme = -qe/me;
qemp = qe/mp;
qema = 2*qe/ma;


% Define red-blue color scale used below
rr = interp1([1 64 128 192 256],[0.0  0.5 1 1.0 0.75],1:256);
gg = interp1([1 64 128 192 256],[0.0  0.5 1 0.5 0.00],1:256);
bb = interp1([1 64 128 192 256],[0.75 1.0 1 0.5 0.00],1:256);
bgrcmap = [rr' gg' bb'];

% Define initial conditions 
na = 1e6; % m^-3
np = 10e6; % m^-3
ne = np + 2*na; % m^-3


Te = 100; % eV
Tp = 3; % eV
Ta = 12; % eV

ve = sqrt(2*qe*Te/me);
va = sqrt(2*qe*Ta/ma);
vp = sqrt(2*qe*Tp/mp);

wpe = sqrt(qe^2*(ne)/me/eps0);
wpp = sqrt(qe^2*(np)/mp/eps0);
wpa = sqrt(4*qe^2*(na)/ma/eps0);
ld = ve/wpe/sqrt(2);

Vbulke =  0*ve;
Vbulka = 100e3;
Vbulkp = 0;

vewidth = [-6 6];
vpwidth = [-10 10];
vawidth = [-10 10];

%% Set up
% Simulation parameters
numpts = 128; % Number of points
ppc = 20000; % particles per cell
numparticles = numpts*ppc;

dt = 0.05/wpe;

Ntimes = 240000;

nn = [1:numparticles 1:numparticles];

dx = ld/2;
xlength = dx*numpts;


KEe = zeros(1,Ntimes);
KEa = zeros(1,Ntimes);
KEp = zeros(1,Ntimes);
KEpa = zeros(1,Ntimes);
EE = zeros(1,Ntimes);
Etot = zeros(1,Ntimes);

Exall = zeros(Ntimes,numpts);
phiall = zeros(Ntimes,numpts);
rhoall = zeros(Ntimes,numpts);

Mome = zeros(1,Ntimes);
Moma = zeros(1,Ntimes);
Momp = zeros(1,Ntimes);
Momtot = zeros(1,Ntimes);

rng('default')
vxe = randn(numparticles,1)*ve/sqrt(2)+Vbulke;
xe = linspace(0,1,numparticles)'*xlength;

rng('default')
vxp = randn(numparticles,1)*vp/sqrt(2)+Vbulkp;
xp = linspace(0,1,numparticles)'*xlength;

rng('default')
vxa = randn(numparticles,1)*va/sqrt(2)+Vbulka;
xa = linspace(0,1,numparticles)'*xlength;

vxeinit = vxe;
vxpinit = vxp;
vxainit = vxa;

% 1D distributions
binnum1ddist = 300;
[Peinit,vxeinit] = hist_dg(vxeinit/ve,'range',vewidth,'nbins',binnum1ddist);
[Ppinit,vxpinit] = hist_dg(vxpinit/vp,'range',vpwidth,'nbins',binnum1ddist);
[Painit,vxainit] = hist_dg(vxainit/va,'range',vawidth,'nbins',binnum1ddist); % Use same binning for protons and alphas
dve1d = median(diff(vxeinit));
dvp1d = median(diff(vxpinit));
dva1d = median(diff(vxainit));

edist1dall = zeros(length(vxeinit),Ntimes);
pdist1dall = zeros(length(vxpinit),Ntimes);
adist1dall = zeros(length(vxainit),Ntimes);

% Set up matrices to store particle trajectories
mididx = floor((ppc*numpts/2));
electronpos = zeros(100,Ntimes);
electronspeed = zeros(100,Ntimes);

protonpos = zeros(100,Ntimes);
protonspeed = zeros(100,Ntimes);

alphapos = zeros(100,Ntimes);
alphaspeed = zeros(100,Ntimes);

xvalsall = (dx/2:dx:xlength)/ld;
timesall = (dt/2:dt:dt*Ntimes)*wpe;

% Set 1D reduced electron distributions
distindices = 0:100:Ntimes;
distindices(1) = 1;
dxx = dx/ld; % in Ld
dvv = 0.1; % 0.1*vth
xrange = [0 xlength/ld];
xpositions = xrange(1)+dxx/2:dxx:xrange(2);
numpositions = length(xpositions);
vxepositions = vewidth(1)+dvv/2:dvv:vewidth(2);
vxppositions = vpwidth(1)+dvv/2:dvv:vpwidth(2);
vxapositions = vawidth(1)+dvv/2:dvv:vawidth(2);
numvxes = length(vxepositions);
numvxps = length(vxppositions);
numvxas = length(vxapositions);
edistmatrix = zeros(length(distindices),length(vxepositions),length(xpositions));
pdistmatrix = zeros(length(distindices),length(vxppositions),length(xpositions));
adistmatrix = zeros(length(distindices),length(vxapositions),length(xpositions));
kk = 1;

%% Main simulation loop
for ii = 1:Ntimes
  % Update position  
  xe = xe + vxe*dt;
  xa = xa + vxa*dt;
  xp = xp + vxp*dt;

  % Apply boundary conditions
  xe(xe < 0) = xe(xe < 0)+xlength; xe(xe > xlength) = xe(xe > xlength)-xlength;
  xa(xa < 0) = xa(xa < 0)+xlength; xa(xa > xlength) = xa(xa > xlength)-xlength; 
  xp(xp < 0) = xp(xp < 0)+xlength; xp(xp > xlength) = xp(xp > xlength)-xlength;

  % Compute rho, phi, and E from particle positions
    deltax = xp/dx;
    positemp = floor(deltax+0.5);
    deltax = deltax-positemp+0.5;
    idx1 = positemp; idx2 = positemp+1;
    idx1(idx1 == 0) = numpts; idx2(idx2 == numpts+1) = 1;
    matp=sparse(nn,[idx1; idx2],[deltax; 1-deltax],numparticles,numpts);
    npx = full(sum(matp,1))/ppc*np;
    
    deltax = xe/dx;
    positemp = floor(deltax+0.5);
    deltax = deltax-positemp+0.5;
    idx1 = positemp; idx2 = positemp+1;
    idx1(idx1 == 0) = numpts; idx2(idx2 == numpts+1) = 1;
    mate=sparse(nn,[idx1; idx2],[deltax; 1-deltax],numparticles,numpts);
    nex = full(sum(mate,1))/ppc*ne;
    
    deltax = xa/dx;
    positemp = floor(deltax+0.5);
    deltax = deltax-positemp+0.5;
    idx1 = positemp; idx2 = positemp+1;
    idx1(idx1 == 0) = numpts; idx2(idx2 == numpts+1) = 1;
    mata=sparse(nn,[idx1; idx2],[deltax; 1-deltax],numparticles,numpts);
    nax = full(sum(mata,1))/ppc*na;

    rhox = qe*(npx + 2*nax - nex);
    
    % Compute E and phi in wave number space
    rhoxfft = fft(rhox);
    
    kmax = pi/dx;
    deltak = 2*kmax/numpts;
    
    kvec = (-kmax:deltak:kmax-deltak);
    kvecs = ifftshift(kvec);

    phik = rhoxfft./kvecs.^2;
    phik(1) = 0;
    Ek = -1i*kvecs.*phik/eps0;
    
    
    phix = real(ifft(phik)/eps0);
    
    Ex = real(ifft(Ek));

  % Apply acceleration to particles
  vxe = vxe + mate*Ex'*dt*qeme;
  vxp = vxp + matp*Ex'*dt*qemp;
  vxa = vxa + mata*Ex'*dt*qema;

  KEe(ii) = sum(0.5*me*vxe.^2)*ne/ppc*dx;
  KEp(ii) = sum(0.5*mp*vxp.^2)*np/ppc*dx;
  KEa(ii) = sum(0.5*ma*vxa.^2)*na/ppc*dx;

  KEpa(ii) = KEe(ii)+KEa(ii)+KEp(ii);
    
  EE(ii) = 0.5*eps0*sum(Ex.^2)*dx;

  Etot(ii) = KEpa(ii) + EE(ii);

  Exall(ii,:) = Ex;
  phiall(ii,:) = phix;
  rhoall(ii,:) = rhox;

    Mome(ii) = sum(me*vxe)*ne/ppc*dx;
    Moma(ii) = sum(ma*vxa)*na/ppc*dx;
    Momp(ii) = sum(mp*vxp)*np/ppc*dx;
    Momtot(ii) = Mome(ii)+Moma(ii)+Momp(ii);

    
    % Store electron distributions
    if ismember(ii,distindices)
      hist2Detemp = calculate_hist2D(xrange,dx/ld,vewidth,dvv,xe/ld,vxe/ve);
      hist2Dptemp = calculate_hist2D(xrange,dx/ld,vpwidth,dvv,xp/ld,vxp/vp);
      hist2Datemp = calculate_hist2D(xrange,dx/ld,vawidth,dvv,xa/ld,vxa/vp);
      
      e1Dreducedtemp = hist2Detemp.n1Dvx;
      e1Dreducedtemp(e1Dreducedtemp == 0) = NaN;
      normfactore = sum(sum(e1Dreducedtemp,'omitnan'),'omitnan')/numpts;
      e1Dreducedtemp = e1Dreducedtemp*ne/normfactore/(dvv*ve);
      edistmatrix(kk,:,:) = e1Dreducedtemp;
      
      p1Dreducedtemp = hist2Dptemp.n1Dvx;
      p1Dreducedtemp(p1Dreducedtemp == 0) = NaN;
      normfactorp = sum(sum(p1Dreducedtemp,'omitnan'),'omitnan')/numpts;
      p1Dreducedtemp = p1Dreducedtemp*np/normfactorp/(dvv*vp);
      pdistmatrix(kk,:,:) = p1Dreducedtemp;
      
      a1Dreducedtemp = hist2Datemp.n1Dvx;
      a1Dreducedtemp(a1Dreducedtemp == 0) = NaN;
      normfactora = sum(sum(a1Dreducedtemp,'omitnan'),'omitnan')/numpts;
      a1Dreducedtemp = a1Dreducedtemp*na/normfactora/(dvv*va);
      adistmatrix(kk,:,:) = a1Dreducedtemp;
      kk = kk+1;
    end

    % Store 1D distributions
    [Petemp,~] = hist_dg(vxe/ve,'range',vewidth,'nbins',binnum1ddist);
    [Pptemp,~] = hist_dg(vxp/vp,'range',vpwidth,'nbins',binnum1ddist);
    [Patemp,~] = hist_dg(vxa/va,'range',vawidth,'nbins',binnum1ddist);
    
    Petemp = Petemp*ne/(dvv*ve)/normfactore/numpts;
    Pptemp = Pptemp*np/(dvv*vp)/normfactorp/numpts;
    Patemp = Patemp*na/(dvv*va)/normfactora/numpts;
    
    edist1dall(:,ii) = Petemp;
    pdist1dall(:,ii) = Pptemp;
    adist1dall(:,ii) = Patemp;

  if mod(ii,10)==0
  % plot
   figure(1)
   set(gcf,'position',[10,10,1200,700])
   
    subplot('Position',[0.6 0.45 0.38 0.50])
    semilogy(vxeinit,Petemp,'k','linewidth',2)
    hold on
    semilogy(vxpinit,Pptemp,'r','linewidth',2)
    semilogy(vxainit,Patemp,'b','linewidth',2)
    hold off
    axis([-10 10 1e-4 2e2])
    xlabel('v/v_e')
    ylabel('f(v) (s m^{-4})')
    
    PEk = Ek.*conj(Ek);
    subplot('Position',[0.6 0.06 0.38 0.32])
    semilogy(kvecs(1:floor(numpts/2))*ld,PEk(1:floor(numpts/2)),'k','linewidth',2)
    xlabel('k \lambda_D')
    ylabel('E_k^2')
   
   subplot('Position',[0.05 0.71 0.5 0.25])
   plot(xa(1:50:end)/ld,vxa(1:50:end)/va,'g.')
   hold on
   plot(xe(1:50:end)/ld,vxe(1:50:end)/ve,'r.')
   plot(xp(1:50:end)/ld,vxp(1:50:end)/vp,'b.')
   hold off
   axis([0 xlength/ld -20 20])
   xlabel('x (\lambda_D)')
   ylabel('v/v_e') 
   title(num2str(ii))
   
   subplot('Position',[0.05 0.38 0.5 0.25])
   plot((dx/2:dx:xlength)/ld,Ex)
   xlabel('x (\lambda_D)')
   ylabel('E (V m^{-1})') 
   axis([0 xlength/ld -1.2*max(abs(Ex)) 1.2*max(abs(Ex))])

   subplot('Position',[0.05 0.06 0.5 0.25])
   semilogy(timesall,EE,'m','linewidth',2)
   hold on;
   plot(timesall,KEe,'r','linewidth',2)
   plot(timesall,KEp,'g','linewidth',2)
   plot(timesall,KEa,'b','linewidth',2)
   plot(timesall,Etot,'k','linewidth',2)
   axis([0 timesall(end) 0 Etot(ii)*2])
   hold off;
   ylabel('En. Density')
   xlabel('t \omega_{pe}') 

    %set(gca, 'Position',[0.6 0.45 0.38 0.50])

  end
end
legend({'W_E','W_e','W_p','W_a','W_{tot}'},'Location','northeast')

% For particle positions remove horizontal lines due to boundary conditions
%%
%electronposdiff = diff(electronpos,1,2);
%electronbposdiff = diff(electronbpos,1,2);
%zerosarr = zeros(100,1);

%electronposdiff = [zerosarr electronposdiff];
%electronbposdiff = [zerosarr electronbposdiff];
%idxe = abs(electronposdiff) > 500;
%idxeb = abs(electronbposdiff) > 500;

%electronpos(idxe) = NaN;
%electronbpos(idxeb) = NaN;

%% Diagnostics
%%  Plot frequency-wave number power spectrum
%pcolor(xvalsall,1:Ntimes,Exall)
%shading flat;

numt = 9000;

Exalls = Exall(end-numt+1:end,:);

Ekalls = fft2(Exalls);

Pk = Ekalls.*conj(Ekalls);

kmax = pi/dx;
deltak = 2*kmax/numpts;

kvec = (-kmax:deltak:kmax-deltak);

Pk2 = fftshift(Pk,2);
Pk2 = fftshift(Pk2,1);

omegamax = pi/dt;
deltaw = 2*omegamax/numt;
omegavec = -(-omegamax:deltaw:omegamax-deltaw);

%wL = sqrt(wpe^2 + 3*ve^2*kvec.^2/2);

figure; pcolor((kvec-deltak/2)*ld,omegavec/wpe,log10(Pk2))
%pcolor((kvec-deltak/2),omegavec/wpe,log10(Pk2))
%hold on
%plot(kvec*ld,wL/wpe,'r')
xlabel('k\lambda_D')
ylabel('\omega/\omega_{pe}')
colorbar;
title('W_E')
shading flat;
caxis([2 10])
axis([-2 2 -6 6])
set(gcf,'color','w')

%hold on; 
%plot(k*ld,real(wc1/wpetot),'k','linewidth',2)
%plot(k*ld,imag(wc1/wpetot),'r','linewidth',2)
%hold off;

%% Field versus position and time
figure; pcolor(xvalsall,timesall,Exall);
shading flat;
xlabel('x/\lambda_D')
ylabel('\omega_{pe} t')
colorbar;
%hold on;
%plot(electronpos(1:100,:)/ld,timesall,'k')
%plot(electronbpos(1:100,:)/ld,timesall,'m')
%hold off;
title('E_x (V m^{-1})')
colormap(bgrcmap);
set(gcf,'color','w')
caxis(max(max(abs(Exall)))*[-1.1 1.1])

%Vphtest = Vbulkeb/2;
%hold on;
%plot([0:1:timesall(end)]*Vphtest/ld,[0:1:timesall(end)]*wpe,'g','linewidth',2)
%hold off;


%% Potential position and time
if 1
figure; pcolor(xvalsall,timesall,phiall);
shading flat;
xlabel('x/\lambda_D')
ylabel('\omega_{pe} t')
colorbar;
hold on;
%plot(electronpos(3,:)/ld,timesall,'r')
%plot(electronbpos(1:10,:)/ld,timesall,'k')
hold off;
title('\phi (V)')
colormap(bgrcmap);
set(gcf,'color','w')
caxis(max(max(abs(phiall)))*[-1.1 1.1])
end
%% Charge density position and time
if 1
figure; pcolor(xvalsall,timesall,rhoall/(np*qe));
shading flat;
xlabel('x/\lambda_D')
ylabel('\omega_{pe} t')
colorbar;
title('\rho/(e n_i)')
colormap(bgrcmap);
set(gcf,'color','w')
caxis(max(max(abs(rhoall/(np*qe))))*[-1.1 1.1])
end
%% Estimate growth rate.
% Rough estimate of the \gamma/\omega_{pe} as a function of time
%gamma = 1/(2*dt)*log(EE(2:end)./EE(1:end-1))/wpe;

dWE = [0 diff(EE)];
gamma = 0.5*dWE./(EE*dt)/wpe;

%gammaest = median(gamma(200:400));

%% Electron distributions final



%%
%dxx = 1; % in Ld
%dvv = 0.1; % 0.1*ve
%hist2De = calculate_hist2D(xrange,1,vxrange,dv,xe/ld,vxe/ve);
%hist2Deb = calculate_hist2D(xrange,1,vxrange,dv,xeb/ld,vxeb/ve);

%e1Dreduced = hist2De.n1Dvx+hist2Deb.n1Dvx*neb/ne;
%e1Dreduced(e1Dreduced == 0) = NaN;
%normfactor = sum(sum(e1Dreduced,'omitnan'),'omitnan')/numpts;
%e1Dreduced = e1Dreduced*(ne + neb)/normfactor/(0.1*ve);
distidxvalue = 401;
figure;
pcolor(xpositions,vxapositions,log10(squeeze(adistmatrix(distidxvalue,:,:)))); shading flat; 
cbar = colorbar; 
ylabel(cbar,'f(v) (s m^{-4})');
hold on;
plot(xvalsall,Exall(distindices(distidxvalue),:)*10,'r')
hold off;
ylabel('v/v_e')
xlabel('x (\lambda_D)')


%% Simple plots of beam and electron speeds
%plot(timesall,electronbspeed)
%plot(timesall,electronspeed(1:100,:))

%% Plot of spatially averaged 1D distributions versus time
figure;
pcolor(timesall,vxpinit,log10(pdist1dall)); shading flat; 
cbar = colorbar;
ylabel(cbar,'f(v) (s m^{-4})');
ylabel('v/v_p')
xlabel('t \omega_{pe}')
%hold on
%plot(timesall,electronbspeed/ve,'b')
%plot(timesall,electronspeed/ve,'r')
%hold off

%%
if 0
vemat = repmat(vxpositions,length(distindices),1,length(xpositions));
fluxmatrix = edistmatrix.*vemat*ve;

pcolor(xpositions,vxpositions,squeeze(fluxmatrix(distidxvalue,:,:))); shading flat; colorbar; 
hold on;
plot(xvalsall,Exall(distindices(distidxvalue),:)*10,'r')
caxis(gca,[-1 1]*max(max(max(abs(fluxmatrix)))))
hold off;
colormap(bgrcmap);
end

%%

% save('ionioninstabilitydata')