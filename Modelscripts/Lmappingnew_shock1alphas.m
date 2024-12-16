%% Init and model conditions
Units = irf_units;
mi = Units.mp;
e = Units.e;

Tiw = 48; % eV. Temperature used for weighted particles. 
Tio = 12; % eV. Temperature of the distribution being modelled. 

ma = mi*4;

vithw = sqrt(2*e*Tiw/ma);
vitho = sqrt(2*e*Tio/ma);

%%

Vsw = -540e3; 
n0 = 2.75*1e6; 
n1 = 8.25*1e6; 
B0 = 19.0*1e-9; 
B1 = 51.0*1e-9; 
nsw = n1-n0;

P0 = 0.088e-9;
P1 = 0.102e-9;

l = 5e3; 

xvec = [-1e3:0.1:1e3]*1e3;
dx = median(diff(xvec));
[xx,yy] = meshgrid(xvec,xvec);
Byxx = -B0*tanh(xx/l)+B1;

Bx = -0.35e-9;
By = -B0*tanh(xvec/l)+B1;
ni = -n0*tanh(xvec/l)+n1;

Units = irf_units;
mu0 = Units.mu0;
e = Units.e;
mi = Units.mp*4;

Jz = -B0*sech(xvec/l).^2/(mu0*l);
Vx = Vsw*nsw./ni;

Ex = -Jz.*By./(e*ni)+P0*sech(xvec/l).^2./(e*ni*l);
Ey = Jz.*Bx./(e*ni);
Ez = -Vsw.*(B1-B0);

phi = cumsum(Ex)*dx;

phi(end)
phitest = 2*B0^2/(e*mu0*n0) + B0*(B1*n0 - B0*n1)/(e*mu0*n0^2)*log((n1+n0)/(n1-n0)) + P0/(e*n0)*log((n1+n0)/(n1-n0))

nd = n1+n0;
nu = n1-n0;
Bd = B1+B0;
Bu = B1-B0;
Pd = P1+P0;
Pu = P1-P0;

phitest2 = (Bd - Bu)^2/(e*mu0*(nd - nu)) + (Bd - Bu)*(Bu*nd - Bd*nu)/(e*mu0*(nd - nu)^2)*log(nd/nu) % + (Pd - Pu)/(e*(nd - nu))*log(nd/nu)

%%

Numparticles = 1e5;
numperloop = 2e3; 
numloops = floor(Numparticles/numperloop)


vxi = randn(Numparticles,1)*vithw/sqrt(2);
vyi = randn(Numparticles,1)*vithw/sqrt(2);
vzi = randn(Numparticles,1)*vithw/sqrt(2);
viabs = sqrt(vxi.^2 + vyi.^2 + vzi.^2);
vxi = vxi + Vsw;

weightsall = vithw^3/vitho^3 * exp(viabs.^2*(vitho^2 - vithw^2)/(vitho^2*vithw^2));

dt = 1e-3;
endtime = 40;
finaltseries = [0:dt:endtime];

opts = odeset('RelTol',1e-7,'AbsTol',1e-9);

% Define grid for calculating ion distributions
xrange = [-2000 500]*1e3;
dxx = 10e3;
vxrange = [-1600 1600]*1e3;
vyrange = [-1600 1600]*1e3;
vzrange = [-1600 1600]*1e3;
dv = 10e3;

xpositions = xrange(1)+dxx/2:dxx:xrange(2);
numpositions = length(xpositions);
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

validpos = -1e15;

for jj = 1:numloops
  
  tempxposarr = zeros(numperloop,length(finaltseries));
  tempvxarr = zeros(numperloop,length(finaltseries));
  tempvyarr = zeros(numperloop,length(finaltseries));
  tempvzarr = zeros(numperloop,length(finaltseries));
  weightstemp = ones(numperloop,1);

  tic
  for ii = 1:numperloop
    x0 = 2000e3;
    y0 = 0;
    z0 = 0;
    initcondition = (jj - 1)*numperloop + ii;
    vx0 = vxi(initcondition);
    vy0 = vyi(initcondition);
    vz0 = vzi(initcondition);
    weights0 = weightsall(initcondition);

    [t1p,d1p] = ode45(@tempp,[0 endtime],[x0; vx0; y0; vy0; z0; vz0],opts);

    posveltemp = interp1(t1p,d1p,finaltseries,'linear');

    tempxposarr(ii,:) = posveltemp(:,1);
    tempvxarr(ii,:) = posveltemp(:,2);
    tempvyarr(ii,:) = posveltemp(:,4);
    tempvzarr(ii,:) = posveltemp(:,6);
    weightstemp(ii) = weights0;
    if max(tempxposarr(:,end)) > validpos
        validpos = max(tempxposarr(:,end));
    end

  end

  histstruct = calculate_histograms(xrange,dxx,vxrange,vyrange,vzrange,dv,tempxposarr,tempvxarr,tempvyarr,tempvzarr,weightstemp);

  n1Dvx = n1Dvx + histstruct.n1Dvx;
  n1Dvy = n1Dvy + histstruct.n1Dvy;
  n1Dvz = n1Dvz + histstruct.n1Dvz;

  n2Dvxvy = n2Dvxvy + histstruct.n2Dvxvy;
  n2Dvxvz = n2Dvxvz + histstruct.n2Dvxvz;
  n2Dvyvz = n2Dvyvz + histstruct.n2Dvyvz;
  
  toc
  jj
end

% Normalize to SI units
nswalpha = 1e6; % Quantities are renormalized by density when calculated
integraldx = sum(n1Dvx(end,:))*dv;
histstruct.n1Dvx = n1Dvx/integraldx*nswalpha;
histstruct.n1Dvy = n1Dvy/integraldx*nswalpha;
histstruct.n1Dvz = n1Dvz/integraldx*nswalpha;

integraldxdy = sum(sum(n2Dvxvy(end,:,:)))*dv^2;
histstruct.n2Dvxvy = n2Dvxvy/integraldxdy*nswalpha;
histstruct.n2Dvxvz = n2Dvxvz/integraldxdy*nswalpha;
histstruct.n2Dvyvz = n2Dvyvz/integraldxdy*nswalpha;

histstruct.Vsw = Vsw;
histstruct.n0 = n0;
histstruct.n1 = n1;
histstruct.B0 = B0;
histstruct.B1 = B1;
histstruct.Bn = Bx;
histstruct.P0 = P0;
histstruct.P1 = P1;
histstruct.l = l;
histstruct.validpos = validpos;
histstruct.dt = dt;
histstruct.dv = dv;
histstruct.Tio = Tio;
histstruct.Tiw = Tiw;
histstruct.ion = 'alpha';
histstruct.nswalpha = nswalpha;

save('histstruct_shock1_alphas.mat','histstruct')


function dydt = tempp(t,q)

l = 5e3;
vsw = -540e3;

n0 = 2.75*1e6; 
n1 = 8.25*1e6; 
B0 = 19.0*1e-9; 
B1 = 51.0*1e-9; 
nsw = n1-n0;

P0 = 0.088e-9;
P1 = 0.102e-9;

mu0 = 1.2566e-06;
eomi = 9.5788e+07/2;
e = 1.6022e-19;

Bx = -0.35e-9;
By = -B0*tanh(q(1)/l)+B1;
ni = -n0*tanh(q(1)/l)+n1;
Jz = -B0*sech(q(1)/l).^2/(mu0*l);
Ex = -Jz.*By./(e*ni)+P0*sech(q(1)/l).^2./(e*ni*l);
Ey = Jz.*Bx./(e*ni);
Ez = -vsw.*(B1-B0);

dydt = [q(2) ; eomi*Ex - eomi*q(6)*By; ...
        q(4) ; eomi*Ey + eomi*q(6)*Bx; ...
        q(6) ; eomi*Ez + eomi*q(2)*By - eomi*q(4)*Bx];
end