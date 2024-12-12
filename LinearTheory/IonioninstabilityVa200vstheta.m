%% Parameters
np = 10.0e6; % m-3
na = 1.0e6;
ne = np+2*na; % m-3
Te = 100;
Ta = 12;% eV
Tp = 3; % eV


Units = irf_units; 
qe = Units.e;
me = Units.me;
mp = Units.mp;
ma = mp*4;
eps = Units.eps0;
kb = Units.kB;

vp = sqrt(2*qe*Tp/mp); % m/s
ve = sqrt(2*qe*Te/me);
va = sqrt(2*qe*Ta/ma);% m/s

Vab = 200e3;

wpe = sqrt(ne*qe^2/(me*eps));
wpa = sqrt(na*2^2*qe^2/(ma*eps));
wpp = sqrt(np*qe^2/(mp*eps));
ld = ve/wpe/sqrt(2);

kvecld = 0.005:0.01:2;
kvec = kvecld/ld;

theta = [0:0.5:90];

cs = ve/sqrt(2)*sqrt(me/mp)*(1+3*Tp/Te);
wr = zeros(1,length(kvec));
wi = zeros(1,length(kvec));

wrmat = zeros(length(theta),length(kvec));
wimat = zeros(length(theta),length(kvec));


for mm = 100:length(theta)
  
  
thetam = theta(mm)
kzvec = kvec*cosd(thetam);

% Define function

for nn = 1:length(kvec)

    xip = @(omega)((omega)/(kvec(nn)*vp));
    xie = @(omega)((omega)/(kvec(nn)*ve));
    xia = @(omega)((omega - kzvec(nn)*Vab)/(kvec(nn)*va));
    disprel = @(w) (1+2*wpp^2/(kvec(nn)^2*vp^2)*(1 + 1i*sqrt(pi)*xip(w).*faddeeva(xip(w),50))+...
                     +2*wpe^2/(kvec(nn)^2*ve^2)*(1 + 1i*sqrt(pi)*xie(w).*faddeeva(xie(w),50))+...
                     +2*wpa^2/(kvec(nn)^2*va^2)*(1 + 1i*sqrt(pi)*xia(w).*faddeeva(xia(w),50)));
    
    if mm < 101                 
      if (nn<5)    
          guessw = kvec(nn)*cs/2+1i*kvec(nn)*cs/10;
      else
          guessw = wr(nn-1)+1i*wi(nn-1);
      end
    %elseif mm < 70
      %guessw = wrmat(mm-1,nn) + 1i*(wimat(mm-1,nn)+2*abs(wimat(mm-1,nn)));
    %elseif mm < 75
    %  guessw = wrmat(mm-1,nn) + 1i*(abs(wimat(mm-1,nn)));
    else 
      guessw = wrmat(mm-1,nn) + 1i*(wimat(mm-1,nn));
    end
    
    options = optimoptions('fsolve','Display','off', 'TolFun', 1e-10, 'TolX', 1e-9,'MaxFunEvals', 1000);
    [x,FVAL,EXITFLAG] = fsolve(@(w) disprel(w), [guessw],options); 
    [~,maxpos] = max(x);
    wr(nn) = real(x(1));
    wi(nn) = imag(x(1));
    
    nn;
end

wrmat(mm,:) = wr;
wimat(mm,:) = wi;

end

for mm = 99:-1:1
  
  
thetam = theta(mm)
kzvec = kvec*cosd(thetam);

% Define function

for nn = 1:length(kvec)

    xip = @(omega)((omega)/(kvec(nn)*vp));
    xie = @(omega)((omega)/(kvec(nn)*ve));
    xia = @(omega)((omega - kzvec(nn)*Vab)/(kvec(nn)*va));
    disprel = @(w) (1+2*wpp^2/(kvec(nn)^2*vp^2)*(1 + 1i*sqrt(pi)*xip(w).*faddeeva(xip(w),50))+...
                     +2*wpe^2/(kvec(nn)^2*ve^2)*(1 + 1i*sqrt(pi)*xie(w).*faddeeva(xie(w),50))+...
                     +2*wpa^2/(kvec(nn)^2*va^2)*(1 + 1i*sqrt(pi)*xia(w).*faddeeva(xia(w),50)));
                     
    if nn < 5
      guessw = wrmat(mm+1,nn) + 1i*(wimat(mm+1,nn));
    else
      guessw = wr(nn-1) + 1i*(wi(nn-1));
    end
    
    options = optimoptions('fsolve','Display','off', 'TolFun', 1e-10, 'TolX', 1e-9,'MaxFunEvals', 1000);
    [x,FVAL,EXITFLAG] = fsolve(@(w) disprel(w), [guessw],options); 
    [~,maxpos] = max(x);
    wr(nn) = real(x(1));
    wi(nn) = imag(x(1));
    
    nn;
end

wrmat(mm,:) = wr;
wimat(mm,:) = wi;

end

%% Save data as a structure
disprelVa200 = struct('ne',ne,'np',np,'na',na,'Te',Te,'Tp',Tp,'Ta',Ta,'Vab',Vab,'ld',ld,...
  'wpe',wpe,'wpp',wpp,'wpa',wpa,'kvec',kvec,'theta',theta,'wrmat',wrmat,'wimat',wimat);
save('disprelVa200.mat','disprelVa200');


%%

rr = interp1([1 64 128 192 256],[0 90 190 250 255]/255,1:256,'pchip');
gg = interp1([1 64 128 192 256],[0 15 55 140 255]/255,1:256,'pchip');
bb = interp1([1 64 128 192 256],[0 110 80 10 0]/255,1:256,'pchip');
cmapinf = [rr' gg' bb'];

rr = interp1([1 64 128 192 256],[0.0  0.5 0.75 1.0 0.75],1:256);
gg = interp1([1 64 128 192 256],[0.0  0.5 0.75 0.5 0.00],1:256);
bb = interp1([1 64 128 192 256],[0.75 1.0 0.75 0.5 0.00],1:256);
cmapbgr = [rr' gg' bb'];

fn=figure;
set(fn,'Position',[10 10 900 400])
h(1)=axes('position',[0.06 0.12 0.4 0.81]); 
h(2)=axes('position',[0.55 0.12 0.4 0.81]); 
ud=get(fn,'userdata');
ud.subplot_handles=h;
set(fn,'userdata',ud);
set(fn,'defaultLineLineWidth',2); 
    
pcolor(h(1),kvecld,theta,wrmat/wpp)
shading(h(1),'flat');
xlabel(h(1),'k \lambda_D','fontsize',14)
ylabel(h(1),'\theta_{k\alpha}','interpreter','tex','fontsize',14)
title(h(1),'V_{\alpha} = 200 km s^{-1}','interpreter','tex','fontsize',14)
c = colorbar('peer',h(1));
c.Label.String = '\omega/\omega_{pp}';
colormap(h(1),cmapinf)
axis(h(1),[0 2 0 90])

pcolor(h(2),kvecld,theta,wimat/wpp)
shading(h(2),'flat');
xlabel(h(2),'k \lambda_D','fontsize',14)
ylabel(h(2),'\theta_{k\alpha}','interpreter','tex','fontsize',14)
title(h(2),'V_{\alpha} = 200 km s^{-1}','interpreter','tex','fontsize',14)
c = colorbar('peer',h(2));
c.Label.String = '\gamma/\omega_{pp}';
colormap(h(2),cmapbgr)
clim(h(2),[-0.15 0.15])
axis(h(2),[0 2 0 90])

set(h(1:2),'fontsize',14);

set(gcf,'color','w')