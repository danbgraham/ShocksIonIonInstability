%% Parameters
np = 10.0e6; % m-3
na = 1.0e6;
ne = np+2*na; % m-3
Te = 100;

na/np


Units = irf_units; 
qe = Units.e;
me = Units.me;
mp = Units.mp;
ma = mp*4;
eps = Units.eps0;
kb = Units.kB;

ve = sqrt(2*qe*Te/me);

Vab = 100e3;

wpe = sqrt(ne*qe^2/(me*eps));
wpa = sqrt(na*2^2*qe^2/(ma*eps));
wpp = sqrt(np*qe^2/(mp*eps));
ld = ve/wpe/sqrt(2);

kvecld = 0.01:0.01:2;
kvec = kvecld/ld;

%wIA = cs*kvec;

wr1 = zeros(1,length(kvec));
wi1 = zeros(1,length(kvec));


Ta = 10.^[0:0.05:3];% eV
Tp = 10.^[0:0.05:3]; % eV

gammamax = zeros(length(Ta),length(Tp));
kmax = zeros(length(Ta),length(Tp));
wmax = zeros(length(Ta),length(Tp));

disprelarrtemp = zeros(length(Tp),length(kvecld));

for ii = 1:length(Ta)
  for jj = 1:length(Tp)

    Tatemp = Ta(ii);
    Tptemp = Tp(jj);
    va = sqrt(2*qe*Tatemp/ma);% m/s
    vp = sqrt(2*qe*Tptemp/mp); % m/s
    cs = ve/sqrt(2)*sqrt(me/mp)*(1+3*Tptemp/Te);
    
    % Define function

    for nn = 1:length(kvec)

      xip = @(omega)((omega)/(kvec(nn)*vp));
      xie = @(omega)((omega)/(kvec(nn)*ve));
      xia = @(omega)((omega - kvec(nn)*Vab)/(kvec(nn)*va));
      disprel = @(w) (1+2*wpp^2/(kvec(nn)^2*vp^2)*(1 + 1i*sqrt(pi)*xip(w).*faddeeva(xip(w),50))+...
        +2*wpe^2/(kvec(nn)^2*ve^2)*(1 + 1i*sqrt(pi)*xie(w).*faddeeva(xie(w),50))+...
        +2*wpa^2/(kvec(nn)^2*va^2)*(1 + 1i*sqrt(pi)*xia(w).*faddeeva(xia(w),50)));

      if ii < 2
        if jj < 2
          if (nn<5)    
            guessw1 = kvec(nn)*cs/5;
          else
            guessw1 = wr1(nn-1)+1i*wi1(nn-1);
          end
        else
          if (nn<5)
            guessw1 = disprelarrtemp(jj-1,nn);
          else
            guessw1 = wr1(nn-1)+1i*wi1(nn-1);
          end
        end
      else
        if (nn<5) 
          guessw1 = disprelarrtemp(jj,nn);
        else
          guessw1 = wr1(nn-1)+1i*wi1(nn-1);
        end
      end

      options = optimoptions('fsolve','Display','off', 'TolFun', 1e-12, 'TolX', 1e-11,'MaxFunEvals', 10000);
      [x,FVAL,EXITFLAG] = fsolve(@(w) disprel(w), guessw1,options); 
      [~,maxpos] = max(x);
      wr1(nn) = real(x(1));
      wi1(nn) = imag(x(1));

    end
    
    disprelarrtemp(jj,:) = wr1+1i*wi1;
    
    [wmaxtemp,idxgam] = max(wi1);
    if wmaxtemp < 0
      gammamax(ii,jj) = NaN;
      kmax(ii,jj) = NaN;
      wmax(ii,jj) = NaN;
    else
      gammamax(ii,jj) = wi1(idxgam);
      kmax(ii,jj) = kvec(idxgam);
      wmax(ii,jj) = wr1(idxgam);
    end
  end
  ii
end

%%
disprelTpTa = struct('ne',ne,'np',np,'na',na,'Te',Te,'Tpvec',Tp,'Tavec',Ta,'Vab',Vab,'ld',ld,...
  'wpe',wpe,'wpp',wpp,'wpa',wpa,'kvec',kvec,'gammamax',gammamax,'kmax',kmax,'wmax',wmax);
%save('disprelTpTa.mat','disprelTpTa');

%% 

rr = interp1([1 64 128 192 256],[0 90 190 250 255]/255,1:256,'pchip');
gg = interp1([1 64 128 192 256],[0 15 55 140 255]/255,1:256,'pchip');
bb = interp1([1 64 128 192 256],[0 110 80 10 0]/255,1:256,'pchip');
cmapinf = [rr' gg' bb'];


Tavec = disprelTpTa.Tavec;
Tpvec = disprelTpTa.Tpvec;
gammamax = disprelTpTa.gammamax;
wmax = disprelTpTa.wmax;
kmax = disprelTpTa.kmax;
wpp = disprelTpTa.wpp;
ld = disprelTpTa.ld;

gammamax(abs(wmax)/wpp > 10) = NaN;
kmax(abs(wmax)/wpp > 10) = NaN;
wmax(abs(wmax)/wpp > 10) = NaN;

Tanom = 12;
Tpnom = 3;


fn=figure;
set(fn,'Position',[10 10 1100 300])
h(1)=axes('position',[0.05 0.18 0.26 0.78]); 
h(2)=axes('position',[0.385 0.18 0.26 0.78]); 
h(3)=axes('position',[0.72 0.18 0.26 0.78]); 
ud=get(fn,'userdata');
ud.subplot_handles=h;
set(fn,'userdata',ud);
set(fn,'defaultLineLineWidth',2); 

%gammamax(gammamax > 1e6) = NaN;
pcolor(h(1),Tpvec,Tavec,gammamax/wpp)
set(h(1),'yscale','log')
set(h(1),'xscale','log')
%hold(h(1),'on')
%plot(h(1),Tpnom,Tanom,'g.')
%hold(h(1),'off')
colormap(h(1),'jet')
c = colorbar('peer',h(1));
c.Label.String = '\gamma_{max}/\omega_{pp}';
shading(h(1),'flat');
xlabel(h(1),'T_p (eV)','fontsize',14)
ylabel(h(1),'T_\alpha (eV)','fontsize',14)
caxis(h(1),[0 0.16]);
colormap(h(1),cmapinf)

pcolor(h(2),Tpvec,Tavec,wmax/wpp)
set(h(2),'yscale','log')
set(h(2),'xscale','log')
colormap(h(2),'jet')
c = colorbar('peer',h(2));
c.Label.String = '\omega_{max}/\omega_{pp}';
shading(h(2),'flat');
xlabel(h(2),'T_p (eV)','fontsize',14)
ylabel(h(2),'T_\alpha (eV)','fontsize',14)
caxis(h(2),[0 1]);
colormap(h(2),cmapinf)

pcolor(h(3),Tpvec,Tavec,kmax*ld)
set(h(3),'yscale','log')
set(h(3),'xscale','log')
colormap(h(3),'jet')
c = colorbar('peer',h(3));
c.Label.String = 'k_{max} \lambda_D';
shading(h(3),'flat');
xlabel(h(3),'T_p (eV)','fontsize',14)
ylabel(h(3),'T_\alpha (eV)','fontsize',14)
caxis(h(3),[0 1.8]);
colormap(h(3),cmapinf)

set(h(1:3),'fontsize',14);
set(h([1:3]),'Color',0.75*[1 1 1]);

set(gcf,'color','w')