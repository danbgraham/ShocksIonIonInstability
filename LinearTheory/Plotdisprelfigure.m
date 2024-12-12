% Plot dispersion relation plots in a single figure
% Takes data stored from the other routines to calculate the dispersion
% properties of the proton-alpha instability and produces a single plot. 
% Written by D. B. Graham
%% Load nominal dispersion relation parameters
load('disprelnominal.mat')
ldnom = disprelnominal.ld;
thetavecnom = disprelnominal.theta;
kvecnom = disprelnominal.kvec;
wrarrnom = disprelnominal.wrmat;
wiarrnom = disprelnominal.wimat;
kvecldnom = kvecnom*ldnom;
wppnom = disprelnominal.wpp;

%% Load dispersion relation for V_\alpha = 200 km s^{-1}
load('disprelVa200.mat')
ld200 = disprelVa200.ld;
thetavec200 = disprelVa200.theta;
kvec200 = disprelVa200.kvec;
wrarr200 = disprelVa200.wrmat;
wiarr200 = disprelVa200.wimat;
kvecld200 = kvec200*ld200;
wpp200 = disprelVa200.wpp;

%% load max. growth rates and other properties versus ion temperatures
load('disprelTpTa.mat');
Tpvec = disprelTpTa.Tpvec;
Tavec = disprelTpTa.Tpvec;
gmaxTpTa = disprelTpTa.gammamax;
wmaxTpTa = disprelTpTa.wmax;
kmaxTpTa = disprelTpTa.kmax;
wppTpTa = disprelTpTa.wpp;
ldTpTa = disprelTpTa.ld;

gmaxTpTa(abs(wmaxTpTa)/wppTpTa > 10) = NaN;
kmaxTpTa(abs(wmaxTpTa)/wppTpTa > 10) = NaN;
wmaxTpTa(abs(wmaxTpTa)/wppTpTa > 10) = NaN;

%% Plot alpha-proton density ratio
load('disprelnanprat.mat')
npvec = disprelnanprat.npvec;
navec = disprelnanprat.navec;
wmaxnanp = disprelnanprat.wmax;
gmaxnanp = disprelnanprat.gammamax;
kmaxnanp = disprelnanprat.kmax;
wppnanp = disprelnanprat.wppvec;
wpananp = disprelnanprat.wpavec;
ldnanp = disprelnanprat.ld;

nanprat = navec./npvec;

%% Plot Dependence on Te and Va
load('disprelTeVa.mat')
wmaxTeVa = disprelTeVa.wmaxarr;
kmaxTeVa = disprelTeVa.kmaxarr;
gmaxTeVa = disprelTeVa.gmaxarr;
vphTeVa = disprelTeVa.vphmaxarr;
csTeVa = disprelTeVa.csarr;
wppTeVa = disprelTeVa.wpp;
Tevec = disprelTeVa.Tearr;
Vavec = disprelTeVa.Vabarr/1e3;
ldTeVa = disprelTeVa.ldarr;
thetaTeVa = disprelTeVa.thetaarr;
gmaxTeVa(gmaxTeVa < 0) = NaN;

Units = irf_units; 
qe = Units.e;
me = Units.me;
mp = Units.mp;

Tecsvec = [3:1:200];
vecsvec = sqrt(2*qe*Tecsvec/me);
csvec = vecsvec./sqrt(2).*sqrt(me/mp).*(1+3*3./Tecsvec);
csvec = csvec/1e3;

%% Define Color maps
rr = interp1([1 64 128 192 256],[0 90 190 250 255]/255,1:256,'pchip');
gg = interp1([1 64 128 192 256],[0 15 55 140 255]/255,1:256,'pchip');
bb = interp1([1 64 128 192 256],[0 110 80 10 0]/255,1:256,'pchip');
cmapinf = [rr' gg' bb'];

rr = interp1([1 64 128 192 256],[0.0  0.5 0.75 1.0 0.75],1:256);
gg = interp1([1 64 128 192 256],[0.0  0.5 0.75 0.5 0.00],1:256);
bb = interp1([1 64 128 192 256],[0.75 1.0 0.75 0.5 0.00],1:256);
cmapbgr = [rr' gg' bb'];

%% Plot Figure
fn=figure;
set(fn,'Position',[10 10 750 800])
h(1)=axes('position',[0.05 0.80 0.25 0.185]); 
h(2)=axes('position',[0.37 0.80 0.25 0.185]); 
h(3)=axes('position',[0.05 0.56 0.25 0.185]); 
h(4)=axes('position',[0.37 0.56 0.25 0.185]); 

h(5)=axes('position',[0.05 0.32 0.25 0.185]); 
h(6)=axes('position',[0.37 0.32 0.25 0.185]); 
h(7)=axes('position',[0.05 0.07 0.25 0.185]); 
h(8)=axes('position',[0.38 0.07 0.26 0.185]); 

h(9)=axes('position',[0.72 0.80 0.27 0.185]); 
h(10)=axes('position',[0.72 0.56 0.27 0.185]); 
h(11)=axes('position',[0.72 0.32 0.27 0.185]); 
h(12)=axes('position',[0.72 0.07 0.27 0.185]); 
ud=get(fn,'userdata');
ud.subplot_handles=h;
set(fn,'userdata',ud);
set(fn,'defaultLineLineWidth',2); 

fontsizeval = 12;

pcolor(h(1),kvecldnom,thetavecnom,wrarrnom/wppnom)
shading(h(1),'flat');
xlabel(h(1),'k \lambda_D','fontsize',fontsizeval)
ylabel(h(1),'\theta (\circ)','interpreter','tex','fontsize',fontsizeval)
irf_legend(h(1),{'(a)'},[0.98, 0.98],'color','b','fontsize',fontsizeval);
irf_legend(h(1),{'V_\alpha = 100 km s^{-1}'},[0.1, 0.98],'color','g','fontsize',fontsizeval);
c = colorbar('peer',h(1));
c.Label.String = '\omega/\omega_{pp}';
clim(h(1),[0 1.5])
colormap(h(1),cmapinf)
axis(h(1),[0 2 0 80])

pcolor(h(2),kvecldnom,thetavecnom,wiarrnom/wppnom)
shading(h(2),'flat');
xlabel(h(2),'k \lambda_D','fontsize',fontsizeval)
ylabel(h(2),'\theta (\circ)','interpreter','tex','fontsize',fontsizeval)
c = colorbar('peer',h(2));
c.Label.String = '\gamma/\omega_{pp}';
irf_legend(h(2),{'(b)'},[0.98, 0.98],'color','r','fontsize',fontsizeval);
colormap(h(2),cmapbgr)
clim(h(2),[-0.12 0.12])
axis(h(2),[0 2 0 80])

pcolor(h(3),kvecld200,thetavec200,wrarr200/wpp200)
shading(h(3),'flat');
xlabel(h(3),'k \lambda_D','fontsize',fontsizeval)
ylabel(h(3),'\theta (\circ)','interpreter','tex','fontsize',fontsizeval)
irf_legend(h(3),{'(c)'},[0.98, 0.98],'color','b','fontsize',fontsizeval);
irf_legend(h(3),{'V_\alpha = 200 km s^{-1}'},[0.2, 0.98],'color','g','fontsize',fontsizeval);
c = colorbar('peer',h(3));
c.Label.String = '\omega/\omega_{pp}';
colormap(h(3),cmapinf)
clim(h(3),[0 1.5])
axis(h(3),[0 2 0 80])

pcolor(h(4),kvecld200,thetavec200,wiarr200/wpp200)
shading(h(4),'flat');
xlabel(h(4),'k \lambda_D','fontsize',fontsizeval)
ylabel(h(4),'\theta (\circ)','interpreter','tex','fontsize',fontsizeval)
c = colorbar('peer',h(4));
c.Label.String = '\gamma/\omega_{pp}';
irf_legend(h(4),{'(d)'},[0.98, 0.98],'color','k','fontsize',fontsizeval);
colormap(h(4),cmapbgr)
clim(h(4),[-0.12 0.12])
axis(h(4),[0 2 0 80])

pcolor(h(5),Tpvec,Tavec,gmaxTpTa/wppTpTa)
set(h(5),'yscale','log')
set(h(5),'xscale','log')
colormap(h(5),'jet')
c = colorbar('peer',h(5));
c.Label.String = '\gamma_{max}/\omega_{pp}';
shading(h(5),'flat');
xlabel(h(5),'T_p (eV)','fontsize',fontsizeval)
ylabel(h(5),'T_\alpha (eV)','fontsize',fontsizeval)
clim(h(5),[0 0.16]);
colormap(h(5),cmapinf)
xticks(h(5),[1e0 1e1 1e2 1e3])
irf_legend(h(5),{'(e)'},[0.98, 0.98],'color','k','fontsize',fontsizeval);
set(h(5),'Color',0.75*[1 1 1]);

pcolor(h(6),Tpvec,Tavec,wmaxTpTa/wppTpTa)
set(h(6),'yscale','log')
set(h(6),'xscale','log')
colormap(h(6),'jet')
c = colorbar('peer',h(6));
c.Label.String = '\omega_{max}/\omega_{pp}';
shading(h(6),'flat');
xlabel(h(6),'T_p (eV)','fontsize',fontsizeval)
ylabel(h(6),'T_\alpha (eV)','fontsize',fontsizeval)
clim(h(6),[0 1]);
colormap(h(6),cmapinf)
xticks(h(6),[1e0 1e1 1e2 1e3])
irf_legend(h(6),{'(f)'},[0.98, 0.98],'color','k','fontsize',fontsizeval);
set(h(6),'Color',0.75*[1 1 1]);

pcolor(h(7),Tpvec,Tavec,kmaxTpTa*ldTpTa)
set(h(7),'yscale','log')
set(h(7),'xscale','log')
c = colorbar('peer',h(7));
c.Label.String = 'k_{max} \lambda_D';
shading(h(7),'flat');
xlabel(h(7),'T_p (eV)','fontsize',fontsizeval)
ylabel(h(7),'T_\alpha (eV)','fontsize',fontsizeval)
clim(h(7),[0 1.8]);
colormap(h(7),cmapinf)
irf_legend(h(7),{'(g)'},[0.98, 0.98],'color','k','fontsize',fontsizeval);
set(h(7),'Color',0.75*[1 1 1]);

semilogx(h(8),nanprat,gmaxnanp./wppnanp,nanprat,wmaxnanp./wppnanp,nanprat,kmaxnanp*ldnanp,'linewidth',2)
xlabel(h(8),'n_\alpha/n_p','fontsize',14)
xticks(h(8),[1e-3 1e-2 1e-1 1e0 1e1])
axis(h(8),[1e-3 1e1 0 1.1])
ylabel(h(8),'\gamma_{max}, \omega_{max}, k_{max}','fontsize',fontsizeval)
legend(h(8),{'\gamma/\omega_{pp}','\omega/\omega_{pp}','k\lambda_D'},'location','northwest','fontsize',fontsizeval)
irf_legend(h(8),{'(h)'},[0.98, 0.88],'color','k','fontsize',fontsizeval);


pcolor(h(9),Tevec,Vavec,gmaxTeVa')
xlabel(h(9),'T_e (eV)','fontsize',fontsizeval)
ylabel(h(9),'V_\alpha (km s^{-1})','interpreter','tex','fontsize',fontsizeval)
c = colorbar('peer',h(9));
c.Label.String = '\gamma_{max}/\omega_{pp}';
shading(h(9),'flat');
irf_legend(h(9),{'(i)'},[0.98, 0.98],'color','k','fontsize',fontsizeval);
%caxis(h(9),[0 0.15])
axis(h(9),[0 200 0 250])
colormap(h(9),cmapinf)
set(h(9),'Color',0.75*[1 1 1]);

pcolor(h(10),Tevec,Vavec,wmaxTeVa')
xlabel(h(10),'T_e (eV)','fontsize',fontsizeval)
ylabel(h(10),'V_\alpha (km s^{-1})','interpreter','tex','fontsize',fontsizeval)
c = colorbar('peer',h(10));
c.Label.String = '\omega_{max}/\omega_{pp}';
shading(h(10),'flat');
irf_legend(h(10),{'(j)'},[0.98, 0.98],'color','k','fontsize',fontsizeval);
clim(h(10),[0 1])
axis(h(10),[0 200 0 250])
colormap(h(10),cmapinf)
set(h(10),'Color',0.75*[1 1 1]);

pcolor(h(11),Tevec,Vavec,thetaTeVa')
xlabel(h(11),'T_e (eV)','fontsize',fontsizeval)
ylabel(h(11),'V_\alpha (km s^{-1})','interpreter','tex','fontsize',fontsizeval)
c = colorbar('peer',h(11));
c.Label.String = '\theta_{max} (\circ)';
shading(h(11),'flat');
irf_legend(h(11),{'(k)'},[0.98, 0.98],'color','k','fontsize',fontsizeval);
clim(h(11),[0 90])
axis(h(11),[0 200 0 250])
colormap(h(11),cmapinf)
set(h(11),'Color',0.75*[1 1 1]);

pcolor(h(12),Tevec,Vavec,kmaxTeVa')
xlabel(h(12),'T_e (eV)','fontsize',fontsizeval)
ylabel(h(12),'V_\alpha (km s^{-1})','interpreter','tex','fontsize',fontsizeval)
c = colorbar('peer',h(12));
c.Label.String = 'k_{max}\lambda_D';
shading(h(12),'flat');
irf_legend(h(12),{'(l)'},[0.98, 0.98],'color','k','fontsize',fontsizeval);
clim(h(12),[0 2])
axis(h(12),[0 200 0 250])
colormap(h(12),cmapinf)
set(h(12),'Color',0.75*[1 1 1]);

set(h(1:12),'fontsize',fontsizeval)

set(gcf,'color','w')