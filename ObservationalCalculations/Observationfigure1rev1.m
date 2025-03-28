%% Load data
ic = 1;
Tint1 = irf.tint('2023-04-24T02:29:00.000Z/2023-04-24T02:29:30.000Z');
Tint5 = irf.tint('2023-04-24T04:20:50.000Z/2023-04-24T04:21:20.000Z');

iPDist1 = mms.get_data('PDi_fpi_brst_l2',Tint1,ic);
iPDist5 = mms.get_data('PDi_fpi_brst_l2',Tint5,ic);

Bxyz1 = mms.get_data('B_dmpa_brst_l2',Tint1,ic);
Bxyz5 = mms.get_data('B_dmpa_brst_l2',Tint5,ic);

Exyz1=mms.get_data('E_dsl_edp_brst_l2',Tint1,ic);
Exyz5=mms.get_data('E_dsl_edp_brst_l2',Tint5,ic);

ne1 = mms.get_data('Ne_fpi_brst_l2',Tint1,ic);
ne5 = mms.get_data('Ne_fpi_brst_l2',Tint5,ic);

iPDist1omni = iPDist1.omni.deflux;
iPDist5omni = iPDist5.omni.deflux;

%% 
ndsl = [0.9203   -0.3886   -0.0445];
t1dsl = [-0.2271   -0.4383   -0.8696];
t2dsl = [0.3184    0.8105   -0.4917];
ndsl1 = ndsl/norm(ndsl);
t1dsl1 = t1dsl/norm(t1dsl);
t2dsl1 = t2dsl/norm(t2dsl);

ndsl = [0.8027   -0.5795    0.1412];
t1dsl = [-0.3078   -0.6052   -0.7341];
t2dsl = [0.5108    0.5458   -0.6642];
ndsl5 = ndsl/norm(ndsl);
t1dsl5 = t1dsl/norm(t1dsl);
t2dsl5 = t2dsl/norm(t2dsl);

%% Rotate data
Bntt1 = irf_newxyz(Bxyz1,ndsl1,t1dsl1,t2dsl1);
Bntt5 = irf_newxyz(Bxyz5,ndsl5,t1dsl5,t2dsl5);

iPDist1.data(:,1:11,:,:) = 0;
iPDist5.data(:,1:11,:,:) = 0;
%iPDist5.data(:,1:11,:,:) = 0;

%% Reduced distributions

load('vpa1.mat');
load('vpa5.mat');

f1Dn1 = vpa1.f1Dn1;
f1Dn5 = vpa5.f1Dn5;

vn1 = vpa1.vn1;
vn5 = vpa5.vn5;

vna1 = vpa1.vna1;
vna5 = vpa5.vna5;

vdiffn1 = vpa1.vdiffn1;
vdiffn5 = vpa5.vdiffn5;

vlimn = [-1800,1400];

%%
dfE = 1/median(diff(Exyz1.time.epochUnix));
fmin = 5; fmax = 4100; %Hz

Exyzhf1 = Exyz1.filt(fmin,0,dfE,5);
Exyzhf5 = Exyz5.filt(fmin,0,dfE,5);

%% Wavelet transforms
nf = 100;
nc = 10;

Ewavelet = irf_wavelet(Exyz1,'nf',nf,'f',[fmin fmax]);

idx = nc/2:nc:length(Ewavelet.t)-nc/2;
Ewavelettimes = Ewavelet.t(idx);
Ewaveletx = zeros(length(idx),nf);
Ewavelety = zeros(length(idx),nf);
Ewaveletz = zeros(length(idx),nf);
for ii = 1:length(idx)
  Ewaveletx(ii,:) = squeeze(irf.nanmean(Ewavelet.p{1,1}(idx(ii)-nc/2+1:idx(ii)+nc/2-1,:),1));
  Ewavelety(ii,:) = squeeze(irf.nanmean(Ewavelet.p{1,2}(idx(ii)-nc/2+1:idx(ii)+nc/2-1,:),1));
  Ewaveletz(ii,:) = squeeze(irf.nanmean(Ewavelet.p{1,3}(idx(ii)-nc/2+1:idx(ii)+nc/2-1,:),1));
end
specE1=struct('t',Ewavelettimes);
specE1.f=Ewavelet.f;
specE1.p=Ewaveletx+Ewavelety+Ewaveletz;
specE1.f_label='';
specE1.p_label={'log_{10} E^2','(mV^2 m^{-2} Hz^{-1})'};

Ewavelet = irf_wavelet(Exyz5,'nf',nf,'f',[fmin fmax]);

idx = nc/2:nc:length(Ewavelet.t)-nc/2;
Ewavelettimes = Ewavelet.t(idx);
Ewaveletx = zeros(length(idx),nf);
Ewavelety = zeros(length(idx),nf);
Ewaveletz = zeros(length(idx),nf);
for ii = 1:length(idx)
  Ewaveletx(ii,:) = squeeze(irf.nanmean(Ewavelet.p{1,1}(idx(ii)-nc/2+1:idx(ii)+nc/2-1,:),1));
  Ewavelety(ii,:) = squeeze(irf.nanmean(Ewavelet.p{1,2}(idx(ii)-nc/2+1:idx(ii)+nc/2-1,:),1));
  Ewaveletz(ii,:) = squeeze(irf.nanmean(Ewavelet.p{1,3}(idx(ii)-nc/2+1:idx(ii)+nc/2-1,:),1));
end
specE5=struct('t',Ewavelettimes);
specE5.f=Ewavelet.f;
specE5.p=Ewaveletx+Ewavelety+Ewaveletz;
specE5.f_label='';
specE5.p_label={'log_{10} E^2','(mV^2 m^{-2} Hz^{-1})'};

%% Compute characteristic frequencies
Units=irf_units; % read in standard units
Me=Units.me;
Mp=Units.mp;
e=Units.e;
epso=Units.eps0;
mu0=Units.mu0;
Mp_Me = Mp/Me;

magB1 = Bxyz1.abs;
magB5 = Bxyz5.abs;

B_SI1=magB1.data*1e-9;
Wpe1 = sqrt(ne1.resample(Bxyz1).data*1e6*e^2/Me/epso);
Wce1 = e*B_SI1/Me;
Wpp1 = sqrt(ne1.resample(Bxyz1).data*1e6*e^2/Mp/epso);
Fce1 = Wce1/2/pi;
Fpe1 = Wpe1/2/pi;
Fcp1 = Fce1/Mp_Me;
Fpp1 = Wpp1/2/pi;
Flh1 = sqrt(Fcp1.*Fce1./(1+Fce1.^2./Fpe1.^2)+Fcp1.^2);
Fce1 = irf.ts_scalar(magB1.time,Fce1);
Flh1 = irf.ts_scalar(magB1.time,Flh1);
Fpp1 = irf.ts_scalar(magB1.time,Fpp1);

B_SI5=magB5.data*1e-9;
Wpe5 = sqrt(ne5.resample(Bxyz5).data*1e6*e^2/Me/epso);
Wce5 = e*B_SI5/Me;
Wpp5 = sqrt(ne5.resample(Bxyz5).data*1e6*e^2/Mp/epso);
Fce5 = Wce5/2/pi;
Fpe5 = Wpe5/2/pi;
Fcp5 = Fce5/Mp_Me;
Fpp5 = Wpp5/2/pi;
Flh5 = sqrt(Fcp5.*Fce5./(1+Fce5.^2./Fpe5.^2)+Fcp5.^2);
Fce5 = irf.ts_scalar(magB5.time,Fce5);
Flh5 = irf.ts_scalar(magB5.time,Flh5);
Fpp5 = irf.ts_scalar(magB5.time,Fpp5);

%% Convert bulk velocities Vn to energy
if 0,
vn1E = 0.5*mp*vn1.data.^2*1e6/e;
vn5E = 0.5*mp*vn5.data.^2*1e6/e;

vna1E = 0.5*mp*vna1.data.^2*1e6/e;
vna5E = 0.5*mp*vna5.data.^2*1e6/e;

vn1E = irf.ts_scalar(vn1.time,vn1E);
vn5E = irf.ts_scalar(vn5.time,vn5E);

vna1E = irf.ts_scalar(vna1.time,vna1E);
vna5E = irf.ts_scalar(vna5.time,vna5E);
end

%% Colormaps

c = [55,137,187;...
  106,193,165;...
  172,220,166;...
  230,244,157;...
  255,254,194;...
  253,223,144;...
  251,173,104;...
  242,109,074;...
  211,064,082]/255;
cmap = interp1(linspace(1,64,size(c,1)),c,1:64);

rr = interp1([1 64 128 192 256],[0 90 190 250 255]/255,1:256,'pchip');
gg = interp1([1 64 128 192 256],[0 15 55 140 255]/255,1:256,'pchip');
bb = interp1([1 64 128 192 256],[0 110 80 10 0]/255,1:256,'pchip');
cmapinf = [rr' gg' bb'];
%% Plot figure

h=irf_plot(12,'newfigure');
%h=irf_figure(540+ic,8);
xSize=1000; ySize=700;
set(gcf,'Position',[10 10 xSize ySize]);

xwidth = 0.40;
ywidth = 0.15;
set(h(1),'position',[0.07 0.97-ywidth xwidth ywidth]);
set(h(2),'position',[0.07 0.97-2*ywidth xwidth ywidth]);
set(h(3),'position',[0.07 0.97-3*ywidth xwidth ywidth]);
set(h(4),'position',[0.07 0.97-4*ywidth xwidth ywidth]);
set(h(5),'position',[0.07 0.97-5*ywidth xwidth ywidth]);
set(h(6),'position',[0.07 0.97-6*ywidth xwidth ywidth]);

set(h(7),'position',[0.570 0.97-ywidth xwidth ywidth]);
set(h(8),'position',[0.570 0.97-2*ywidth xwidth ywidth]);
set(h(9),'position',[0.570 0.97-3*ywidth xwidth ywidth]);
set(h(10),'position',[0.570 0.97-4*ywidth xwidth ywidth]);
set(h(11),'position',[0.570 0.97-5*ywidth xwidth ywidth]);
set(h(12),'position',[0.570 0.97-6*ywidth xwidth ywidth]);

h(1)=irf_panel('Bntt1');
irf_plot(h(1),Bntt1);
ylabel(h(1),{'B (nT)'},'Interpreter','tex');
irf_zoom(h(1),'y',[-15 75]);
irf_legend(h(1),{'B_{n}','B_{t1}','B_{t2}'},[0.98 0.35],'fontsize',14)
irf_legend(h(1),'(a)',[0.99 0.8],'color','k','fontsize',14)
title(h(1),'Shock 1','fontsize',14);

h(2)=irf_panel('idist1');
irf_spectrogram(h(2),iPDist1omni.specrec,'log','donotfitcolorbarlabel');
%hold(h(2),'on')
%irf_plot(h(2),vn1E,'k')
%irf_plot(h(2),vn1E,'k.')
%irf_plot(h(2),vna1E,'r')
%irf_plot(h(2),vna1E,'r.')
%hold(h(2),'off')
irf_legend(h(2),'(b)',[0.99 0.98],'color','k','fontsize',14);
set(h(2),'yscale','log');
set(h(2),'ytick',[1e1 1e2 1e3 1e4]);
ylabel(h(2),'E_{i} (eV)','fontsize',14,'Interpreter','tex');

specrecn = f1Dn1.specrec('1D_velocity');
specrecn.p_label={'log_{10}f_{i,n}','(s m^{-4})'};
specrecn.p(specrecn.p < 1e-3) = NaN;
h(3)=irf_panel('f1Dn1');
irf_spectrogram(h(3),specrecn);
hold(h(3),'on')
irf_plot(h(3),vn1,'k')
irf_plot(h(3),vn1,'k.')
irf_plot(h(3),vna1,'r')
irf_plot(h(3),vna1,'r.')
hold(h(3),'off')
ylabel(h(3),{'v_n (km s^{-1})'},'interpreter','tex')
irf_legend(h(3),{'V_p',' ','V_\alpha'},[0.1 0.98],'fontsize',14)
colormap(h(3),cmap)
irf_zoom(h(3),'y',vlimn)
clim(h(3),[-3 2])
irf_legend(h(3),'(c)',[0.99 0.98],'color','k','fontsize',14)

h(4)=irf_panel('dVn1');
irf_plot(h(4),vdiffn1);
hold(h(4),'on')
irf_plot(h(4),vdiffn1,'k.')
hold(h(4),'off')
ylabel(h(4),{'\Delta V_n (km s^{-1})'},'Interpreter','tex');
irf_legend(h(4),'(d)',[0.99 0.98],'color','k','fontsize',14)

h(5)=irf_panel('Ehf1');
irf_plot(h(5),Exyzhf1);
hold(h(5),'on')
irf_plot(h(5),Exyzhf1.x,'k');
irf_plot(h(5),Exyzhf1.y,'b');
hold(h(5),'off')
ylabel(h(5),{'\delta E (mV m^{-1})'},'Interpreter','tex');
irf_legend(h(5),{'E_{x}','E_{y}','E_{z}'},[0.98 0.12],'fontsize',14)
irf_legend(h(5),'(e)',[0.99 0.94],'color','k','fontsize',14)
irf_legend(h(5),sprintf('f > %.1f Hz',fmin),[0.01 0.1],'color','k','fontsize',14)

h(6)=irf_panel('Espec1');
irf_spectrogram(h(6),specE1,'log');
hold(h(6),'on');
irf_plot(h(6),Flh1,'color','b','LineWidth',1.5)
irf_plot(h(6),Fce1,'color','r','LineWidth',1.5)
irf_plot(h(6),Fpp1,'color',[0.8 0.8 0],'LineWidth',1.5)
hold(h(6),'off');
irf_legend(h(6),'(f)',[0.99 0.65],'color','w','fontsize',14)
irf_legend(h(6),'f_{LH}',[0.10 0.55],'color','b','fontsize',14,'backgroundcolor','w')
irf_legend(h(6),'f_{ce}',[0.02 0.55],'color','r','fontsize',14,'backgroundcolor','w')
irf_legend(h(6),'f_{pp}',[0.18 0.55],'color',[0.8 0.8 0],'fontsize',14,'backgroundcolor','w')
clim(h(6),[-6 2]);
set(h(6),'yscale','log');
set(h(6),'ytick',[1e1 1e2 1e3 1e4]);
ylabel(h(6),{'f (Hz)'},'fontsize',14,'Interpreter','tex');

colormap(h(2),cmap);
colormap(h(6),cmapinf);

irf_plot_axis_align(h(1:6));
irf_zoom(h(1:6),'x',Tint1);


h(7)=irf_panel('Bntt5');
irf_plot(h(7),Bntt5);
ylabel(h(7),{'B (nT)'},'Interpreter','tex');
irf_zoom(h(7),'y',[-15 75]);
irf_legend(h(7),'(g)',[0.99 0.8],'color','k','fontsize',14)
title(h(7),'Shock 5','fontsize',14);

h(8)=irf_panel('idist5');
irf_spectrogram(h(8),iPDist5omni .specrec,'log','donotfitcolorbarlabel');
%hold(h(8),'on')
%irf_plot(h(8),vn5E,'k')
%irf_plot(h(8),vn5E,'k.')
%irf_plot(h(8),vna5E,'r')
%irf_plot(h(8),vna5E,'r.')
%hold(h(8),'off')
irf_legend(h(8),'(h)',[0.99 0.98],'color','k','fontsize',14);
set(h(8),'yscale','log');
set(h(8),'ytick',[1e1 1e2 1e3 1e4]);
ylabel(h(8),'E_{i} (eV)','fontsize',14,'Interpreter','tex');

specrecn = f1Dn5.specrec('1D_velocity');
specrecn.p_label={'log_{10}f_{i,n}','(s m^{-4})'};
specrecn.p(specrecn.p < 1e-3) = NaN;
h(9)=irf_panel('f1Dn5');
irf_spectrogram(h(9),specrecn);
hold(h(9),'on')
irf_plot(h(9),vn5,'k')
irf_plot(h(9),vn5,'k.')
irf_plot(h(9),vna5,'r')
irf_plot(h(9),vna5,'r.')
hold(h(9),'off')
ylabel(h(9),{'v_n (km s^{-1})'},'interpreter','tex')
colormap(h(9),cmap)
irf_zoom(h(9),'y',vlimn)
clim(h(9),[-3 2])
irf_legend(h(9),'(i)',[0.99 0.98],'color','k','fontsize',14)

h(10)=irf_panel('dVn5');
irf_plot(h(10),vdiffn5);
hold(h(10),'on')
irf_plot(h(10),vdiffn5,'k.')
hold(h(10),'off')
ylabel(h(10),{'\Delta V_n (km s^{-1})'},'Interpreter','tex');
irf_legend(h(10),'(j)',[0.99 0.98],'color','k','fontsize',14)

h(11)=irf_panel('Ehf5');
irf_plot(h(11),Exyzhf5);
hold(h(11),'on')
irf_plot(h(11),Exyzhf5.x,'k');
irf_plot(h(11),Exyzhf5.y,'b');
hold(h(11),'off')
ylabel(h(11),{'\delta E (mV m^{-1})'},'Interpreter','tex');
%irf_legend(h(11),{'E_{x}','E_{y}','E_{z}'},[0.98 0.12],'fontsize',14)
irf_legend(h(11),'(k)',[0.99 0.94],'color','k','fontsize',14)
irf_legend(h(11),sprintf('f > %.1f Hz',fmin),[0.01 0.1],'color','k','fontsize',14)

h(12)=irf_panel('Espec5');
irf_spectrogram(h(12),specE5,'log');
hold(h(12),'on');
irf_plot(h(12),Flh5,'color','b','LineWidth',1.5)
irf_plot(h(12),Fce5,'color','r','LineWidth',1.5)
irf_plot(h(12),Fpp5,'color',[0.8 0.8 0],'LineWidth',1.5)
hold(h(12),'off');
irf_legend(h(12),'(l)',[0.99 0.65],'color','w','fontsize',14)
clim(h(12),[-6 2]);
set(h(12),'yscale','log');
set(h(12),'ytick',[1e1 1e2 1e3 1e4]);
ylabel(h(12),{'f (Hz)'},'fontsize',14,'Interpreter','tex');

colormap(h(8),cmap);
colormap(h(12),cmapinf);

set(h([3 9]),'Color',0.75*[1 1 1]);
grid(h(3),'off')
grid(h(9),'off')

irf_plot_axis_align(h(7:12));
irf_zoom(h(7:12),'x',Tint5);

set(h(1:12),'fontsize',14);
xtickangle(h(1:12),0)