%% Shock times
Tint1 = irf.tint('2023-04-24T02:28:13.00Z/2023-04-24T02:31:03.00Z');
Tint2 = irf.tint('2023-04-24T03:49:33.00Z/2023-04-24T03:51:03.00Z');
Tint3 = irf.tint('2023-04-24T04:00:03.00Z/2023-04-24T04:03:33.00Z');
Tint4 = irf.tint('2023-04-24T04:14:03.00Z/2023-04-24T04:16:33.00Z');
Tint5 = irf.tint('2023-04-24T04:20:03.00Z/2023-04-24T04:22:23.00Z');

Tints1 = irf.tint('2023-04-24T02:29:11.85Z/2023-04-24T02:29:11.85Z');
Tints2 = irf.tint('2023-04-24T03:50:12.3Z/2023-04-24T03:50:12.3Z');
Tints3 = irf.tint('2023-04-24T04:02:03.384Z/2023-04-24T04:02:04.384Z');
Tints4 = irf.tint('2023-04-24T04:15:25.0Z/2023-04-24T04:15:25.0Z');
Tints5 = irf.tint('2023-04-24T04:21:00.00Z/2023-04-24T04:21:00.00Z');

Vsn1 = 50;
Vsn2 = -140;
Vsn3 = 60;
Vsn4 = -50;
Vsn5 = 40;

%% Shock coordinate systems
ndsl = [0.9203   -0.3886   -0.0445];
t1dsl = [-0.2271   -0.4383   -0.8696];
t2dsl = [0.3184    0.8105   -0.4917];
ndsl1 = ndsl/norm(ndsl);
t1dsl1 = t1dsl/norm(t1dsl);
t2dsl1 = t2dsl/norm(t2dsl);

ndsl = [0.9015    0.1799   -0.3935];
t1dsl = [-0.1538   -0.7168   -0.6801];
t2dsl = [-0.4045    0.6736   -0.6186];
ndsl2 = ndsl/norm(ndsl);
t1dsl2 = t1dsl/norm(t1dsl);
t2dsl2 = t2dsl/norm(t2dsl);

ndsl = [0.8871   -0.4602   -0.0356];
t1dsl = [-0.3382   -0.5956   -0.7286];
t2dsl = [0.3142    0.6584   -0.6840];
ndsl3 = ndsl/norm(ndsl);
t1dsl3 = t1dsl/norm(t1dsl);
t2dsl3 = t2dsl/norm(t2dsl);

ndsl = [0.9489   -0.2611   -0.1774];
t1dsl = [-0.3127   -0.7001   -0.6419];
t2dsl = [0.0434    0.6646   -0.7460];
ndsl4 = ndsl/norm(ndsl);
t1dsl4 = t1dsl/norm(t1dsl);
t2dsl4 = t2dsl/norm(t2dsl);

ndsl = [0.8027   -0.5795    0.1412];
t1dsl = [-0.3078   -0.6052   -0.7341];
t2dsl = [0.5108    0.5458   -0.6642];
ndsl5 = ndsl/norm(ndsl);
t1dsl5 = t1dsl/norm(t1dsl);
t2dsl5 = t2dsl/norm(t2dsl);

%% Load data
Bxyz1 = mms.get_data('B_dmpa_brst_l2',Tint1,1);
Bxyz2 = mms.get_data('B_dmpa_brst_l2',Tint2,1);
Bxyz3 = mms.get_data('B_dmpa_brst_l2',Tint3,1);
Bxyz4 = mms.get_data('B_dmpa_brst_l2',Tint4,1);
Bxyz5 = mms.get_data('B_dmpa_brst_l2',Tint5,1);

Exyz1 = mms.get_data('E_dsl_edp_brst_l2',Tint1,1);
Exyz2 = mms.get_data('E_dsl_edp_brst_l2',Tint2,1);
Exyz3 = mms.get_data('E_dsl_edp_brst_l2',Tint3,1);
Exyz4 = mms.get_data('E_dsl_edp_brst_l2',Tint4,1);
Exyz5 = mms.get_data('E_dsl_edp_brst_l2',Tint5,1);

Te1 = mms.get_data('Te_dbcs_fpi_brst_l2',Tint1,1);
Te2 = mms.get_data('Te_dbcs_fpi_brst_l2',Tint2,1);
Te3 = mms.get_data('Te_dbcs_fpi_brst_l2',Tint3,1);
Te4 = mms.get_data('Te_dbcs_fpi_brst_l2',Tint4,1);
Te5 = mms.get_data('Te_dbcs_fpi_brst_l2',Tint5,1);

ne1 = mms.get_data('Ne_fpi_brst_l2',Tint1,1);
ne2 = mms.get_data('Ne_fpi_brst_l2',Tint2,1);
ne3 = mms.get_data('Ne_fpi_brst_l2',Tint3,1);
ne4 = mms.get_data('Ne_fpi_brst_l2',Tint4,1);
ne5 = mms.get_data('Ne_fpi_brst_l2',Tint5,1);

%% Scalar temperature
Tes1 = Te1.trace/3;
Tes2 = Te2.trace/3;
Tes3 = Te3.trace/3;
Tes4 = Te4.trace/3;
Tes5 = Te5.trace/3;

%% Rotate of vectors

Bntt1 = irf_newxyz(Bxyz1,ndsl1,t1dsl1,t2dsl1);
Bntt2 = irf_newxyz(Bxyz2,ndsl2,t1dsl2,t2dsl2);
Bntt3 = irf_newxyz(Bxyz3,ndsl3,t1dsl3,t2dsl3);
Bntt4 = irf_newxyz(Bxyz4,ndsl4,t1dsl4,t2dsl4);
Bntt5 = irf_newxyz(Bxyz5,ndsl5,t1dsl5,t2dsl5);

Bntt1 = Bntt1.resample(Exyz1);
Bntt2 = Bntt2.resample(Exyz2);
Bntt3 = Bntt3.resample(Exyz3);
Bntt4 = Bntt4.resample(Exyz4);
Bntt5 = Bntt5.resample(Exyz5);

Tes1r = Tes1.resample(Exyz1);
Tes2r = Tes2.resample(Exyz2);
Tes3r = Tes3.resample(Exyz3);
Tes4r = Tes4.resample(Exyz4);
Tes5r = Tes5.resample(Exyz5);

dExyz1 = Exyz1.filt(200,0,8192,5);
dExyz2 = Exyz2.filt(200,0,8192,5);
dExyz3 = Exyz3.filt(200,0,8192,5);
dExyz4 = Exyz4.filt(200,0,8192,5);
dExyz5 = Exyz5.filt(200,0,8192,5);

%% Delta V 
load('vpa1.mat');
load('vpa2.mat');
load('vpa3.mat');
load('vpa4.mat');
load('vpa5.mat');

vpn1 = vpa1.vn1;
vpn2 = vpa2.vn2;
vpn3 = vpa3.vn3;
vpn4 = vpa4.vn4;
vpn5 = vpa5.vn5;

vdiffn1 = vpa1.vdiffn1;
vdiffn2 = vpa2.vdiffn2;
vdiffn3 = vpa3.vdiffn3;
vdiffn4 = vpa4.vdiffn4;
vdiffn5 = vpa5.vdiffn5;

%% Convert position to time
Times1 = Bntt1.time.epochUnix - Tints1(1).epochUnix;
Times2 = Bntt2.time.epochUnix - Tints2(1).epochUnix;
Times3 = Bntt3.time.epochUnix - Tints3(1).epochUnix;
Times4 = Bntt4.time.epochUnix - Tints4(1).epochUnix;
Times5 = Bntt5.time.epochUnix - Tints5(1).epochUnix;

Timesi1 = vpn1.time.epochUnix - Tints1(1).epochUnix;
Timesi2 = vpn2.time.epochUnix - Tints2(1).epochUnix;
Timesi3 = vpn3.time.epochUnix - Tints3(1).epochUnix;
Timesi4 = vpn4.time.epochUnix - Tints4(1).epochUnix;
Timesi5 = vpn5.time.epochUnix - Tints5(1).epochUnix;

npos1 = -Times1*Vsn1;
npos2 = -Times2*Vsn2;
npos3 = -Times3*Vsn3;
npos4 = -Times4*Vsn4;
npos5 = -Times5*Vsn5;

nposi1 = -Timesi1*Vsn1;
nposi2 = -Timesi2*Vsn2;
nposi3 = -Timesi3*Vsn3;
nposi4 = -Timesi4*Vsn4;
nposi5 = -Timesi5*Vsn5;

%% Load electrostatic wave properties
load('SH_2023_04_24_T02_27_MMS_1_stat_results.mat');
FPL1 = FPL;
FSC1 = FSC;
K_lFAC1 = K_lFAC;
K_nt1t21 = K_nt1t2;
l3d1 = l3d;
TKB1 = TKB;
Vphpl1 = Vphpl;
Vphsc1 = Vphsc;

load('SH_2023_04_24_T03_49_MMS_1_stat_results.mat');
FPL2 = FPL;
FSC2 = FSC;
K_lFAC2 = K_lFAC;
K_nt1t22 = K_nt1t2;
l3d2 = l3d;
TKB2 = TKB;
Vphpl2 = Vphpl;
Vphsc2 = Vphsc;

load('SH_2023_04_24_T04_00_MMS_1_stat_results.mat');
FPL3 = FPL;
FSC3 = FSC;
K_lFAC3 = K_lFAC;
K_nt1t23 = K_nt1t2;
l3d3 = l3d;
TKB3 = TKB;
Vphpl3 = Vphpl;
Vphsc3 = Vphsc;

load('SH_2023_04_24_T04_14_MMS_1_stat_results.mat');
FPL4 = FPL;
FSC4 = FSC;
K_lFAC4 = K_lFAC;
K_nt1t24 = K_nt1t2;
l3d4 = l3d;
TKB4 = TKB;
Vphpl4 = Vphpl;
Vphsc4 = Vphsc;

load('SH_2023_04_24_T04_19_MMS_1_stat_results.mat');
FPL5 = FPL;
FSC5 = FSC;
K_lFAC5 = K_lFAC;
K_nt1t25 = K_nt1t2;
l3d5 = l3d;
TKB5 = TKB;
Vphpl5 = Vphpl;
Vphsc5 = Vphsc;

clear FPL FSC K_lFAC K_nt1t2 l3d TKB Vphpl Vphsc;


Units = irf_units;
me = Units.me;
qe = Units.e;
eps = Units.eps0;
mp = Units.mp;

c_eval('fpi? = sqrt(qe^2*ne?.data*1e6/(mp*eps))/2/pi;',1:5);
c_eval('fpi? = irf.ts_scalar(ne?.time,fpi?);',1:5);
c_eval('fpi? = fpi?.resample(l3d?);',1:5);
c_eval('freqfpp? = irf.ts_scalar(l3d?.time,FPL?.abs.data./fpi?.data);',1:5);
c_eval('freqfppsc? = irf.ts_scalar(l3d?.time,FSC?.abs.data./fpi?.data);',1:5);

c_eval('cs? = sqrt(qe*Tes?.data/mp);',1:5);
c_eval('cs? = irf.ts_scalar(ne?.time,cs?);',1:5);
c_eval('cs? = cs?.resample(l3d?);',1:5);

c_eval('lambdaD? = sqrt(eps*Tes?.data./(ne?.data*1e6*qe));',1:5);
c_eval('lambdaD? = irf.ts_scalar(ne?.time,lambdaD?);',1:5);
c_eval('lambdaD? = lambdaD?.resample(l3d?);',1:5);

c_eval('kval? = 2*pi./l3d?.data.*lambdaD?.data;',1:5);
c_eval('kval? = irf.ts_scalar(l3d?.time,kval?);',1:5);

c_eval('thetakB? = TKB?.data;',1:5);
c_eval('idx? = thetakB? > 90;',1:5);
c_eval('thetakB?(idx?) = 180-thetakB?(idx?);',1:5);
c_eval('thetakB? = irf.ts_scalar(TKB?.time,thetakB?);',1:5);

c_eval('vphplcs? = Vphpl?.abs.data./cs?.data;',1:5);
c_eval('vphplcs? = irf.ts_scalar(l3d?.time,vphplcs?)',1:5);

c_eval('vphplcssc? = Vphsc?.abs.data./cs?.data;',1:5);
c_eval('vphplcssc? = irf.ts_scalar(l3d?.time,vphplcssc?)',1:5);

%% Get position
Timesw1 = l3d1.time.epochUnix - Tints1(1).epochUnix;
Timesw2 = l3d2.time.epochUnix - Tints2(1).epochUnix;
Timesw3 = l3d3.time.epochUnix - Tints3(1).epochUnix;
Timesw4 = l3d4.time.epochUnix - Tints4(1).epochUnix;
Timesw5 = l3d5.time.epochUnix - Tints5(1).epochUnix;

nposw1 = -Timesw1*Vsn1;
nposw2 = -Timesw2*Vsn2;
nposw3 = -Timesw3*Vsn3;
nposw4 = -Timesw4*Vsn4;
nposw5 = -Timesw5*Vsn5;

%% Average properties
c_eval('idx? = nposw? < 50 & nposw? > -100;',1:5);
c_eval('median(vphplcssc?.abs.data(idx?))',1:5);

%% Figure 

h=irf_plot(8,'newfigure');
%h=irf_figure(540+ic,8);
xSize=700; ySize=800;
set(gcf,'Position',[10 10 xSize ySize]);

xwidth = 0.90;
ywidth = 0.113;
set(h(1),'position',[0.080 0.995-ywidth xwidth ywidth]);
set(h(2),'position',[0.080 0.990-2*ywidth xwidth ywidth]);
set(h(3),'position',[0.080 0.985-3*ywidth xwidth ywidth]);
set(h(4),'position',[0.080 0.980-4*ywidth xwidth ywidth]);
set(h(5),'position',[0.080 0.975-5*ywidth xwidth ywidth]);
set(h(6),'position',[0.080 0.970-6*ywidth xwidth ywidth]);
set(h(7),'position',[0.080 0.965-7*ywidth xwidth ywidth]);
set(h(8),'position',[0.080 0.960-8*ywidth xwidth ywidth]);


h(1)=irf_panel('Bn');
plot(h(1),npos1,Bntt1.y.data,'k','linewidth',1.5)
hold(h(1),'on')
plot(h(1),npos2,Bntt2.y.data,'color',[0.3 0.3 1],'linewidth',1.5)
plot(h(1),npos3,Bntt3.y.data,'r','linewidth',1.5)
plot(h(1),npos4,Bntt4.y.data,'color',[0 0.6 0],'linewidth',1.5)
plot(h(1),npos5,Bntt5.y.data,'m','linewidth',1.5)
hold(h(1),'off')
axis(h(1),[-600 600 0 75])
ylabel(h(1),'B_{t1} (nT)')
set(h(1),'xticklabel',[])
irf_legend(h(1),'(a)',[0.99 0.97])

h(2)=irf_panel('Te');
plot(h(2),npos1,Tes1r.data,'k','linewidth',1.5)
hold(h(2),'on')
plot(h(2),npos2,Tes2r.data,'color',[0.3 0.3 1],'linewidth',1.5)
plot(h(2),npos3,Tes3r.data,'r','linewidth',1.5)
plot(h(2),npos4,Tes4r.data,'color',[0 0.6 0],'linewidth',1.5)
plot(h(2),npos5,Tes5r.data,'m','linewidth',1.5)
hold(h(2),'off')
axis(h(2),[-600 600 0 165])
ylabel(h(2),'T_{e} (eV)')
set(h(2),'xticklabel',[])
irf_legend(h(2),'(b)',[0.99 0.97])

h(3)=irf_panel('Vpn');
plot(h(3),nposi1,vpn1.data-Vsn1,'k','linewidth',1.5)
hold(h(3),'on')
plot(h(3),nposi2,vpn2.data-Vsn2,'color',[0.3 0.3 1],'linewidth',1.5)
plot(h(3),nposi3,vpn3.data-Vsn3,'r','linewidth',1.5)
plot(h(3),nposi4,vpn4.data-Vsn4,'color',[0 0.6 0],'linewidth',1.5)
plot(h(3),nposi5,vpn5.data-Vsn5,'m','linewidth',1.5)
hold(h(3),'off')
axis(h(3),[-600 600 -670 60])
ylabel(h(3),'V_{p,n} (km s^{-1})')
set(h(3),'xticklabel',[])
irf_legend(h(3),'(c)',[0.99 0.97])

h(4)=irf_panel('Vdiffn');
plot(h(4),nposi1,vdiffn1.data,'k','linewidth',1.5)
hold(h(4),'on')
plot(h(4),nposi2,vdiffn2.data,'color',[0.3 0.3 1],'linewidth',1.5)
plot(h(4),nposi3,vdiffn3.data,'r','linewidth',1.5)
plot(h(4),nposi4,vdiffn4.data,'color',[0 0.6 0],'linewidth',1.5)
plot(h(4),nposi5,vdiffn5.data,'m','linewidth',1.5)
plot(h(4),nposi1,vdiffn1.data,'k.','linewidth',1.5,'MarkerSize',10)
plot(h(4),nposi2,vdiffn2.data,'.','color',[0.3 0.3 1],'linewidth',1.5,'MarkerSize',10)
plot(h(4),nposi3,vdiffn3.data,'r.','linewidth',1.5,'MarkerSize',10)
plot(h(4),nposi4,vdiffn4.data,'.','color',[0 0.6 0],'linewidth',1.5,'MarkerSize',10)
plot(h(4),nposi5,vdiffn5.data,'m.','linewidth',1.5,'MarkerSize',10)
plot(h(4),[-600 600],[0 0],'k--','linewidth',0.5)
hold(h(4),'off')
legend(h(4),{'Shock 1','Shock 2','Shock 3','Shock 4','Shock 5'},'location','northwest')
axis(h(4),[-600 600 -60 300])
ylabel(h(4),'\DeltaV_n (km s^{-1})')
set(h(4),'xticklabel',[])
irf_legend(h(4),'(d)',[0.99 0.97])

h(5)=irf_panel('Eabs');
plot(h(5),npos1,dExyz1.abs.data,'k','linewidth',1.5)
hold(h(5),'on')
plot(h(5),npos2,dExyz2.abs.data,'color',[0.3 0.3 1],'linewidth',1.5)
plot(h(5),npos3,dExyz3.abs.data,'r','linewidth',1.5)
plot(h(5),npos4,dExyz4.abs.data,'color',[0 0.6 0],'linewidth',1.5)
plot(h(5),npos5,dExyz5.abs.data,'m','linewidth',1.5)
plot(h(5),npos2,dExyz2.abs.data,'color',[0.3 0.3 1],'linewidth',1.5)
hold(h(5),'off')
axis(h(5),[-600 600 0 1100])
ylabel(h(5),'|\delta E| (mV m^{-1})','fontsize',13)
set(h(5),'xticklabel',[])
irf_legend(h(5),'(e)',[0.99 0.97])
%xlabel(h(5),'n (km)')

h(6)=irf_panel('klambda');
plot(h(6),nposw1,kval1.data,'k.','MarkerSize',10)
hold(h(6),'on')
plot(h(6),nposw2,kval2.data,'.','color',[0.3 0.3 1],'MarkerSize',10)
plot(h(6),nposw3,kval3.data,'r.','MarkerSize',10)
plot(h(6),nposw4,kval4.data,'.','color',[0 0.6 0],'MarkerSize',10)
plot(h(6),nposw5,kval5.data,'m.','MarkerSize',10)
hold(h(6),'off')
axis(h(6),[-600 600 0 1.4])
ylabel(h(6),'k\lambda_D')
set(h(6),'xticklabel',[])
irf_legend(h(6),'(f)',[0.99 0.97])

h(7)=irf_panel('kn');
plot(h(7),nposw1,K_nt1t21.x.data,'k.','MarkerSize',10)
hold(h(7),'on')
plot(h(7),nposw2,K_nt1t22.x.data,'.','color',[0.3 0.3 1],'MarkerSize',10)
plot(h(7),nposw3,K_nt1t23.x.data,'r.','MarkerSize',10)
plot(h(7),nposw4,K_nt1t24.x.data,'.','color',[0 0.6 0],'MarkerSize',10)
plot(h(7),nposw5,K_nt1t25.x.data,'m.','MarkerSize',10)
hold(h(7),'off')
axis(h(7),[-600 600 -1.2 1.2])
ylabel(h(7),'k_n/k')
irf_legend(h(7),'(g)',[0.99 0.97])
set(h(7),'xticklabel',[])

h(8)=irf_panel('thetakB');
plot(h(8),nposw1,thetakB1.data,'k.','MarkerSize',10)
hold(h(8),'on')
plot(h(8),nposw2,thetakB2.data,'.','color',[0.3 0.3 1],'MarkerSize',10)
plot(h(8),nposw3,thetakB3.data,'r.','MarkerSize',10)
plot(h(8),nposw4,thetakB4.data,'.','color',[0 0.6 0],'MarkerSize',10)
plot(h(8),nposw5,thetakB5.data,'m.','MarkerSize',10)
hold(h(8),'off')
axis(h(8),[-600 600 0 100])
ylabel(h(8),'\theta_{kB} (\circ)')
irf_legend(h(8),'(h)',[0.99 0.97])
yticks(h(8),[0 45 90])
%set(h(8),'xticklabel',[])
xlabel(h(8),'n (km)')


set(h(1:8),'fontsize',13);
