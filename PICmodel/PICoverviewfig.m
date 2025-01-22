%% Load data
load('ionioninstabilitydata.mat');
clear vxe vxp vxa xe xp xa mate matp mata;

%% Data analysis
Units = irf_units;
me = Units.me;
mp = Units.mp;

% Add correction to normalization
edist1dallc = edist1dall*dvv/dve1d;
pdist1dallc = pdist1dall*dvv/dvp1d;
adist1dallc = adist1dall*dvv/dva1d;

% Change time to proton angular plasma frequency
timesallwpp = timesall*wpp/wpe;

%% Average for fields data
% Time average fields data
ntotal = length(timesall);
nptsav = 120;
ntotals = floor(ntotal/nptsav);
Exaveraged = zeros(ntotals,128);
phiaveraged = zeros(ntotals,128);
rhoaveraged = zeros(ntotals,128);
timesav = timesallwpp(floor(nptsav/2):nptsav:ntotal);

for ii = 1:ntotals
  Exaveraged(ii,:) = mean(Exall((ii-1)*nptsav+1:ii*nptsav,:),1); 
  phiaveraged(ii,:) = mean(phiall((ii-1)*nptsav+1:ii*nptsav,:),1); 
  rhoaveraged(ii,:) = mean(rhoall((ii-1)*nptsav+1:ii*nptsav,:),1); 
end


%% Average for 1D particle data
% Apply time average to 1D distributions
edist1dallav = zeros(ntotals,300);
pdist1dallav = zeros(ntotals,300);
adist1dallav = zeros(ntotals,300);

for ii = 1:ntotals
  edist1dallav(ii,:) = mean(edist1dallc(:,(ii-1)*nptsav+1:ii*nptsav),2); 
  pdist1dallav(ii,:) = mean(pdist1dallc(:,(ii-1)*nptsav+1:ii*nptsav),2); 
  adist1dallav(ii,:) = mean(adist1dallc(:,(ii-1)*nptsav+1:ii*nptsav),2); 
end

%% Averaged distributions at the start and end of simulation
edist1Dinit = mean(edist1dallc(:,1:500),2,'omitnan');
pdist1Dinit = mean(pdist1dallc(:,1:500),2,'omitnan');
adist1Dinit = mean(adist1dallc(:,1:500),2,'omitnan');

edist1Dfinal = mean(edist1dallc(:,end-500:end),2,'omitnan');
pdist1Dfinal = mean(pdist1dallc(:,end-500:end),2,'omitnan');
adist1Dfinal = mean(adist1dallc(:,end-500:end),2,'omitnan');

[~,idxt100] = min(abs(timesallwpp-100));

edist1Dt100 = mean(edist1dallc(:,idxt100-250:idxt100-249),2,'omitnan');
pdist1Dt100 = mean(pdist1dallc(:,idxt100-250:idxt100-249),2,'omitnan');
adist1Dt100 = mean(adist1dallc(:,idxt100-250:idxt100-249),2,'omitnan');

%% Renormalize energies to energy densities in units of eV/m^3
Wec = KEe/numpts/dx/qe;
Wpc = KEp/numpts/dx/qe;
Wac = KEa/numpts/dx/qe;
WEc = EE/numpts/dx/qe;
WEc = smooth(WEc,120,'moving');

Wtot = Wec+Wpc+Wac+WEc';


ttmp = [0.007:0.00001:0.014];
gamma = 4.16e3*0.116; % Prediction from linear growth rate
E2 = 1.5e1*exp(2*gamma*ttmp);

%% Calculate spatially-averaged particle moments
%electrons
vxeinitmat = vxeinit'*ones(1,length(edist1dallc));

nemoms = sum(edist1dallc,1)*dve1d*ve;
Vemoms = sum(vxeinitmat*ve.*edist1dallc,1)*dve1d*ve;
Vemoms = Vemoms./nemoms;
Vemomsmat = repmat(Vemoms,300,1);
Pemoms = me*sum(edist1dallc.*(vxeinitmat*ve-Vemomsmat).^2*dve1d*ve,1);
Temoms = Pemoms./nemoms/qe;
vethsim = sqrt(2*qe*Temoms/me);

vxpinitmat = vxpinit'*ones(1,length(pdist1dallc));
npmoms = sum(pdist1dallc,1)*dvp1d*vp;
Vpmoms = sum(vxpinitmat*vp.*pdist1dallc,1)*dvp1d*vp;
Vpmoms = Vpmoms./npmoms;
Vpmomsmat = repmat(Vpmoms,300,1);
Ppmoms = mp*sum(pdist1dallc.*(vxpinitmat*vp-Vpmomsmat).^2*dvp1d*vp,1);
Tpmoms = Ppmoms./npmoms/qe;
vpthsim = sqrt(2*qe*Tpmoms/mp);

vxainitmat = vxainit'*ones(1,length(adist1dallc));
namoms = sum(adist1dallc,1)*dva1d*va;
Vamoms = sum(vxainitmat*va.*adist1dallc,1)*dva1d*va;
Vamoms = Vamoms./namoms;
Vamomsmat = repmat(Vamoms,300,1);
Pamoms = ma*sum(adist1dallc.*(vxainitmat*va-Vamomsmat).^2*dva1d*va,1);
Tamoms = Pamoms./namoms/qe;
vathsim = sqrt(2*qe*Tamoms/ma);

%% Plot Figure

idx1 = 420;
idx2 = 590;
idx3 = 1200;
idx4 = 2360;

tvalidx1 = timesallwpp(distindices(idx1));
tvalidx2 = timesallwpp(distindices(idx2));
tvalidx3 = timesallwpp(distindices(idx3));
tvalidx4 = timesallwpp(distindices(idx4));


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

fn=figure;
set(fn,'Position',[10 10 900 800])
h(1)=axes('position',[0.06 0.51 0.18 0.46]); 
h(2)=axes('position',[0.30 0.51 0.18 0.46]); 
h(3)=axes('position',[0.53 0.87 0.25 0.12]); 
h(4)=axes('position',[0.53 0.75 0.25 0.12]); 
h(5)=axes('position',[0.53 0.63 0.25 0.12]); 
h(6)=axes('position',[0.53 0.51 0.25 0.12]); 

h(7)=axes('position',[0.89 0.87 0.105 0.12]); 
h(8)=axes('position',[0.89 0.69 0.105 0.12]); 
h(9)=axes('position',[0.89 0.51 0.105 0.12]); 

h(10)=axes('position',[0.045 0.35 0.44 0.095]);
h(11)=axes('position',[0.045 0.255 0.44 0.095]);
h(12)=axes('position',[0.045 0.155 0.44 0.095]);
h(13)=axes('position',[0.045 0.06 0.44 0.095]);
h(14)=axes('position',[0.49 0.35 0.44 0.095]);
h(15)=axes('position',[0.49 0.255 0.44 0.095]);
h(16)=axes('position',[0.49 0.155 0.44 0.095]);
h(17)=axes('position',[0.49 0.06 0.44 0.095]);

ud=get(fn,'userdata');
ud.subplot_handles=h;
set(fn,'userdata',ud);
set(fn,'defaultLineLineWidth',2); 

% Fields
pcolor(h(1),xvalsall,timesav,Exaveraged);
shading(h(1),'flat');
hold(h(1),'on')
plot(h(1),[0 64],[tvalidx1 tvalidx1],'--','color',[0 0.7 0],'linewidth',1)
plot(h(1),[0 64],[tvalidx2 tvalidx2],'--','color',[0 0.7 0],'linewidth',1)
plot(h(1),[0 64],[tvalidx3 tvalidx3],'--','color',[0 0.7 0],'linewidth',1)
plot(h(1),[0 64],[tvalidx4 tvalidx4],'--','color',[0 0.7 0],'linewidth',1)
hold(h(1),'off')
xlabel(h(1),'n/\lambda_D','fontsize',13)
ylabel(h(1),'\omega_{pp} t','fontsize',13)
cbar = colorbar('peer',h(1));
title(h(1),'E (V m^{-1})','fontsize',12)
colormap(h(1),bgrcmap);
clim(h(1),max(max(abs(Exaveraged)))*[-1.0 1.0])
axis(h(1),[0 64 0 timesallwpp(end)])
xticks(h(1),[0 20 40 60 80])
irf_legend(h(1),'(a)',[0.98 0.995],'fontsize',13,'color','k')

pcolor(h(2),xvalsall,timesav,phiaveraged);
shading(h(2),'flat');
hold(h(2),'on')
plot(h(2),[0 64],[tvalidx1 tvalidx1],'--','color',[0 0.7 0],'linewidth',1)
plot(h(2),[0 64],[tvalidx2 tvalidx2],'--','color',[0 0.7 0],'linewidth',1)
plot(h(2),[0 64],[tvalidx3 tvalidx3],'--','color',[0 0.7 0],'linewidth',1)
plot(h(2),[0 64],[tvalidx4 tvalidx4],'--','color',[0 0.7 0],'linewidth',1)
hold(h(2),'off')
xlabel(h(2),'n/\lambda_D','fontsize',13)
ylabel(h(2),'\omega_{pp} t','fontsize',13)
cbar = colorbar('peer',h(2));
title(h(2),'\phi (V)','fontsize',12)
colormap(h(2),bgrcmap);
clim(h(2),max(max(abs(phiaveraged)))*[-1.0 1.0])
axis(h(2),[0 64 0 timesallwpp(end)])
xticks(h(2),[0 20 40 60 80])
irf_legend(h(2),'(b)',[0.98 0.995],'fontsize',13,'color','k')

% Energy densities
semilogy(h(3),timesallwpp,Wec,'color','k','linewidth',2)
hold(h(3),'on')
semilogy(h(3),timesallwpp,Wpc,'color','r','linewidth',2)
semilogy(h(3),timesallwpp,Wac,'color','b','linewidth',2)
semilogy(h(3),timesallwpp,WEc,'color',[0 0.5 0],'linewidth',2)
semilogy(h(3),timesallwpp,Wtot,'color','m','linewidth',2)
semilogy(h(3),ttmp*wpp,E2,'k--')
hold(h(3),'off')
set(h(3),'xticklabel',[])
ylabel(h(3),'W (eV m^{-3})')
axis(h(3),[0 timesallwpp(end) 5e3 1e9])
yticks(h(3),[1e3 1e5 1e7 1e9])
irf_legend(h(3),'(c)',[0.99 0.65],'fontsize',13,'color','k')
irf_legend(h(3),'W_e',[0.3 0.12],'fontsize',13,'color','k')
irf_legend(h(3),'W_p',[0.4 0.12],'fontsize',13,'color','r')
irf_legend(h(3),'W_\alpha',[0.57 0.12],'fontsize',13,'color','b')
irf_legend(h(3),'W_E',[0.67 0.12],'fontsize',13,'color',[0 0.5 0])
irf_legend(h(3),'W_{tot}',[0.79 0.12],'fontsize',13,'color','m')

% Particle distributions
edist1dallav(edist1dallav == 0) = NaN;
pcolor(h(4),timesav,vxeinit,log10(edist1dallav'));
shading(h(4),'flat');
hold(h(4),'on')
plot(h(4),timesallwpp,Vemoms/ve,'color','k','linewidth',1.0)
plot(h(4),timesallwpp,(Vemoms+vethsim)/ve,'k--','linewidth',1.0)
plot(h(4),timesallwpp,(Vemoms-vethsim)/ve,'k--','linewidth',1.0)
hold(h(4),'off')
xlabel(h(4),' ','fontsize',13)
ylabel(h(4),'v/v_{e,init}','fontsize',13)
cbar = colorbar(h(4),'position',[0.782 0.755 0.01 0.11]);
cbar.Label.String = 'log_{10}\langlef_e\rangle (s m^{-4})';
colormap(h(4),cmap);
set(h(4),'xticklabel',[])
clim(h(4),[-4.5 0.2]);
irf_legend(h(4),'Electrons',[0.01 0.98],'fontsize',13,'color','k')
axis(h(4),[0 timesallwpp(end) -4.9 4.9])
irf_legend(h(4),'(d)',[0.99 0.98],'fontsize',13,'color','k')

adist1dallav(adist1dallav == 0) = NaN;
pcolor(h(5),timesav,vxainit,log10(adist1dallav'));
shading(h(5),'flat');
hold(h(5),'on')
plot(h(5),timesallwpp,Vamoms/va,'color','k','linewidth',1.0)
plot(h(5),timesallwpp,(Vamoms+vathsim)/va,'k--','linewidth',1.0)
plot(h(5),timesallwpp,(Vamoms-vathsim)/va,'k--','linewidth',1.0)
hold(h(5),'off')
xlabel(h(5),' ','fontsize',13)
ylabel(h(5),'v/v_{\alpha,init}','fontsize',13)
cbar = colorbar(h(5),'position',[0.782 0.635 0.01 0.11]);
cbar.Label.String = 'log_{10}\langlef_\alpha\rangle (s m^{-4})';
colormap(h(5),cmap);
set(h(5),'xticklabel',[])
clim(h(5),[-3.5 1.5]);
irf_legend(h(5),'Alphas',[0.01 0.08],'fontsize',13,'color','k')
axis(h(5),[0 timesallwpp(end) -3.9 9.5])
irf_legend(h(5),'(e)',[0.99 0.98],'fontsize',13,'color','k')

pdist1dallav(pdist1dallav == 0) = NaN;
pcolor(h(6),timesav,vxpinit,log10(pdist1dallav'));
shading(h(6),'flat');
hold(h(6),'on')
plot(h(6),timesallwpp,Vpmoms/vp,'color','k','linewidth',1.0)
plot(h(6),timesallwpp,(Vpmoms+vpthsim)/vp,'k--','linewidth',1.0)
plot(h(6),timesallwpp,(Vpmoms-vpthsim)/vp,'k--','linewidth',1.0)
hold(h(6),'off')
xlabel(h(6),' ','fontsize',13)
ylabel(h(6),'v/v_{p,init}','fontsize',13)
cbar = colorbar(h(6),'position',[0.782 0.515 0.01 0.11]);
cbar.Label.String = 'log_{10}\langlef_p\rangle (s m^{-4})';
colormap(h(6),cmap);
clim(h(6),[-2.5 2.2]);
irf_legend(h(6),'Protons',[0.01 0.98],'fontsize',13,'color','k')
xlabel(h(6),'\omega_{pp} t','fontsize',13)
axis(h(6),[0 timesallwpp(end) -3.9 9.5])
irf_legend(h(6),'(f)',[0.99 0.98],'fontsize',13,'color','k')

set(h([4:6]),'Color',0.75*[1 1 1]);

% 1D distributions at beginning and end of run
plot(h(7),vxeinit,edist1Dinit,'color','k','linewidth',2)
hold(h(7),'on')
plot(h(7),vxeinit,edist1Dt100,'color',[0 0.5 0],'linewidth',2)
plot(h(7),vxeinit,edist1Dfinal,'color','r','linewidth',2)
hold(h(7),'off')
set(h(7),'yscale','log')
xlabel(h(7),'v/v_e','fontsize',13)
ylabel(h(7),'\langlef_e\rangle (s m^{-4})','fontsize',13)
axis(h(7),[-2.4 2.4 5e-3 2e0])
yticks(h(7),[1e-4 1e-3 1e-2 1e-1 1e0])
irf_legend(h(7),'Initial',[0.4 0.35],'fontsize',13,'color','k')
irf_legend(h(7),'\omega_{pp} t = 100',[0.15 0.20],'fontsize',13,'color',[0 0.5 0])
irf_legend(h(7),'Final',[0.4 0.02],'fontsize',13,'color','r')
irf_legend(h(7),'(g)',[0.98 0.98],'fontsize',13,'color','k')

plot(h(8),vxainit,adist1Dinit,'color','k','linewidth',2)
hold(h(8),'on')
plot(h(8),vxainit,adist1Dt100,'color',[0 0.5 0],'linewidth',2)
plot(h(8),vxainit,adist1Dfinal,'color','r','linewidth',2)
hold(h(8),'off')
set(h(8),'yscale','log')
xlabel(h(8),'v/v_\alpha','fontsize',13)
ylabel(h(8),'\langlef_\alpha\rangle (s m^{-4})','fontsize',13)
axis(h(8),[0 7 5e-2 5e1])
yticks(h(8),[1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2])
xticks(h(8),[0 2 4 6 8 10])
irf_legend(h(8),'(h)',[0.98 0.98],'fontsize',13,'color','k')

plot(h(9),vxpinit,pdist1Dinit,'color','k','linewidth',2)
hold(h(9),'on')
plot(h(9),vxpinit,pdist1Dt100,'color',[0 0.5 0],'linewidth',2)
plot(h(9),vxpinit,pdist1Dfinal,'color','r','linewidth',2)
hold(h(9),'off')
set(h(9),'yscale','log')
xlabel(h(9),'v/v_p','fontsize',13)
ylabel(h(9),'\langlef_p\rangle (s m^{-4})','fontsize',13)
axis(h(9),[-3.5 7.5 5e-2 5e2])
yticks(h(9),[1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2])
irf_legend(h(9),'(i)',[0.98 0.98],'fontsize',13,'color','k')

% 1D ion distrbutions as a function of position
pcolor(h(10),xpositions,vxapositions,log10(squeeze(adistmatrix(idx1,:,:)))); 
shading(h(10),'flat')
hold(h(10),'on');
plot(h(10),xvalsall,mean(Exall(distindices(idx1)-60:distindices(idx1)+60,:),1)*10,'b','linewidth',2)
hold(h(10),'off');
clim(h(10),[-2 2.5]);
ylabel(h(10),'v/v_{\alpha,init}','fontsize',14)
set(h(10),'xticklabel',[])
colormap(h(10),cmap)
irf_legend(h(10),'(j)',[0.99 0.98],'fontsize',13,'color','k')
irf_legend(h(10),'10E (V m^{-1})',[0.98 0.1],'fontsize',13,'color','b')
irf_legend(h(10),['\omega_{pp}t = ' num2str(tvalidx1,2)],[0.02 0.14],'fontsize',13,'color',[0 0.7 0])
axis(h(10),[0 64 -8 8])

pcolor(h(11),xpositions,vxppositions,log10(squeeze(pdistmatrix(idx1,:,:)))); 
shading(h(11),'flat')
hold(h(11),'on');
plot(h(11),xvalsall,mean(Exall(distindices(idx1)-60:distindices(idx1)+60,:),1)*10,'b','linewidth',2)
hold(h(11),'off');
ylabel(h(11),'v/v_{p,init}','fontsize',14)
set(h(11),'xticklabel',[])
clim(h(11),[-2 2.5]);
colormap(h(11),cmap)
irf_legend(h(11),'(k)',[0.99 0.98],'fontsize',13,'color','k')
axis(h(11),[0 64 -8 8])

pcolor(h(12),xpositions,vxapositions,log10(squeeze(adistmatrix(idx2,:,:)))); 
shading(h(12),'flat')
hold(h(12),'on');
plot(h(12),xvalsall,mean(Exall(distindices(idx2)-60:distindices(idx2)+60,:),1)*10,'b','linewidth',2)
hold(h(12),'off');
ylabel(h(12),'v/v_{\alpha,init}','fontsize',14)
set(h(12),'xticklabel',[])
clim(h(12),[-2 2.5]);
colormap(h(12),cmap)
irf_legend(h(12),'(l)',[0.99 0.98],'fontsize',13,'color','k')
irf_legend(h(12),['\omega_{pp}t = ' num2str(tvalidx2,2)],[0.02 0.14],'fontsize',13,'color',[0 0.7 0])
axis(h(12),[0 64 -8 8])

pcolor(h(13),xpositions,vxppositions,log10(squeeze(pdistmatrix(idx2,:,:)))); 
shading(h(13),'flat')
hold(h(13),'on');
plot(h(13),xvalsall,mean(Exall(distindices(idx2)-60:distindices(idx2)+60,:),1)*10,'b','linewidth',2)
hold(h(13),'off');
ylabel(h(13),'v/v_{p,init}','fontsize',14)
xlabel(h(13),'n/\lambda_D','fontsize',14)
clim(h(13),[-2 2.5]);
colormap(h(13),cmap)
irf_legend(h(13),'(m)',[0.99 0.98],'fontsize',13,'color','k')
axis(h(13),[0 64 -8 8])

pcolor(h(14),xpositions,vxapositions,log10(squeeze(adistmatrix(idx3,:,:)))); 
shading(h(14),'flat')
hold(h(14),'on');
plot(h(14),xvalsall,mean(Exall(distindices(idx3)-60:distindices(idx3)+60,:),1)*10,'b','linewidth',2)
hold(h(14),'off');
ylabel(h(14),' ','fontsize',14)
set(h(14),'xticklabel',[])
set(h(14),'yticklabel',[])
clim(h(14),[-2 2.5]);
colormap(h(14),cmap)
c = colorbar(h(14),'position',[0.935 0.07 0.010 0.36]);
irf_legend(h(14),'(n)',[0.99 0.98],'fontsize',13,'color','k')
irf_legend(h(14),['\omega_{pp}t = ' num2str(tvalidx3,3)],[0.02 0.14],'fontsize',13,'color',[0 0.7 0])
c.Label.String = 'log_{10}f_i (s m^{-4})';
axis(h(14),[0 64 -8 8])


pcolor(h(15),xpositions,vxppositions,log10(squeeze(pdistmatrix(idx3,:,:))));
shading(h(15),'flat')
hold(h(15),'on');
plot(h(15),xvalsall,mean(Exall(distindices(idx3)-60:distindices(idx3)+60,:),1)*10,'b','linewidth',2)
hold(h(15),'off');
ylabel(h(15),' ','fontsize',14)
set(h(15),'xticklabel',[])
set(h(15),'yticklabel',[])
clim(h(15),[-2 2.5]);
colormap(h(15),cmap)
irf_legend(h(15),'(o)',[0.99 0.98],'fontsize',13,'color','k')
axis(h(15),[0 64 -8 8])

pcolor(h(16),xpositions,vxapositions,log10(squeeze(adistmatrix(idx4,:,:)))); 
shading(h(16),'flat')
hold(h(16),'on');
plot(h(16),xvalsall,mean(Exall(distindices(idx4)-60:distindices(idx4)+60,:),1)*10,'b','linewidth',2)
hold(h(16),'off');
ylabel(h(16),' ','fontsize',14)
set(h(16),'xticklabel',[])
set(h(16),'yticklabel',[])
clim(h(16),[-2 2.5]);
colormap(h(16),cmap)
irf_legend(h(16),'(p)',[0.99 0.98],'fontsize',13,'color','k')
irf_legend(h(16),['\omega_{pp}t = ' num2str(tvalidx4,3)],[0.02 0.14],'fontsize',13,'color',[0 0.7 0])
axis(h(16),[0 64 -8 8])

pcolor(h(17),xpositions,vxppositions,log10(squeeze(pdistmatrix(idx4,:,:)))); 
shading(h(17),'flat')
hold(h(17),'on');
plot(h(17),xvalsall,mean(Exall(distindices(idx4)-60:distindices(idx4)+60,:),1)*10,'b','linewidth',2)
hold(h(17),'off');
ylabel(h(17),' ','fontsize',14)
xlabel(h(17),'n/\lambda_D','fontsize',14)
set(h(17),'yticklabel',[])
clim(h(17),[-2 2.5]);
colormap(h(17),cmap)
irf_legend(h(17),'(q)',[0.99 0.98],'fontsize',13,'color','k')
axis(h(17),[0 64 -8 8])

set(h([10:17]),'Color',0.75*[1 1 1]);

set(h(1:17),'fontsize',13);
set(gcf,'color','w')
