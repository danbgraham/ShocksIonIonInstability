%% Parameters
ne = 12e6; % m-3
Te = 100;

Ta = 12; % eV
Tp = 3; % eV


Units = irf_units; 
qe = Units.e;
me = Units.me;
mp = Units.mp;
ma = mp*4;
eps = Units.eps0;
kb = Units.kB;

ve = sqrt(2*qe*Te/me);
va = sqrt(2*qe*Ta/ma);% m/s
vp = sqrt(2*qe*Tp/mp); % m/s
cs = ve/sqrt(2)*sqrt(me/mp)*(1+3*Tp/Te);

Vab = 100e3;

wpe = sqrt(ne*qe^2/(me*eps));
ld = ve/wpe/sqrt(2);

kvecld = 0.01:0.001:1.5;
kvec = kvecld/ld;

%wIA = cs*kvec;

wr1 = zeros(1,length(kvec));
wi1 = zeros(1,length(kvec));

dispreltemp = zeros(size(wr1));


nanprat = 10.^[-3:0.05:1];

gammamax = zeros(size(nanprat));
kmax = zeros(size(nanprat));
wmax = zeros(size(nanprat));
wpparr = zeros(size(nanprat));
wpaarr = zeros(size(nanprat));
naarr = zeros(size(nanprat));
nparr = zeros(size(nanprat));

for ii = 1:length(nanprat)
    
    % Define function
    np = 12e6/(2*nanprat(ii) + 1); % m-3
    na = (12e6-np)/2; % m-3
    wpa = sqrt(na*2^2*qe^2/(ma*eps));
    wpp = sqrt(np*qe^2/(mp*eps));
    wpparr(ii) = wpp;
    wpaarr(ii) = wpa;
    nparr(ii) = np;
    naarr(ii) = na;

    for nn = 1:length(kvec)

      xip = @(omega)((omega)/(kvec(nn)*vp));
      xie = @(omega)((omega)/(kvec(nn)*ve));
      xia = @(omega)((omega - kvec(nn)*Vab)/(kvec(nn)*va));
      disprel = @(w) (1+2*wpp^2/(kvec(nn)^2*vp^2)*(1 + 1i*sqrt(pi)*xip(w).*faddeeva(xip(w),50))+...
        +2*wpe^2/(kvec(nn)^2*ve^2)*(1 + 1i*sqrt(pi)*xie(w).*faddeeva(xie(w),50))+...
        +2*wpa^2/(kvec(nn)^2*va^2)*(1 + 1i*sqrt(pi)*xia(w).*faddeeva(xia(w),50)));

      if ii < 2
        if (nn<5) 
          guessw1 = cs*kvec(nn)/5;
        else
          guessw1 = wr1(nn-1)+1i*wi1(nn-1);
        end
      else
        if (nn<5) 
          guessw1 = dispreltemp(nn);
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
    
    dispreltemp = wr1+1i*wi1;
    
    [wmaxtemp,idxgam] = max(wi1);
    if wmaxtemp < 0
      gammamax(ii) = NaN;
      kmax(ii) = NaN;
      wmax(ii) = NaN;
    else
      gammamax(ii) = wi1(idxgam);
      kmax(ii) = kvec(idxgam);
      wmax(ii) = wr1(idxgam);
    end
    ii
end

%%
disprelnanprat = struct('ne',ne,'npvec',nparr,'navec',naarr,'Te',Te,'Tpvec',Tp,'Tavec',Ta,'Vab',Vab,'ld',ld,...
  'wpe',wpe,'wppvec',wpparr,'wpavec',wpaarr,'kvec',kvec,'gammamax',gammamax,'kmax',kmax,'wmax',wmax);
save('disprelnanprat.mat','disprelnanprat');

%%

semilogx(nanprat,gammamax./wpparr,nanprat,wmax./wpparr,nanprat,kmax*ld,'linewidth',2)
xlabel('n_\alpha/n_p','fontsize',14)
ylabel('\gamma_{max}/\omega_{pp}, \omega_{max}/\omega_{pp}, k_{max} \lambda_D','fontsize',14)
legend({'\gamma','\omega','k'},'location','northwest')
set(gca,'fontsize',14)
set(gcf,'color','w')