function [wmaxwpp,gmaxwpp,thetakamax,kldmax,vphmax,wppval,csval,ldval] = getmaxparameters(Teval,Vabval)
%GETMAXPARAMETERS Function to compute instability properties at peak growth
%rate for a given T_e and alpha speed. Written by D. B. Graham.

np = 10.0e6; % m-3
na = 1.0e6;
ne = np+2*na; % m-3
Te = Teval;
Ta = 12;% eV
Tp = 3; % eV


Units = irf_units; 
qe = Units.e;
me = Units.me;
mp = Units.mp;
ma = mp*4;
eps = Units.eps0;

vp = sqrt(2*qe*Tp/mp); % m/s
ve = sqrt(2*qe*Te/me);
va = sqrt(2*qe*Ta/ma);% m/s

Vab = Vabval;

wpe = sqrt(ne*qe^2/(me*eps));
wpa = sqrt(na*2^2*qe^2/(ma*eps));
wpp = sqrt(np*qe^2/(mp*eps));
ld = ve/wpe/sqrt(2);

kvecld = 0.005:0.01:2.5;
kvec = kvecld/ld;

theta = [0:1:90];

cs = ve/sqrt(2)*sqrt(me/mp)*(1+3*Tp/Te);
wr = zeros(1,length(kvec));
wi = zeros(1,length(kvec));

wrmat = zeros(length(theta),length(kvec));
wimat = zeros(length(theta),length(kvec));


for mm = 1:length(theta)
  
  
thetam = theta(mm);
kzvec = kvec*cosd(thetam);

for nn = 1:length(kvec)

    xip = @(omega)((omega)/(kvec(nn)*vp));
    xie = @(omega)((omega)/(kvec(nn)*ve));
    xia = @(omega)((omega - kzvec(nn)*Vab)/(kvec(nn)*va));
    disprel = @(w) (1+2*wpp^2/(kvec(nn)^2*vp^2)*(1 + 1i*sqrt(pi)*xip(w).*faddeeva(xip(w),50))+...
                     +2*wpe^2/(kvec(nn)^2*ve^2)*(1 + 1i*sqrt(pi)*xie(w).*faddeeva(xie(w),50))+...
                     +2*wpa^2/(kvec(nn)^2*va^2)*(1 + 1i*sqrt(pi)*xia(w).*faddeeva(xia(w),50)));
    
    if mm < 1000                 
      if (nn<10)    
          guessw = kvec(nn)*cs/2;
      else
          guessw = wr(nn-1)+1i*wi(nn-1)+1i*0.01*wpp;
      end
    else 
      guessw = wrmat(mm-1,nn) + 1i*(wimat(mm-1,nn))+1i*0.01*wpp;
    end
    
    options = optimoptions('fsolve','Display','off', 'TolFun', 1e-10, 'TolX', 1e-9,'MaxFunEvals', 1000);
    [x,FVAL,EXITFLAG] = fsolve(@(w) disprel(w), [guessw],options); 
    [~,maxpos] = max(x);
    wr(nn) = real(x(1));
    wi(nn) = imag(x(1));
    
end

wrmat(mm,:) = wr;
wimat(mm,:) = wi;

end


idxmax = find(wimat == max(max(wimat)));
[idxtheta,idxk] = ind2sub(size(wimat),idxmax);
thetakamax = theta(idxtheta);
kldmax = kvecld(idxk);
gmaxwpp = max(max(wimat))/wpp;
wmaxwpp = wrmat(idxtheta,idxk)/wpp;
vphmax = wmaxwpp*wpp*ld/kldmax;
wppval = wpp;
csval = cs;
ldval = ld;

if gmaxwpp < 0
  thetakamax = NaN;
  kldmax = NaN;
  wmaxwpp = NaN;
  vphmax = NaN;
end

 
end


