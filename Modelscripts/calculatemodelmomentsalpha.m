function momentsstruct = calculatemodelmomentsalpha(histstruct)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

xposp = histstruct.xpositions;
vxvec = histstruct.vxpositions;
vyvec = histstruct.vypositions;
vzvec = histstruct.vzpositions;

dv = histstruct.dv;

vxmat = ones(size(xposp'))*vxvec;
vymat = ones(size(xposp'))*vyvec;
vzmat = ones(size(xposp'))*vzvec;

fallvx = histstruct.n1Dvx;
fallvy = histstruct.n1Dvy;
fallvz = histstruct.n1Dvz;

Units = irf_units;
e = Units.e;
mi = Units.mp*4;

np = squeeze(sum(fallvx,2))*dv;

vx = sum(fallvx.*vxmat*dv,2)./np;
vy = sum(fallvy.*vymat*dv,2)./np;
vz = sum(fallvz.*vzmat*dv,2)./np;

vxbulkmat = vx*ones(size(vxvec));
vybulkmat = vy*ones(size(vyvec));
vzbulkmat = vz*ones(size(vzvec));

Pxx = mi*sum(fallvx.*(vxmat - vxbulkmat).^2*dv,2);
Pyy = mi*sum(fallvy.*(vymat - vybulkmat).^2*dv,2);
Pzz = mi*sum(fallvz.*(vzmat - vzbulkmat).^2*dv,2);

Ptotxx = mi*sum(fallvx.*(vxmat).^2*dv,2);

Txx = Pxx./np/e;
Tyy = Pyy./np/e;
Tzz = Pzz./np/e;

Pdynx = mi.*np.*vx.^2;

Ts = (Txx + Tyy + Tzz)/3;



momentsstruct = struct('np',np,'vx',vx,'vy',vy,'vz',vz,'Ts',Ts,'xpos',xposp,'Txx',Txx,'Tyy',Tyy,'Tzz',Tzz,'Pxx',Pxx,'Pyy',Pyy,'Pzz',Pzz,'Pdynx',Pdynx,'Ptotxx',Ptotxx);

end

