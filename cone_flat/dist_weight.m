function g2=dist_weight(pf,DSD,rDu,rDv,delt_u,delt_v)

dist=weight(DSD,rDu,rDv,delt_u,delt_v);
g2=pf.*dist;
end


function dist=weight(DSD,rDu,rDv,delt_u,delt_v)
% pf_weight=gpuArray(pf);
u=[-rDu:delt_u:rDu];
v=[-rDv:delt_v:rDv];
[uu,vv]=meshgrid(u,v);
uu=uu';
vv=vv';
dist=DSD./sqrt(uu.^2+vv.^2+DSD.^2);
end