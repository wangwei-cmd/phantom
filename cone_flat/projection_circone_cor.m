function pf=projection_circone_cor(img,theta,DSD,DSO,u_cor,v_cor,x_cor,y_cor,z_cor,delt_t)
nproj=length(theta);
[x,y,z]=size(img);
[cor_x,cor_y,cor_z]=meshgrid(x_cor,y_cor,z_cor);
pf=zeros(length(u_cor),length(v_cor),length(theta));
pf=gpuArray(pf);
[uu,vv]=meshgrid(u_cor,v_cor);
uu=uu';
vv=vv';
rr=ceil(sqrt(x.^2+y.^2)/2);
tt=[DSO-rr:delt_t:DSO+rr];
x0=DSO*cos(theta);
y0=DSO*sin(theta);

tt=repmat(tt',[1,size(uu,1),size(uu,2)]);
tt=permute(tt,[2,3,1]);
tt=gpuArray(tt);
x0=gpuArray(x0);
y0=gpuArray(y0);
uu=gpuArray(uu);
vv=gpuArray(vv);
cor_x=gpuArray(cor_x);cor_y=gpuArray(cor_y);cor_z=gpuArray(cor_z);
for i=1:nproj
    eu=[-sin(theta(i)),cos(theta(i)),0];
    ew=[cos(theta(i)),sin(theta(i)),0];
    xx=eu(1).*uu-DSD.*ew(1);
    yy=eu(2).*uu-DSD.*ew(2);
    zz=vv;
    LL=sqrt(uu.^2+vv.^2+DSD.^2);
    xq=x0(i)+tt.*xx./LL;
    yq=y0(i)+tt.*yy./LL;
    zq=tt.*zz./LL;
    tmp=interp3(cor_x,cor_y,cor_z,permute(img,[2,1,3]),xq,yq,zq);
    tmp(isnan(tmp))=0;
    pf(:,:,i)=sum(tmp,3)*delt_t;
end
% pf=gather(pf);
pf=permute(pf,[3,1,2]);
end
