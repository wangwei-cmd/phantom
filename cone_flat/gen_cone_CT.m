clear;
doc='E:\CT_image\AMP\circle_conebeam\train';
type='.IMA';
ct_name=glob(doc,type);
L=50;
num=50;
assert(L*num<=length(ct_name));

f=(zeros(512,512,num,L));

fx=512;fy=512;fz=50;
rx=floor(fx/2);
ry=floor(fy/2);
rz=floor(fz/2);
rr=sqrt(2)*max(rx,ry);
DSO=1000;
DSD=ceil(1000+sqrt(2)*256);

tan1=rr/DSO;
rDu=ceil(DSD*tan1);
tan2=rz/(DSO-rr);
rDv=ceil(DSD*tan2);

delt_u=1;
delt_v=1;
delt_t=1;
x_cor=[-rx:1:fx-rx-1];y_cor=[-ry:1:fy-ry-1];z_cor=[-rz:1:fz-rz-1];
u_cor=[rot90(0:-delt_u:-rDu,2),delt_u:delt_u:rDu];
v_cor=[rot90(0:-delt_v:-rDv,2),delt_v:delt_v:rDv];

% theta=[0:1:180+2*ceil(atand(tan1))]*pi/180;
% post='short_scan';

theta=[0:1:180]*pi/180;
post='super_short_scan';

sin=zeros(length(u_cor),length(v_cor),length(theta),L);

for i=1:1
    for j=1:num
    f(:,:,j,i)=(double(dicomread(ct_name{i*num+j-num})));
    end
    tmp=projection_circone_cor(f(:,:,:,i),theta,DSD,DSO,...
      u_cor,v_cor,x_cor,y_cor,z_cor,delt_t);
    sin(:,:,:,i)=gather(tmp);
end
save(['cone_ct_',post,'.mat'],'f','sin');