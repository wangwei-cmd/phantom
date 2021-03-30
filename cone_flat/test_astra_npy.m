clear;
img=load('./img_512x512x25.mat');
img=double(img.img);
img=img(:,:,1:25);

% addpath(genpath('G:\GitLab\circle_cone_beam_walnut/astra-1.9.0.dev11-matlab-win-x64'));
[fx,fy,fz]=size(img);

rx=floor(fx/2);
ry=floor(fy/2);
rz=floor(fz/2);
rr=sqrt(2)*max(rx,ry);
% DSO=1000;
DSO=500;
DSD=ceil(DSO+sqrt(2)*256);

delt_u=1;
delt_v=1;
delt_t=1;

tan1=rr/sqrt(DSO^2-rr^2);
rDu=ceil(DSD*tan1)/delt_u;
tan2=rz/(DSO-rr);
rDv=ceil(DSD*tan2)/delt_v;

beta=[0:1:180+2*ceil(atand(tan1))]*pi/180;
% beta=[0:1:359]*pi/180;

x_cor=[-rx:1:fx-rx-1]+0.5;y_cor=[-ry:1:fy-ry-1]+0.5;z_cor=[-rz:1:fz-rz-1];
u_cor=[-rDu:rDu]*delt_u;
v_cor=[-rDv:rDv]*delt_v;


% pf=projection_circone_cor(img,beta,DSD,DSO,u_cor,v_cor,x_cor,y_cor,z_cor,delt_t);

det_row_count=length(v_cor);
det_col_count=length(u_cor);
theta=beta;

% source_origin=DSO;
% origin_det=DSD-DSO;
% vol_geom = astra_create_vol_geom(fy,fx,fz);
% proj_geom = astra_create_proj_geom('cone',delt_v,delt_u, det_row_count,det_col_count,theta,source_origin,origin_det);
% [proj_id, sin] = astra_create_sino3d_cuda(img, proj_geom, vol_geom);

rf1=readNPY('ct.npy');
pf=readNPY('sin.npy');
% pf=gpuArray(readNPY('sin.npy'));



[u_star,v_star,w_star]=compute_star_uvw_cor(beta,DSD,DSO,x_cor,y_cor,z_cor);
NN=300;
extend_L=[u_cor(1)-NN*delt_u:delt_u:u_cor(1)-delt_u];
extend_R=[u_cor(end)+delt_u:delt_u:u_cor(end)+NN*delt_u];
u_cor_extend=[extend_L,u_cor,extend_R];%%%%Hilbert filter extend


d=8*pi/180;
W=weight_ace(beta,u_cor,DSD,d);
pf1=pf.*W;
g1=compute_grad_s_cor(pf1,DSD,u_cor,v_cor,beta);
g2=weight_cor(g1,DSD,u_cor,v_cor);
g2=padarray(g2,[0 NN],0,'both');
g4=Htransform_cor(g2,u_cor_extend);
g4=gather(g4);
rf2=backproject_cor(g4,u_star,v_star,w_star,beta,u_cor_extend,v_cor)*2;
rf2=gather(rf2);
rf2=rot90(rf2,-1);
% figure;imshow(rf2(:,:,8),[])

[cor_x,cor_y]=meshgrid(x_cor,y_cor);
cor_x=cor_x';cor_y=cor_y';
theta1=solve_chord(beta(1),cor_x,cor_y,DSO);
[index1]=index_theta(beta(1)*ones(size(theta1)),theta1,beta);
theta2=solve_chord(beta(end),cor_x,cor_y,DSO);
[index2]=index_theta_1(theta2,beta(end)*ones(size(theta2)),beta);
index=(index1+index2)/2;

g1=compute_grad_s_cor(pf,DSD,u_cor,v_cor,beta);
g2=weight_cor(g1,DSD,u_cor,v_cor);
g2=padarray(g2,[0 NN],0,'both');
g4=Htransform_cor(g2,u_cor_extend);
g4=gather(g4);
rf3=backproject_index_cor(g4,index,u_star,v_star,w_star,beta,u_cor_extend,v_cor);
rf3=gather(rf3);
rf3=rot90(rf3,-1);


% figure;imshow(img(:,:,8),[])
% figure;imshow(rf1(:,:,8),[])
% figure;imshow(rf2(:,:,8),[])
% figure;imshow(rf3(:,:,8),[])
figure;imshow(ls1(rf1(:,:,8)),[0,0.6])
figure;imshow(ls1(rf2(:,:,8)),[0,0.6])
figure;imshow(ls1(rf3(:,:,8)),[0,0.6])
figure;imshow(ls1(abs(img(:,:,8)-rf1(:,:,8))),[0.,0.4])
figure;imshow(ls1(abs(img(:,:,8)-rf2(:,:,8))),[0.,0.4])
figure;imshow(ls1(abs(img(:,:,8)-rf3(:,:,8))),[0.,0.4])
for i=1:size(img,3)
    M=max(max(img(:,:,i)));
pp1(i)=psnr(rf1(:,:,i),img(:,:,i),M);
qq1(i)=ssim(rf1(:,:,i)/M,img(:,:,i)/M);
e1(i)=calrmse(rf1(:,:,i),img(:,:,i));
pp2(i)=psnr(rf2(:,:,i),img(:,:,i),M);
qq2(i)=ssim(rf2(:,:,i)/M,img(:,:,i)/M);
e2(i)=calrmse(rf2(:,:,i),img(:,:,i));
pp3(i)=psnr(rf3(:,:,i),img(:,:,i),M);
qq3(i)=ssim(rf3(:,:,i)/M,img(:,:,i)/M);
e3(i)=calrmse(rf3(:,:,i),img(:,:,i));
end
pp1
pp2
pp3
sum(e1)/length(e1)
sum(e2)/length(e2)
sum(e3)/length(e3)