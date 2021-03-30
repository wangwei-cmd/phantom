
clear;
addpath(genpath('../astra-1.9.0.dev11-matlab-win-x64'));
% f=dicomread('/home/wangwei/sdc/AMP/L067/full_1mm_sharp/L067_FD_1_SHARP_1.CT.0002.0523.2016.01.21.18.11.40.977560.404643926.IMA');
% f=double(f);
% f=f/max(f(:));

f=dicomread('1.IMA');
f=double(f);

[fx,fy]=size(f);
rx=floor(fx/2);
ry=floor(fy/2);
x_cor=[-rx:fx-rx-1]+0.5;
y_cor=[-ry:fy-ry-1]+0.5;

delt_beta=1;
R=1000;
DSD=ceil(R+sqrt(2)*256);
rr=sqrt(2)*(fx/2);
tan1=rr/R;
delt_u=1;
rDu=ceil(DSD*tan1)/delt_u;
 
u_cor=[-rDu:rDu]*delt_u;

beta=[0:1:180+2*ceil(atand(tan1))]*pi/180;

N=512;
vol_geom = astra_create_vol_geom(N, N,x_cor(1),x_cor(end),y_cor(1),y_cor(end));
proj_geom = astra_create_proj_geom('fanflat', 1.0, length(u_cor), beta, R, DSD-R);
proj_id =astra_create_projector('line_fanflat', proj_geom, vol_geom);
% [sinogram_id, sin] = astra_create_sino(f, proj_id);
[sinogram_id, sin] = astra_create_sino_cuda(f, proj_geom,vol_geom);
% sin=projection_fan(f,u_cor,x_cor,y_cor,beta,R,DSD);

parker=readNPY('parker_weighting.npy');
data=sin.*parker;
% data=sin;
sinogram_id = astra_mex_data2d('create', '-sino', proj_geom,data);
rec_id = astra_mex_data2d('create', '-vol', vol_geom,0);
cfg = astra_struct('FBP_CUDA');
cfg.ReconstructionDataId = rec_id;
cfg.ProjectionDataId = sinogram_id;
cfg.ProjectorId = proj_id;
% cfg.option.FilterType = 'Ram-Lak';
% cfg.option.ShortScan='yes';
alg_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('run', alg_id);
rf1 = astra_mex_data2d('get', rec_id);
astra_mex_algorithm('delete', alg_id);

d=5*pi/180;
W=weight_ace(beta,u_cor,DSD,d);
s_fan=sin.*W;
rf2=back_ace2(s_fan,u_cor,x_cor,y_cor,R,DSD,beta)*2;


theta1=solve_chord(beta(1),x_cor,y_cor,R);
[index1]=index_theta(beta(1)*ones(size(theta1)),theta1,beta);
theta2=solve_chord(beta(end),x_cor,y_cor,R);
[index2]=index_theta_1(theta2,beta(end)*ones(size(theta2)),beta);
rf3=back_ace2_index(sin,(index1+index2)/2,u_cor,x_cor,y_cor,R,DSD,beta);


figure;imshow(f,[]);
figure;imshow(rf1,[])
figure;imshow(rf2,[]);
figure;imshow(rf3,[]);
% psnr(rf1,f,max(f(:)))
% psnr(rf2,f,max(f(:)))
% psnr(rf3,f,max(f(:)))
calrmse(rf1,f)
calrmse(rf2,f)
calrmse(rf3,f)
ssim(uint8(rf1/max(f(:))*255),uint8(f/max(f(:))*255))
ssim(uint8(rf2/max(f(:))*255),uint8(f/max(f(:))*255))
ssim(uint8(rf3/max(f(:))*255),uint8(f/max(f(:))*255))

% figure;imshow(ls1(f),[0,0.6]);
figure;imshow(ls1(rf1),[0,0.6])
figure;imshow(ls1(rf2),[0,0.6]);
figure;imshow(ls1(rf3),[0,0.6]);
figure;imshow(ls1(abs(f-rf1)),[0.,0.4])
figure;imshow(ls1(abs(f-rf2)),[0.,0.4])
figure;imshow(ls1(abs(f-rf3)),[0.,0.4])

rmpath(genpath('G:\GitLab\circle_cone_beam_walnut\astra-1.9.0.dev11-matlab-win-x64'));