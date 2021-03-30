function [u_star,v_star,w_star]=compute_star_uvw_cor(theta,DSD,DSO,cor_x,cor_y,cor_z)
fx=length(cor_x);fy=length(cor_y);fz=length(cor_z);
u_star=zeros(length(theta),fx,fy,fz);
% u_star=gpuArray(u_star);
v_star=u_star;
w_star=u_star;


% cor_x=gpuArray(cor_x);cor_y=gpuArray(cor_y);cor_z=gpuArray(cor_z);
[x_cor,y_cor,z_cor]=meshgrid(cor_x,cor_y,cor_z);
x_cor=permute(x_cor,[2,1,3]);y_cor=permute(y_cor,[2,1,3]);
for i=1:length(theta)
    eu=[-sin(theta(i)),cos(theta(i)),0];
    ew=[cos(theta(i)),sin(theta(i)),0];
%     eu=gpuArray(eu);
%     ev=gpuArray([0,0,1]);
%     ew=gpuArray(ew);
    w_star(i,:,:,:)=gather(DSO-(x_cor.*ew(1)+y_cor.*ew(2)));
    u_star(i,:,:,:)=gather(DSD.*(x_cor.*eu(1)+y_cor.*eu(2))./squeeze(w_star(i,:,:,:)));
    v_star(i,:,:,:)=gather(DSD.*z_cor./squeeze(w_star(i,:,:,:)));
%     i
end