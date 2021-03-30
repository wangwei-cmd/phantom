function [u_star,v_star]=compute_u(beta,x_cor,y_cor,R,DSD)
% x_cor=-floor(fx/2):fx-floor(fx/2)-1;
% y_cor=-floor(fy/2):fy-floor(fy/2)-1;
fx=length(x_cor);fy=length(y_cor);
[x_cor,y_cor]=meshgrid(x_cor,y_cor);
x_cor=x_cor';y_cor=y_cor';
u_star=zeros(length(beta),fx,fy);
v_star=u_star;
for i=1:length(beta)
    v_star(i,:,:)=R-x_cor.*cos(beta(i))-y_cor.*sin(beta(i));
    u_star(i,:,:)=DSD*(-x_cor.*sin(beta(i))+y_cor.*cos(beta(i)))./squeeze(v_star(i,:,:));
end
    