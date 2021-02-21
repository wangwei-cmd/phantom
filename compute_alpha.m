function [alpha_star,v_star]=compute_alpha(beta,x_cor,y_cor,R)
fx=length(x_cor);fy=length(y_cor);
x_cor=gpuArray(x_cor);
y_cor=gpuArray(y_cor);
[x_cor,y_cor]=meshgrid(x_cor,y_cor);
x_cor=x_cor';y_cor=y_cor';
alpha_star=zeros(length(beta),fx,fy);
alpha_star=gpuArray(alpha_star);
v_star=alpha_star;
for i=1:length(beta)
    v_star(i,:,:)=R-x_cor.*cos(beta(i))-y_cor.*sin(beta(i));
    alpha_star(i,:,:)=atan((-x_cor.*sin(beta(i))+y_cor.*cos(beta(i)))./squeeze(v_star(i,:,:)));
end