function pf=projection_fan_circle(f,cor_alpha,x_cor,y_cor,beta,R)
% f=gpuArray(f);
% cor_alpha=gpuArray(cor_alpha);
% beta=gpuArray(beta);
% [fx,fy]=size(f);
% alpha=[-rDalpha:rDalpha]*delt_alpha;

% alpha=atan(u./DSD);
Lb=length(beta);
La=length(cor_alpha);
rvk=ceil(sqrt(2)*max(max(abs(x_cor),abs(y_cor))));

delt_t=min(x_cor(2)-x_cor(1),y_cor(2)-y_cor(1));
delt_t=abs(delt_t);
pf=zeros(Lb,La);
pf=gpuArray(pf);
[s_cor,alpha_cor]=meshgrid(beta,cor_alpha);
s_cor=s_cor';alpha_cor=alpha_cor';
for l=R-rvk-1:delt_t:R+rvk+1  
    cor_x=R*cos(s_cor)-l*cos(alpha_cor-s_cor);
    cor_y=R*sin(s_cor)+l*sin(alpha_cor-s_cor);
    tmp=interp2(x_cor,y_cor,f',cor_x,cor_y);
    tmp(isnan(tmp))=0;
    pf=pf+tmp*delt_t;
end
% pf=gather(pf);       