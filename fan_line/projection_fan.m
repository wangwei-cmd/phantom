function pf=projection_fan(f,u_cor,x_cor,y_cor,beta,R,DSD)
[fx,fy]=size(f);
% u_cor=[-rDu:rDu]*delt_u;
% alpha=atan(u./DSD);
Lb=length(beta);
La=length(u_cor);
vk=ceil(sqrt(2)*max(fx,fy));
rvk=floor(vk/2);

pf=zeros(Lb,La);

[s_cor,u_cor]=meshgrid(beta,u_cor);
s_cor=s_cor';u_cor=u_cor';
LL=sqrt(DSD^2+u_cor.^2);
for l=R-rvk-1:R+rvk+1  
    cor_x=R*cos(s_cor)-l*(u_cor.*sin(s_cor)+DSD.*cos(s_cor))./LL;
    cor_y=R*sin(s_cor)+l*(u_cor.*cos(s_cor)-DSD.*sin(s_cor))./LL;
    tmp=interp2(x_cor,y_cor,f',cor_x,cor_y);
    tmp(isnan(tmp))=0;
    pf=pf+tmp;
end
       