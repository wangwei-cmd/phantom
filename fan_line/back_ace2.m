function rf=back_ace2(pf,u_cor,x_cor,y_cor,R,DSD,beta)
% u_cor=[-rDu:rDu]*delt_u;
delt_u=u_cor(2)-u_cor(1);
dist=DSD./sqrt(u_cor.^2+DSD^2);
dist1=(u_cor.^2+DSD^2)./DSD;
% alpha=atan(u./DSD);
pf=pf.*dist;
d_beta=zeros(size(pf));
d_u=zeros(size(pf));
d_beta(2:end-1,:,:)=(pf(3:end,:,:)-pf(1:end-2,:,:))/2/(beta(2)-beta(1));
d_u(:,2:end-1)=(pf(:,3:end)-pf(:,1:end-2))/(2*delt_u).*dist1(2:end-1);
pf_grad=d_beta+d_u;

NN=30;
extend_L=[u_cor(1)-NN*delt_u:delt_u:u_cor(1)-delt_u];
extend_R=[u_cor(end)+delt_u:delt_u:u_cor(end)+NN*delt_u];
u_cor_extend=[extend_L,u_cor,extend_R];%%%%Hilbert filter extend
pf_grad=padarray(pf_grad,[0 NN],0,'both');
g4=Htransform(pf_grad,u_cor_extend);
g4=g4(:,length(extend_L)+1:length(extend_L)+length(u_cor));

[u_star,v_star]=compute_u(beta,x_cor,y_cor,R,DSD);
g5=zeros(size(u_star));
for i=1:size(u_star,1)
    g5(i,:,:)=interp1(u_cor,squeeze(g4(i,:)),u_star(i,:,:));
end
g5=g5./v_star;
g5(isnan(g5))=0;
rf=squeeze(sum(g5,1))/2/pi/2*(beta(2)-beta(1));


        

