function [g_s,g_q,g_u,g_v]=compute_grad_s_cor(pf,D,cor_u,cor_v,theta)
delt_theta=theta(2)-theta(1);
delt_u=cor_u(2)-cor_u(1);
delt_v=cor_v(2)-cor_v(1);
[ns,nu,nv]=size(pf);
% ru=floor(nu/2);
% rv=floor(nv/2);
[u_cor,v_cor]=meshgrid(cor_u,cor_v);
u_cor=u_cor';
v_cor=v_cor';
% pf=permute(pf,[3,1,2]);

% [g_q,g_u,g_v]=gradient(pf);

g_q=zeros(size(pf));g_u=zeros(size(pf));g_v=zeros(size(pf));
g_q=gpuArray(g_q);g_u=gpuArray(g_u);g_v=gpuArray(g_v);
g_q(2:end-1,:,:)=(pf(3:end,:,:)-pf(1:end-2,:,:))./2/delt_theta;
g_u(:,2:end-1,:)=(pf(:,3:end,:)-pf(:,1:end-2,:))./2/delt_u;
g_v(:,:,2:end-1)=(pf(:,:,3:end)-pf(:,:,1:end-2))./2/delt_v;


% g_q=zeros(size(pf));g_u=zeros(size(pf));g_v=zeros(size(pf));
% g_q(1:end-1,:,:)=(pf(2:end,:,:)-pf(1:end-1,:,:))./delt_theta;
% g_u(:,1:end-1,:)=(pf(:,2:end,:)-pf(:,1:end-1,:))./delt_u;
% g_v(:,:,1:end-1)=(pf(:,:,2:end)-pf(:,:,1:end-1))./delt_v;


g_s=pf;
for i=1:ns
g_s(i,:,:)=squeeze(g_q(i,:,:))+(u_cor.^2+D.^2).*squeeze(g_u(i,:,:))./D...
+u_cor.*v_cor.*squeeze(g_v(i,:,:))./D;
end
