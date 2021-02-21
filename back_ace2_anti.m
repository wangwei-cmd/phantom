function rf=back_ace2_anti(pf,alpha_cor,x_cor,y_cor,R,beta)
% alpha_cor=[-rDalpha:rDalpha]*delt_alpha;
delt_alpha=alpha_cor(2)-alpha_cor(1);
% alpha=atan(u./DSD);
d_beta=zeros(size(pf));
d_u=zeros(size(pf));
d_beta=gpuArray(d_beta);
d_u=gpuArray(d_u);
d_beta(2:end-1,:,:)=(pf(3:end,:,:)-pf(1:end-2,:,:))/2/(beta(2)-beta(1));
d_u(:,2:end-1)=(pf(:,3:end)-pf(:,1:end-2))/(2*delt_alpha);
pf_grad=d_beta+d_u;
% alpha_cor_shift=alpha_cor+delt_alpha*0.5;
NN=300;
extend_L=[alpha_cor(1)-NN*delt_alpha:delt_alpha:alpha_cor(1)-delt_alpha];
extend_R=[alpha_cor(end)+delt_alpha:delt_alpha:alpha_cor(end)+NN*delt_alpha];
alpha_cor_extend=[extend_L,alpha_cor,extend_R];%%%%Hilbert filter extend
pf_grad=padarray(pf_grad,[0 NN],0,'both');
g4=Htransform(pf_grad,alpha_cor_extend);
g4=g4.*cos(alpha_cor_extend);
[alpha_star,v_star]=compute_alpha(beta,x_cor,y_cor,R);
g5=zeros(size(alpha_star));
g5=gpuArray(g5);
for i=1:size(alpha_star,1)
    g5(i,:,:)=interp1(alpha_cor_extend,squeeze(g4(i,:)),alpha_star(i,:,:));
end
g5=g5./v_star;
g5(isnan(g5))=0;
rf=squeeze(sum(g5,1))/2/pi/2*(beta(2)-beta(1));
rf=rot90(rf,-1);
rf=fliplr(rf);