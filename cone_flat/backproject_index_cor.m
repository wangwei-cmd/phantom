function rf=backproject_index_cor(g4,index,u_star,v_star,w_star,theta,u_cor,v_cor)

[ntheta,nx,ny,nz]=size(v_star);

theta_q=repmat(theta',1,nx,ny,nz);
[theta_cor,u_cor,v_cor]=meshgrid(theta,u_cor,v_cor);
g5=interp3(theta_cor,u_cor,v_cor,permute(g4,[2,1,3]),theta_q,u_star,v_star);

% g5=zeros(ntheta,nx,ny,nz);
% [u_cor,v_cor]=meshgrid(u_cor,v_cor);
% for i=1:ntheta
%     tmp=interp2(u_cor,v_cor,squeeze(g4(i,:,:))',gpuArray(u_star(i,:,:,:)),gpuArray(v_star(i,:,:,:)));
%     g5(i,:,:,:)=gather(tmp);
% end


g5(isnan(g5))=0;
f=index.*g5./w_star;
rf=squeeze(sum(f,1))/2/pi;
end