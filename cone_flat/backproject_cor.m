function rf=backproject_cor(g4,u_star,v_star,w_star,theta,u_cor,v_cor)

[ntheta,nx,ny,nz]=size(v_star);

theta_q=repmat(theta',1,nx,ny,nz);
[theta_cor,u_cor,v_cor]=meshgrid(theta,u_cor,v_cor);
g5=interp3(theta_cor,u_cor,v_cor,permute(g4,[2,1,3]),theta_q,u_star,v_star);
g5(isnan(g5))=0;
f=g5./w_star;
rf=squeeze(sum(f,1))/2/pi*(theta(2)-theta(1))/2;
end