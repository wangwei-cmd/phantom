function g2=weight_cor(g1,DSD,u_cor,v_cor)
% u_cor=[-rDu:1:rDu]*delt_u;
% v_cor=[-rDv:1:rDv]*delt_v;
[cor_u,cor_v]=meshgrid(u_cor,v_cor);
cor_u=cor_u';cor_v=cor_v';
dist=DSD./sqrt(DSD.^2+cor_u.^2+cor_v.^2);
g2=permute(g1,[2,3,1]).*dist;
% g2=permute(g1,[2,3,1]).*dist;
g2=permute(g2,[3,1,2]); 