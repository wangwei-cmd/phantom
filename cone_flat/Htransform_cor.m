function g4=Htransform_cor(g3,u_cor)
delt_u=u_cor(2)-u_cor(1);
% k_h=filter(nu);
h_k=makefilter(u_cor);

g4=imfilter(permute(g3,[2,1,3]),h_k','conv');
g4=permute(g4,[2,1,3]);

g4=g4*delt_u;


function h_k=makefilter(u_cor)
% u_cor=[-rDu:1:rDu]*delt_u;
delt_u=u_cor(2)-u_cor(1);
b=1/2/delt_u;
h_k=hilbert_fun(b,u_cor);
% h_k(rDu+1)=(h_k(rDu)+h_k(rDu+2))/2;

function h_t=hilbert_fun(b,t)
h_t=(1-cos(2*pi*b*t))./(pi*t);
h_t(t==0)=0;

% function h_t=hilbert_fun(b,t)
% h_t=1./(pi*t);
% h_t(t==0)=0;