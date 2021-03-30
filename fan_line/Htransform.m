function g4=Htransform(g3,u_cor)
g4=g3;
[nbeta,nu]=size(g3);
h_k=makefilter(u_cor);

delt_u=u_cor(2)-u_cor(1);

g4=imfilter(g3',h_k','conv')*delt_u;
g4=g4';


% for i=1:nbeta
%         tmp=conv(squeeze(g3(i,:)),h_k,'same');  
%         g4(i,:)=tmp;
% end
% g4=g4*delt_u;


function h_k=makefilter(u_cor)
% u_cor=[-rDu:1:rDu]*delt_u;
delt_u=u_cor(2)-u_cor(1);
b=1/2/delt_u;
h_k=hilbert_fun(b,u_cor);
% h_k(rDalpha+1)=(h_k(rDalpha)+h_k(rDalpha+2))/2;

function h_t=hilbert_fun(b,t)
h_t=(1-cos(2*pi*b*t))./(pi*t);
h_t(t==0)=0;

