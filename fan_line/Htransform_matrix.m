function g4=Htransform_matrix(g3,u_cor)
g4=g3;
[nbeta,nu]=size(g3);
h_matrix=makefilter_matrix(u_cor);


delt_u=u_cor(2)-u_cor(1);
% g4=imfilter(g3',h_k','conv')*delt_u;
g4=h_matrix*g3'*delt_u;
g4=g4';


% for i=1:nbeta
%         tmp=conv(squeeze(g3(i,:)),h_k,'same');  
%         g4(i,:)=tmp;
% end
% g4=g4*delt_u;


function h_matrix=makefilter_matrix(u_cor)
% u_cor=[-rDu:1:rDu]*delt_u;
delt_u=u_cor(2)-u_cor(1);
L=length(u_cor);
uu=[-L+1:L-1]*delt_u;
b=1/2/delt_u;
% sin_u=sin(uu);
h_k=hilbert_fun(b,uu);
% h_k=u_cor./sin_u.*hilbert_fun(b,u_cor);
h_k(isnan(h_k))=0;
h_matrix=zeros(L,L);
for i=1:L
    h_matrix(i,:)=h_k(i+L-1:-1:i);
end



function h_t=hilbert_fun(b,t)
h_t=(1-cos(2*pi*b*t))./(pi*t);
h_t(t==0)=0;