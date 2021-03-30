function h_k=makefilter(rDu,delt_u)
u=[-rDu:1:rDu]*delt_u;
b=1/2/delt_u;
h_k=hilbert_fun(b,u);
% h_k(rDalpha+1)=(h_k(rDalpha)+h_k(rDalpha+2))/2;
end

function h_t=hilbert_fun(b,t)
h_t=(1-cos(2*pi*b*t))./(pi*t);
h_t(t==0)=0;
end