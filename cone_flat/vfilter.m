function v=vfilter(g,DSD,theta,rDu,delt_u,rDv,delt_v)
% L=[-rDu:delt_u:rDu];
[cor_u,cor_v]=meshgrid([-rDu:1:rDu]*delt_u,[-rDv:1:rDv]*delt_v);
cor_u=cor_u';cor_v=cor_v';
w_u=w_function(pi,[-rDu:1:rDu]*delt_u);



v1=imfilter(permute(g,[2,1,3]),w_u','conv');
v1=permute(v1,[2,1,3]);

v=permute(v1,[2,3,1]).*delt_u.*(cor_u.^2+DSD.^2)/DSD;
v=permute(v*2*pi,[3,1,2]);
end


function u=u_function(s)
u=zeros(size(s));
    u(s==0)=1/2;
    v=s(s~=0);
u(s~=0)=(cos(v)-1)./(v.^2)+sin(v)./v;
end


function w_b=w_function(b,s)
w_b=u_function(b*s)*b^2/(4*pi^2);
end