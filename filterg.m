function v=filterg(g,alpha)

cos_alpha=cos(alpha);
sin_alpha=sin(alpha);
h=alpha(2)-alpha(1);
w_b=w_function(pi/h,alpha);
t1=(alpha./sin_alpha).^2;
w_c_sin=t1.*w_b;
w_c_sin(alpha==0)=w_b(alpha==0);

[gx,~]=size(g);
v=g;
for i=1:gx 
    v(i,:)=conv(g(i,:).*cos_alpha,w_c_sin,'same');
end
v=h*v;
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
