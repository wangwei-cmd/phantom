function w=weight_ace(beta,u_cor,DSD,d)
beta_s=beta(1);
beta_e=beta(end);
w=zeros(length(beta),length(u_cor));
t1=c_function(beta,beta_s,beta_e,d);
for i=1:length(beta)
    theta=mod(beta(i)+pi-2*atan(u_cor./DSD)',2*pi);
    t2=c_function(theta,beta_s,beta_e,d);
    w(i,:)=t1(i)./(t1(i)+t2);
end
end

function c=c_function(beta,beta_s,beta_e,d)
beta_1=beta_s+d;
beta_2=beta_e-d;
id1=(beta>=beta_s)&(beta<beta_1);
id2=(beta>=beta_1)&(beta<beta_2);
id3=(beta>=beta_2)&(beta<=beta_e);
id4=(beta>beta_e)|(beta<beta_s);
assert(sum(id1+id2+id3+id4==1)==length(id1));
c=zeros(1,length(beta));
tmp1=(cos(pi*(beta-beta_s-d)/2/d)).^2;
c(id1)=tmp1(id1);
c(id2)=1;
tmp3=(cos(pi*(beta-beta_s+d)/2/d)).^2;
c(id3)=tmp3(id3);
c(id4)=0;
end