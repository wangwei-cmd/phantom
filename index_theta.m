function [index]=index_theta(s_b,s_t,theta)
s_t(s_b>s_t)=s_t(s_b>s_t)+2*pi;
delt_theta=theta(2)-theta(1);
[nx,ny,nz]=size(s_b);
ntheta=length(theta);
index=zeros(ntheta,nx,ny,nz);
index=gpuArray(index);
index_R=(s_b-theta(1))/delt_theta;
s_b_R=ceil(index_R)*delt_theta+theta(1);
index_L=(s_t-theta(1))/delt_theta;
s_t_L=floor(index_L)*delt_theta+theta(1);
theta(theta<s_b(1,1))=theta(theta<s_b(1,1))+2*pi;
for i=1:ntheta
id=(s_b_R<=theta(i))&(theta(i)<=s_t_L);
tmp=zeros(nx,ny,nz);
tmp=gpuArray(tmp);
tmp(id)=delt_theta;
id1=(s_b<=theta(i))&(theta(i)<s_b_R);
tmp(id1)=s_b_R(id1)-s_b(id1);
id2=(s_t_L<theta(i))&(theta(i)<=s_t);
tmp(id2)=s_t(id2)-s_t_L(id2);
assert(sum(id+id1+id2>1,'all')==0);
index(i,:,:,:)=tmp;
end
bad=(s_t>theta(end));
LL=s_t-s_b;
index(bad)=index(bad)./LL(bad);
index=gather(index);