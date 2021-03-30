function [index,p1]=index_theta_2(s_start,chord,theta)
ntheta=length(theta);
s_b=theta(1)*ones(size(chord));
s_t=theta(end)*ones(size(chord));
p1=(chord>=s_b)&(chord<=s_t);
p2=logical(1-p1);
s_start=s_start*ones(size(chord));


p11=logical((s_start>chord).*p1);
p12=logical((s_start<=chord).*p1);
index11=index_theta(s_start,chord,theta);
index12=index_theta_1(chord,s_start,theta);
p2=repmat(p2,[1,1,ntheta]);
p11=repmat(p11,[1,1,ntheta]);
p12=repmat(p12,[1,1,ntheta]);
p2=permute(p2,[3,1,2]);p11=permute(p11,[3,1,2]);p12=permute(p12,[3,1,2]);

index11(p11)=0;
index12(p12)=0;
index=index11+index12;
index(p2)=0;