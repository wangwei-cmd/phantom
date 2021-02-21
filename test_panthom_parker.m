% test_phantom_parker
clear
[phantom] = analytical_phantom(1, 1);
x_cor = ((0:511)-255.5)*0.055;
y_cor = ((0:511)-255.5)*0.055;
xcoord = ones(512,1)*x_cor;
ycoord = transpose(y_cor)*ones(1,512);
[f] = discrete_phantom(xcoord, ycoord, phantom);

R=500*0.055;
alpha_cor = [-35:0.1:35]*pi/180;
beta = [0:0.25:180+70]*pi/180;
M=length(beta);N=length(alpha_cor);
alphac = ones(M,1)*alpha_cor;
lambdac = transpose(beta)*ones(1,N);
% [alphac,lambdac]=meshgrid(alpha,lambda);
scoord = R*sin(alphac);
theta = lambdac+pi/2-alphac;
[s_fan] = line_integrals(scoord, theta, phantom);

alpha_thres=(beta(end)-pi)/2;
% alpha_thres=40*pi/180;
weight=parker_weight(beta,alpha_cor,alpha_thres);
rf1=backproj(s_fan,weight,beta,alpha_cor,R,x_cor,y_cor);
rf1=gather(rf1);
f=gather(f);
figure;imshow(rf1,[]);
figure;imshow(f,[]);
psnr(rf1,f,max(f(:)))
ssim(uint8(rf1/max(f(:))*255),uint8(f/max(f(:))*255))
imwrite(uint8(rf1/max(f(:))*255),'./rf/Parker.png')

function [rf]=backproj(g,weight,beta,alpha_cor,R,x_cor,y_cor)
xcor=gpuArray(x_cor);
ycor=gpuArray(y_cor);
[x_cor,y_cor]=meshgrid(xcor,ycor);
x_cor=x_cor';y_cor=y_cor';
% [alpha_star,~]=compute_alpha_antirotation(beta,xcor,ycor,R);
[alpha_star,~]=compute_alpha(beta,xcor,ycor,R);
g1=filterg(g.*weight,alpha_cor);
a1=R*cos(beta);
fx=length(xcor);fy=length(ycor);
a1=repmat(a1,[fx,1,fy]);
a1=permute(a1,[1,3,2]);
a2=R*sin(beta);
a2=repmat(a2,[fx,1,fy]);
a2=permute(a2,[1,3,2]);
dist=R./((a1-x_cor).^2+(a2-y_cor).^2);
dist=permute(dist,[3,1,2]);
g_interp=zeros(length(beta),fx,fy);
g_interp=gpuArray(g_interp);
for i=1:length(beta)
g_interp(i,:,:)=interp1(alpha_cor,g1(i,:),alpha_star(i,:,:));
end
g_interp(isnan(g_interp))=0;
tmp=g_interp.*dist;
rf=sum(tmp,1)*(beta(2)-beta(1));
rf=2*squeeze(rf);
rf=rot90(rf,-1);
rf=fliplr(rf);
end

function weight=parker_weight(beta,alpha_cor,alpha_thre)
[cor_beta,cor_alpha]=meshgrid(beta,alpha_cor);
cor_beta=cor_beta';cor_alpha=cor_alpha';
alpha_1=2*(alpha_thre+cor_alpha);
alpha_2=pi+2*cor_alpha;
alpha_3=pi+2*alpha_thre;
id1=(cor_beta>=0)&(cor_beta<alpha_1);
id2=(cor_beta>=alpha_1)&(cor_beta<alpha_2);
id3=(cor_beta>=alpha_2)&(cor_beta<alpha_3);
id4=(cor_beta>=alpha_3)&(cor_beta<=2*pi);
weight=zeros(size(id1));
weight=gpuArray(weight);
tmp1=(sin(pi/4*cor_beta./(alpha_thre+cor_alpha))).^2;
weight(id1)=tmp1(id1);
weight(id2)=1;
tmp3=(sin(pi/4*(pi+2*alpha_thre-cor_beta)./(alpha_thre-cor_alpha))).^2; 
weight(id3)=tmp3(id3);
weight(id4)=0;
end

function v=filterg(g,alpha)

cos_alpha=cos(alpha);
sin_alpha=sin(alpha);
h=alpha(2)-alpha(1);

% c=(R+rvk)*pi;
% w_c_sin=w_function(c,sin_alpha);

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
u=gpuArray(u);
    u(s==0)=1/2;
    v=s(s~=0);
    u(s~=0)=(cos(v)-1)./(v.^2)+sin(v)./v;

end


function w_b=w_function(b,s)
w_b=u_function(b*s)*b^2/(4*pi^2);
end