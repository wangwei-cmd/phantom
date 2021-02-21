% test_panthom_Noo
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

d=6*pi/180;
W=weight_ace(beta,alpha_cor,d);
s_fan=s_fan.*W;
rf=back_ace2_anti(s_fan,alpha_cor,x_cor,y_cor,R,beta)*2;
rf=gather(rf);
f=gather(f);

figure;imshow(rf,[]);
figure;imshow(f,[]);
psnr(rf,f,max(f(:)))
ssim(uint8(rf/max(f(:))*255),uint8(f/max(f(:))*255))
imwrite(uint8(rf/max(f(:))*255),'./rf/ace.png')