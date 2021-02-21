% test_panthom_index
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
% figure;imshow(image,[]);
% figure;imshow(sino,[]);

theta1=solve_chord(beta(1),x_cor,y_cor,R);
[index1]=index_theta(beta(1)*ones(size(theta1)),theta1,beta);

theta2=solve_chord(beta(end),x_cor,y_cor,R);
[index2]=index_theta_1(theta2,beta(end)*ones(size(theta2)),beta);

rf3=back_ace2_index_anti(s_fan,(index1+index2)/2,alpha_cor,x_cor,y_cor,R,beta);
rf3=gather(rf3);
f=gather(f);
figure;imshow(rf3,[]);
figure;imshow(f,[]);
psnr(rf3,f,max(f(:)))
ssim(uint8(rf3/max(f(:))*255),uint8(f/max(f(:))*255))
imwrite(uint8(rf3/max(f(:))*255),'./rf/our.png')
imwrite(uint8(f/max(f(:))*255),'./rf/orig.png')