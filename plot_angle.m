R=DSO;
% [x,y] = pol2cart([beta(40):0.1*pi/180:beta(end)],R);
[x,y] = pol2cart([beta(1):0.1*pi/180:beta(end)],R);
Mx=ceil(max(x));mx=floor(min(x));My=ceil(max(y));my=floor(min(y));
f=zeros([2*R+1,2*R+1]);
% [x_cor,y_cor]=meshgrid([-256:255],[-256,255]);
center=ceil(size(f)/2);
% f([-256:255]+center(1),[-256:255]+center(2))=gather(ls1(rf3));
f([-256:255]+center(1),[-256:255]+center(2))=gather(ls1(rf1));
plot(round(x),round(y));
xx=round(x)+center(1);
yy=round(y)+center(2);
f=repmat(f,[1,1,3]);
for i=1:length(xx)
f(xx(i),yy(i),1)=10000;
end
figure;imshow(f,[])