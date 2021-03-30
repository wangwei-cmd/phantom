function theta1=solve_chord(theta,x,y,R)
% y_s=[R*cos(theta),R*sin(theta)];

cos_theta1=(R*(-R^2 + x.^2 + y.^2).*cos(theta) - ...
 2*x.*(-R^2 + x.*R*cos(theta) + y.*R*sin(theta)))./(R*(R^2 + x.^2 + y.^2 - ...
   2*x.*R*cos(theta) - 2*y.*R*sin(theta)));
sin_theta1=(-2*y.*(-R^2 + x.*R*cos(theta) + y.*R*sin(theta)) + ...
 R.*(-R^2 + x.^2 + y.^2).*sin(theta))./(R*(R^2 + x.^2 + y.^2 - ...
   2*x.*R*cos(theta) - 2*y.*R*sin(theta)));

% t1=(x.^2+y.^2+R^2-2*R*(x.*cos(theta)+y.*sin(theta))).^(1/2);
% cos_theta1=(x-R*cos(theta))./t1;
% sin_theta1=(y-R*sin(theta))./t1;


theta1=acos(cos_theta1);
theta1=real(theta1);
theta1(sin_theta1<0)=2*pi-theta1(sin_theta1<0);
theta1=mod(theta1,2*pi);

