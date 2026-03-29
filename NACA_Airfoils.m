function [x_b, y_b, x, yc] = NACA_Airfoils(m,p,t,c,N)

%cosine spacing 
theta=linspace(-pi,pi,N+1);
x=(c/2)*(1-cos(theta));

%Thicknesss distribution 
yt= (c*t/0.2)*((0.2969*sqrt(x/c))- (0.126*(x/c))-(0.3516*(x/c).^2)+(0.2843*(x/c).^3)-(0.1036*(x/c).^4));

% camber line and slope as vectors
if p==0
yc= zeros(size(x));
dycdx= zeros(size(x));
else
    x_cam= x/c;
    yc=(m/p^2)*(2*p*x_cam-x_cam.^2).*(x_cam <=p)+ (m/(1-p)^2)*((1-2*p)+ 2*p*x_cam- x_cam.^2) .*(x_cam>p);
    dycdx= (2*m/p^2)*(p-x_cam).*(x_cam <= p)+(2*m/(1-p)^2)*(p-x_cam).*(x_cam >p);
end

% Angle
Xi= atan(dycdx);

% Upper Surface
xU=x-yt .*sin(Xi);
yU=yc+yt .*cos(Xi);

%Lower Surface
xL= x+yt .*sin(Xi);
yL=yc-yt.*cos(Xi);

%% Combine surfaces into one (clockwise)
% Upper surface Ttailing edge to leading edge
xU=flip(xU);
yU=flip(yU);

% Lower surface leading edge to trailing edge
xL= xL(2:end);   %skip 1st point of the lower surface so doesn't duplicate Leaing edge
yL= yL(2:end);

% Final boundary points
x_b= [xU,xL];
y_b= [yU,yL];
end

